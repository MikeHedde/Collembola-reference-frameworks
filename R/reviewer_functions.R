# Reviewer robustness helper functions ----------------------------------------

required_packages(c("dplyr", "tidyr", "tibble", "ggplot2", "stringr"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(stringr)
})

winsorise <- function(x, probs = c(0.01, 0.99)) {
  q <- stats::quantile(x, probs = probs, na.rm = TRUE, type = 7)
  pmax(pmin(x, q[[2]]), q[[1]])
}

extreme_iqr_flag <- function(x, k = 3) {
  q <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE, type = 7)
  iqr <- q[[2]] - q[[1]]
  x < (q[[1]] - k * iqr) | x > (q[[2]] + k * iqr)
}

balanced_subsampling_breaks <- function(ref_data, n, B, seed,
                                        frameworks = FRAMEWORKS,
                                        indicators = MANUSCRIPT_INDICATORS) {
  set.seed(seed)
  n_by_framework <- table(as.character(ref_data$Projet))
  if (any(n_by_framework[frameworks] < n)) {
    stop("At least one framework has fewer than n = ", n, " eligible sites.")
  }

  result <- vector("list", B * length(frameworks) * length(indicators))
  k <- 1L
  for (b in seq_len(B)) {
    sampled_rows <- lapply(frameworks, function(fw) {
      idx <- which(as.character(ref_data$Projet) == fw)
      sample(idx, size = n, replace = FALSE)
    })
    names(sampled_rows) <- frameworks

    for (fw in frameworks) {
      d <- ref_data[sampled_rows[[fw]], , drop = FALSE]
      for (ind in indicators) {
        tmp <- breaks_to_table(score_breaks(d[[ind]]), fw, ind, "balanced_n58")
        tmp$replicate <- b
        result[[k]] <- tmp
        k <- k + 1L
      }
    }
  }
  bind_rows(result)
}

pairwise_score_agreement <- function(score_long) {
  frameworks <- unique(score_long$framework)
  pairs <- combn(frameworks, 2, simplify = FALSE)
  out <- list()
  k <- 1L

  for (ind in unique(score_long$indicator)) {
    dat_ind <- score_long %>% filter(.data$indicator == ind)
    for (pair in pairs) {
      wide <- dat_ind %>%
        filter(.data$framework %in% pair) %>%
        select(site_id, framework, score) %>%
        distinct() %>%
        pivot_wider(names_from = "framework", values_from = "score")

      a <- wide[[pair[1]]]
      b <- wide[[pair[2]]]
      ok <- complete.cases(a, b)
      a <- a[ok]
      b <- b[ok]

      broad_a <- cut(a, breaks = c(-Inf, 2, 4, Inf), labels = c("0-2", "3-4", "5-6"))
      broad_b <- cut(b, breaks = c(-Inf, 2, 4, Inf), labels = c("0-2", "3-4", "5-6"))

      out[[k]] <- tibble(
        indicator = ind,
        framework_1 = pair[1],
        framework_2 = pair[2],
        n = length(a),
        spearman_rho = if (length(a) >= 3L) suppressWarnings(stats::cor(a, b, method = "spearman")) else NA_real_,
        exact_agreement = if (length(a)) mean(a == b) else NA_real_,
        within_one_class = if (length(a)) mean(abs(a - b) <= 1) else NA_real_,
        broad_class_agreement = if (length(a)) mean(broad_a == broad_b) else NA_real_
      )
      k <- k + 1L
    }
  }
  bind_rows(out)
}

build_balanced_pooled_breaks <- function(ref_data, n, B, seed,
                                         indicators = MANUSCRIPT_INDICATORS) {
  set.seed(seed)
  boot <- vector("list", B * length(indicators))
  k <- 1L
  for (b in seq_len(B)) {
    idx <- unlist(lapply(FRAMEWORKS, function(fw) {
      candidate <- which(as.character(ref_data$Projet) == fw)
      sample(candidate, n, replace = FALSE)
    }))
    d <- ref_data[idx, , drop = FALSE]
    for (ind in indicators) {
      tmp <- breaks_to_table(score_breaks(d[[ind]]), "POOLED_BALANCED", ind, "pooled_balanced")
      tmp$replicate <- b
      boot[[k]] <- tmp
      k <- k + 1L
    }
  }
  boot <- bind_rows(boot)
  summary <- boot %>%
    group_by(.data$framework, .data$indicator, .data$boundary, .data$probability, .data$scenario) %>%
    summarise(
      value = median(.data$value, na.rm = TRUE),
      q025 = quantile(.data$value, 0.025, na.rm = TRUE),
      q975 = quantile(.data$value, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  list(draws = boot, summary = summary)
}

# ---- Explicit Bioindicateur2 control mapping ---------------------------------

read_control_mapping <- function(path) {
  if (!file.exists(path)) return(tibble())
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("Nom_site", "control_Site", "control_Traitement", "control_Contamination", "include")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Control mapping is missing columns: ", paste(missing, collapse = ", "))
  x$include <- tolower(trimws(as.character(x$include))) %in% c("true", "t", "1", "yes", "y", "oui")
  x
}

control_mapping_is_ready <- function(mapping) {
  if (!nrow(mapping)) return(FALSE)
  any(mapping$include & nzchar(trimws(as.character(mapping$Nom_site))))
}

write_observed_control_levels <- function(validation_data, path) {
  needed <- intersect(c("Nom_site", "Traitement", "Contamination", "Usage", "Site"), names(validation_data))
  if (!length(needed)) return(invisible(NULL))
  x <- validation_data %>%
    select(all_of(needed)) %>%
    distinct() %>%
    arrange(across(any_of(c("Nom_site", "Traitement", "Contamination", "Site"))))
  write_csv_safe(x, path)
}

match_mapping_value <- function(x, target) {
  target <- trimws(as.character(target))
  if (is.na(target) || !nzchar(target) || target == "*") return(rep(TRUE, length(x)))
  trimws(as.character(x)) == target
}

make_explicit_lrr_table <- function(validation_data, score_long, mapping,
                                    indicators = MANUSCRIPT_INDICATORS) {
  if (!control_mapping_is_ready(mapping)) {
    return(list(lrr = tibble(), diagnostics = tibble(status = "mapping_not_ready")))
  }
  if (!all(c("Nom_site", "Traitement", "Site") %in% names(validation_data))) {
    stop("Bioindicateur2 must contain Nom_site, Traitement and Site for explicit control mapping.")
  }

  out <- list()
  diag <- list()
  k <- 1L
  kd <- 1L

  for (i in seq_len(nrow(mapping))) {
    m <- mapping[i, , drop = FALSE]
    if (!isTRUE(m$include)) next
    block <- validation_data[as.character(validation_data$Nom_site) == as.character(m$Nom_site), , drop = FALSE]

    if (!nrow(block)) {
      diag[[kd]] <- tibble(Nom_site = m$Nom_site, status = "site_not_found", n_block = 0L, n_controls = 0L)
      kd <- kd + 1L
      next
    }

    # Identify the control explicitly by Site whenever supplied. Treatment and
    # contamination are retained as cross-checks. This avoids ambiguous labels
    # such as Yvetot, where two rows are both labelled "blé".
    control_flag <- match_mapping_value(block$Site, m$control_Site)
    control_flag <- control_flag & match_mapping_value(block$Traitement, m$control_Traitement)
    if ("Contamination" %in% names(block)) {
      control_flag <- control_flag & match_mapping_value(block$Contamination, m$control_Contamination)
    }
    controls <- block[control_flag, , drop = FALSE]
    treatments <- block[!control_flag, , drop = FALSE]

    if (!nrow(controls) || !nrow(treatments)) {
      diag[[kd]] <- tibble(
        Nom_site = m$Nom_site,
        status = if (!nrow(controls)) "no_control_match" else "no_treatment_rows",
        n_block = nrow(block),
        n_controls = nrow(controls)
      )
      kd <- kd + 1L
      next
    }

    for (ind in indicators) {
      raw_control <- mean(controls[[ind]], na.rm = TRUE)
      control_keys <- as.character(controls$Site)

      for (j in seq_len(nrow(treatments))) {
        raw_treatment <- treatments[[ind]][j]
        lrr <- if (is.finite(raw_treatment) && is.finite(raw_control) && raw_treatment > 0 && raw_control > 0) {
          log(raw_treatment / raw_control)
        } else NA_real_

        treatment_key <- as.character(treatments$Site[j])
        for (fw in FRAMEWORKS) {
          ts <- score_long %>%
            filter(.data$site_id == treatment_key, .data$indicator == ind, .data$framework == fw) %>%
            pull(score)
          cs <- score_long %>%
            filter(.data$site_id %in% control_keys, .data$indicator == ind, .data$framework == fw) %>%
            summarise(v = mean(.data$score, na.rm = TRUE)) %>%
            pull(v)

          out[[k]] <- tibble(
            Nom_site = as.character(m$Nom_site),
            treatment_site = treatment_key,
            treatment = as.character(treatments$Traitement[j]),
            contamination = if ("Contamination" %in% names(treatments)) as.character(treatments$Contamination[j]) else NA_character_,
            control_site = as.character(m$control_Site),
            control_treatment = as.character(m$control_Traitement),
            control_contamination = as.character(m$control_Contamination),
            indicator = ind,
            raw_treatment = raw_treatment,
            raw_control = raw_control,
            LRR = lrr,
            framework = fw,
            score_treatment = if (length(ts) == 1L) ts else NA_real_,
            score_control = if (length(cs) == 1L) cs else NA_real_,
            delta_score = if (length(ts) == 1L && length(cs) == 1L) ts - cs else NA_real_
          )
          k <- k + 1L
        }
      }
    }

    diag[[kd]] <- tibble(
      Nom_site = as.character(m$Nom_site),
      status = "used_explicit_mapping",
      n_block = nrow(block),
      n_controls = nrow(controls),
      control_sites = paste(as.character(controls$Site), collapse = ";")
    )
    kd <- kd + 1L
  }

  list(
    lrr = if (length(out)) bind_rows(out) else tibble(),
    diagnostics = if (length(diag)) bind_rows(diag) else tibble()
  )
}
