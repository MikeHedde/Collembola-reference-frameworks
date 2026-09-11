# Legacy indicator and scoring definitions ------------------------------------
# These functions reproduce the definitions used in H. Suarez's 2025 scripts.
# Reviewer robustness checks must not silently redefine them.

required_packages(c("dplyr", "tidyr", "tibble", "vegan"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(vegan)
})

FRAMEWORKS <- c(
  "RMQS_Biodiversite",
  "RMQS_BioDiv_Bretagne",
  "ANDRA",
  "TIGA_rural"
)

FRAMEWORK_LABELS <- c(
  RMQS_Biodiversite = "RMQS Biodiversite",
  RMQS_BioDiv_Bretagne = "RMQS BioDiv Bretagne",
  ANDRA = "ANDRA",
  TIGA_rural = "TIGA rural"
)

# All metrics computed in the historical analysis are retained for provenance.
LEGACY_INDICATORS <- c("RSr", "Shannon", "Pielou", "ab")
# The submitted manuscript analyses/discusses rarefied richness and density.
MANUSCRIPT_INDICATORS <- c("RSr", "ab")

N_SCORE_CLASSES <- 7L
RAREFACTION_N <- 500L
AGRICULTURAL_CLC1 <- 2

get_taxon_columns <- function(df) {
  anchor <- match("CLC_niveau3", names(df))
  if (is.na(anchor) || anchor == ncol(df)) {
    stop("Could not identify taxon columns: 'CLC_niveau3' is missing or is the last column.")
  }
  names(df)[seq.int(anchor + 1L, ncol(df))]
}

prepare_taxon_matrix <- function(df, taxon_cols = get_taxon_columns(df)) {
  x <- as.data.frame(df[, taxon_cols, drop = FALSE])
  x[] <- lapply(x, function(v) {
    v <- suppressWarnings(as.numeric(v))
    v[is.na(v)] <- 0
    round(v)
  })
  x
}

compute_legacy_indicators <- function(df, rarefaction_n = RAREFACTION_N) {
  taxon_cols <- get_taxon_columns(df)
  comm_all <- prepare_taxon_matrix(df, taxon_cols)

  # Historical rule from Analyse statistique.R: a space in the taxon column name
  # is used as the operational filter for genus/species-level identifications.
  strict_cols <- taxon_cols[grepl(" ", taxon_cols, fixed = TRUE)]
  comm_strict <- comm_all[, strict_cols, drop = FALSE]

  if (!length(strict_cols)) stop("No strict taxon columns identified by historical name rule.")

  total_strict <- rowSums(comm_strict)
  below_rarefaction <- total_strict < rarefaction_n

  # vegan::rarefy returns NaN/NA where the requested sample exceeds abundance;
  # preserve that behaviour rather than silently changing the estimator.
  rsr <- suppressWarnings(vegan::rarefy(comm_strict, sample = rarefaction_n))
  shannon <- vegan::diversity(comm_strict, index = "shannon")

  df %>%
    mutate(
      RSr = as.numeric(rsr),
      ab = rowSums(comm_all),
      Shannon = as.numeric(shannon),
      Pielou = ifelse(RSr <= 1 | is.na(RSr), NA_real_, Shannon / log(RSr)),
      .strict_total = total_strict,
      .below_rarefaction_n = below_rarefaction
    )
}

score_breaks <- function(x, nb_class = N_SCORE_CLASSES) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) stop("At least two finite reference values are required to build scores.")
  if (nb_class < 3L) stop("nb_class must be >= 3.")

  nb_internal <- nb_class - 2L
  internal_breaks <- as.numeric(stats::quantile(
    x,
    probs = seq(0, 1, length.out = nb_internal + 1L),
    na.rm = TRUE,
    names = FALSE,
    type = 7
  ))
  c(-Inf, internal_breaks, Inf)
}

score_values <- function(x, breaks) {
  # Match cut(..., include.lowest=TRUE, right=FALSE) used by GiveScores.R.
  out <- rep(NA_integer_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  out[ok] <- findInterval(x[ok], breaks, rightmost.closed = FALSE, all.inside = TRUE) - 1L
  out
}

breaks_to_table <- function(breaks, framework, indicator, scenario) {
  internal <- breaks[is.finite(breaks)]
  probs <- seq(0, 1, length.out = length(internal))
  tibble(
    scenario = scenario,
    framework = framework,
    indicator = indicator,
    boundary = paste0("q", sprintf("%02d", round(100 * probs))),
    probability = probs,
    value = internal
  )
}

build_break_table <- function(data, frameworks = FRAMEWORKS,
                              indicators = MANUSCRIPT_INDICATORS,
                              scenario = "full_untrimmed") {
  out <- list()
  k <- 1L
  for (fw in frameworks) {
    ref_fw <- data[data$Projet == fw, , drop = FALSE]
    for (ind in indicators) {
      out[[k]] <- breaks_to_table(score_breaks(ref_fw[[ind]]), fw, ind, scenario)
      k <- k + 1L
    }
  }
  bind_rows(out)
}

get_breaks_from_table <- function(tbl, framework, indicator) {
  # IMPORTANT: framework and indicator are also column names in tbl. Using
  # .env$... prevents dplyr's data mask from resolving the RHS to the column
  # itself (which would make every row TRUE and mix all frameworks together).
  vals <- tbl %>%
    filter(
      .data$framework == .env$framework,
      .data$indicator == .env$indicator
    ) %>%
    arrange(.data$probability) %>%
    pull("value")

  if (!length(vals)) stop("No breaks found for ", framework, " / ", indicator)
  if (anyNA(vals) || any(!is.finite(vals))) {
    stop("Non-finite internal breaks found for ", framework, " / ", indicator)
  }
  if (is.unsorted(vals, strictly = FALSE)) {
    stop(
      "Internal breaks are not sorted for ", framework, " / ", indicator,
      ". Values: ", paste(vals, collapse = ", ")
    )
  }
  c(-Inf, vals, Inf)
}

score_validation_data <- function(validation_data, break_table,
                                  frameworks = unique(break_table$framework),
                                  indicators = unique(break_table$indicator),
                                  scenario_name = "full_untrimmed") {
  site_id <- if ("Site" %in% names(validation_data)) as.character(validation_data$Site) else as.character(seq_len(nrow(validation_data)))
  out <- list()
  k <- 1L
  for (fw in frameworks) {
    for (ind in indicators) {
      br <- get_breaks_from_table(break_table, fw, ind)
      out[[k]] <- tibble(
        site_id = site_id,
        indicator = ind,
        framework = fw,
        scenario = scenario_name,
        raw_value = validation_data[[ind]],
        score = score_values(validation_data[[ind]], br)
      )
      k <- k + 1L
    }
  }
  bind_rows(out)
}

legacy_tiga_ab_breaks <- function(reference_agri, threshold = 75000) {
  # Historical Analyse statistique.R applied the TIGA_rural abundance filter only
  # immediately before *density* scoring. Richness scores had already been built
  # from the unfiltered TIGA reference. Reproduce that sequence exactly.
  full <- build_break_table(
    reference_agri,
    indicators = MANUSCRIPT_INDICATORS,
    scenario = "legacy_tiga_ab_lt_75000"
  )

  tiga_ab <- reference_agri[
    reference_agri$Projet == "TIGA_rural" &
      !(is.finite(reference_agri$ab) & reference_agri$ab >= threshold),
    ,
    drop = FALSE
  ]
  replacement <- breaks_to_table(
    score_breaks(tiga_ab$ab),
    framework = "TIGA_rural",
    indicator = "ab",
    scenario = "legacy_tiga_ab_lt_75000"
  )

  full %>%
    filter(!(.data$framework == "TIGA_rural" & .data$indicator == "ab")) %>%
    bind_rows(replacement) %>%
    arrange(match(.data$framework, FRAMEWORKS), match(.data$indicator, MANUSCRIPT_INDICATORS), .data$probability)
}
