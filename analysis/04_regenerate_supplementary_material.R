#!/usr/bin/env Rscript

# ==============================================================================
# 04 - Regenerate historical and reviewer-driven Supplementary Material
# ==============================================================================
# Purpose
#   A. Rebuild the Supplementary Materials (SM1-SM11) that accompanied the
#      submitted manuscript, using the cleaned current workflow as source of
#      truth rather than opaque historical intermediate files.
#   B. Generate candidate supplementary figures/tables for the robustness
#      analyses requested by reviewers.
#
# Historical supplement structure in the submitted manuscript
#   SM1  Site counts by land-cover category
#   SM2  Four indicator distributions - all sites
#   SM3  Descriptive statistics - all sites
#   SM4  Pairwise Cramer-von Mises tests - all sites
#   SM5  Four indicator distributions - agricultural sites
#   SM6  Descriptive statistics - agricultural sites
#   SM7  Pairwise Cramer-von Mises tests - agricultural sites
#   SM8  Illustration of the score-class construction
#   SM9  Framework-specific scoring schemes - agricultural sites
#   SM10 Bioindicateur2 score distributions - four indicators
#   SM11 Friedman + Nemenyi tests on Bioindicateur2 scores
#
# Revision rule
#   The current reproducible pipeline is authoritative when it differs from an
#   historical intermediate. The ad hoc TIGA abundance < 75000 filter is NOT
#   used in revised Supplementary Materials. Its effect remains documented in
#   output/baseline and output/reviewer.
# ============================================================================== 

# Locate project root -----------------------------------------------------------
.project_paths_candidates <- c(file.path(getwd(), "R", "project_paths.R"))
.args_all <- commandArgs(trailingOnly = FALSE)
.file_arg <- grep("^--file=", .args_all, value = TRUE)
if (length(.file_arg)) {
  .script_file <- normalizePath(sub("^--file=", "", .file_arg[1]), mustWork = FALSE)
  .project_paths_candidates <- c(
    .project_paths_candidates,
    file.path(dirname(.script_file), "..", "R", "project_paths.R")
  )
}
.ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (!is.null(.ofile)) {
  .project_paths_candidates <- c(
    .project_paths_candidates,
    file.path(dirname(normalizePath(.ofile, mustWork = FALSE)), "..", "R", "project_paths.R")
  )
}
.project_paths_candidates <- unique(normalizePath(.project_paths_candidates, mustWork = FALSE))
.project_paths_file <- .project_paths_candidates[file.exists(.project_paths_candidates)][1]
if (is.na(.project_paths_file)) {
  stop("Could not locate R/project_paths.R. Open the .Rproj or run inside the repository.")
}
source(.project_paths_file)
rm(.project_paths_candidates, .args_all, .file_arg, .ofile, .project_paths_file)
if (exists(".script_file")) rm(.script_file)

source_repo("R/io_helpers.R")
source_repo("R/data_preparation.R")
source_repo("R/indicator_scoring.R")
source_repo("R/reviewer_functions.R")

required_packages(c(
  "dplyr", "tidyr", "tibble", "ggplot2", "cowplot", "scales", "vegan",
  "e1071", "cramer", "PMCMRplus"
))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(cowplot)
})

dir.create(SUPPLEMENTARY_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(SUPPLEMENTARY_OUTPUT_DIR, "reviewer_candidates"), recursive = TRUE, showWarnings = FALSE)

message("[04] Regenerating Supplementary Material...")

# Load/rebuild analysis-ready data ---------------------------------------------
ref_rds <- file.path(DERIVED_DIR, "reference_frameworks_merged.rds")
bio_rds <- file.path(DERIVED_DIR, "bioindicateur2_curated.rds")
if (!file.exists(ref_rds) || !file.exists(bio_rds)) {
  message("[04] Derived datasets absent; running analysis/01_build_analysis_dataset.R first.")
  source_repo("analysis/01_build_analysis_dataset.R")
}

reference_raw <- readRDS(ref_rds)
validation_raw <- readRDS(bio_rds)
reference <- compute_legacy_indicators(reference_raw)
validation <- compute_legacy_indicators(validation_raw)

reference_all <- reference %>%
  filter(.data$Projet %in% FRAMEWORKS) %>%
  mutate(Projet = factor(as.character(.data$Projet), levels = rev(FRAMEWORKS)))

reference_agri <- reference %>%
  filter(
    .data$Projet %in% FRAMEWORKS,
    as.character(.data$CLC_niveau1) == as.character(AGRICULTURAL_CLC1)
  ) %>%
  mutate(Projet = factor(as.character(.data$Projet), levels = rev(FRAMEWORKS)))

validation_agri <- validation %>%
  filter(as.character(.data$CLC_niveau1) == as.character(AGRICULTURAL_CLC1))

framework_labels <- c(
  RMQS_Biodiversite = "RMQS Biodiversity",
  RMQS_BioDiv_Bretagne = "RMQS BioDiv Bretagne",
  ANDRA = "ANDRA",
  TIGA_rural = "TIGA rural"
)

framework_colors <- c(
  RMQS_Biodiversite = "#56B4E9",
  RMQS_BioDiv_Bretagne = "#009E73",
  ANDRA = "#F0E442",
  TIGA_rural = "#D55E00"
)

score_colors <- c(
  "0" = "#8B4513",
  "1" = "#D73027",
  "2" = "#FC8D59",
  "3" = "#FEE08B",
  "4" = "#D9EF8B",
  "5" = "#91BFDB",
  "6" = "#4575B4"
)

indicator_labels <- c(
  RSr = "Species richness",
  Shannon = "Shannon index",
  Pielou = "Pielou index",
  ab = "Density per m²"
)

# Export helpers ---------------------------------------------------------------
save_plot_pair <- function(plot, stem, width, height, dpi = 300, directory = SUPPLEMENTARY_OUTPUT_DIR) {
  pdf_path <- file.path(directory, paste0(stem, ".pdf"))
  png_path <- file.path(directory, paste0(stem, ".png"))
  ggsave(pdf_path, plot = plot, width = width, height = height, units = "in")
  ggsave(png_path, plot = plot, width = width, height = height, units = "in", dpi = dpi)
  invisible(c(pdf = pdf_path, png = png_path))
}

# Descriptive statistics -------------------------------------------------------
legacy_descriptive_stats <- function(data, indicators = LEGACY_INDICATORS) {
  bind_rows(lapply(FRAMEWORKS, function(fw) {
    d <- data %>% filter(as.character(.data$Projet) == fw)
    bind_rows(lapply(indicators, function(ind) {
      x <- d[[ind]]
      tibble(
        framework = fw,
        indicator = ind,
        n = sum(is.finite(x)),
        mean = mean(x, na.rm = TRUE),
        sd = stats::sd(x, na.rm = TRUE),
        skewness = e1071::skewness(x, na.rm = TRUE),
        # e1071::kurtosis() already returns excess kurtosis (normal = 0).
        # The historical script subtracted 3 a second time. We retain that
        # legacy value as an audit field but report the corrected statistic.
        kurtosis = e1071::kurtosis(x, na.rm = TRUE),
        kurtosis_historical_double_subtracted = e1071::kurtosis(x, na.rm = TRUE) - 3
      )
    }))
  }))
}

# Pairwise Cramer-von Mises + Bonferroni --------------------------------------
pairwise_cramer <- function(data, indicator, seed = 1234, replicates = 10000) {
  groups <- FRAMEWORKS
  pairs <- utils::combn(groups, 2, simplify = FALSE)
  set.seed(seed)
  raw <- bind_rows(lapply(seq_along(pairs), function(k) {
    pair <- pairs[[k]]
    d1 <- data %>% filter(as.character(.data$Projet) == pair[1])
    d2 <- data %>% filter(as.character(.data$Projet) == pair[2])
    x1 <- d1[[indicator]]
    x2 <- d2[[indicator]]
    x1 <- x1[is.finite(x1)]
    x2 <- x2[is.finite(x2)]
    tst <- cramer::cramer.test(x1, x2, replicates = replicates, sim = "ordinary")
    tibble(
      indicator = indicator,
      framework_1 = pair[1],
      framework_2 = pair[2],
      W = as.numeric(tst$statistic),
      p_raw = as.numeric(tst$p.value)
    )
  }))
  raw %>%
    mutate(
      p_bonferroni = pmin(1, stats::p.adjust(.data$p_raw, method = "bonferroni")),
      significant_0_05 = .data$p_bonferroni < 0.05
    )
}

pairwise_cramer_all <- function(data, indicators = LEGACY_INDICATORS, seed = 1234, replicates = 10000) {
  bind_rows(lapply(indicators, function(ind) pairwise_cramer(data, ind, seed = seed, replicates = replicates)))
}

# Compact-letter display without an additional R dependency -------------------
# We enumerate maximal cliques of the graph of non-significant pairwise
# comparisons. Groups sharing a clique share a letter. Significant pairs never
# share a letter. For four frameworks this is transparent and deterministic.
compact_letters_from_pairs <- function(pair_tbl, alpha = 0.05, group_order = FRAMEWORKS,
                                       p_col = "p_bonferroni") {
  g <- group_order
  n <- length(g)
  nonsig <- diag(TRUE, n)
  dimnames(nonsig) <- list(g, g)

  for (i in seq_len(nrow(pair_tbl))) {
    a <- as.character(pair_tbl$framework_1[i])
    b <- as.character(pair_tbl$framework_2[i])
    p <- pair_tbl[[p_col]][i]
    keep <- is.finite(p) && p >= alpha
    nonsig[a, b] <- keep
    nonsig[b, a] <- keep
  }

  subsets <- unlist(lapply(seq_len(n), function(k) utils::combn(g, k, simplify = FALSE)), recursive = FALSE)
  is_clique <- vapply(subsets, function(s) {
    if (length(s) <= 1) return(TRUE)
    all(nonsig[s, s])
  }, logical(1))
  cliques <- subsets[is_clique]

  is_maximal <- vapply(seq_along(cliques), function(i) {
    s <- cliques[[i]]
    !any(vapply(seq_along(cliques), function(j) {
      if (i == j) return(FALSE)
      t <- cliques[[j]]
      length(t) > length(s) && all(s %in% t)
    }, logical(1)))
  }, logical(1))
  cliques <- cliques[is_maximal]

  clique_key <- vapply(cliques, function(s) {
    idx <- match(s, g)
    paste(sprintf("%02d", sort(idx)), collapse = "-")
  }, character(1))
  cliques <- cliques[order(clique_key)]

  letter_pool <- c(letters, paste0("a", letters), paste0("b", letters))
  if (length(cliques) > length(letter_pool)) stop("Too many CLD cliques for available labels.")
  clique_letters <- letter_pool[seq_along(cliques)]

  out <- setNames(rep("", n), g)
  for (i in seq_along(cliques)) {
    out[cliques[[i]]] <- paste0(out[cliques[[i]]], clique_letters[i])
  }
  out
}

# Distribution figures --------------------------------------------------------
make_distribution_panel <- function(data, indicator, cramer_tbl, seed = 123) {
  # Use .env$indicator below: the function argument has the same name as a data column.
  d <- data %>%
    transmute(
      framework = factor(as.character(.data$Projet), levels = rev(FRAMEWORKS)),
      value = .data[[indicator]]
    )

  letters_this <- compact_letters_from_pairs(
    cramer_tbl %>% filter(.data$indicator == .env$indicator),
    group_order = FRAMEWORKS
  )
  ann <- tibble(
    framework = factor(names(letters_this), levels = rev(FRAMEWORKS)),
    letter = unname(letters_this)
  )
  rng <- range(d$value[is.finite(d$value)], na.rm = TRUE)
  x_ann <- rng[2] + 0.06 * diff(rng)
  if (!is.finite(x_ann) || diff(rng) == 0) x_ann <- rng[2]
  ann$x <- x_ann

  set.seed(seed)
  p <- ggplot(d, aes(x = .data$value, y = .data$framework, fill = .data$framework)) +
    geom_violin(trim = TRUE, alpha = 0.70, na.rm = TRUE) +
    geom_boxplot(width = 0.10, fill = "white", outlier.shape = NA, na.rm = TRUE) +
    geom_point(
      position = position_jitter(width = 0, height = 0.16, seed = seed),
      alpha = 0.50, size = 0.90, na.rm = TRUE
    ) +
    geom_text(
      data = ann,
      aes(x = .data$x, y = .data$framework, label = .data$letter),
      inherit.aes = FALSE, size = 3.5
    ) +
    scale_fill_manual(values = framework_colors, guide = "none") +
    scale_y_discrete(labels = framework_labels) +
    labs(x = unname(indicator_labels[indicator]), y = NULL) +
    theme_minimal(base_size = 9.5) +
    theme(panel.grid.minor = element_blank())

  if (indicator == "ab") {
    p <- p + scale_x_continuous(labels = scales::label_number(big.mark = " "), expand = expansion(mult = c(0.02, 0.10)))
  } else {
    p <- p + scale_x_continuous(expand = expansion(mult = c(0.02, 0.10)))
  }
  p
}

make_distribution_figure <- function(data, cramer_tbl, stem) {
  p1 <- make_distribution_panel(data, "RSr", cramer_tbl)
  p2 <- make_distribution_panel(data, "Shannon", cramer_tbl)
  p3 <- make_distribution_panel(data, "Pielou", cramer_tbl)
  p4 <- make_distribution_panel(data, "ab", cramer_tbl)
  combined <- cowplot::plot_grid(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), ncol = 2)
  save_plot_pair(combined, stem, width = 11.2, height = 8.4)
}

# ==============================================================================
# SM1 - site/observation counts by land-cover category
# ==============================================================================
clc_labels <- c(
  `1` = "Artificialized",
  `2` = "Agricultural",
  `3` = "Forests and semi-natural environments",
  `4` = "Wetlands"
)

sm1_long <- reference %>%
  filter(.data$Projet %in% FRAMEWORKS) %>%
  mutate(
    framework = as.character(.data$Projet),
    CLC_niveau1_chr = as.character(.data$CLC_niveau1),
    land_cover = unname(clc_labels[.data$CLC_niveau1_chr])
  ) %>%
  filter(.data$CLC_niveau1_chr %in% names(clc_labels)) %>%
  count(.data$land_cover, .data$framework, name = "n") %>%
  tidyr::complete(
    land_cover = unname(clc_labels),
    framework = FRAMEWORKS,
    fill = list(n = 0L)
  )

sm1_wide <- sm1_long %>%
  mutate(
    land_cover = factor(.data$land_cover, levels = unname(clc_labels)),
    framework = factor(.data$framework, levels = FRAMEWORKS)
  ) %>%
  arrange(.data$land_cover, .data$framework) %>%
  pivot_wider(names_from = "framework", values_from = "n")

write_csv_safe(sm1_long, file.path(SUPPLEMENTARY_OUTPUT_DIR, "SM01_site_counts_by_land_cover_long.csv"))
write_csv_safe(sm1_wide, file.path(SUPPLEMENTARY_OUTPUT_DIR, "SM01_site_counts_by_land_cover.csv"))

sm1_missing <- reference %>%
  filter(.data$Projet %in% FRAMEWORKS) %>%
  mutate(framework = as.character(.data$Projet)) %>%
  group_by(.data$framework) %>%
  summarise(
    n_total = n(),
    n_missing_or_other_CLC1 = sum(is.na(.data$CLC_niveau1) | !as.character(.data$CLC_niveau1) %in% names(clc_labels)),
    .groups = "drop"
  )
write_csv_safe(sm1_missing, file.path(SUPPLEMENTARY_OUTPUT_DIR, "SM01_audit_missing_land_cover.csv"))

# ==============================================================================
# SM2-SM4 - all-site distributions, descriptive stats, Cramer tests
# ==============================================================================
message("[04] SM2-SM4: all-site distributions and tests...")
sm4_cramer_all <- pairwise_cramer_all(reference_all, LEGACY_INDICATORS, seed = 1234, replicates = 10000)
sm3_stats_all <- legacy_descriptive_stats(reference_all)

write_csv_safe(sm3_stats_all, file.path(SUPPLEMENTARY_OUTPUT_DIR, "SM03_descriptive_statistics_all_sites.csv"))
write_csv_safe(sm4_cramer_all, file.path(SUPPLEMENTARY_OUTPUT_DIR, "SM04_cramer_von_mises_all_sites.csv"))
make_distribution_figure(reference_all, sm4_cramer_all, "SM02_indicator_distributions_all_sites")

# ==============================================================================
# SM5-SM7 - agricultural-site distributions, descriptive stats, Cramer tests
# ==============================================================================
message("[04] SM5-SM7: agricultural distributions and tests...")
sm7_cramer_agri <- pairwise_cramer_all(reference_agri, LEGACY_INDICATORS, seed = 1234, replicates = 10000)
sm6_stats_agri <- legacy_descriptive_stats(reference_agri)

write_csv_safe(sm6_stats_agri, file.path(SUPPLEMENTARY_OUTPUT_DIR, "SM06_descriptive_statistics_agricultural_sites.csv"))
write_csv_safe(sm7_cramer_agri, file.path(SUPPLEMENTARY_OUTPUT_DIR, "SM07_cramer_von_mises_agricultural_sites.csv"))
make_distribution_figure(reference_agri, sm7_cramer_agri, "SM05_indicator_distributions_agricultural_sites")

# ==============================================================================
# SM8 - method used to divide a distribution into seven score classes
# ==============================================================================
message("[04] SM8: scoring-method illustration...")
x_extension <- 2
x_vals <- seq(-4 - x_extension, 4 + x_extension, length.out = 1000)
dens_vals <- stats::dnorm(x_vals)
quants <- stats::qnorm(c(0.00001, 0.2, 0.4, 0.6, 0.8, 0.99999))
bounds <- c(min(x_vals), quants, max(x_vals))
bar_df <- tibble(
  xmin = bounds[-length(bounds)],
  xmax = bounds[-1],
  score = factor(0:6, levels = 0:6)
)

area_layers <- lapply(seq_len(nrow(bar_df)), function(i) {
  x_seg <- seq(bar_df$xmin[i], bar_df$xmax[i], length.out = 300)
  y_seg <- stats::dnorm(x_seg)
  geom_area(
    data = data.frame(x = x_seg, y = y_seg),
    aes(x = .data$x, y = .data$y),
    fill = score_colors[as.character(bar_df$score[i])], alpha = 0.20
  )
})

p_method <- ggplot() +
  area_layers +
  geom_line(aes(x = x_vals, y = dens_vals), colour = "black", linewidth = 0.8) +
  geom_rect(
    data = bar_df,
    aes(xmin = .data$xmin, xmax = .data$xmax, ymin = -0.020, ymax = -0.005, fill = .data$score),
    colour = "black", linewidth = 0.25
  ) +
  geom_text(
    data = bar_df,
    aes(x = (.data$xmin + .data$xmax) / 2, y = -0.0125, label = .data$score),
    size = 3.4
  ) +
  scale_fill_manual(values = score_colors, guide = "none") +
  coord_cartesian(xlim = range(x_vals), ylim = c(-0.05, 0.45)) +
  labs(x = "Indicator value", y = "Distribution in the reference framework") +
  theme_minimal(base_size = 10) +
  theme(axis.title.y = element_text(margin = margin(r = 10)))

save_plot_pair(p_method, "SM08_scoring_class_construction", width = 11.0, height = 6.5)
write_csv_safe(bar_df, file.path(SUPPLEMENTARY_OUTPUT_DIR, "SM08_scoring_class_illustration_bounds.csv"))

# ==============================================================================
# SM9 - scoring schemes for each framework and all four historical indicators
# ==============================================================================
message("[04] SM9: framework-specific scoring schemes...")
sm9_breaks <- build_break_table(
  reference_agri,
  frameworks = FRAMEWORKS,
  indicators = LEGACY_INDICATORS,
  scenario = "revised_full_untrimmed"
)
write_csv_safe(sm9_breaks, file.path(SUPPLEMENTARY_OUTPUT_DIR, "SM09_scoring_thresholds_agricultural_sites.csv"))

score_display_caps <- c(RSr = 25, Shannon = 3, Pielou = 1, ab = 110000)

breaks_to_rectangles <- function(break_table, indicator) {
  bind_rows(lapply(FRAMEWORKS, function(fw) {
    vals <- break_table %>%
      filter(.data$framework == .env$fw, .data$indicator == .env$indicator) %>%
      arrange(.data$probability) %>%
      pull(value)
    br <- c(-Inf, vals, Inf)
    tibble(
      framework = fw,
      score = factor(0:6, levels = 0:6),
      xmin = head(br, -1),
      xmax = tail(br, -1)
    ) %>%
      mutate(
        xmin = ifelse(is.infinite(.data$xmin), 0, .data$xmin),
        xmax = ifelse(is.infinite(.data$xmax), score_display_caps[[indicator]], .data$xmax),
        xmax = pmin(.data$xmax, score_display_caps[[indicator]])
      ) %>%
      filter(.data$xmax > .data$xmin)
  }))
}

sm9_rects <- bind_rows(lapply(LEGACY_INDICATORS, function(ind) {
  breaks_to_rectangles(sm9_breaks, ind) %>% mutate(indicator = ind)
}))
write_csv_safe(sm9_rects, file.path(SUPPLEMENTARY_OUTPUT_DIR, "SM09_scoring_scheme_rectangles.csv"))

make_scoring_scheme_panel <- function(indicator, show_legend = FALSE) {
  d <- sm9_rects %>%
    filter(.data$indicator == .env$indicator) %>%
    mutate(framework = factor(.data$framework, levels = rev(FRAMEWORKS)))

  p <- ggplot(d) +
    geom_rect(
      aes(
        xmin = .data$xmin, xmax = .data$xmax,
        ymin = as.numeric(.data$framework) - 0.38,
        ymax = as.numeric(.data$framework) + 0.38,
        fill = .data$score
      ),
      colour = "black", linewidth = 0.20
    ) +
    scale_y_continuous(
      breaks = seq_along(rev(FRAMEWORKS)),
      labels = framework_labels[rev(FRAMEWORKS)]
    ) +
    scale_fill_manual(values = score_colors, name = "Assigned score", drop = FALSE) +
    labs(x = unname(indicator_labels[indicator]), y = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = if (show_legend) "right" else "none"
    )

  if (indicator == "ab") p <- p + scale_x_continuous(labels = scales::label_number(big.mark = " "))
  p
}

p_s9_1 <- make_scoring_scheme_panel("RSr")
p_s9_2 <- make_scoring_scheme_panel("Shannon")
p_s9_3 <- make_scoring_scheme_panel("Pielou")
p_s9_4 <- make_scoring_scheme_panel("ab", show_legend = TRUE)
p_s9 <- cowplot::plot_grid(p_s9_1, p_s9_2, p_s9_3, p_s9_4, labels = c("A", "B", "C", "D"), ncol = 1, rel_heights = c(1, 1, 1, 1.05))
save_plot_pair(p_s9, "SM09_framework_specific_scoring_schemes", width = 10.5, height = 10.0)

# ==============================================================================
# SM10-SM11 - Bioindicateur2 scores and Friedman/Nemenyi tests
# ==============================================================================
message("[04] SM10-SM11: Bioindicateur2 score distributions and tests...")
sm10_scores <- score_validation_data(
  validation_agri,
  sm9_breaks,
  frameworks = FRAMEWORKS,
  indicators = LEGACY_INDICATORS,
  scenario_name = "revised_full_untrimmed"
)
write_csv_safe(sm10_scores, file.path(SUPPLEMENTARY_OUTPUT_DIR, "SM10_Bioindicateur2_scores_long.csv"))

score_counts <- sm10_scores %>%
  filter(!is.na(.data$score)) %>%
  mutate(
    framework = factor(.data$framework, levels = FRAMEWORKS),
    score_factor = factor(.data$score, levels = 6:0)
  ) %>%
  count(.data$indicator, .data$framework, .data$score_factor, name = "n") %>%
  tidyr::complete(
    indicator = LEGACY_INDICATORS,
    framework = factor(FRAMEWORKS, levels = FRAMEWORKS),
    score_factor = factor(6:0, levels = 6:0),
    fill = list(n = 0L)
  )
write_csv_safe(score_counts, file.path(SUPPLEMENTARY_OUTPUT_DIR, "SM10_Bioindicateur2_score_counts.csv"))

extract_nemenyi_pairs <- function(pmat, frameworks = FRAMEWORKS) {
  pmat <- as.matrix(pmat)
  out <- list()
  k <- 1L
  for (i in seq_len(length(frameworks) - 1L)) {
    for (j in (i + 1L):length(frameworks)) {
      a <- frameworks[i]
      b <- frameworks[j]
      p <- NA_real_
      if (a %in% rownames(pmat) && b %in% colnames(pmat)) p <- pmat[a, b]
      if ((!is.finite(p)) && b %in% rownames(pmat) && a %in% colnames(pmat)) p <- pmat[b, a]
      out[[k]] <- tibble(framework_1 = a, framework_2 = b, p_nemenyi = as.numeric(p))
      k <- k + 1L
    }
  }
  bind_rows(out)
}

run_score_tests <- function(score_long, indicator) {
  wide <- score_long %>%
    filter(.data$indicator == .env$indicator) %>%
    select("site_id", "framework", "score") %>%
    pivot_wider(names_from = "framework", values_from = "score") %>%
    drop_na(all_of(FRAMEWORKS))

  mat <- as.matrix(wide[, FRAMEWORKS, drop = FALSE])
  storage.mode(mat) <- "numeric"
  if (nrow(mat) < 2L) stop("Too few complete Bioindicateur2 sites for ", indicator)

  # PMCMRplus::frdAllPairsNemenyiTest() delegates matrix handling to
  # frdRanks(), which uses matrix row names as block labels. A matrix
  # created from a tibble has NULL row names, which can trigger
  # `levels<-.factor`: "number of levels differs" inside PMCMRplus.
  # Give each complete Bioindicateur2 site an explicit, unique block label.
  rownames(mat) <- make.unique(as.character(wide$site_id))

  fr <- stats::friedman.test(mat)
  nm <- PMCMRplus::frdAllPairsNemenyiTest(mat)
  pairs <- extract_nemenyi_pairs(nm$p.value, frameworks = FRAMEWORKS) %>%
    mutate(indicator = indicator, .before = 1)

  list(
    friedman = tibble(
      indicator = indicator,
      n_complete_sites = nrow(mat),
      statistic = as.numeric(fr$statistic),
      df = as.numeric(fr$parameter),
      p_value = as.numeric(fr$p.value)
    ),
    pairs = pairs
  )
}

score_tests <- lapply(LEGACY_INDICATORS, function(ind) run_score_tests(sm10_scores, ind))
sm11_friedman <- bind_rows(lapply(score_tests, `[[`, "friedman"))
sm11_nemenyi <- bind_rows(lapply(score_tests, `[[`, "pairs"))
write_csv_safe(sm11_friedman, file.path(SUPPLEMENTARY_OUTPUT_DIR, "SM11_Friedman_tests.csv"))
write_csv_safe(sm11_nemenyi, file.path(SUPPLEMENTARY_OUTPUT_DIR, "SM11_Nemenyi_pairwise_tests.csv"))

make_score_panel <- function(indicator, show_legend = FALSE) {
  d <- score_counts %>% filter(.data$indicator == .env$indicator)
  pair_tbl <- sm11_nemenyi %>%
    filter(.data$indicator == .env$indicator) %>%
    transmute(
      framework_1 = .data$framework_1,
      framework_2 = .data$framework_2,
      p_nemenyi = .data$p_nemenyi
    )
  cld <- compact_letters_from_pairs(pair_tbl, p_col = "p_nemenyi", group_order = FRAMEWORKS)
  top <- d %>%
    group_by(.data$framework) %>%
    summarise(y = sum(.data$n), .groups = "drop") %>%
    mutate(letter = unname(cld[as.character(.data$framework)]))

  ggplot(d, aes(x = .data$framework, y = .data$n, fill = .data$score_factor)) +
    geom_col(colour = "black", width = 0.70, linewidth = 0.25) +
    geom_text(
      data = top,
      aes(x = .data$framework, y = .data$y + 0.7, label = .data$letter),
      inherit.aes = FALSE, size = 3.4
    ) +
    scale_fill_manual(values = score_colors, name = "Score", drop = FALSE) +
    scale_x_discrete(labels = framework_labels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = "Scoring scheme used", y = "Number of sites", title = unname(indicator_labels[indicator])) +
    theme_minimal(base_size = 9.5) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(size = 7.5),
      legend.position = if (show_legend) "right" else "none"
    )
}

p_s10_1 <- make_score_panel("RSr")
p_s10_2 <- make_score_panel("Shannon")
p_s10_3 <- make_score_panel("Pielou")
p_s10_4 <- make_score_panel("ab", show_legend = TRUE)
p_s10 <- cowplot::plot_grid(p_s10_1, p_s10_2, p_s10_3, p_s10_4, labels = c("A", "B", "C", "D"), ncol = 2, rel_widths = c(1, 1.08))
save_plot_pair(p_s10, "SM10_Bioindicateur2_scores_four_indicators", width = 11.0, height = 8.2)

# ==============================================================================
# Reviewer-driven candidate Supplementary Material
# ==============================================================================
# Run reviewer analysis only if the expected output is absent. This keeps 04
# standalone while avoiding duplicate computation in analysis/run_all.R.
reviewer_key <- file.path(REVIEWER_OUTPUT_DIR, "Table_02_balanced_n58_threshold_stability.csv")
if (!file.exists(reviewer_key)) {
  message("[04] Reviewer outputs absent; running analysis/reviewer_robustness_analysis.R...")
  source_repo("analysis/reviewer_robustness_analysis.R")
}

candidate_dir <- file.path(SUPPLEMENTARY_OUTPUT_DIR, "reviewer_candidates")

# Candidate S12 / robustness: balanced n = 58 threshold stability
for (ext in c("pdf", "png")) {
  src <- file.path(REVIEWER_OUTPUT_DIR, paste0("Figure_S1_balanced_n58_threshold_stability.", ext))
  if (file.exists(src)) {
    file.copy(src, file.path(candidate_dir, paste0("Candidate_balanced_n58_threshold_stability.", ext)), overwrite = TRUE)
  }
}

# Candidate outlier figure: original versus 1-99% winsorised q20-q80.
full_thr_path <- file.path(REVIEWER_OUTPUT_DIR, "Table_01_full_untrimmed_thresholds.csv")
win_thr_path <- file.path(REVIEWER_OUTPUT_DIR, "Table_S2_winsorised_thresholds.csv")
if (file.exists(full_thr_path) && file.exists(win_thr_path)) {
  full_thr <- utils::read.csv(full_thr_path, stringsAsFactors = FALSE)
  win_thr <- utils::read.csv(win_thr_path, stringsAsFactors = FALSE)
  outlier_plot_data <- full_thr %>%
    select(framework, indicator, boundary, probability, full = value) %>%
    left_join(
      win_thr %>% select(framework, indicator, boundary, winsorised = value),
      by = c("framework", "indicator", "boundary")
    ) %>%
    filter(.data$boundary %in% c("q20", "q40", "q60", "q80")) %>%
    mutate(
      framework = factor(.data$framework, levels = FRAMEWORKS),
      indicator_label = recode(.data$indicator, RSr = "Species richness", ab = "Density")
    )
  write_csv_safe(outlier_plot_data, file.path(candidate_dir, "Candidate_outlier_q20_q80_plot_data.csv"))

  p_out <- ggplot(outlier_plot_data, aes(x = .data$full, y = .data$winsorised, colour = .data$framework, shape = .data$boundary)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.45, linetype = 2) +
    geom_point(size = 2.6) +
    facet_wrap(~indicator_label, scales = "free") +
    scale_colour_manual(values = framework_colors, labels = framework_labels, name = "Framework") +
    labs(x = "Original q20-q80 threshold", y = "Winsorised 1-99% q20-q80 threshold", shape = "Quantile") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank())
  save_plot_pair(p_out, "Candidate_outlier_internal_threshold_stability", 8.5, 4.6, directory = candidate_dir)
}

# Candidate cross-framework concordance figure
for (ext in c("pdf", "png")) {
  src <- file.path(REVIEWER_OUTPUT_DIR, paste0("Figure_01_cross_framework_score_agreement.", ext))
  if (file.exists(src)) {
    file.copy(src, file.path(candidate_dir, paste0("Candidate_cross_framework_score_agreement.", ext)), overwrite = TRUE)
  }
}

# Candidate pooled raw vs balanced central thresholds
pool_path <- file.path(REVIEWER_OUTPUT_DIR, "Table_05_pooled_reference_thresholds.csv")
if (file.exists(pool_path)) {
  pool <- utils::read.csv(pool_path, stringsAsFactors = FALSE) %>%
    filter(.data$boundary %in% c("q20", "q40", "q60", "q80")) %>%
    mutate(
      scenario_label = recode(.data$scenario, pooled_raw = "Pooled raw", pooled_balanced = "Pooled balanced"),
      indicator_label = recode(.data$indicator, RSr = "Species richness", ab = "Density")
    )
  write_csv_safe(pool, file.path(candidate_dir, "Candidate_pooled_thresholds_q20_q80.csv"))

  p_pool <- ggplot(pool, aes(x = .data$boundary, y = .data$value, group = .data$scenario_label, colour = .data$scenario_label)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.3) +
    facet_wrap(~indicator_label, scales = "free_y") +
    labs(x = "Quantile threshold", y = "Indicator value", colour = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = "top")
  save_plot_pair(p_pool, "Candidate_pooled_raw_vs_balanced_thresholds", 8.0, 4.6, directory = candidate_dir)
}

# Candidate raw LRR vs Delta-score figure
for (ext in c("pdf", "png")) {
  src <- file.path(REVIEWER_OUTPUT_DIR, paste0("Figure_02_raw_LRR_vs_delta_score.", ext))
  if (file.exists(src)) {
    file.copy(src, file.path(candidate_dir, paste0("Candidate_raw_LRR_vs_delta_score.", ext)), overwrite = TRUE)
  }
}

# Copy key reviewer tables into the candidate supplementary folder. Keep their
# original names after a clear prefix so manuscript placement can be decided later.
reviewer_tables_to_copy <- c(
  "Table_02_balanced_n58_threshold_stability.csv",
  "Table_S1b_balanced_n58_pairwise_threshold_differences.csv",
  "Table_S3_outlier_threshold_sensitivity.csv",
  "Table_S3b_outlier_score_sensitivity.csv",
  "Table_S4_extreme_observation_diagnostics.csv",
  "Table_05_pooled_reference_thresholds.csv",
  "Table_07_pooled_raw_vs_balanced_agreement.csv",
  "Table_08_raw_LRR_and_delta_score.csv",
  "Table_09_LRR_delta_score_summary.csv"
)
for (fn in reviewer_tables_to_copy) {
  src <- file.path(REVIEWER_OUTPUT_DIR, fn)
  if (file.exists(src)) {
    file.copy(src, file.path(candidate_dir, paste0("Candidate_", fn)), overwrite = TRUE)
  }
}

# Manifest --------------------------------------------------------------------
manifest <- c(
  paste0("Generated: ", format(Sys.time(), tz = "UTC"), " UTC"),
  "Supplementary Material regeneration: COMPLETE",
  "",
  "Historical supplementary items regenerated from the current reproducible pipeline:",
  " SM1  Site counts by land-cover category",
  " SM2  Indicator distributions, all sites (RSr, Shannon, Pielou, density)",
  " SM3  Descriptive statistics, all sites",
  " SM4  Pairwise Cramer-von Mises tests, all sites, Bonferroni-adjusted",
  " SM5  Indicator distributions, agricultural sites",
  " SM6  Descriptive statistics, agricultural sites",
  " SM7  Pairwise Cramer-von Mises tests, agricultural sites, Bonferroni-adjusted",
  " SM8  Seven-class scoring-method illustration",
  " SM9  Framework-specific scoring schemes, agricultural sites",
  " SM10 Bioindicateur2 score distributions for four indicators",
  " SM11 Friedman and Nemenyi tests on Bioindicateur2 scores",
  "",
  "Revision decisions:",
  " - current reproducible inputs/indicators/scores are authoritative;",
  " - the historical ad hoc TIGA abundance <75000 filter is not used here;",
  " - Shannon and Pielou are retained as historical secondary/supplementary indicators;",
  " - kurtosis is corrected because e1071::kurtosis() already returns excess kurtosis; the historical double-subtracted value is retained as an audit column;",
  " - Cramer-test p-values use 10,000 bootstrap replicates with fixed seed 1234 to reduce Monte-Carlo instability;",
  " - the SM9 density display extends to 110,000 ind. m^-2 so the complete TIGA score range is visible;",
  " - reviewer robustness analyses remain focused on manuscript indicators RSr and density.",
  "",
  "Reviewer-driven candidate supplementary items are stored under:",
  " output/supplementary/reviewer_candidates/",
  "Their final numbering/placement should be decided during manuscript revision.",
  "",
  "Original submitted supplementary rasters are retained under:",
  " legacy/submitted_supplementary_figures/"
)
write_text_lines(manifest, file.path(SUPPLEMENTARY_OUTPUT_DIR, "SUPPLEMENTARY_REPRODUCTION_MANIFEST.txt"))

capture.output(sessionInfo(), file = file.path(SUPPLEMENTARY_OUTPUT_DIR, "sessionInfo_supplementary.txt"))
message("[04] Supplementary Material written to: ", SUPPLEMENTARY_OUTPUT_DIR)
