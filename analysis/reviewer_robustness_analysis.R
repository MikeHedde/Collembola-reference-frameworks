#!/usr/bin/env Rscript

# ==============================================================================
# Reviewer robustness analyses for Collembola empirical reference frameworks
# ==============================================================================
# Runs the five reviewer-driven checks on the analysis dataset rebuilt from the
# separate curated workbooks. The historical merged Excel file is never used as
# an analysis input.
#
# Main manuscript metrics: rarefied richness (RSr) and density (ab).
# Shannon and Pielou are still recomputed in the baseline reproduction to retain
# provenance with H. Suarez's original scripts, but are not expanded into every
# reviewer sensitivity unless the manuscript scope changes.
# ============================================================================== 

# Locate the project even when Rscript is launched outside the repository root.
.project_paths_candidates <- c(file.path(getwd(), "R", "project_paths.R"))
.args_all <- commandArgs(trailingOnly = FALSE)
.file_arg <- grep("^--file=", .args_all, value = TRUE)
if (length(.file_arg)) {
  .script_file <- normalizePath(sub("^--file=", "", .file_arg[1]), mustWork = FALSE)
  .project_paths_candidates <- c(.project_paths_candidates, file.path(dirname(.script_file), "..", "R", "project_paths.R"))
}
.ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (!is.null(.ofile)) {
  .project_paths_candidates <- c(.project_paths_candidates, file.path(dirname(normalizePath(.ofile, mustWork = FALSE)), "..", "R", "project_paths.R"))
}
.project_paths_candidates <- unique(normalizePath(.project_paths_candidates, mustWork = FALSE))
.project_paths_file <- .project_paths_candidates[file.exists(.project_paths_candidates)][1]
if (is.na(.project_paths_file)) stop("Could not locate R/project_paths.R. Open the .Rproj or run inside the repository.")
source(.project_paths_file)
rm(.project_paths_candidates, .args_all, .file_arg, .ofile, .project_paths_file)
if (exists(".script_file")) rm(.script_file)
source_repo("R/io_helpers.R")
source_repo("R/data_preparation.R")
source_repo("R/indicator_scoring.R")
source_repo("R/reviewer_functions.R")

required_packages(c("dplyr", "tidyr", "tibble", "ggplot2", "vegan"))

SEED <- 20260818L
N_REPEATS <- 1000L
N_BALANCED <- 58L
WINSOR_PROBS <- c(0.01, 0.99)
CONTROL_MAPPING_FILE <- file.path(CONFIG_DIR, "bioindicateur2_controls.csv")

ref_rds <- file.path(DERIVED_DIR, "reference_frameworks_merged.rds")
bio_rds <- file.path(DERIVED_DIR, "bioindicateur2_curated.rds")
if (!file.exists(ref_rds) || !file.exists(bio_rds)) {
  message("[reviewer] Derived datasets absent; running analysis/01_build_analysis_dataset.R first.")
  source_repo("analysis/01_build_analysis_dataset.R")
}

message("[reviewer] Recomputing manuscript indicators...")
reference <- compute_legacy_indicators(readRDS(ref_rds))
validation <- compute_legacy_indicators(readRDS(bio_rds))

reference_agri <- reference %>%
  filter(.data$Projet %in% FRAMEWORKS, as.character(.data$CLC_niveau1) == as.character(AGRICULTURAL_CLC1)) %>%
  mutate(Projet = factor(as.character(.data$Projet), levels = FRAMEWORKS))
validation_agri <- validation %>%
  filter(as.character(.data$CLC_niveau1) == as.character(AGRICULTURAL_CLC1))

sample_sizes <- reference_agri %>% count(.data$Projet, name = "n_sites")
write_csv_safe(sample_sizes, file.path(REVIEWER_OUTPUT_DIR, "Table_00_reference_sample_sizes_agricultural.csv"))
if (min(sample_sizes$n_sites) != N_BALANCED) {
  stop("Expected smallest agricultural framework n = 58, observed n = ", min(sample_sizes$n_sites), ". Run audit/build scripts first.")
}

# ---- Baseline thresholds used by robustness comparisons ----------------------
full_breaks <- build_break_table(
  reference_agri,
  indicators = MANUSCRIPT_INDICATORS,
  scenario = "full_untrimmed"
)
write_csv_safe(full_breaks, file.path(REVIEWER_OUTPUT_DIR, "Table_01_full_untrimmed_thresholds.csv"))

# Historical ad hoc TIGA abundance filter is made explicit, not silently kept or
# silently removed. This table quantifies the consequence of that decision.
legacy_breaks <- legacy_tiga_ab_breaks(reference_agri, threshold = 75000)
legacy_compare <- full_breaks %>%
  select(framework, indicator, boundary, untrimmed = value) %>%
  left_join(
    legacy_breaks %>% select(framework, indicator, boundary, legacy_tiga_filter = value),
    by = c("framework", "indicator", "boundary")
  ) %>%
  mutate(shift = .data$legacy_tiga_filter - .data$untrimmed)
write_csv_safe(legacy_compare, file.path(REVIEWER_OUTPUT_DIR, "Table_S0_historical_TIGA_filter_sensitivity.csv"))

# ==============================================================================
# Q1. Unequal sample sizes: repeated balanced n = 58 subsampling
# ==============================================================================
message("[reviewer] Q1: balanced n=58 repeated subsampling (", N_REPEATS, " repeats)...")
balanced_breaks <- balanced_subsampling_breaks(
  reference_agri,
  n = N_BALANCED,
  B = N_REPEATS,
  seed = SEED,
  indicators = MANUSCRIPT_INDICATORS
)
write_csv_safe(balanced_breaks, file.path(REVIEWER_OUTPUT_DIR, "Table_S1_balanced_n58_all_repeats.csv"))

balanced_summary <- balanced_breaks %>%
  group_by(.data$framework, .data$indicator, .data$boundary, .data$probability) %>%
  summarise(
    median = median(.data$value, na.rm = TRUE),
    q025 = quantile(.data$value, 0.025, na.rm = TRUE),
    q975 = quantile(.data$value, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    full_breaks %>% select(framework, indicator, boundary, full_value = value),
    by = c("framework", "indicator", "boundary")
  )
write_csv_safe(balanced_summary, file.path(REVIEWER_OUTPUT_DIR, "Table_02_balanced_n58_threshold_stability.csv"))

framework_pairs <- combn(FRAMEWORKS, 2, simplify = FALSE)
balanced_pairwise <- bind_rows(lapply(framework_pairs, function(pair) {
  w <- balanced_breaks %>%
    filter(.data$framework %in% pair) %>%
    select(replicate, indicator, boundary, framework, value) %>%
    pivot_wider(names_from = "framework", values_from = "value")
  w$diff <- w[[pair[1]]] - w[[pair[2]]]
  w %>%
    group_by(.data$indicator, .data$boundary) %>%
    summarise(
      framework_1 = pair[1],
      framework_2 = pair[2],
      median_difference = median(.data$diff, na.rm = TRUE),
      q025 = quantile(.data$diff, 0.025, na.rm = TRUE),
      q975 = quantile(.data$diff, 0.975, na.rm = TRUE),
      probability_difference_gt_0 = mean(.data$diff > 0, na.rm = TRUE),
      .groups = "drop"
    )
}))
write_csv_safe(balanced_pairwise, file.path(REVIEWER_OUTPUT_DIR, "Table_S1b_balanced_n58_pairwise_threshold_differences.csv"))

p_balanced <- balanced_summary %>%
  filter(!.data$boundary %in% c("q00", "q100")) %>%
  ggplot(aes(x = .data$framework, y = .data$median, ymin = .data$q025, ymax = .data$q975)) +
  geom_pointrange() +
  geom_point(aes(y = .data$full_value), shape = 1, size = 2.5) +
  facet_grid(indicator ~ boundary, scales = "free_y") +
  scale_x_discrete(labels = FRAMEWORK_LABELS) +
  labs(
    x = NULL,
    y = "Empirical score boundary",
    caption = "Filled: median of repeated balanced n=58 subsamples; bars: 95% interval; open: full untrimmed reference."
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(REVIEWER_OUTPUT_DIR, "Figure_S1_balanced_n58_threshold_stability.pdf"), p_balanced, width = 11, height = 6.5)

# ==============================================================================
# Q2. Outlier sensitivity: untrimmed vs 1-99% winsorisation + >3 IQR diagnostics
# ==============================================================================
message("[reviewer] Q2: outlier sensitivity...")
winsorised_reference <- reference_agri
for (fw in FRAMEWORKS) {
  idx <- which(as.character(winsorised_reference$Projet) == fw)
  for (ind in MANUSCRIPT_INDICATORS) {
    winsorised_reference[[ind]][idx] <- winsorise(winsorised_reference[[ind]][idx], WINSOR_PROBS)
  }
}
winsor_breaks <- build_break_table(winsorised_reference, indicators = MANUSCRIPT_INDICATORS, scenario = "winsorised_01_99")
write_csv_safe(winsor_breaks, file.path(REVIEWER_OUTPUT_DIR, "Table_S2_winsorised_thresholds.csv"))

outlier_threshold_compare <- full_breaks %>%
  select(framework, indicator, boundary, full_value = value) %>%
  left_join(
    winsor_breaks %>% select(framework, indicator, boundary, winsor_value = value),
    by = c("framework", "indicator", "boundary")
  ) %>%
  mutate(
    absolute_shift = .data$winsor_value - .data$full_value,
    relative_shift = ifelse(.data$full_value == 0, NA_real_, .data$absolute_shift / abs(.data$full_value))
  )
write_csv_safe(outlier_threshold_compare, file.path(REVIEWER_OUTPUT_DIR, "Table_S3_outlier_threshold_sensitivity.csv"))

extreme_counts <- bind_rows(lapply(FRAMEWORKS, function(fw) {
  d <- reference_agri %>% filter(as.character(.data$Projet) == fw)
  bind_rows(lapply(MANUSCRIPT_INDICATORS, function(ind) {
    flag <- extreme_iqr_flag(d[[ind]], k = 3)
    tibble(
      framework = fw,
      indicator = ind,
      n = sum(is.finite(d[[ind]])),
      n_extreme_3IQR = sum(flag, na.rm = TRUE),
      proportion_extreme_3IQR = mean(flag, na.rm = TRUE)
    )
  }))
}))
write_csv_safe(extreme_counts, file.path(REVIEWER_OUTPUT_DIR, "Table_S4_extreme_observation_diagnostics.csv"))

score_full <- score_validation_data(validation_agri, full_breaks, frameworks = FRAMEWORKS, indicators = MANUSCRIPT_INDICATORS, scenario_name = "full_untrimmed")
score_winsor <- score_validation_data(validation_agri, winsor_breaks, frameworks = FRAMEWORKS, indicators = MANUSCRIPT_INDICATORS, scenario_name = "winsorised_01_99")
outlier_score_sensitivity <- score_full %>%
  select(site_id, indicator, framework, score_full = score) %>%
  left_join(
    score_winsor %>% select(site_id, indicator, framework, score_winsor = score),
    by = c("site_id", "indicator", "framework")
  ) %>%
  group_by(.data$indicator, .data$framework) %>%
  summarise(
    n = sum(complete.cases(.data$score_full, .data$score_winsor)),
    spearman_rho = suppressWarnings(cor(.data$score_full, .data$score_winsor, method = "spearman", use = "complete.obs")),
    exact_agreement = mean(.data$score_full == .data$score_winsor, na.rm = TRUE),
    within_one_class = mean(abs(.data$score_full - .data$score_winsor) <= 1, na.rm = TRUE),
    .groups = "drop"
  )
write_csv_safe(outlier_score_sensitivity, file.path(REVIEWER_OUTPUT_DIR, "Table_S3b_outlier_score_sensitivity.csv"))

# ==============================================================================
# Q3. Cross-framework concordance on the same Bioindicateur2 observations
# ==============================================================================
message("[reviewer] Q3: cross-framework concordance...")
write_csv_safe(score_full, file.path(REVIEWER_OUTPUT_DIR, "Table_03_validation_scores_by_framework.csv"))
agreement <- pairwise_score_agreement(score_full)
write_csv_safe(agreement, file.path(REVIEWER_OUTPUT_DIR, "Table_04_cross_framework_score_agreement.csv"))

agreement_long <- agreement %>%
  select(indicator, framework_1, framework_2, exact_agreement, within_one_class, broad_class_agreement) %>%
  pivot_longer(
    cols = c("exact_agreement", "within_one_class", "broad_class_agreement"),
    names_to = "agreement_metric", values_to = "agreement"
  )
p_agreement <- ggplot(agreement_long, aes(x = .data$framework_1, y = .data$framework_2, fill = .data$agreement)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", .data$agreement)), size = 3) +
  facet_grid(indicator ~ agreement_metric) +
  scale_fill_gradient(low = "white", high = "black", limits = c(0, 1)) +
  scale_x_discrete(labels = FRAMEWORK_LABELS) +
  scale_y_discrete(labels = FRAMEWORK_LABELS) +
  labs(x = NULL, y = NULL, fill = "Agreement") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(file.path(REVIEWER_OUTPUT_DIR, "Figure_01_cross_framework_score_agreement.pdf"), p_agreement, width = 11, height = 5.5)

# ==============================================================================
# Q4. Pooled raw (site-weighted) vs pooled balanced (equal framework weight)
# ==============================================================================
message("[reviewer] Q4: pooled raw versus pooled balanced...")
pooled_raw <- bind_rows(lapply(MANUSCRIPT_INDICATORS, function(ind) {
  breaks_to_table(score_breaks(reference_agri[[ind]]), "POOLED_RAW", ind, "pooled_raw")
}))
pooled_balanced_result <- build_balanced_pooled_breaks(reference_agri, N_BALANCED, N_REPEATS, SEED + 1L, MANUSCRIPT_INDICATORS)
pooled_balanced <- pooled_balanced_result$summary %>% select(scenario, framework, indicator, boundary, probability, value)
write_csv_safe(bind_rows(pooled_raw, pooled_balanced), file.path(REVIEWER_OUTPUT_DIR, "Table_05_pooled_reference_thresholds.csv"))
write_csv_safe(pooled_balanced_result$draws, file.path(REVIEWER_OUTPUT_DIR, "Table_S5_pooled_balanced_all_repeats.csv"))

score_pooled_raw <- score_validation_data(validation_agri, pooled_raw, frameworks = "POOLED_RAW", indicators = MANUSCRIPT_INDICATORS, scenario_name = "pooled_raw")
score_pooled_balanced <- score_validation_data(validation_agri, pooled_balanced, frameworks = "POOLED_BALANCED", indicators = MANUSCRIPT_INDICATORS, scenario_name = "pooled_balanced")
score_pooled <- bind_rows(score_pooled_raw, score_pooled_balanced)
write_csv_safe(score_pooled, file.path(REVIEWER_OUTPUT_DIR, "Table_06_validation_scores_pooled_references.csv"))

pooled_agreement <- score_pooled %>%
  select(site_id, indicator, framework, score) %>%
  pivot_wider(names_from = "framework", values_from = "score") %>%
  group_by(.data$indicator) %>%
  summarise(
    n = sum(complete.cases(.data$POOLED_RAW, .data$POOLED_BALANCED)),
    spearman_rho = suppressWarnings(cor(.data$POOLED_RAW, .data$POOLED_BALANCED, method = "spearman", use = "complete.obs")),
    exact_agreement = mean(.data$POOLED_RAW == .data$POOLED_BALANCED, na.rm = TRUE),
    within_one_class = mean(abs(.data$POOLED_RAW - .data$POOLED_BALANCED) <= 1, na.rm = TRUE),
    .groups = "drop"
  )
write_csv_safe(pooled_agreement, file.path(REVIEWER_OUTPUT_DIR, "Table_07_pooled_raw_vs_balanced_agreement.csv"))

# ==============================================================================
# Q5. Raw LRR vs framework-dependent Delta-score -- explicit controls only
# ==============================================================================
message("[reviewer] Q5: preparing LRR / Delta-score analysis...")
write_observed_control_levels(validation_agri, file.path(REVIEWER_OUTPUT_DIR, "diagnostic_Bioindicateur2_treatment_levels.csv"))
control_mapping <- read_control_mapping(CONTROL_MAPPING_FILE)

if (control_mapping_is_ready(control_mapping)) {
  lrr_result <- make_explicit_lrr_table(validation_agri, score_full, control_mapping, MANUSCRIPT_INDICATORS)
  write_csv_safe(lrr_result$diagnostics, file.path(REVIEWER_OUTPUT_DIR, "diagnostic_control_mapping.csv"))
  if (nrow(lrr_result$lrr)) {
    write_csv_safe(lrr_result$lrr, file.path(REVIEWER_OUTPUT_DIR, "Table_08_raw_LRR_and_delta_score.csv"))

    lrr_summary <- lrr_result$lrr %>%
      group_by(.data$Nom_site, .data$indicator, .data$framework) %>%
      summarise(
        n = n(),
        median_LRR = median(.data$LRR, na.rm = TRUE),
        median_delta_score = median(.data$delta_score, na.rm = TRUE),
        .groups = "drop"
      )
    write_csv_safe(lrr_summary, file.path(REVIEWER_OUTPUT_DIR, "Table_09_LRR_delta_score_summary.csv"))

    p_lrr <- lrr_result$lrr %>%
      filter(is.finite(.data$LRR), is.finite(.data$delta_score)) %>%
      ggplot(aes(x = .data$LRR, y = .data$delta_score)) +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_vline(xintercept = 0, linetype = 2) +
      geom_point(alpha = 0.65) +
      facet_grid(indicator ~ framework) +
      labs(x = "Raw-value log response ratio: ln(treatment / control)", y = "Framework-dependent score contrast (Delta-score)") +
      theme_bw(base_size = 10)
    ggsave(file.path(REVIEWER_OUTPUT_DIR, "Figure_02_raw_LRR_vs_delta_score.pdf"), p_lrr, width = 11, height = 6)
  }
} else {
  message(
    "[reviewer] LRR intentionally NOT calculated: config/bioindicateur2_controls.csv has no enabled explicit mapping.\n",
    "           Inspect output/reviewer/diagnostic_Bioindicateur2_treatment_levels.csv, fill the mapping, and rerun."
  )
  write_csv_safe(
    tibble(
      status = "LRR_NOT_RUN_EXPLICIT_CONTROL_MAPPING_REQUIRED",
      mapping_file = CONTROL_MAPPING_FILE,
      diagnostic_levels = file.path(REVIEWER_OUTPUT_DIR, "diagnostic_Bioindicateur2_treatment_levels.csv")
    ),
    file.path(REVIEWER_OUTPUT_DIR, "diagnostic_control_mapping.csv")
  )
}

# ---- Reproducibility record ---------------------------------------------------
run_metadata <- tibble(
  parameter = c("seed", "repeats", "balanced_n", "score_classes", "rarefaction_n", "agricultural_CLC1", "winsor_lower", "winsor_upper", "manuscript_indicators"),
  value = c(SEED, N_REPEATS, N_BALANCED, N_SCORE_CLASSES, RAREFACTION_N, AGRICULTURAL_CLC1, WINSOR_PROBS[1], WINSOR_PROBS[2], paste(MANUSCRIPT_INDICATORS, collapse = ","))
)
write_csv_safe(run_metadata, file.path(REVIEWER_OUTPUT_DIR, "run_metadata.csv"))
capture.output(sessionInfo(), file = file.path(REVIEWER_OUTPUT_DIR, "sessionInfo.txt"))
message("[reviewer] Done. Results written to: ", REVIEWER_OUTPUT_DIR)
