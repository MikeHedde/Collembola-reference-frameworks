#!/usr/bin/env Rscript

# ==============================================================================
# 02 - Reproduce the legacy indicator/scoring core before reviewer sensitivities
# ==============================================================================
# This script is deliberately separate from the robustness analyses. It exposes
# exactly which parts of the historical analysis are reproduced and highlights
# the undocumented TIGA abundance filter (< 75000) as a distinct legacy scenario.
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

required_packages(c("dplyr", "tidyr", "tibble", "vegan"))

ref_rds <- file.path(DERIVED_DIR, "reference_frameworks_merged.rds")
bio_rds <- file.path(DERIVED_DIR, "bioindicateur2_curated.rds")
if (!file.exists(ref_rds) || !file.exists(bio_rds)) {
  message("[02] Derived datasets absent; running analysis/01_build_analysis_dataset.R first.")
  source_repo("analysis/01_build_analysis_dataset.R")
}

message("[02] Recomputing legacy indicators and score boundaries...")
reference_raw <- readRDS(ref_rds)
validation_raw <- readRDS(bio_rds)

reference <- compute_legacy_indicators(reference_raw)
validation <- compute_legacy_indicators(validation_raw)

reference_agri <- reference %>%
  filter(.data$Projet %in% FRAMEWORKS, as.character(.data$CLC_niveau1) == as.character(AGRICULTURAL_CLC1)) %>%
  mutate(Projet = factor(as.character(.data$Projet), levels = FRAMEWORKS))
validation_agri <- validation %>%
  filter(as.character(.data$CLC_niveau1) == as.character(AGRICULTURAL_CLC1))

sample_sizes <- reference_agri %>% count(.data$Projet, name = "n_sites")
write_csv_safe(sample_sizes, file.path(BASELINE_OUTPUT_DIR, "baseline_reference_sample_sizes_agricultural.csv"))

indicator_summary <- reference_agri %>%
  group_by(.data$Projet) %>%
  summarise(
    across(
      all_of(LEGACY_INDICATORS),
      list(
        n = ~ sum(is.finite(.x)),
        mean = ~ mean(.x, na.rm = TRUE),
        sd = ~ stats::sd(.x, na.rm = TRUE),
        median = ~ stats::median(.x, na.rm = TRUE),
        min = ~ min(.x, na.rm = TRUE),
        max = ~ max(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )
write_csv_safe(indicator_summary, file.path(BASELINE_OUTPUT_DIR, "baseline_indicator_summary.csv"))

# Main revised baseline: all observations retained.
full_untrimmed <- build_break_table(
  reference_agri,
  indicators = MANUSCRIPT_INDICATORS,
  scenario = "full_untrimmed"
)
write_csv_safe(full_untrimmed, file.path(BASELINE_OUTPUT_DIR, "baseline_thresholds_full_untrimmed.csv"))

# Historical scoring exception found in Analyse statistique.R.
legacy_breaks <- legacy_tiga_ab_breaks(reference_agri, threshold = 75000)
write_csv_safe(legacy_breaks, file.path(BASELINE_OUTPUT_DIR, "baseline_thresholds_legacy_tiga_ab_lt_75000.csv"))

legacy_filter_diagnostic <- reference_agri %>%
  filter(.data$Projet == "TIGA_rural") %>%
  summarise(
    n_total = n(),
    n_ab_ge_75000 = sum(is.finite(.data$ab) & .data$ab >= 75000),
    max_ab = max(.data$ab, na.rm = TRUE)
  )
write_csv_safe(legacy_filter_diagnostic, file.path(BASELINE_OUTPUT_DIR, "diagnostic_legacy_tiga_ab_filter.csv"))

threshold_scenario_compare <- full_untrimmed %>%
  select(framework, indicator, boundary, full_untrimmed = value) %>%
  left_join(
    legacy_breaks %>% select(framework, indicator, boundary, legacy = value),
    by = c("framework", "indicator", "boundary")
  ) %>%
  mutate(difference_legacy_minus_untrimmed = .data$legacy - .data$full_untrimmed)
write_csv_safe(threshold_scenario_compare, file.path(BASELINE_OUTPUT_DIR, "baseline_legacy_vs_untrimmed_thresholds.csv"))

# Scores on the independent Bioindicateur2 agricultural subset under both
# scenarios. This makes any consequence of the historical ad hoc rule explicit.
score_untrimmed <- score_validation_data(
  validation_agri, full_untrimmed,
  frameworks = FRAMEWORKS,
  indicators = MANUSCRIPT_INDICATORS,
  scenario_name = "full_untrimmed"
)
score_legacy <- score_validation_data(
  validation_agri, legacy_breaks,
  frameworks = FRAMEWORKS,
  indicators = MANUSCRIPT_INDICATORS,
  scenario_name = "legacy_tiga_ab_lt_75000"
)
write_csv_safe(bind_rows(score_untrimmed, score_legacy), file.path(BASELINE_OUTPUT_DIR, "baseline_validation_scores_scenarios.csv"))

score_scenario_change <- score_untrimmed %>%
  select(site_id, indicator, framework, score_untrimmed = score) %>%
  left_join(
    score_legacy %>% select(site_id, indicator, framework, score_legacy = score),
    by = c("site_id", "indicator", "framework")
  ) %>%
  mutate(score_difference = .data$score_legacy - .data$score_untrimmed)
write_csv_safe(score_scenario_change, file.path(BASELINE_OUTPUT_DIR, "baseline_score_change_legacy_filter.csv"))

rarefaction_diagnostic <- bind_rows(
  reference_agri %>%
    transmute(dataset = "reference", Projet = as.character(.data$Projet), below_500 = .data$.below_rarefaction_n) %>%
    count(.data$dataset, .data$Projet, .data$below_500, name = "n"),
  validation_agri %>%
    transmute(dataset = "validation", Projet = "Bioindicateur2", below_500 = .data$.below_rarefaction_n) %>%
    count(.data$dataset, .data$Projet, .data$below_500, name = "n")
)
write_csv_safe(rarefaction_diagnostic, file.path(BASELINE_OUTPUT_DIR, "diagnostic_rarefaction_support.csv"))

capture.output(sessionInfo(), file = file.path(BASELINE_OUTPUT_DIR, "sessionInfo_baseline.txt"))
message("[02] Baseline reproduction outputs written to: ", BASELINE_OUTPUT_DIR)
