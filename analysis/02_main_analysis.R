#!/usr/bin/env Rscript

# ==============================================================================
# 02 - Main analysis
# ==============================================================================
#
# Starting from site-level Collembola indicator values, this script:
#   1. selects the agricultural reference populations;
#   2. describes their RSr and density distributions;
#   3. derives framework-specific empirical quantile boundaries;
#   4. assigns framework-specific scores to the same 26 independent
#      Bioindicateur2 observations.
#
# ==============================================================================

source("R/project_paths.R")
source_repo("R/io_helpers.R")
source_repo("R/indicator_scoring.R")

required_packages(
  c(
    "dplyr",
    "tidyr",
    "tibble"
  )
)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})

message(
  "[02] Running main analysis..."
)

# ------------------------------------------------------------------------------
# Load analysis-ready data
# ------------------------------------------------------------------------------

reference_file <- file.path(
  ANALYSIS_DATA_DIR,
  "reference_indicator_values.csv"
)

bioindicateur2_file <- file.path(
  ANALYSIS_DATA_DIR,
  "bioindicateur2_indicator_values.csv"
)

assert_files_exist(
  c(
    reference_file,
    bioindicateur2_file
  ),
  context = "Analysis-ready input"
)

reference <- read.csv(
  reference_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

bioindicateur2 <- read.csv(
  bioindicateur2_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------------------------
# Agricultural reference populations
# ------------------------------------------------------------------------------

reference_agri <- reference %>%
  filter(
    .data$framework %in%
      FRAMEWORKS,
    .data$CLC_niveau1 ==
      AGRICULTURAL_CLC1
  ) %>%
  mutate(
    framework = factor(
      .data$framework,
      levels = FRAMEWORKS
    )
  )

# ------------------------------------------------------------------------------
# Reference sample sizes
# ------------------------------------------------------------------------------

sample_sizes <- reference_agri %>%
  count(
    .data$framework,
    name = "n_sites"
  )

write_csv_safe(
  sample_sizes,
  file.path(
    ANALYSIS_OUTPUT_DIR,
    "reference_sample_sizes_agricultural.csv"
  )
)

# ------------------------------------------------------------------------------
# Descriptive statistics
# ------------------------------------------------------------------------------

indicator_summary <- bind_rows(
  lapply(
    FRAMEWORKS,
    function(fw) {

      d <- reference_agri %>%
        filter(
          as.character(
            .data$framework
          ) == fw
        )

      bind_rows(
        lapply(
          MANUSCRIPT_INDICATORS,
          function(ind) {

            x <- d[[ind]]
            x <- x[
              is.finite(x)
            ]

            tibble(
              framework = fw,
              indicator = ind,
              n = length(x),
              mean = mean(x),
              sd = stats::sd(x),
              median =
                stats::median(x),
              min = min(x),
              max = max(x)
            )
          }
        )
      )
    }
  )
)

write_csv_safe(
  indicator_summary,
  file.path(
    ANALYSIS_OUTPUT_DIR,
    "reference_indicator_summary.csv"
  )
)

# ------------------------------------------------------------------------------
# Framework-specific empirical score boundaries
# ------------------------------------------------------------------------------

reference_thresholds <-
  build_break_table(
    reference_agri,
    frameworks = FRAMEWORKS,
    indicators =
      MANUSCRIPT_INDICATORS,
    scenario =
      "reference_framework"
  )

write_csv_safe(
  reference_thresholds,
  file.path(
    ANALYSIS_OUTPUT_DIR,
    "reference_thresholds.csv"
  )
)

# ------------------------------------------------------------------------------
# Scores assigned to the same Bioindicateur2 observations
# ------------------------------------------------------------------------------

bioindicateur2_scores <-
  score_validation_data(
    bioindicateur2,
    reference_thresholds,
    frameworks = FRAMEWORKS,
    indicators =
      MANUSCRIPT_INDICATORS,
    scenario_name =
      "reference_framework"
  )

write_csv_safe(
  bioindicateur2_scores,
  file.path(
    ANALYSIS_OUTPUT_DIR,
    "bioindicateur2_scores.csv"
  )
)

score_distribution <-
  bioindicateur2_scores %>%
  count(
    .data$indicator,
    .data$framework,
    .data$score,
    name = "n"
  ) %>%
  group_by(
    .data$indicator,
    .data$framework
  ) %>%
  mutate(
    proportion =
      .data$n /
      sum(.data$n)
  ) %>%
  ungroup()

write_csv_safe(
  score_distribution,
  file.path(
    ANALYSIS_OUTPUT_DIR,
    "bioindicateur2_score_distribution.csv"
  )
)

# ------------------------------------------------------------------------------
# Rarefaction-support metadata retained from the original processing
# ------------------------------------------------------------------------------

rarefaction_support <-
  reference %>%
  mutate(
    below_rarefaction_n =
      as.logical(
        .data$below_rarefaction_n
      )
  ) %>%
  count(
    .data$framework,
    .data$below_rarefaction_n,
    name = "n"
  )

write_csv_safe(
  rarefaction_support,
  file.path(
    ANALYSIS_OUTPUT_DIR,
    "rarefaction_support_reference.csv"
  )
)

# ------------------------------------------------------------------------------
# Reproducibility record
# ------------------------------------------------------------------------------

dir.create(
  file.path(
    ROOT_DIR,
    "reproducibility"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    ROOT_DIR,
    "reproducibility",
    "sessionInfo.txt"
  )
)

message(
  "[02] Main analysis complete."
)

message(
  "[02] Results written to: ",
  ANALYSIS_OUTPUT_DIR
)
