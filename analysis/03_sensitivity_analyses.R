#!/usr/bin/env Rscript

# ==============================================================================
# 03 - Sensitivity and contrast-translation analyses
# ==============================================================================
#
# Starting from analysis-ready site-level indicator values, this script tests:
#   1. sensitivity of empirical thresholds to unequal reference sample size;
#   2. sensitivity to extreme observations;
#   3. cross-framework concordance of scores;
#   4. sensitivity of a pooled reference to dataset weighting;
#   5. translation of raw treatment-control contrasts into ordinal score space.
#
# ==============================================================================

source("R/project_paths.R")
source_repo("R/io_helpers.R")
source_repo("R/indicator_scoring.R")
source_repo("R/analysis_functions.R")

required_packages(c(
  "dplyr",
  "tidyr",
  "tibble"
))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})

SEED <- 20260818L
N_REPEATS <- 1000L
N_BALANCED <- 58L
WINSOR_PROBS <- c(0.01, 0.99)

# Numerical tolerance used only when assigning a direction to LRR.
# It is not a biological effect-size threshold.
LRR_TOLERANCE <- 1e-8

message(
  "[03] Running sensitivity and contrast-translation analyses..."
)

# ------------------------------------------------------------------------------
# Load analysis-ready inputs
# ------------------------------------------------------------------------------

reference <- read.csv(
  file.path(
    ANALYSIS_DATA_DIR,
    "reference_indicator_values.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

bioindicateur2 <- read.csv(
  file.path(
    ANALYSIS_DATA_DIR,
    "bioindicateur2_indicator_values.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

contrasts <- read.csv(
  file.path(
    ANALYSIS_DATA_DIR,
    "bioindicateur2_contrasts.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

reference_agri <- reference %>%
  filter(
    .data$framework %in% FRAMEWORKS,
    .data$CLC_niveau1 ==
      AGRICULTURAL_CLC1
  ) %>%
  mutate(
    framework = factor(
      .data$framework,
      levels = FRAMEWORKS
    )
  )

sample_sizes <- reference_agri %>%
  count(
    .data$framework,
    name = "n_sites"
  )

if (
  min(sample_sizes$n_sites) !=
    N_BALANCED
) {
  stop(
    "Expected smallest agricultural reference population n = 58; found ",
    min(sample_sizes$n_sites),
    "."
  )
}

# ------------------------------------------------------------------------------
# Full reference thresholds and scores
# ------------------------------------------------------------------------------

full_breaks <- build_break_table(
  reference_agri,
  frameworks = FRAMEWORKS,
  indicators = MANUSCRIPT_INDICATORS,
  scenario = "full_reference"
)

score_full <- score_validation_data(
  bioindicateur2,
  full_breaks,
  frameworks = FRAMEWORKS,
  indicators = MANUSCRIPT_INDICATORS,
  scenario_name = "full_reference"
)

write_csv_safe(
  full_breaks,
  file.path(
    SENSITIVITY_OUTPUT_DIR,
    "reference_thresholds_full.csv"
  )
)

# ==============================================================================
# 1. Equal-size sensitivity: n = 58 in each reference framework
# ==============================================================================

message(
  "[03] Equal-size resampling: n = ",
  N_BALANCED,
  ", B = ",
  N_REPEATS,
  "..."
)

balanced_breaks <-
  balanced_subsampling_breaks(
    reference_agri,
    n = N_BALANCED,
    B = N_REPEATS,
    seed = SEED,
    indicators =
      MANUSCRIPT_INDICATORS
  )

write_csv_safe(
  balanced_breaks,
  file.path(
    SENSITIVITY_OUTPUT_DIR,
    "balanced_n58_all_repeats.csv"
  )
)

balanced_summary <-
  balanced_breaks %>%
  group_by(
    .data$framework,
    .data$indicator,
    .data$boundary,
    .data$probability
  ) %>%
  summarise(
    median =
      median(
        .data$value,
        na.rm = TRUE
      ),
    q025 =
      quantile(
        .data$value,
        0.025,
        na.rm = TRUE
      ),
    q975 =
      quantile(
        .data$value,
        0.975,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  left_join(
    full_breaks %>%
      select(
        framework,
        indicator,
        boundary,
        full_value = value
      ),
    by = c(
      "framework",
      "indicator",
      "boundary"
    )
  )

write_csv_safe(
  balanced_summary,
  file.path(
    SENSITIVITY_OUTPUT_DIR,
    "balanced_n58_threshold_stability.csv"
  )
)

framework_pairs <- combn(
  FRAMEWORKS,
  2,
  simplify = FALSE
)

balanced_pairwise <- bind_rows(
  lapply(
    framework_pairs,
    function(pair) {

      w <- balanced_breaks %>%
        filter(
          .data$framework %in%
            pair
        ) %>%
        select(
          replicate,
          indicator,
          boundary,
          framework,
          value
        ) %>%
        pivot_wider(
          names_from = "framework",
          values_from = "value"
        )

      w$diff <-
        w[[pair[1]]] -
        w[[pair[2]]]

      w %>%
        group_by(
          .data$indicator,
          .data$boundary
        ) %>%
        summarise(
          framework_1 =
            pair[1],
          framework_2 =
            pair[2],
          median_difference =
            median(
              .data$diff,
              na.rm = TRUE
            ),
          q025 =
            quantile(
              .data$diff,
              0.025,
              na.rm = TRUE
            ),
          q975 =
            quantile(
              .data$diff,
              0.975,
              na.rm = TRUE
            ),
          probability_difference_gt_0 =
            mean(
              .data$diff > 0,
              na.rm = TRUE
            ),
          .groups = "drop"
        )
    }
  )
)

write_csv_safe(
  balanced_pairwise,
  file.path(
    SENSITIVITY_OUTPUT_DIR,
    "balanced_n58_pairwise_threshold_differences.csv"
  )
)

internal_persistence <-
  balanced_pairwise %>%
  filter(
    .data$boundary %in%
      c(
        "q20",
        "q40",
        "q60",
        "q80"
      )
  ) %>%
  mutate(
    interval_excludes_zero =
      .data$q025 > 0 |
      .data$q975 < 0
  ) %>%
  group_by(
    .data$indicator
  ) %>%
  summarise(
    n_pairwise_boundaries = n(),
    n_intervals_excluding_zero =
      sum(
        .data$interval_excludes_zero
      ),
    .groups = "drop"
  )

write_csv_safe(
  internal_persistence,
  file.path(
    SENSITIVITY_OUTPUT_DIR,
    "balanced_n58_internal_threshold_persistence.csv"
  )
)

# ==============================================================================
# 2. Extreme-value sensitivity
# ==============================================================================

message(
  "[03] Extreme-value sensitivity..."
)

winsorised_reference <-
  reference_agri

for (fw in FRAMEWORKS) {

  idx <- which(
    as.character(
      winsorised_reference$framework
    ) == fw
  )

  for (ind in MANUSCRIPT_INDICATORS) {

    winsorised_reference[[ind]][idx] <-
      winsorise(
        winsorised_reference[[ind]][idx],
        WINSOR_PROBS
      )
  }
}

winsor_breaks <- build_break_table(
  winsorised_reference,
  frameworks = FRAMEWORKS,
  indicators =
    MANUSCRIPT_INDICATORS,
  scenario =
    "winsorised_01_99"
)

outlier_threshold_compare <-
  full_breaks %>%
  select(
    framework,
    indicator,
    boundary,
    full_value = value
  ) %>%
  left_join(
    winsor_breaks %>%
      select(
        framework,
        indicator,
        boundary,
        winsor_value = value
      ),
    by = c(
      "framework",
      "indicator",
      "boundary"
    )
  ) %>%
  mutate(
    absolute_shift =
      .data$winsor_value -
      .data$full_value
  )

write_csv_safe(
  outlier_threshold_compare,
  file.path(
    SENSITIVITY_OUTPUT_DIR,
    "outlier_threshold_sensitivity.csv"
  )
)

extreme_counts <- bind_rows(
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

            flag <-
              extreme_iqr_flag(
                d[[ind]],
                k = 3
              )

            tibble(
              framework = fw,
              indicator = ind,
              n =
                sum(
                  is.finite(
                    d[[ind]]
                  )
                ),
              n_extreme_3IQR =
                sum(
                  flag,
                  na.rm = TRUE
                ),
              proportion_extreme_3IQR =
                mean(
                  flag,
                  na.rm = TRUE
                )
            )
          }
        )
      )
    }
  )
)

write_csv_safe(
  extreme_counts,
  file.path(
    SENSITIVITY_OUTPUT_DIR,
    "extreme_observation_diagnostics.csv"
  )
)

score_winsor <- score_validation_data(
  bioindicateur2,
  winsor_breaks,
  frameworks = FRAMEWORKS,
  indicators =
    MANUSCRIPT_INDICATORS,
  scenario_name =
    "winsorised_01_99"
)

outlier_score_sensitivity <-
  score_full %>%
  select(
    site_id,
    indicator,
    framework,
    score_full = score
  ) %>%
  left_join(
    score_winsor %>%
      select(
        site_id,
        indicator,
        framework,
        score_winsor = score
      ),
    by = c(
      "site_id",
      "indicator",
      "framework"
    )
  ) %>%
  group_by(
    .data$indicator,
    .data$framework
  ) %>%
  summarise(
    n = n(),
    exact_agreement =
      mean(
        .data$score_full ==
          .data$score_winsor
      ),
    within_one_class =
      mean(
        abs(
          .data$score_full -
            .data$score_winsor
        ) <= 1
      ),
    .groups = "drop"
  )

write_csv_safe(
  outlier_score_sensitivity,
  file.path(
    SENSITIVITY_OUTPUT_DIR,
    "outlier_score_sensitivity.csv"
  )
)

# ==============================================================================
# 3. Cross-framework score concordance
# ==============================================================================

message(
  "[03] Cross-framework score concordance..."
)

agreement <-
  pairwise_score_agreement(
    score_full
  )

write_csv_safe(
  agreement,
  file.path(
    SENSITIVITY_OUTPUT_DIR,
    "cross_framework_score_agreement.csv"
  )
)

# ==============================================================================
# 4. Site-weighted versus equal-weight pooled references
# ==============================================================================

message(
  "[03] Pooled-reference weighting sensitivity..."
)

pooled_raw <- bind_rows(
  lapply(
    MANUSCRIPT_INDICATORS,
    function(ind) {

      breaks_to_table(
        score_breaks(
          reference_agri[[ind]]
        ),
        "POOLED_RAW",
        ind,
        "pooled_raw"
      )
    }
  )
)

pooled_balanced_result <-
  build_balanced_pooled_breaks(
    reference_agri,
    n = N_BALANCED,
    B = N_REPEATS,
    seed = SEED + 1L,
    indicators =
      MANUSCRIPT_INDICATORS
  )

pooled_balanced <-
  pooled_balanced_result$summary %>%
  select(
    scenario,
    framework,
    indicator,
    boundary,
    probability,
    value
  )

write_csv_safe(
  bind_rows(
    pooled_raw,
    pooled_balanced
  ),
  file.path(
    SENSITIVITY_OUTPUT_DIR,
    "pooled_reference_thresholds.csv"
  )
)

score_pooled_raw <-
  score_validation_data(
    bioindicateur2,
    pooled_raw,
    frameworks =
      "POOLED_RAW",
    indicators =
      MANUSCRIPT_INDICATORS,
    scenario_name =
      "pooled_raw"
  )

score_pooled_balanced <-
  score_validation_data(
    bioindicateur2,
    pooled_balanced,
    frameworks =
      "POOLED_BALANCED",
    indicators =
      MANUSCRIPT_INDICATORS,
    scenario_name =
      "pooled_balanced"
  )

score_pooled <- bind_rows(
  score_pooled_raw,
  score_pooled_balanced
)

pooled_agreement <-
  score_pooled %>%
  select(
    site_id,
    indicator,
    framework,
    score
  ) %>%
  pivot_wider(
    names_from = "framework",
    values_from = "score"
  ) %>%
  group_by(
    .data$indicator
  ) %>%
  summarise(
    n =
      sum(
        complete.cases(
          .data$POOLED_RAW,
          .data$POOLED_BALANCED
        )
      ),
    spearman_rho =
      suppressWarnings(
        cor(
          .data$POOLED_RAW,
          .data$POOLED_BALANCED,
          method = "spearman",
          use = "complete.obs"
        )
      ),
    exact_agreement =
      mean(
        .data$POOLED_RAW ==
          .data$POOLED_BALANCED,
        na.rm = TRUE
      ),
    within_one_class =
      mean(
        abs(
          .data$POOLED_RAW -
            .data$POOLED_BALANCED
        ) <= 1,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

write_csv_safe(
  pooled_agreement,
  file.path(
    SENSITIVITY_OUTPUT_DIR,
    "pooled_raw_vs_balanced_agreement.csv"
  )
)

# ==============================================================================
# 5. Raw treatment-control contrasts translated into score space
# ==============================================================================

message(
  "[03] Treatment-control contrast translation..."
)

score_lookup <-
  score_full %>%
  select(
    site_id,
    indicator,
    framework,
    score
  )

contrast_scores <-
  tidyr::crossing(
    contrasts,
    framework =
      FRAMEWORKS
  ) %>%
  left_join(
    score_lookup %>%
      rename(
        treatment_site =
          site_id,
        score_treatment =
          score
      ),
    by = c(
      "treatment_site",
      "indicator",
      "framework"
    )
  ) %>%
  left_join(
    score_lookup %>%
      rename(
        control_site =
          site_id,
        score_control =
          score
      ),
    by = c(
      "control_site",
      "indicator",
      "framework"
    )
  ) %>%
  mutate(
    delta_score =
      .data$score_treatment -
      .data$score_control,
    LRR_status = case_when(
      !is.finite(.data$LRR) ~
        "undefined",
      abs(.data$LRR) <
        LRR_TOLERANCE ~
        "numerical_zero",
      .data$LRR > 0 ~
        "positive",
      .data$LRR < 0 ~
        "negative",
      TRUE ~
        "zero"
    ),
    directional_LRR =
      is.finite(.data$LRR) &
      abs(.data$LRR) >=
        LRR_TOLERANCE,
    translation = case_when(
      !.data$directional_LRR ~
        NA_character_,
      .data$delta_score == 0 ~
        "compressed_to_zero",
      sign(.data$delta_score) ==
        sign(.data$LRR) ~
        "same_direction",
      TRUE ~
        "opposite_direction"
    )
  )

write_csv_safe(
  contrast_scores,
  file.path(
    SENSITIVITY_OUTPUT_DIR,
    "raw_LRR_and_delta_score.csv"
  )
)

directional_by_framework <-
  contrast_scores %>%
  filter(
    .data$directional_LRR
  ) %>%
  group_by(
    .data$framework
  ) %>%
  summarise(
    n_contrasts_evaluated = n(),
    nonzero_delta_same_direction =
      sum(
        .data$translation ==
          "same_direction"
      ),
    compressed_to_delta_zero =
      sum(
        .data$translation ==
          "compressed_to_zero"
      ),
    delta_opposite_to_LRR =
      sum(
        .data$translation ==
          "opposite_direction"
      ),
    same_direction_percent =
      100 *
      .data$nonzero_delta_same_direction /
      .data$n_contrasts_evaluated,
    compressed_to_zero_percent =
      100 *
      .data$compressed_to_delta_zero /
      .data$n_contrasts_evaluated,
    .groups = "drop"
  )

directional_all <-
  contrast_scores %>%
  filter(
    .data$directional_LRR
  ) %>%
  summarise(
    framework = "All",
    n_contrasts_evaluated = n(),
    nonzero_delta_same_direction =
      sum(
        .data$translation ==
          "same_direction"
      ),
    compressed_to_delta_zero =
      sum(
        .data$translation ==
          "compressed_to_zero"
      ),
    delta_opposite_to_LRR =
      sum(
        .data$translation ==
          "opposite_direction"
      ),
    same_direction_percent =
      100 *
      .data$nonzero_delta_same_direction /
      .data$n_contrasts_evaluated,
    compressed_to_zero_percent =
      100 *
      .data$compressed_to_delta_zero /
      .data$n_contrasts_evaluated
  )

directional_summary <-
  bind_rows(
    directional_by_framework,
    directional_all
  )

write_csv_safe(
  directional_summary,
  file.path(
    SENSITIVITY_OUTPUT_DIR,
    "directional_contrast_summary.csv"
  )
)

# ------------------------------------------------------------------------------
# Metadata
# ------------------------------------------------------------------------------

run_metadata <- tibble(
  parameter = c(
    "seed",
    "repeats",
    "balanced_n",
    "winsor_lower",
    "winsor_upper",
    "LRR_numerical_zero_tolerance",
    "manuscript_indicators"
  ),
  value = c(
    SEED,
    N_REPEATS,
    N_BALANCED,
    WINSOR_PROBS[1],
    WINSOR_PROBS[2],
    LRR_TOLERANCE,
    paste(
      MANUSCRIPT_INDICATORS,
      collapse = ","
    )
  )
)

write_csv_safe(
  run_metadata,
  file.path(
    SENSITIVITY_OUTPUT_DIR,
    "run_metadata.csv"
  )
)

message(
  "[03] Sensitivity analyses complete."
)
message(
  "[03] Results written to: ",
  SENSITIVITY_OUTPUT_DIR
)
