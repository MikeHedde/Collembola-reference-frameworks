#!/usr/bin/env Rscript

# ==============================================================================
# 05 - Supplementary material
# ==============================================================================
#
# Generates the final supplementary figures and tables from:
#   - analysis-ready site-level indicator values;
#   - outputs of the main analysis (02);
#   - outputs of the sensitivity analyses (03).
#
# Supplementary structure:
#
#   Figure S1  Indicator distributions, all land-cover categories
#   Table  S1  Land-cover composition and sample sizes
#   Table  S2  Descriptive statistics, all sites
#   Table  S3  Cramer-von Mises tests, all sites
#
#   Figure S2  Indicator distributions, agricultural sites
#   Table  S4  Descriptive statistics, agricultural sites
#   Table  S5  Cramer-von Mises tests, agricultural sites
#
#   Figure S3  Construction of the seven-level scoring system
#   Figure S4  Framework-specific class boundaries
#   Figure S5  Bioindicateur2 scores for four indicators
#   Table  S6  Friedman and Nemenyi tests
#
#   Figure S6  Equal-size n = 58 threshold sensitivity
#   Table  S7  Pairwise threshold differences at equal sample size
#   Table  S8  Extreme-value / winsorisation sensitivity
#
#   Figure S7  Cross-framework score agreement
#   Table  S9  Cross-framework ranking and classification concordance
#
#   Table S10  Pooled raw versus equal-weight reference
#   Table S11  Translation of directional raw contrasts into score space
#
# ==============================================================================

source("R/project_paths.R")
source_repo("R/io_helpers.R")
source_repo("R/indicator_scoring.R")

required_packages(c(
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "cowplot",
  "scales",
  "e1071",
  "cramer",
  "PMCMRplus"
))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(cowplot)
})

message("[05] Generating supplementary material...")

dir.create(
  SUPPLEMENTARY_OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------------------------
# Inputs
# ------------------------------------------------------------------------------

reference_file <- file.path(
  ANALYSIS_DATA_DIR,
  "reference_indicator_values.csv"
)

bio_file <- file.path(
  ANALYSIS_DATA_DIR,
  "bioindicateur2_indicator_values.csv"
)

sensitivity_files <- c(
  balanced_summary = file.path(
    SENSITIVITY_OUTPUT_DIR,
    "balanced_n58_threshold_stability.csv"
  ),
  balanced_pairwise = file.path(
    SENSITIVITY_OUTPUT_DIR,
    "balanced_n58_pairwise_threshold_differences.csv"
  ),
  extreme = file.path(
    SENSITIVITY_OUTPUT_DIR,
    "extreme_observation_diagnostics.csv"
  ),
  outlier_scores = file.path(
    SENSITIVITY_OUTPUT_DIR,
    "outlier_score_sensitivity.csv"
  ),
  outlier_thresholds = file.path(
    SENSITIVITY_OUTPUT_DIR,
    "outlier_threshold_sensitivity.csv"
  ),
  agreement = file.path(
    SENSITIVITY_OUTPUT_DIR,
    "cross_framework_score_agreement.csv"
  ),
  pooled_thresholds = file.path(
    SENSITIVITY_OUTPUT_DIR,
    "pooled_reference_thresholds.csv"
  ),
  pooled_agreement = file.path(
    SENSITIVITY_OUTPUT_DIR,
    "pooled_raw_vs_balanced_agreement.csv"
  ),
  directional = file.path(
    SENSITIVITY_OUTPUT_DIR,
    "directional_contrast_summary.csv"
  )
)

assert_files_exist(
  c(
    reference_file,
    bio_file,
    sensitivity_files
  ),
  context = "Supplementary-material input"
)

reference <- read.csv(
  reference_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

bioindicateur2 <- read.csv(
  bio_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

reference_all <- reference %>%
  filter(
    .data$framework %in%
      FRAMEWORKS
  ) %>%
  mutate(
    framework = factor(
      .data$framework,
      levels = FRAMEWORKS
    )
  )

reference_agri <- reference_all %>%
  filter(
    .data$CLC_niveau1 ==
      AGRICULTURAL_CLC1
  )

# ------------------------------------------------------------------------------
# Labels
# ------------------------------------------------------------------------------

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
  RSr = "Rarefied species richness",
  Shannon = "Shannon diversity",
  Pielou = "Pielou evenness",
  ab = "Density per m²"
)

ALL_INDICATORS <- c(
  "RSr",
  "Shannon",
  "Pielou",
  "ab"
)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

save_plot_pair <- function(
  plot,
  stem,
  width,
  height,
  dpi = 300
) {

  ggsave(
    file.path(
      SUPPLEMENTARY_OUTPUT_DIR,
      paste0(stem, ".pdf")
    ),
    plot = plot,
    width = width,
    height = height,
    units = "in"
  )

  ggsave(
    file.path(
      SUPPLEMENTARY_OUTPUT_DIR,
      paste0(stem, ".png")
    ),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi
  )
}


descriptive_statistics <- function(
  data,
  indicators = ALL_INDICATORS
) {

  bind_rows(
    lapply(
      FRAMEWORKS,
      function(fw) {

        d <- data %>%
          filter(
            as.character(
              .data$framework
            ) == fw
          )

        bind_rows(
          lapply(
            indicators,
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
                skewness =
                  e1071::skewness(
                    x,
                    na.rm = TRUE
                  ),
                kurtosis_excess =
                  e1071::kurtosis(
                    x,
                    na.rm = TRUE
                  )
              )
            }
          )
        )
      }
    )
  )
}


pairwise_cramer <- function(
  data,
  indicator,
  seed = 1234L,
  replicates = 10000L
) {

  pairs <- combn(
    FRAMEWORKS,
    2,
    simplify = FALSE
  )

  set.seed(seed)

  out <- bind_rows(
    lapply(
      pairs,
      function(pair) {

        x1 <- data %>%
          filter(
            as.character(
              .data$framework
            ) == pair[1]
          ) %>%
          pull(
            all_of(indicator)
          )

        x2 <- data %>%
          filter(
            as.character(
              .data$framework
            ) == pair[2]
          ) %>%
          pull(
            all_of(indicator)
          )

        x1 <- x1[
          is.finite(x1)
        ]

        x2 <- x2[
          is.finite(x2)
        ]

        tst <- cramer::cramer.test(
          x1,
          x2,
          replicates = replicates,
          sim = "ordinary"
        )

        tibble(
          indicator = indicator,
          framework_1 = pair[1],
          framework_2 = pair[2],
          W =
            as.numeric(
              tst$statistic
            ),
          p_raw =
            as.numeric(
              tst$p.value
            )
        )
      }
    )
  )

  out %>%
    mutate(
      p_bonferroni =
        stats::p.adjust(
          .data$p_raw,
          method = "bonferroni"
        ),
      significant_0_05 =
        .data$p_bonferroni < 0.05
    )
}


pairwise_cramer_all <- function(
  data,
  indicators = ALL_INDICATORS
) {

  bind_rows(
    lapply(
      indicators,
      function(ind) {

        pairwise_cramer(
          data,
          ind,
          seed = 1234L,
          replicates = 10000L
        )
      }
    )
  )
}


# Compact-letter display derived from non-significant framework pairs.
compact_letters_from_pairs <- function(
  pair_tbl,
  p_col,
  alpha = 0.05,
  group_order = FRAMEWORKS
) {

  g <- group_order
  n <- length(g)

  nonsig <- diag(
    TRUE,
    n
  )

  dimnames(nonsig) <- list(
    g,
    g
  )

  for (i in seq_len(
    nrow(pair_tbl)
  )) {

    a <-
      as.character(
        pair_tbl$framework_1[i]
      )

    b <-
      as.character(
        pair_tbl$framework_2[i]
      )

    p <-
      pair_tbl[[p_col]][i]

    keep <-
      is.finite(p) &&
      p >= alpha

    nonsig[a, b] <- keep
    nonsig[b, a] <- keep
  }

  subsets <- unlist(
    lapply(
      seq_len(n),
      function(k) {
        combn(
          g,
          k,
          simplify = FALSE
        )
      }
    ),
    recursive = FALSE
  )

  is_clique <- vapply(
    subsets,
    function(s) {

      if (length(s) <= 1L) {
        return(TRUE)
      }

      all(
        nonsig[
          s,
          s,
          drop = FALSE
        ]
      )
    },
    logical(1)
  )

  cliques <-
    subsets[
      is_clique
    ]

  is_maximal <- vapply(
    seq_along(cliques),
    function(i) {

      s <- cliques[[i]]

      !any(
        vapply(
          seq_along(cliques),
          function(j) {

            if (i == j) {
              return(FALSE)
            }

            z <- cliques[[j]]

            length(z) >
              length(s) &&
              all(
                s %in% z
              )
          },
          logical(1)
        )
      )
    },
    logical(1)
  )

  cliques <-
    cliques[
      is_maximal
    ]

  clique_key <- vapply(
    cliques,
    function(s) {

      paste(
        sprintf(
          "%02d",
          sort(
            match(
              s,
              g
            )
          )
        ),
        collapse = "-"
      )
    },
    character(1)
  )

  cliques <-
    cliques[
      order(
        clique_key
      )
    ]

  letter_pool <- c(
    letters,
    paste0(
      "a",
      letters
    ),
    paste0(
      "b",
      letters
    )
  )

  out <- setNames(
    rep(
      "",
      n
    ),
    g
  )

  for (i in seq_along(
    cliques
  )) {

    out[
      cliques[[i]]
    ] <- paste0(
      out[
        cliques[[i]]
      ],
      letter_pool[i]
    )
  }

  out
}


make_distribution_panel <- function(
  data,
  indicator,
  cramer_tbl,
  seed = 123L
) {

  d <- data %>%
    transmute(
      framework = factor(
        as.character(
          .data$framework
        ),
        levels =
          rev(
            FRAMEWORKS
          )
      ),
      value =
        .data[[indicator]]
    )

  pair_tbl <- cramer_tbl %>%
    filter(
      .data$indicator ==
        .env$indicator
    )

  letters_this <-
    compact_letters_from_pairs(
      pair_tbl,
      p_col =
        "p_bonferroni"
    )

  ann <- tibble(
    framework = factor(
      names(
        letters_this
      ),
      levels =
        rev(
          FRAMEWORKS
        )
    ),
    letter =
      unname(
        letters_this
      )
  )

  rng <- range(
    d$value[
      is.finite(
        d$value
      )
    ],
    na.rm = TRUE
  )

  span <- diff(rng)

  if (
    !is.finite(span) ||
    span == 0
  ) {
    span <- 1
  }

  ann$x <-
    rng[2] +
    0.06 * span

  set.seed(seed)

  p <- ggplot(
    d,
    aes(
      x = .data$value,
      y = .data$framework,
      fill = .data$framework
    )
  ) +
    geom_violin(
      trim = TRUE,
      alpha = 0.70,
      na.rm = TRUE
    ) +
    geom_boxplot(
      width = 0.10,
      fill = "white",
      outlier.shape = NA,
      na.rm = TRUE
    ) +
    geom_point(
      position =
        position_jitter(
          width = 0,
          height = 0.16,
          seed = seed
        ),
      alpha = 0.50,
      size = 0.90,
      na.rm = TRUE
    ) +
    geom_text(
      data = ann,
      aes(
        x = .data$x,
        y = .data$framework,
        label = .data$letter
      ),
      inherit.aes = FALSE,
      size = 3.5
    ) +
    scale_fill_manual(
      values =
        framework_colors,
      guide = "none"
    ) +
    scale_y_discrete(
      labels =
        FRAMEWORK_LABELS
    ) +
    labs(
      x =
        unname(
          indicator_labels[
            indicator
          ]
        ),
      y = NULL
    ) +
    theme_minimal(
      base_size = 9.5
    ) +
    theme(
      panel.grid.minor =
        element_blank()
    )

  if (indicator == "ab") {

    p <- p +
      scale_x_continuous(
        labels =
          scales::label_number(
            big.mark = " "
          ),
        expand =
          expansion(
            mult = c(
              0.02,
              0.10
            )
          )
      )

  } else {

    p <- p +
      scale_x_continuous(
        expand =
          expansion(
            mult = c(
              0.02,
              0.10
            )
          )
      )
  }

  p
}


make_distribution_figure <- function(
  data,
  cramer_tbl,
  stem
) {

  panels <- lapply(
    ALL_INDICATORS,
    function(ind) {
      make_distribution_panel(
        data,
        ind,
        cramer_tbl
      )
    }
  )

  combined <-
    cowplot::plot_grid(
      plotlist = panels,
      labels = c(
        "A",
        "B",
        "C",
        "D"
      ),
      ncol = 2
    )

  save_plot_pair(
    combined,
    stem,
    width = 11.2,
    height = 8.4
  )
}


extract_nemenyi_pairs <- function(
  pmat,
  frameworks = FRAMEWORKS
) {

  pmat <- as.matrix(
    pmat
  )

  out <- list()
  k <- 1L

  for (
    i in seq_len(
      length(
        frameworks
      ) - 1L
    )
  ) {

    for (
      j in (
        i + 1L
      ):length(
        frameworks
      )
    ) {

      a <- frameworks[i]
      b <- frameworks[j]

      p <- NA_real_

      if (
        a %in%
          rownames(
            pmat
          ) &&
        b %in%
          colnames(
            pmat
          )
      ) {
        p <- pmat[a, b]
      }

      if (
        !is.finite(p) &&
        b %in%
          rownames(
            pmat
          ) &&
        a %in%
          colnames(
            pmat
          )
      ) {
        p <- pmat[b, a]
      }

      out[[k]] <- tibble(
        framework_1 = a,
        framework_2 = b,
        p_nemenyi =
          as.numeric(p)
      )

      k <- k + 1L
    }
  }

  bind_rows(out)
}


run_score_tests <- function(
  score_long,
  indicator
) {

  wide <- score_long %>%
    filter(
      .data$indicator ==
        .env$indicator
    ) %>%
    select(
      site_id,
      framework,
      score
    ) %>%
    pivot_wider(
      names_from =
        "framework",
      values_from =
        "score"
    ) %>%
    drop_na(
      all_of(
        FRAMEWORKS
      )
    )

  mat <- as.matrix(
    wide[
      ,
      FRAMEWORKS,
      drop = FALSE
    ]
  )

  storage.mode(mat) <-
    "numeric"

  rownames(mat) <-
    make.unique(
      as.character(
        wide$site_id
      )
    )

  fr <-
    stats::friedman.test(
      mat
    )

  nm <-
    PMCMRplus::frdAllPairsNemenyiTest(
      mat
    )

  pairs <-
    extract_nemenyi_pairs(
      nm$p.value,
      frameworks =
        FRAMEWORKS
    ) %>%
    mutate(
      indicator = indicator,
      .before = 1
    )

  list(
    friedman = tibble(
      indicator = indicator,
      n_complete_sites =
        nrow(mat),
      statistic =
        as.numeric(
          fr$statistic
        ),
      df =
        as.numeric(
          fr$parameter
        ),
      p_value =
        as.numeric(
          fr$p.value
        )
    ),
    pairs = pairs
  )
}

# ==============================================================================
# Table S1 - land-cover composition
# ==============================================================================

message("[05] Table S1...")

land_cover_levels <- c(
  "Artificialized",
  "Agricultural",
  "Forests and semi-natural environments",
  "Wetlands",
  "Missing / not assigned"
)

table_s1_long <- reference_all %>%
  mutate(
    land_cover =
      case_when(
        .data$CLC_niveau1 == 1 ~
          "Artificialized",
        .data$CLC_niveau1 == 2 ~
          "Agricultural",
        .data$CLC_niveau1 == 3 ~
          "Forests and semi-natural environments",
        .data$CLC_niveau1 == 4 ~
          "Wetlands",
        TRUE ~
          "Missing / not assigned"
      )
  ) %>%
  count(
    .data$land_cover,
    .data$framework,
    name = "n"
  ) %>%
  tidyr::complete(
    land_cover =
      land_cover_levels,
    framework =
      factor(
        FRAMEWORKS,
        levels = FRAMEWORKS
      ),
    fill = list(
      n = 0L
    )
  )

table_s1_wide <- table_s1_long %>%
  mutate(
    land_cover = factor(
      .data$land_cover,
      levels =
        land_cover_levels
    )
  ) %>%
  arrange(
    .data$land_cover
  ) %>%
  pivot_wider(
    names_from =
      "framework",
    values_from =
      "n"
  ) %>%
  mutate(
    land_cover =
      as.character(
        .data$land_cover
      )
  )

total_row <- tibble(
  land_cover = "Total",
  RMQS_Biodiversite =
    sum(
      reference_all$framework ==
        "RMQS_Biodiversite"
    ),
  RMQS_BioDiv_Bretagne =
    sum(
      reference_all$framework ==
        "RMQS_BioDiv_Bretagne"
    ),
  ANDRA =
    sum(
      reference_all$framework ==
        "ANDRA"
    ),
  TIGA_rural =
    sum(
      reference_all$framework ==
        "TIGA_rural"
    )
)

table_s1_wide <- bind_rows(
  table_s1_wide,
  total_row
)

write_csv_safe(
  table_s1_wide,
  file.path(
    SUPPLEMENTARY_OUTPUT_DIR,
    "Table_S1_land_cover_composition.csv"
  )
)

# ==============================================================================
# Figure S1 + Tables S2-S3 - all sites
# ==============================================================================

message("[05] Figure S1 and Tables S2-S3...")

table_s2 <-
  descriptive_statistics(
    reference_all
  )

table_s3 <-
  pairwise_cramer_all(
    reference_all
  )

write_csv_safe(
  table_s2,
  file.path(
    SUPPLEMENTARY_OUTPUT_DIR,
    "Table_S2_descriptive_statistics_all_sites.csv"
  )
)

write_csv_safe(
  table_s3,
  file.path(
    SUPPLEMENTARY_OUTPUT_DIR,
    "Table_S3_cramer_von_mises_all_sites.csv"
  )
)

make_distribution_figure(
  reference_all,
  table_s3,
  "Figure_S1_distributions_all_sites"
)

# ==============================================================================
# Figure S2 + Tables S4-S5 - agricultural reference populations
# ==============================================================================

message("[05] Figure S2 and Tables S4-S5...")

table_s4 <-
  descriptive_statistics(
    reference_agri
  )

table_s5 <-
  pairwise_cramer_all(
    reference_agri
  )

write_csv_safe(
  table_s4,
  file.path(
    SUPPLEMENTARY_OUTPUT_DIR,
    "Table_S4_descriptive_statistics_agricultural.csv"
  )
)

write_csv_safe(
  table_s5,
  file.path(
    SUPPLEMENTARY_OUTPUT_DIR,
    "Table_S5_cramer_von_mises_agricultural.csv"
  )
)

make_distribution_figure(
  reference_agri,
  table_s5,
  "Figure_S2_distributions_agricultural"
)

# ==============================================================================
# Figure S3 - construction of the seven-level scoring system
# ==============================================================================

message("[05] Figure S3...")

illustration_values <- reference_agri %>%
  filter(
    as.character(
      .data$framework
    ) ==
      "RMQS_Biodiversite"
  ) %>%
  pull(
    RSr
  )

illustration_values <-
  illustration_values[
    is.finite(
      illustration_values
    )
  ]

illustration_quantiles <-
  as.numeric(
    stats::quantile(
      illustration_values,
      probs = c(
        0,
        0.2,
        0.4,
        0.6,
        0.8,
        1
      ),
      type = 7,
      names = FALSE
    )
  )

illustration_span <-
  diff(
    range(
      illustration_values
    )
  )

display_min <-
  illustration_quantiles[1] -
  0.20 *
    illustration_span

display_max <-
  illustration_quantiles[6] +
  0.20 *
    illustration_span

display_bounds <- c(
  display_min,
  illustration_quantiles,
  display_max
)

score_bar <- tibble(
  score = factor(
    0:6,
    levels = 0:6
  ),
  xmin =
    head(
      display_bounds,
      -1
    ),
  xmax =
    tail(
      display_bounds,
      -1
    )
)

density_obj <-
  stats::density(
    illustration_values,
    from = display_min,
    to = display_max,
    n = 800
  )

density_df <- tibble(
  x = density_obj$x,
  y = density_obj$y
)

quantile_lines <- tibble(
  x =
    illustration_quantiles,
  label = c(
    "q00",
    "q20",
    "q40",
    "q60",
    "q80",
    "q100"
  )
)

p_s3 <- ggplot() +
  geom_line(
    data = density_df,
    aes(
      x = .data$x,
      y = .data$y
    ),
    linewidth = 0.8
  ) +
  geom_vline(
    data = quantile_lines,
    aes(
      xintercept =
        .data$x
    ),
    linetype = 2,
    linewidth = 0.35
  ) +
  geom_rect(
    data = score_bar,
    aes(
      xmin = .data$xmin,
      xmax = .data$xmax,
      ymin = -0.020,
      ymax = -0.005,
      fill = .data$score
    ),
    colour = "black",
    linewidth = 0.25
  ) +
  geom_text(
    data = score_bar,
    aes(
      x =
        (
          .data$xmin +
            .data$xmax
        ) / 2,
      y = -0.0125,
      label = .data$score
    ),
    size = 3.4
  ) +
  geom_text(
    data = quantile_lines,
    aes(
      x = .data$x,
      y =
        max(
          density_df$y
        ) *
        1.04,
      label = .data$label
    ),
    size = 3.0
  ) +
  scale_fill_manual(
    values =
      score_colors,
    guide = "none"
  ) +
  coord_cartesian(
    xlim = c(
      display_min,
      display_max
    ),
    ylim = c(
      -0.04,
      max(
        density_df$y
      ) *
        1.13
    )
  ) +
  labs(
    x = "Indicator value",
    y =
      "Distribution in the empirical reference population"
  ) +
  theme_minimal(
    base_size = 10
  )

save_plot_pair(
  p_s3,
  "Figure_S3_score_construction",
  width = 10.5,
  height = 5.5
)

write_csv_safe(
  score_bar,
  file.path(
    SUPPLEMENTARY_OUTPUT_DIR,
    "Figure_S3_score_intervals.csv"
  )
)

# ==============================================================================
# Figure S4 - framework-specific score boundaries
# ==============================================================================

message("[05] Figure S4...")

all_breaks <- build_break_table(
  reference_agri,
  frameworks =
    FRAMEWORKS,
  indicators =
    ALL_INDICATORS,
  scenario =
    "reference_framework"
)

write_csv_safe(
  all_breaks,
  file.path(
    SUPPLEMENTARY_OUTPUT_DIR,
    "Figure_S4_framework_thresholds.csv"
  )
)

make_rectangles <- function(
  indicator
) {

  x_all <-
    reference_agri[[indicator]]

  x_all <-
    x_all[
      is.finite(
        x_all
      )
    ]

  global_min <-
    min(
      0,
      min(
        x_all
      )
    )

  global_max <-
    max(
      x_all
    )

  global_range <-
    global_max -
    global_min

  display_max <-
    global_max +
    max(
      0.05 *
        global_range,
      .Machine$double.eps
    )

  if (indicator == "Pielou") {
    display_max <- max(
      1,
      global_max
    )
  }

  bind_rows(
    lapply(
      FRAMEWORKS,
      function(fw) {

        q <- all_breaks %>%
          filter(
            .data$framework ==
              fw,
            .data$indicator ==
              .env$indicator
          ) %>%
          arrange(
            .data$probability
          ) %>%
          pull(
            value
          )

        boundaries <- c(
          global_min,
          q,
          display_max
        )

        tibble(
          framework = fw,
          indicator = indicator,
          score = factor(
            0:6,
            levels = 0:6
          ),
          xmin =
            head(
              boundaries,
              -1
            ),
          xmax =
            tail(
              boundaries,
              -1
            )
        ) %>%
          filter(
            .data$xmax >
              .data$xmin
          )
      }
    )
  )
}

s4_rectangles <- bind_rows(
  lapply(
    ALL_INDICATORS,
    make_rectangles
  )
)

write_csv_safe(
  s4_rectangles,
  file.path(
    SUPPLEMENTARY_OUTPUT_DIR,
    "Figure_S4_score_rectangles.csv"
  )
)

make_s4_panel <- function(
  indicator,
  show_legend = FALSE
) {

  d <- s4_rectangles %>%
    filter(
      .data$indicator ==
        .env$indicator
    ) %>%
    mutate(
      framework = factor(
        .data$framework,
        levels =
          rev(
            FRAMEWORKS
          )
      )
    )

  p <- ggplot(
    d
  ) +
    geom_rect(
      aes(
        xmin = .data$xmin,
        xmax = .data$xmax,
        ymin =
          as.numeric(
            .data$framework
          ) - 0.38,
        ymax =
          as.numeric(
            .data$framework
          ) + 0.38,
        fill = .data$score
      ),
      colour = "black",
      linewidth = 0.20
    ) +
    scale_y_continuous(
      breaks =
        seq_along(
          rev(
            FRAMEWORKS
          )
        ),
      labels =
        FRAMEWORK_LABELS[
          rev(
            FRAMEWORKS
          )
        ]
    ) +
    scale_fill_manual(
      values =
        score_colors,
      name =
        "Assigned score",
      drop = FALSE
    ) +
    labs(
      x =
        unname(
          indicator_labels[
            indicator
          ]
        ),
      y = NULL
    ) +
    theme_minimal(
      base_size = 9
    ) +
    theme(
      panel.grid.major.y =
        element_blank(),
      panel.grid.minor =
        element_blank(),
      legend.position =
        if (
          show_legend
        ) {
          "right"
        } else {
          "none"
        }
    )

  if (indicator == "ab") {
    p <- p +
      scale_x_continuous(
        labels =
          scales::label_number(
            big.mark = " "
          )
      )
  }

  p
}

p_s4 <- cowplot::plot_grid(
  make_s4_panel("RSr"),
  make_s4_panel("Shannon"),
  make_s4_panel("Pielou"),
  make_s4_panel(
    "ab",
    show_legend = TRUE
  ),
  labels = c(
    "A",
    "B",
    "C",
    "D"
  ),
  ncol = 1,
  rel_heights = c(
    1,
    1,
    1,
    1.05
  )
)

save_plot_pair(
  p_s4,
  "Figure_S4_framework_specific_score_boundaries",
  width = 10.5,
  height = 10
)

# ==============================================================================
# Figure S5 + Table S6 - Bioindicateur2 scores and statistical comparison
# ==============================================================================

message("[05] Figure S5 and Table S6...")

scores_four <- score_validation_data(
  bioindicateur2,
  all_breaks,
  frameworks =
    FRAMEWORKS,
  indicators =
    ALL_INDICATORS,
  scenario_name =
    "reference_framework"
)

write_csv_safe(
  scores_four,
  file.path(
    SUPPLEMENTARY_OUTPUT_DIR,
    "Figure_S5_Bioindicateur2_scores_long.csv"
  )
)

score_tests <- lapply(
  ALL_INDICATORS,
  function(ind) {
    run_score_tests(
      scores_four,
      ind
    )
  }
)

table_s6_friedman <- bind_rows(
  lapply(
    score_tests,
    `[[`,
    "friedman"
  )
)

table_s6_nemenyi <- bind_rows(
  lapply(
    score_tests,
    `[[`,
    "pairs"
  )
)

write_csv_safe(
  table_s6_friedman,
  file.path(
    SUPPLEMENTARY_OUTPUT_DIR,
    "Table_S6_Friedman_tests.csv"
  )
)

write_csv_safe(
  table_s6_nemenyi,
  file.path(
    SUPPLEMENTARY_OUTPUT_DIR,
    "Table_S6_Nemenyi_pairwise_tests.csv"
  )
)

score_counts <- scores_four %>%
  filter(
    !is.na(
      .data$score
    )
  ) %>%
  mutate(
    framework = factor(
      .data$framework,
      levels =
        FRAMEWORKS
    ),
    score_factor = factor(
      .data$score,
      levels = 6:0
    )
  ) %>%
  count(
    .data$indicator,
    .data$framework,
    .data$score_factor,
    name = "n"
  ) %>%
  complete(
    indicator =
      ALL_INDICATORS,
    framework =
      factor(
        FRAMEWORKS,
        levels =
          FRAMEWORKS
      ),
    score_factor =
      factor(
        6:0,
        levels = 6:0
      ),
    fill = list(
      n = 0L
    )
  )

make_s5_panel <- function(
  indicator,
  show_legend = FALSE
) {

  d <- score_counts %>%
    filter(
      .data$indicator ==
        .env$indicator
    )

  pair_tbl <-
    table_s6_nemenyi %>%
    filter(
      .data$indicator ==
        .env$indicator
    )

  cld <-
    compact_letters_from_pairs(
      pair_tbl,
      p_col =
        "p_nemenyi"
    )

  top <- d %>%
    group_by(
      .data$framework
    ) %>%
    summarise(
      y = sum(
        .data$n
      ),
      .groups = "drop"
    ) %>%
    mutate(
      letter =
        unname(
          cld[
            as.character(
              .data$framework
            )
          ]
        )
    )

  ggplot(
    d,
    aes(
      x = .data$framework,
      y = .data$n,
      fill =
        .data$score_factor
    )
  ) +
    geom_col(
      colour = "black",
      width = 0.70,
      linewidth = 0.25
    ) +
    geom_text(
      data = top,
      aes(
        x = .data$framework,
        y = .data$y + 0.7,
        label = .data$letter
      ),
      inherit.aes = FALSE,
      size = 3.4
    ) +
    scale_fill_manual(
      values =
        score_colors,
      name = "Score",
      drop = FALSE
    ) +
    scale_x_discrete(
      labels =
        FRAMEWORK_LABELS
    ) +
    scale_y_continuous(
      expand =
        expansion(
          mult = c(
            0,
            0.12
          )
        )
    ) +
    labs(
      x =
        "Empirical reference framework",
      y =
        "Number of observations",
      title =
        unname(
          indicator_labels[
            indicator
          ]
        )
    ) +
    theme_minimal(
      base_size = 9.5
    ) +
    theme(
      panel.grid.major.x =
        element_blank(),
      panel.grid.minor =
        element_blank(),
      plot.title =
        element_text(
          hjust = 0.5
        ),
      axis.text.x =
        element_text(
          size = 7.5
        ),
      legend.position =
        if (
          show_legend
        ) {
          "right"
        } else {
          "none"
        }
    )
}

p_s5 <- cowplot::plot_grid(
  make_s5_panel("RSr"),
  make_s5_panel("Shannon"),
  make_s5_panel("Pielou"),
  make_s5_panel(
    "ab",
    show_legend = TRUE
  ),
  labels = c(
    "A",
    "B",
    "C",
    "D"
  ),
  ncol = 2,
  rel_widths = c(
    1,
    1.08
  )
)

save_plot_pair(
  p_s5,
  "Figure_S5_Bioindicateur2_scores_four_indicators",
  width = 11,
  height = 8.2
)

# ==============================================================================
# Figure S6 + Table S7 - equal-size sensitivity
# ==============================================================================

message("[05] Figure S6 and Table S7...")

balanced_summary <- read.csv(
  sensitivity_files[
    "balanced_summary"
  ],
  stringsAsFactors = FALSE
)

balanced_pairwise <- read.csv(
  sensitivity_files[
    "balanced_pairwise"
  ],
  stringsAsFactors = FALSE
)

s6_plot_data <- balanced_summary %>%
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
    framework = factor(
      .data$framework,
      levels =
        FRAMEWORKS
    ),
    indicator_label =
      recode(
        .data$indicator,
        RSr =
          "Rarefied species richness",
        ab =
          "Density"
      )
  )

write_csv_safe(
  s6_plot_data,
  file.path(
    SUPPLEMENTARY_OUTPUT_DIR,
    "Figure_S6_balanced_n58_plot_data.csv"
  )
)

p_s6 <- ggplot(
  s6_plot_data,
  aes(
    x = .data$framework,
    y = .data$median
  )
) +
  geom_errorbar(
    data = s6_plot_data %>%
      filter(
        as.character(
          .data$framework
        ) !=
          "RMQS_Biodiversite"
      ),
    aes(
      ymin = .data$q025,
      ymax = .data$q975
    ),
    width = 0.15
  ) +
  geom_point(
    shape = 16,
    size = 2.4
  ) +
  geom_point(
    aes(
      y =
        .data$full_value
    ),
    shape = 1,
    size = 3
  ) +
  facet_grid(
    indicator_label ~ boundary,
    scales = "free_y"
  ) +
  scale_x_discrete(
    labels =
      FRAMEWORK_LABELS
  ) +
  labs(
    x = NULL,
    y =
      "Empirical score boundary"
  ) +
  theme_bw(
    base_size = 9.5
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 35,
        hjust = 1
      ),
    panel.grid.minor =
      element_blank()
  )

save_plot_pair(
  p_s6,
  "Figure_S6_balanced_n58_threshold_stability",
  width = 11,
  height = 6.5
)

table_s7 <- balanced_pairwise %>%
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
  )

write_csv_safe(
  table_s7,
  file.path(
    SUPPLEMENTARY_OUTPUT_DIR,
    "Table_S7_balanced_n58_pairwise_threshold_differences.csv"
  )
)

# ==============================================================================
# Table S8 - extreme observations and winsorisation
# ==============================================================================

message("[05] Table S8...")

extreme <- read.csv(
  sensitivity_files[
    "extreme"
  ],
  stringsAsFactors = FALSE
)

outlier_scores <- read.csv(
  sensitivity_files[
    "outlier_scores"
  ],
  stringsAsFactors = FALSE
)

outlier_thresholds <- read.csv(
  sensitivity_files[
    "outlier_thresholds"
  ],
  stringsAsFactors = FALSE
)

internal_shift <- outlier_thresholds %>%
  filter(
    .data$boundary %in%
      c(
        "q20",
        "q40",
        "q60",
        "q80"
      )
  ) %>%
  group_by(
    .data$framework,
    .data$indicator
  ) %>%
  summarise(
    max_internal_threshold_shift =
      max(
        abs(
          .data$absolute_shift
        ),
        na.rm = TRUE
      ),
    .groups = "drop"
  )

table_s8 <- extreme %>%
  left_join(
    outlier_scores,
    by = c(
      "framework",
      "indicator"
    )
  ) %>%
  left_join(
    internal_shift,
    by = c(
      "framework",
      "indicator"
    )
  ) %>%
  transmute(
    framework =
      .data$framework,
    indicator =
      .data$indicator,
    reference_sample_size =
      .data$n.x,
    n_extreme_3IQR =
      .data$n_extreme_3IQR,
    extreme_percent =
      100 *
      .data$proportion_extreme_3IQR,
    exact_score_agreement_percent =
      100 *
      .data$exact_agreement,
    within_one_class_percent =
      100 *
      .data$within_one_class,
    max_internal_threshold_shift =
      .data$max_internal_threshold_shift
  )

write_csv_safe(
  table_s8,
  file.path(
    SUPPLEMENTARY_OUTPUT_DIR,
    "Table_S8_extreme_value_winsorisation_sensitivity.csv"
  )
)

# ==============================================================================
# Figure S7 + Table S9 - cross-framework concordance
# ==============================================================================

message("[05] Figure S7 and Table S9...")

agreement <- read.csv(
  sensitivity_files[
    "agreement"
  ],
  stringsAsFactors = FALSE
)

table_s9 <- agreement %>%
  mutate(
    exact_agreement_percent =
      100 *
      .data$exact_agreement,
    within_one_class_percent =
      100 *
      .data$within_one_class,
    broad_class_agreement_percent =
      100 *
      .data$broad_class_agreement
  )

write_csv_safe(
  table_s9,
  file.path(
    SUPPLEMENTARY_OUTPUT_DIR,
    "Table_S9_cross_framework_concordance.csv"
  )
)

agreement_long <- agreement %>%
  select(
    indicator,
    framework_1,
    framework_2,
    exact_agreement,
    within_one_class,
    broad_class_agreement
  ) %>%
  pivot_longer(
    cols = c(
      "exact_agreement",
      "within_one_class",
      "broad_class_agreement"
    ),
    names_to =
      "agreement_metric",
    values_to =
      "agreement"
  ) %>%
  mutate(
    agreement_metric =
      recode(
        .data$agreement_metric,
        exact_agreement =
          "Exact score",
        within_one_class =
          "Within ±1 class",
        broad_class_agreement =
          "Broad class"
      ),
    indicator =
      recode(
        .data$indicator,
        RSr =
          "Rarefied species richness",
        ab =
          "Density"
      ),
    framework_1 = factor(
      .data$framework_1,
      levels =
        FRAMEWORKS
    ),
    framework_2 = factor(
      .data$framework_2,
      levels =
        rev(
          FRAMEWORKS
        )
    )
  )

p_s7 <- ggplot(
  agreement_long,
  aes(
    x = .data$framework_1,
    y = .data$framework_2,
    fill = .data$agreement
  )
) +
  geom_tile(
    colour = "grey75"
  ) +
  geom_text(
    aes(
      label =
        sprintf(
          "%.1f%%",
          100 *
            .data$agreement
        )
    ),
    size = 3
  ) +
  facet_grid(
    indicator ~ agreement_metric
  ) +
  scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(
      0,
      1
    ),
    labels =
      scales::label_percent(
        accuracy = 1
      )
  ) +
  scale_x_discrete(
    labels =
      FRAMEWORK_LABELS
  ) +
  scale_y_discrete(
    labels =
      FRAMEWORK_LABELS
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill =
      "Agreement"
  ) +
  theme_bw(
    base_size = 9
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      ),
    panel.grid =
      element_blank()
  )

save_plot_pair(
  p_s7,
  "Figure_S7_cross_framework_score_agreement",
  width = 11,
  height = 5.5
)

# ==============================================================================
# Table S10 - pooled reference weighting
# ==============================================================================

message("[05] Table S10...")

pooled_thresholds <- read.csv(
  sensitivity_files[
    "pooled_thresholds"
  ],
  stringsAsFactors = FALSE
)

pooled_agreement <- read.csv(
  sensitivity_files[
    "pooled_agreement"
  ],
  stringsAsFactors = FALSE
)

table_s10_thresholds <- pooled_thresholds %>%
  select(
    indicator,
    boundary,
    scenario,
    value
  ) %>%
  pivot_wider(
    names_from =
      "scenario",
    values_from =
      "value"
  ) %>%
  mutate(
    difference_balanced_minus_raw =
      .data$pooled_balanced -
      .data$pooled_raw
  ) %>%
  arrange(
    match(
      .data$indicator,
      MANUSCRIPT_INDICATORS
    ),
    match(
      .data$boundary,
      c(
        "q00",
        "q20",
        "q40",
        "q60",
        "q80",
        "q100"
      )
    )
  )

table_s10_agreement <- pooled_agreement %>%
  mutate(
    exact_agreement_percent =
      100 *
      .data$exact_agreement,
    within_one_class_percent =
      100 *
      .data$within_one_class
  )

write_csv_safe(
  table_s10_thresholds,
  file.path(
    SUPPLEMENTARY_OUTPUT_DIR,
    "Table_S10_pooled_reference_thresholds.csv"
  )
)

write_csv_safe(
  table_s10_agreement,
  file.path(
    SUPPLEMENTARY_OUTPUT_DIR,
    "Table_S10_pooled_score_agreement.csv"
  )
)

# ==============================================================================
# Table S11 - directional raw contrasts translated into score space
# ==============================================================================

message("[05] Table S11...")

directional <- read.csv(
  sensitivity_files[
    "directional"
  ],
  stringsAsFactors = FALSE
)

table_s11 <- directional %>%
  transmute(
    framework =
      .data$framework,
    n_contrasts_evaluated =
      .data$n_contrasts_evaluated,
    nonzero_delta_same_direction =
      .data$nonzero_delta_same_direction,
    compressed_to_delta_zero =
      .data$compressed_to_delta_zero,
    delta_opposite_to_LRR =
      .data$delta_opposite_to_LRR,
    same_direction_percent =
      .data$same_direction_percent,
    compressed_to_zero_percent =
      .data$compressed_to_zero_percent
  )

write_csv_safe(
  table_s11,
  file.path(
    SUPPLEMENTARY_OUTPUT_DIR,
    "Table_S11_directional_contrast_translation.csv"
  )
)

# ==============================================================================
# Internal checks
# ==============================================================================

message("[05] Running consistency checks...")

if (
  sum(
    table_s1_wide[
      table_s1_wide$land_cover ==
        "Total",
      FRAMEWORKS
    ]
  ) != 756
) {
  stop(
    "Table S1 totals do not sum to 756."
  )
}

rmqs_missing <- table_s1_wide %>%
  filter(
    .data$land_cover ==
      "Missing / not assigned"
  ) %>%
  pull(
    RMQS_Biodiversite
  )

if (
  length(
    rmqs_missing
  ) != 1L ||
  rmqs_missing != 3L
) {
  stop(
    "Expected 3 RMQS-Biodiversite observations without CLC level 1."
  )
}

s7_check <- table_s7 %>%
  group_by(
    .data$indicator
  ) %>%
  summarise(
    n_excluding_zero =
      sum(
        .data$interval_excludes_zero
      ),
    .groups = "drop"
  )

if (
  s7_check$n_excluding_zero[
    s7_check$indicator ==
      "RSr"
  ] != 17L ||
  s7_check$n_excluding_zero[
    s7_check$indicator ==
      "ab"
  ] != 15L
) {
  stop(
    "Equal-n sensitivity results do not match the validated 17/24 and 15/24 results."
  )
}

all_row <- table_s11 %>%
  filter(
    .data$framework ==
      "All"
  )

if (
  nrow(all_row) != 1L ||
  all_row$n_contrasts_evaluated != 104L ||
  all_row$nonzero_delta_same_direction != 75L ||
  all_row$compressed_to_delta_zero != 29L ||
  all_row$delta_opposite_to_LRR != 0L
) {
  stop(
    "Table S11 does not reproduce the validated 104 / 75 / 29 / 0 result."
  )
}

# ------------------------------------------------------------------------------
# Manifest and session information
# ------------------------------------------------------------------------------

manifest <- c(
  paste0(
    "Generated: ",
    format(
      Sys.time(),
      tz = "UTC"
    ),
    " UTC"
  ),
  "",
  "Final supplementary material: COMPLETE",
  "",
  "Figure S1  Indicator distributions across all sampled land-cover categories",
  "Table S1   Land-cover composition and sample-size distribution",
  "Table S2   Descriptive statistics across all sampled sites",
  "Table S3   Pairwise Cramer-von Mises tests across all sampled sites",
  "",
  "Figure S2  Indicator distributions in agricultural reference populations",
  "Table S4   Descriptive statistics in agricultural reference populations",
  "Table S5   Pairwise Cramer-von Mises tests in agricultural reference populations",
  "",
  "Figure S3  Construction of the seven-level relative scoring system",
  "Figure S4  Framework-specific class boundaries",
  "Figure S5  Bioindicateur2 scores for four Collembola metrics",
  "Table S6   Friedman and Nemenyi score comparisons",
  "",
  "Figure S6  Equal-size n=58 threshold sensitivity",
  "Table S7   Pairwise threshold differences under equal sample size",
  "Table S8   Extreme-value and winsorisation sensitivity",
  "",
  "Figure S7  Cross-framework score agreement",
  "Table S9   Ranking and classification concordance",
  "",
  "Table S10  Raw pooled versus equal-weight pooled reference",
  "Table S11  Translation of directional raw contrasts into ordinal score differences",
  "",
  "Cramer-von Mises tests use 10,000 ordinary resampling replicates and Bonferroni adjustment.",
  "Equal-size sensitivity uses 1,000 repeated samples of n=58 per reference framework.",
  "LRR numerical-zero tolerance is defined in analysis/03_sensitivity_analyses.R."
)

write_text_lines(
  manifest,
  file.path(
    SUPPLEMENTARY_OUTPUT_DIR,
    "SUPPLEMENTARY_MATERIAL_MANIFEST.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    SUPPLEMENTARY_OUTPUT_DIR,
    "sessionInfo_supplementary.txt"
  )
)

message(
  "[05] Supplementary material written to: ",
  SUPPLEMENTARY_OUTPUT_DIR
)

message(
  "[05] All consistency checks PASSED."
)
