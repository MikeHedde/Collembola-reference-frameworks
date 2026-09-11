#!/usr/bin/env Rscript

# ==============================================================================
# 04 - Manuscript figures
# ==============================================================================
#
# Regenerates the four figures used in the revised manuscript from the
# analysis-ready data and the validated outputs of scripts 02 and 03.
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
  "maps",
  "scales"
))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(cowplot)
})

message("[04] Generating manuscript figures...")

dir.create(
  MANUSCRIPT_OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------------------------
# Required inputs
# ------------------------------------------------------------------------------

reference_file <- file.path(
  ANALYSIS_DATA_DIR,
  "reference_indicator_values.csv"
)

map_file <- file.path(
  ANALYSIS_DATA_DIR,
  "map_points.csv"
)

scores_file <- file.path(
  ANALYSIS_OUTPUT_DIR,
  "bioindicateur2_scores.csv"
)

contrast_file <- file.path(
  SENSITIVITY_OUTPUT_DIR,
  "raw_LRR_and_delta_score.csv"
)

assert_files_exist(
  c(
    reference_file,
    map_file,
    scores_file,
    contrast_file
  ),
  context = "Figure-generation input"
)

reference <- read.csv(
  reference_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

map_points <- read.csv(
  map_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

bio_scores <- read.csv(
  scores_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

contrast_scores <- read.csv(
  contrast_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

reference_agri <- reference %>%
  filter(
    .data$framework %in% FRAMEWORKS,
    .data$CLC_niveau1 == AGRICULTURAL_CLC1
  ) %>%
  mutate(
    framework = factor(
      .data$framework,
      levels = FRAMEWORKS
    )
  )

# ------------------------------------------------------------------------------
# Stable labels and palettes
# ------------------------------------------------------------------------------

framework_labels_fig <- FRAMEWORK_LABELS

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

# Retains the palette of the final treatment-control contrast figure.
contrast_framework_colors <- c(
  ANDRA = "#F8766D",
  RMQS_Biodiversite = "#7CAE00",
  RMQS_BioDiv_Bretagne = "#00BFC4",
  TIGA_rural = "#C77CFF"
)

contrast_framework_order <- c(
  "ANDRA",
  "RMQS_Biodiversite",
  "RMQS_BioDiv_Bretagne",
  "TIGA_rural"
)

# ------------------------------------------------------------------------------
# Export helper
# ------------------------------------------------------------------------------

save_plot_pair <- function(
  plot,
  stem,
  width,
  height,
  dpi = 300
) {

  pdf_path <- file.path(
    MANUSCRIPT_OUTPUT_DIR,
    paste0(stem, ".pdf")
  )

  png_path <- file.path(
    MANUSCRIPT_OUTPUT_DIR,
    paste0(stem, ".png")
  )

  ggsave(
    pdf_path,
    plot = plot,
    width = width,
    height = height,
    units = "in"
  )

  ggsave(
    png_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi
  )

  invisible(
    c(
      pdf = pdf_path,
      png = png_path
    )
  )
}

# ==============================================================================
# Figure 1 - spatial distribution of datasets
# ==============================================================================

message("[04] Figure 1...")

map_order <- c(
  FRAMEWORKS,
  "Bioindicateur2"
)

map_colors <- c(
  framework_colors,
  Bioindicateur2 = "#BEBEBE"
)

map_labels <- c(
  framework_labels_fig,
  Bioindicateur2 = "Bioindicateur2"
)

map_points <- map_points %>%
  filter(
    is.finite(.data$longitude),
    is.finite(.data$latitude),
    between(.data$longitude, -7, 12),
    between(.data$latitude, 40, 52.5)
  ) %>%
  mutate(
    framework = factor(
      .data$framework,
      levels = map_order
    )
  )

map_counts <- map_points %>%
  count(
    .data$framework,
    name = "n_mapped_sites",
    .drop = FALSE
  )

write_csv_safe(
  map_counts,
  file.path(
    MANUSCRIPT_OUTPUT_DIR,
    "Figure_1_map_counts.csv"
  )
)

france_map <- ggplot2::map_data(
  "world"
) %>%
  filter(
    .data$region == "France"
  )

p_map <- ggplot() +
  geom_polygon(
    data = france_map,
    aes(
      x = .data$long,
      y = .data$lat,
      group = .data$group
    ),
    fill = "grey94",
    colour = "grey55",
    linewidth = 0.25
  ) +
  geom_point(
    data = map_points,
    aes(
      x = .data$longitude,
      y = .data$latitude,
      colour = .data$framework
    ),
    size = 1.8,
    alpha = 0.85
  ) +
  annotate(
    "segment",
    x = 3.9,
    xend = 6.35,
    y = 42.1,
    yend = 42.1,
    linewidth = 1.1
  ) +
  annotate(
    "segment",
    x = 3.9,
    xend = 3.9,
    y = 42.02,
    yend = 42.18,
    linewidth = 0.8
  ) +
  annotate(
    "segment",
    x = 6.35,
    xend = 6.35,
    y = 42.02,
    yend = 42.18,
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = 5.12,
    y = 42.35,
    label = "200 km",
    size = 3.2
  ) +
  annotate(
    "text",
    x = 8.15,
    y = 49.9,
    label = "N",
    size = 5
  ) +
  annotate(
    "segment",
    x = 8.15,
    xend = 8.15,
    y = 48.8,
    yend = 49.65,
    arrow = grid::arrow(
      length = grid::unit(
        0.18,
        "inches"
      )
    ),
    linewidth = 0.7
  ) +
  scale_colour_manual(
    values = map_colors,
    labels = map_labels,
    drop = FALSE,
    name = "Dataset"
  ) +
  coord_quickmap(
    xlim = c(-5.8, 10.0),
    ylim = c(41.0, 51.2),
    expand = FALSE
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_void(
    base_size = 11
  ) +
  theme(
    legend.position = c(
      0.79,
      0.32
    ),
    legend.background =
      element_rect(
        fill = scales::alpha(
          "white",
          0.82
        ),
        colour = NA
      ),
    legend.title =
      element_text(
        face = "bold"
      ),
    plot.margin =
      margin(
        4,
        4,
        4,
        4
      )
  )

save_plot_pair(
  p_map,
  "Figure_1_reference_framework_map",
  width = 7.0,
  height = 6.0
)

# ==============================================================================
# Figure 2 - empirical reference distributions
# ==============================================================================

message("[04] Figure 2...")

fig2_data <- reference_agri %>%
  transmute(
    framework = factor(
      as.character(
        .data$framework
      ),
      levels = rev(FRAMEWORKS)
    ),
    RSr = .data$RSr,
    ab = .data$ab
  )

write_csv_safe(
  fig2_data,
  file.path(
    MANUSCRIPT_OUTPUT_DIR,
    "Figure_2_distribution_values.csv"
  )
)

# Compact-letter annotations obtained from the final pairwise
# Cramér-von Mises comparisons with Bonferroni adjustment.
fig2_letters_rsr <- tibble(
  framework = factor(
    FRAMEWORKS,
    levels = rev(FRAMEWORKS)
  ),
  letter = c(
    "a",
    "b",
    "b",
    "a"
  )
)

fig2_letters_ab <- tibble(
  framework = factor(
    FRAMEWORKS,
    levels = rev(FRAMEWORKS)
  ),
  letter = c(
    "a",
    "b",
    "a",
    "a"
  )
)

write_csv_safe(
  bind_rows(
    fig2_letters_rsr %>%
      mutate(
        indicator = "RSr"
      ),
    fig2_letters_ab %>%
      mutate(
        indicator = "ab"
      )
  ),
  file.path(
    MANUSCRIPT_OUTPUT_DIR,
    "Figure_2_significance_letters.csv"
  )
)

set.seed(123)

p_rsr <- ggplot(
  fig2_data,
  aes(
    x = .data$RSr,
    y = .data$framework,
    fill = .data$framework
  )
) +
  geom_violin(
    trim = TRUE,
    alpha = 0.70
  ) +
  geom_boxplot(
    width = 0.10,
    fill = "white",
    outlier.shape = NA
  ) +
  geom_point(
    position =
      position_jitter(
        width = 0,
        height = 0.16,
        seed = 123
      ),
    alpha = 0.55,
    size = 1.15
  ) +
  geom_text(
    data = fig2_letters_rsr,
    aes(
      x =
        max(
          fig2_data$RSr,
          na.rm = TRUE
        ) * 1.04,
      y = .data$framework,
      label = .data$letter
    ),
    inherit.aes = FALSE,
    size = 3.6
  ) +
  scale_fill_manual(
    values = framework_colors,
    labels = framework_labels_fig,
    guide = "none"
  ) +
  scale_y_discrete(
    labels = framework_labels_fig
  ) +
  scale_x_continuous(
    expand =
      expansion(
        mult = c(
          0.03,
          0.10
        )
      )
  ) +
  labs(
    x = "Rarefied species richness",
    y = NULL
  ) +
  theme_minimal(
    base_size = 10
  ) +
  theme(
    panel.grid.minor =
      element_blank()
  )

set.seed(123)

p_ab <- ggplot(
  fig2_data,
  aes(
    x = .data$ab,
    y = .data$framework,
    fill = .data$framework
  )
) +
  geom_violin(
    trim = TRUE,
    alpha = 0.70
  ) +
  geom_boxplot(
    width = 0.10,
    fill = "white",
    outlier.shape = NA
  ) +
  geom_point(
    position =
      position_jitter(
        width = 0,
        height = 0.16,
        seed = 123
      ),
    alpha = 0.55,
    size = 1.05
  ) +
  geom_text(
    data = fig2_letters_ab,
    aes(
      x =
        max(
          fig2_data$ab,
          na.rm = TRUE
        ) * 1.04,
      y = .data$framework,
      label = .data$letter
    ),
    inherit.aes = FALSE,
    size = 3.6
  ) +
  scale_fill_manual(
    values = framework_colors,
    labels = framework_labels_fig,
    guide = "none"
  ) +
  scale_y_discrete(
    labels = framework_labels_fig
  ) +
  scale_x_continuous(
    labels =
      scales::label_number(
        big.mark = " "
      ),
    expand =
      expansion(
        mult = c(
          0.03,
          0.10
        )
      )
  ) +
  labs(
    x = expression(
      "Density per " * m^2
    ),
    y = NULL
  ) +
  theme_minimal(
    base_size = 10
  ) +
  theme(
    panel.grid.minor =
      element_blank()
  )

p_fig2 <- cowplot::plot_grid(
  p_rsr,
  p_ab,
  labels = c(
    "A",
    "B"
  ),
  ncol = 2,
  rel_widths = c(
    1,
    1
  )
)

save_plot_pair(
  p_fig2,
  "Figure_2_reference_distributions",
  width = 10.2,
  height = 5.0
)

# ==============================================================================
# Figure 3 - Bioindicateur2 score distributions
# ==============================================================================

message("[04] Figure 3...")

score_counts <- bio_scores %>%
  filter(
    !is.na(.data$score)
  ) %>%
  mutate(
    framework = factor(
      .data$framework,
      levels = FRAMEWORKS
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
  tidyr::complete(
    indicator =
      MANUSCRIPT_INDICATORS,
    framework =
      factor(
        FRAMEWORKS,
        levels = FRAMEWORKS
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

write_csv_safe(
  score_counts,
  file.path(
    MANUSCRIPT_OUTPUT_DIR,
    "Figure_3_score_counts.csv"
  )
)

# Compact-letter annotations from the final Friedman + Nemenyi analyses.
fig3_letters_rsr <- c(
  RMQS_Biodiversite = "a",
  RMQS_BioDiv_Bretagne = "bc",
  ANDRA = "c",
  TIGA_rural = "ab"
)

fig3_letters_ab <- c(
  RMQS_Biodiversite = "a",
  RMQS_BioDiv_Bretagne = "b",
  ANDRA = "a",
  TIGA_rural = "a"
)

fig3_letter_table <- bind_rows(
  tibble(
    indicator = "RSr",
    framework =
      names(
        fig3_letters_rsr
      ),
    letter =
      unname(
        fig3_letters_rsr
      )
  ),
  tibble(
    indicator = "ab",
    framework =
      names(
        fig3_letters_ab
      ),
    letter =
      unname(
        fig3_letters_ab
      )
  )
)

write_csv_safe(
  fig3_letter_table,
  file.path(
    MANUSCRIPT_OUTPUT_DIR,
    "Figure_3_significance_letters.csv"
  )
)

make_score_panel <- function(
  indicator_name,
  title_text,
  letter_vector,
  show_legend = FALSE
) {

  d <- score_counts %>%
    filter(
      .data$indicator ==
        indicator_name
    )

  top <- d %>%
    group_by(
      .data$framework
    ) %>%
    summarise(
      y = sum(.data$n),
      .groups = "drop"
    )

  top$letter <-
    unname(
      letter_vector[
        as.character(
          top$framework
        )
      ]
    )

  ggplot(
    d,
    aes(
      x = .data$framework,
      y = .data$n,
      fill = .data$score_factor
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
      size = 3.5
    ) +
    scale_fill_manual(
      values = score_colors,
      name = "Score",
      drop = FALSE
    ) +
    scale_x_discrete(
      labels = framework_labels_fig
    ) +
    scale_y_continuous(
      expand =
        expansion(
          mult = c(
            0,
            0.10
          )
        )
    ) +
    labs(
      x = "Empirical reference framework",
      y = "Number of observations",
      title = title_text
    ) +
    theme_minimal(
      base_size = 10
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
          size = 8
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

p_score_rsr <- make_score_panel(
  "RSr",
  "Rarefied species richness",
  fig3_letters_rsr,
  show_legend = FALSE
)

p_score_ab <- make_score_panel(
  "ab",
  expression(
    "Density per " * m^2
  ),
  fig3_letters_ab,
  show_legend = TRUE
)

p_fig3 <- cowplot::plot_grid(
  p_score_rsr,
  p_score_ab,
  labels = c(
    "A",
    "B"
  ),
  ncol = 2,
  rel_widths = c(
    1,
    1.12
  )
)

save_plot_pair(
  p_fig3,
  "Figure_3_Bioindicateur2_scores",
  width = 10.2,
  height = 5.2
)

# ==============================================================================
# Figure 4 - signed treatment-control score contrasts
# ==============================================================================

message("[04] Figure 4...")

required_sites <- c(
  "QualiAgro",
  "Thil",
  "Yvetot",
  "MetalEurope"
)

fig4_data <- contrast_scores %>%
  filter(
    .data$Nom_site %in%
      required_sites,
    is.finite(
      .data$delta_score
    )
  ) %>%
  mutate(
    Nom_site = factor(
      .data$Nom_site,
      levels =
        sort(
          required_sites
        )
    ),
    indicator_label =
      ifelse(
        .data$indicator == "ab",
        "Density",
        "Species richness"
      ),
    indicator_label =
      factor(
        .data$indicator_label,
        levels = c(
          "Density",
          "Species richness"
        )
      ),
    framework = factor(
      .data$framework,
      levels =
        contrast_framework_order
    )
  )

treatment_levels <- fig4_data %>%
  distinct(
    .data$Nom_site,
    .data$treatment
  ) %>%
  arrange(
    .data$Nom_site,
    .data$treatment
  ) %>%
  pull(
    "treatment"
  ) %>%
  unique()

fig4_data$treatment <- factor(
  fig4_data$treatment,
  levels = treatment_levels
)

write_csv_safe(
  fig4_data,
  file.path(
    MANUSCRIPT_OUTPUT_DIR,
    "Figure_4_treatment_control_score_contrasts.csv"
  )
)

pd <- position_dodge2(
  width = 0.70,
  preserve = "single"
)

p_fig4 <- ggplot(
  fig4_data,
  aes(
    x = .data$delta_score,
    y = .data$treatment,
    group = .data$framework
  )
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.35
  ) +
  geom_errorbar(
    aes(
      xmin =
        pmin(
          0,
          .data$delta_score
        ),
      xmax =
        pmax(
          0,
          .data$delta_score
        ),
      colour =
        .data$framework
    ),
    orientation = "y",
    linewidth = 0.60,
    position = pd
  ) +
  geom_point(
    aes(
      colour =
        .data$framework
    ),
    position = pd,
    size = 2.3
  ) +
  facet_grid(
    Nom_site ~ indicator_label,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_colour_manual(
    values =
      contrast_framework_colors,
    labels =
      framework_labels_fig,
    name =
      "Empirical reference framework"
  ) +
  coord_cartesian(
    xlim = c(
      -5,
      5
    )
  ) +
  labs(
    x = expression(
      Delta *
        " score (treatment - local control)"
    ),
    y = NULL
  ) +
  theme_bw(
    base_size = 9
  ) +
  theme(
    strip.text =
      element_text(
        size = 9
      ),
    axis.text.y =
      element_text(
        size = 6.5
      ),
    panel.grid.minor =
      element_blank(),
    legend.position =
      "right"
  )

save_plot_pair(
  p_fig4,
  "Figure_4_treatment_control_score_contrasts",
  width = 11.0,
  height = 8.0
)

# ==============================================================================
# Figure-generation manifest
# ==============================================================================

manifest <- c(
  paste0(
    "Generated: ",
    format(
      Sys.time(),
      tz = "UTC"
    ),
    " UTC"
  ),
  "Final manuscript figures: COMPLETE",
  "",
  "Figure 1:",
  " - site locations from data/analysis/map_points.csv;",
  " - static mainland-France outline; no external basemap dependency.",
  "",
  "Figure 2:",
  " - agricultural empirical reference populations;",
  " - rarefied species richness and individual density;",
  " - letters represent final Bonferroni-adjusted Cramer-von Mises comparisons.",
  "",
  "Figure 3:",
  " - same 26 Bioindicateur2 observations scored against each framework;",
  " - letters represent final Nemenyi post-hoc comparisons following Friedman tests.",
  "",
  "Figure 4:",
  " - QualiAgro, Thil, Yvetot and MetalEurope;",
  " - signed Delta score retained;",
  " - full empirical reference populations used throughout."
)

write_text_lines(
  manifest,
  file.path(
    MANUSCRIPT_OUTPUT_DIR,
    "FIGURE_GENERATION_MANIFEST.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    MANUSCRIPT_OUTPUT_DIR,
    "sessionInfo_manuscript_figures.txt"
  )
)

message(
  "[04] Manuscript figures written to: ",
  MANUSCRIPT_OUTPUT_DIR
)
