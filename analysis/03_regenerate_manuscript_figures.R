#!/usr/bin/env Rscript

# ==============================================================================
# 03 - Regenerate the central manuscript figures
# ==============================================================================
# Purpose
#   1. Recreate Figures 1-4 from the cleaned, analysis-ready workflow.
#   2. Preserve a "submitted-style" rendering for reproducibility checks.
#   3. Produce revised candidates where reviewer-driven methodological decisions
#      differ from the submitted analysis (untrimmed TIGA abundance; signed Delta).
#
# Important provenance notes
#   - Figures 2-4 can be regenerated from the supplied analytical data/scripts.
#   - The original Figure 1 used a satellite/GIS basemap, but the corresponding
#     GIS project / tile source was not supplied. We therefore regenerate the
#     sampling locations on a static France outline (no web dependency) and keep
#     the submitted raster under legacy/submitted_figures/ for visual comparison.
#   - The submitted Figure 4 analysis used FOUR Bioindicateur2 sites:
#     QualiAgro, Thil, Yvetot and MetalEurope. This corrects the three-site
#     description accidentally given in the submitted Methods.
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
  "dplyr", "tidyr", "tibble", "ggplot2", "cowplot", "maps", "scales", "vegan"
))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(cowplot)
})

dir.create(MANUSCRIPT_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("[03] Regenerating central manuscript figures...")

# Load/rebuild analysis-ready inputs -------------------------------------------
ref_rds <- file.path(DERIVED_DIR, "reference_frameworks_merged.rds")
bio_rds <- file.path(DERIVED_DIR, "bioindicateur2_curated.rds")
if (!file.exists(ref_rds) || !file.exists(bio_rds)) {
  message("[03] Derived datasets absent; running analysis/01_build_analysis_dataset.R first.")
  source_repo("analysis/01_build_analysis_dataset.R")
}

reference_raw <- readRDS(ref_rds)
validation_raw <- readRDS(bio_rds)
reference <- compute_legacy_indicators(reference_raw)
validation <- compute_legacy_indicators(validation_raw)

reference_agri <- reference %>%
  filter(
    .data$Projet %in% FRAMEWORKS,
    as.character(.data$CLC_niveau1) == as.character(AGRICULTURAL_CLC1)
  ) %>%
  mutate(Projet = factor(as.character(.data$Projet), levels = FRAMEWORKS))

validation_agri <- validation %>%
  filter(as.character(.data$CLC_niveau1) == as.character(AGRICULTURAL_CLC1))

# Stable labels/palettes copied from the submitted plotting logic --------------
framework_labels_fig <- c(
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

# Original ggplot2 default hue palette for four framework levels used in the
# historical four-site Delta-score script. Explicit values make rendering stable.
delta_framework_colors <- c(
  ANDRA = "#F8766D",
  RMQS = "#7CAE00",
  RMQS_Bretagne = "#00BFC4",
  TIGA_rural = "#C77CFF"
)

# Utility exporters ------------------------------------------------------------
save_plot_pair <- function(plot, stem, width, height, dpi = 300) {
  pdf_path <- file.path(MANUSCRIPT_OUTPUT_DIR, paste0(stem, ".pdf"))
  png_path <- file.path(MANUSCRIPT_OUTPUT_DIR, paste0(stem, ".png"))
  ggsave(pdf_path, plot = plot, width = width, height = height, units = "in")
  ggsave(png_path, plot = plot, width = width, height = height, units = "in", dpi = dpi)
  invisible(c(pdf = pdf_path, png = png_path))
}

# ==============================================================================
# Figure 1 - sampling locations
# ==============================================================================
# The submitted caption/figure used the following effective plotting population:
# RMQS Biodiversity: all 96 sites; Bretagne: all 98; ANDRA: agricultural sites;
# TIGA: all 430; Bioindicateur2: agricultural sites (26). We reproduce that
# selection here for a submitted-figure audit, while making the rule explicit.

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

map_ref <- reference %>%
  mutate(
    X_WGS84_num = as_num(.data$X_WGS84),
    Y_WGS84_num = as_num(.data$Y_WGS84),
    Projet_chr = as.character(.data$Projet)
  ) %>%
  filter(
    .data$Projet_chr %in% FRAMEWORKS,
    .data$Projet_chr != "ANDRA" |
      as.character(.data$CLC_niveau1) == as.character(AGRICULTURAL_CLC1),
    is.finite(.data$X_WGS84_num),
    is.finite(.data$Y_WGS84_num)
  ) %>%
  distinct(.data$Projet_chr, .data$Site, .data$X_WGS84_num, .data$Y_WGS84_num) %>%
  transmute(
    framework = .data$Projet_chr,
    Site = as.character(.data$Site),
    longitude = .data$X_WGS84_num,
    latitude = .data$Y_WGS84_num
  )

map_bio <- validation_agri %>%
  mutate(
    longitude = as_num(.data$X_WGS84),
    latitude = as_num(.data$Y_WGS84)
  ) %>%
  filter(is.finite(.data$longitude), is.finite(.data$latitude)) %>%
  distinct(.data$Site, .data$longitude, .data$latitude) %>%
  transmute(
    framework = "Bioindicateur2",
    Site = as.character(.data$Site),
    longitude = .data$longitude,
    latitude = .data$latitude
  )

map_points <- bind_rows(map_ref, map_bio) %>%
  filter(
    between(.data$longitude, -7, 12),
    between(.data$latitude, 40, 52.5)
  )

map_counts <- map_points %>%
  count(.data$framework, name = "n_unique_mapped_sites")
write_csv_safe(map_counts, file.path(MANUSCRIPT_OUTPUT_DIR, "Figure_01_map_counts.csv"))
write_csv_safe(map_points, file.path(MANUSCRIPT_OUTPUT_DIR, "Figure_01_map_points.csv"))

france_map <- ggplot2::map_data("world") %>%
  filter(.data$region == "France")

map_colors <- c(
  framework_colors,
  Bioindicateur2 = "#BEBEBE"
)
map_labels <- c(
  framework_labels_fig,
  Bioindicateur2 = "Bioindicateur2"
)
map_order <- c(FRAMEWORKS, "Bioindicateur2")
map_points$framework <- factor(map_points$framework, levels = map_order)

p_map <- ggplot() +
  geom_polygon(
    data = france_map,
    aes(x = .data$long, y = .data$lat, group = .data$group),
    fill = "grey94", colour = "grey55", linewidth = 0.25
  ) +
  geom_point(
    data = map_points,
    aes(x = .data$longitude, y = .data$latitude, colour = .data$framework),
    size = 1.8, alpha = 0.85
  ) +
  # Approximate 200-km scale bar near 43 N: 2.45 longitude degrees ~ 200 km.
  annotate("segment", x = 3.9, xend = 6.35, y = 42.1, yend = 42.1, linewidth = 1.1) +
  annotate("segment", x = 3.9, xend = 3.9, y = 42.02, yend = 42.18, linewidth = 0.8) +
  annotate("segment", x = 6.35, xend = 6.35, y = 42.02, yend = 42.18, linewidth = 0.8) +
  annotate("text", x = 5.12, y = 42.35, label = "200 km", size = 3.2) +
  annotate("text", x = 8.15, y = 49.9, label = "N", size = 5) +
  annotate(
    "segment", x = 8.15, xend = 8.15, y = 48.8, yend = 49.65,
    arrow = grid::arrow(length = grid::unit(0.18, "inches")), linewidth = 0.7
  ) +
  scale_colour_manual(values = map_colors, labels = map_labels, drop = FALSE, name = "Sites") +
  coord_quickmap(xlim = c(-5.8, 10.0), ylim = c(41.0, 51.2), expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 11) +
  theme(
    legend.position = c(0.79, 0.32),
    legend.background = element_rect(fill = scales::alpha("white", 0.82), colour = NA),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(4, 4, 4, 4)
  )

save_plot_pair(p_map, "Figure_1_reference_framework_map_reconstructed", width = 7.0, height = 6.0)

# ==============================================================================
# Figure 2 - distributions of rarefied richness and density
# ==============================================================================
fig2_data <- reference_agri %>%
  transmute(
    framework = factor(as.character(.data$Projet), levels = rev(FRAMEWORKS)),
    RSr = .data$RSr,
    ab = .data$ab
  )
write_csv_safe(fig2_data, file.path(MANUSCRIPT_OUTPUT_DIR, "Figure_02_distribution_values.csv"))

# Letters are the submitted multiple-comparison annotations. They are preserved
# here as visual-reproduction metadata rather than silently recomputed with a
# potentially different package/version.
fig2_letters_rsr <- tibble(
  framework = factor(FRAMEWORKS, levels = rev(FRAMEWORKS)),
  letter = c("a", "b", "b", "a")
)
fig2_letters_ab <- tibble(
  framework = factor(FRAMEWORKS, levels = rev(FRAMEWORKS)),
  letter = c("a", "b", "a", "a")
)

set.seed(123)
p_rsr <- ggplot(fig2_data, aes(x = .data$RSr, y = .data$framework, fill = .data$framework)) +
  geom_violin(trim = TRUE, alpha = 0.70) +
  geom_boxplot(width = 0.10, fill = "white", outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0, height = 0.16, seed = 123), alpha = 0.55, size = 1.15) +
  geom_text(
    data = fig2_letters_rsr,
    aes(x = max(fig2_data$RSr, na.rm = TRUE) * 1.04, y = .data$framework, label = .data$letter),
    inherit.aes = FALSE, size = 3.6
  ) +
  scale_fill_manual(values = framework_colors, labels = framework_labels_fig, guide = "none") +
  scale_y_discrete(labels = framework_labels_fig) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.10))) +
  labs(x = "Species richness", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank())

set.seed(123)
p_ab <- ggplot(fig2_data, aes(x = .data$ab, y = .data$framework, fill = .data$framework)) +
  geom_violin(trim = TRUE, alpha = 0.70) +
  geom_boxplot(width = 0.10, fill = "white", outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0, height = 0.16, seed = 123), alpha = 0.55, size = 1.05) +
  geom_text(
    data = fig2_letters_ab,
    aes(x = max(fig2_data$ab, na.rm = TRUE) * 1.04, y = .data$framework, label = .data$letter),
    inherit.aes = FALSE, size = 3.6
  ) +
  scale_fill_manual(values = framework_colors, labels = framework_labels_fig, guide = "none") +
  scale_y_discrete(labels = framework_labels_fig) +
  scale_x_continuous(labels = scales::label_number(big.mark = " "), expand = expansion(mult = c(0.03, 0.10))) +
  labs(x = expression("Density per " * m^2), y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank())

p_fig2 <- cowplot::plot_grid(p_rsr, p_ab, labels = c("A", "B"), ncol = 2, rel_widths = c(1, 1))
save_plot_pair(p_fig2, "Figure_2_reference_distributions_submitted_style", width = 10.2, height = 5.0)

# ==============================================================================
# Figure 3 - Bioindicateur2 score distributions across frameworks
# ==============================================================================
full_breaks <- build_break_table(
  reference_agri,
  indicators = MANUSCRIPT_INDICATORS,
  scenario = "full_untrimmed"
)
submitted_breaks <- legacy_tiga_ab_breaks(reference_agri, threshold = 75000)

make_score_bar_figure <- function(break_table, scenario_name, stem, submitted_letters = FALSE) {
  score_long <- score_validation_data(
    validation_agri,
    break_table,
    frameworks = FRAMEWORKS,
    indicators = MANUSCRIPT_INDICATORS,
    scenario_name = scenario_name
  )

  counts <- score_long %>%
    filter(!is.na(.data$score)) %>%
    mutate(
      framework = factor(.data$framework, levels = FRAMEWORKS),
      score_factor = factor(.data$score, levels = 6:0)
    ) %>%
    count(.data$indicator, .data$framework, .data$score_factor, name = "n") %>%
    tidyr::complete(
      indicator = MANUSCRIPT_INDICATORS,
      framework = factor(FRAMEWORKS, levels = FRAMEWORKS),
      score_factor = factor(6:0, levels = 6:0),
      fill = list(n = 0L)
    )

  write_csv_safe(
    counts,
    file.path(MANUSCRIPT_OUTPUT_DIR, paste0(stem, "_score_counts.csv"))
  )

  make_panel <- function(ind, title_text, letter_vec = NULL, show_legend = FALSE) {
    d <- counts %>% filter(.data$indicator == ind)
    p <- ggplot(d, aes(x = .data$framework, y = .data$n, fill = .data$score_factor)) +
      geom_col(colour = "black", width = 0.70, linewidth = 0.25) +
      scale_fill_manual(values = score_colors, name = "Score", drop = FALSE) +
      scale_x_discrete(labels = framework_labels_fig) +
      labs(x = "Scoring scheme used", y = "Number of sites", title = title_text) +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 8),
        legend.position = if (show_legend) "right" else "none"
      )

    if (!is.null(letter_vec)) {
      top <- d %>% group_by(.data$framework) %>% summarise(y = sum(.data$n), .groups = "drop")
      top$letter <- unname(letter_vec[as.character(top$framework)])
      p <- p +
        geom_text(
          data = top,
          aes(x = .data$framework, y = .data$y + 0.7, label = .data$letter),
          inherit.aes = FALSE, size = 3.5
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.10)))
    }
    p
  }

  rsr_letters <- if (submitted_letters) c(
    RMQS_Biodiversite = "a", RMQS_BioDiv_Bretagne = "bc", ANDRA = "c", TIGA_rural = "ab"
  ) else NULL
  ab_letters <- if (submitted_letters) c(
    RMQS_Biodiversite = "a", RMQS_BioDiv_Bretagne = "b", ANDRA = "a", TIGA_rural = "a"
  ) else NULL

  p1 <- make_panel("RSr", "Species richness", rsr_letters, show_legend = FALSE)
  p2 <- make_panel("ab", expression("Density per " * m^2), ab_letters, show_legend = TRUE)
  combined <- cowplot::plot_grid(p1, p2, labels = c("A", "B"), ncol = 2, rel_widths = c(1, 1.12))
  save_plot_pair(combined, stem, width = 10.2, height = 5.2)
  invisible(score_long)
}

score_submitted <- make_score_bar_figure(
  submitted_breaks,
  "submitted_legacy_tiga_ab_lt_75000",
  "Figure_3_Bioindicateur2_scores_submitted_legacy",
  submitted_letters = TRUE
)
score_revised <- make_score_bar_figure(
  full_breaks,
  "revised_full_untrimmed",
  "Figure_3_Bioindicateur2_scores_revised_untrimmed",
  submitted_letters = FALSE
)

# ==============================================================================
# Figure 4 - four-site local contrasts in score space
# ==============================================================================
control_file <- file.path(CONFIG_DIR, "bioindicateur2_controls.csv")
mapping <- read_control_mapping(control_file)
required_sites <- c("QualiAgro", "Thil", "Yvetot", "MetalEurope")
ready_sites <- as.character(mapping$Nom_site[mapping$include])
missing_controls <- setdiff(required_sites, ready_sites)
if (length(missing_controls)) {
  stop(
    "Figure 4 requires explicit control mappings for four sites. Missing: ",
    paste(missing_controls, collapse = ", "),
    "\nEdit config/bioindicateur2_controls.csv."
  )
}

make_delta_data <- function(score_long, scenario_label) {
  res <- make_explicit_lrr_table(
    validation_agri,
    score_long,
    mapping,
    indicators = MANUSCRIPT_INDICATORS
  )
  d <- res$lrr %>%
    filter(.data$Nom_site %in% required_sites, !is.na(.data$delta_score)) %>%
    mutate(
      scenario = scenario_label,
      indicator_label = ifelse(.data$indicator == "ab", "Density", "Species richness"),
      abs_delta = abs(.data$delta_score),
      framework_plot = recode(
        .data$framework,
        RMQS_Biodiversite = "RMQS",
        RMQS_BioDiv_Bretagne = "RMQS_Bretagne",
        ANDRA = "ANDRA",
        TIGA_rural = "TIGA_rural"
      )
    )
  list(data = d, diagnostics = res$diagnostics)
}

delta_submitted_obj <- make_delta_data(score_submitted, "submitted_legacy")
delta_revised_obj <- make_delta_data(score_revised, "revised_untrimmed")
delta_submitted <- delta_submitted_obj$data
delta_revised <- delta_revised_obj$data

write_csv_safe(delta_submitted, file.path(MANUSCRIPT_OUTPUT_DIR, "Figure_04_delta_scores_submitted_legacy.csv"))
write_csv_safe(delta_revised, file.path(MANUSCRIPT_OUTPUT_DIR, "Figure_04_delta_scores_revised_untrimmed.csv"))
write_csv_safe(delta_submitted_obj$diagnostics, file.path(MANUSCRIPT_OUTPUT_DIR, "Figure_04_control_diagnostics.csv"))

rank_delta <- function(d) {
  eps <- 0.01
  d %>%
    group_by(.data$Nom_site, .data$treatment_site, .data$treatment, .data$indicator) %>%
    mutate(
      max_abs = max(.data$abs_delta, na.rm = TRUE),
      is_best = abs(.data$abs_delta - .data$max_abs) <= eps,
      n_best = sum(.data$is_best, na.rm = TRUE),
      best_type = case_when(
        .data$is_best & .data$n_best == 1 ~ "best_solo",
        .data$is_best & .data$n_best > 1 ~ "best_co",
        TRUE ~ NA_character_
      )
    ) %>%
    ungroup()
}

rank_summary <- function(rank_table) {
  n_cases <- rank_table %>%
    distinct(.data$Nom_site, .data$treatment_site, .data$indicator) %>%
    nrow()
  rank_table %>%
    filter(!is.na(.data$best_type)) %>%
    count(.data$framework, .data$best_type, name = "n") %>%
    mutate(n_cases = n_cases, pct = 100 * .data$n / .data$n_cases) %>%
    arrange(.data$best_type, desc(.data$n))
}

rank_submitted <- rank_delta(delta_submitted)
rank_revised <- rank_delta(delta_revised)
write_csv_safe(rank_summary(rank_submitted), file.path(MANUSCRIPT_OUTPUT_DIR, "Figure_04_rank_summary_submitted_legacy.csv"))
write_csv_safe(rank_summary(rank_revised), file.path(MANUSCRIPT_OUTPUT_DIR, "Figure_04_rank_summary_revised_untrimmed.csv"))

# Historical submitted-style Figure 4: absolute Delta score and opacity emphasis
# for best/co-best frameworks, matching scores.R logic.
plot_delta_abs <- function(rank_table, stem) {
  d <- rank_table %>%
    mutate(
      is_best_any = !is.na(.data$best_type),
      Nom_site = factor(.data$Nom_site, levels = sort(required_sites)),
      indicator_label = factor(.data$indicator_label, levels = c("Density", "Species richness"))
    )

  # Keep treatment labels stable within each site. facet_grid with free_y will
  # show only treatments present in the corresponding site.
  d$treatment <- factor(d$treatment, levels = unique(d$treatment[order(d$Nom_site, d$treatment)]))
  pd <- position_dodge2(width = 0.70, preserve = "single")

  p <- ggplot(d, aes(x = .data$abs_delta, y = .data$treatment, group = .data$framework_plot)) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_errorbar(
      aes(
        xmin = 0,
        xmax = pmax(0, .data$abs_delta),
        colour = .data$framework_plot,
        alpha = .data$is_best_any
      ),
      linewidth = 0.60,
      position = pd
    ) +
    geom_point(
      aes(
        colour = .data$framework_plot,
        size = .data$is_best_any,
        alpha = .data$is_best_any
      ),
      position = pd,
      stroke = 0.9
    ) +
    facet_grid(Nom_site ~ indicator_label, scales = "free_y", space = "free_y") +
    scale_colour_manual(values = delta_framework_colors, name = "Framework") +
    scale_size_manual(values = c(`FALSE` = 2.0, `TRUE` = 3.0), guide = "none") +
    scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 1.0), guide = "none") +
    coord_cartesian(xlim = c(0, 5)) +
    labs(x = expression("absolute " * Delta * " score (treatment - local control)"), y = NULL) +
    theme_bw(base_size = 9) +
    theme(
      strip.text = element_text(size = 9),
      axis.text.y = element_text(size = 6.5),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  save_plot_pair(p, stem, width = 11.0, height = 8.0)
}

plot_delta_abs(rank_submitted, "Figure_4_local_contrast_scores_submitted_legacy")

# Revised candidate: keep the sign of Delta scores and use the untrimmed reference.
plot_delta_signed <- function(d, stem) {
  d <- d %>%
    mutate(
      Nom_site = factor(.data$Nom_site, levels = sort(required_sites)),
      indicator_label = factor(.data$indicator_label, levels = c("Density", "Species richness")),
      treatment = factor(.data$treatment, levels = unique(.data$treatment[order(.data$Nom_site, .data$treatment)]))
    )
  pd <- position_dodge2(width = 0.70, preserve = "single")

  p <- ggplot(d, aes(x = .data$delta_score, y = .data$treatment, group = .data$framework_plot)) +
    geom_vline(xintercept = 0, linewidth = 0.35) +
    geom_errorbar(
      aes(
        xmin = pmin(0, .data$delta_score),
        xmax = pmax(0, .data$delta_score),
        colour = .data$framework_plot
      ),
      linewidth = 0.60,
      position = pd
    ) +
    geom_point(aes(colour = .data$framework_plot), position = pd, size = 2.3) +
    facet_grid(Nom_site ~ indicator_label, scales = "free_y", space = "free_y") +
    scale_colour_manual(values = delta_framework_colors, name = "Framework") +
    coord_cartesian(xlim = c(-5, 5)) +
    labs(x = expression(Delta * " score (treatment - local control)"), y = NULL) +
    theme_bw(base_size = 9) +
    theme(
      strip.text = element_text(size = 9),
      axis.text.y = element_text(size = 6.5),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  save_plot_pair(p, stem, width = 11.0, height = 8.0)
}

plot_delta_signed(delta_revised, "Figure_4_local_contrast_scores_revised_signed")

# Reproduction manifest --------------------------------------------------------
manifest <- c(
  paste0("Generated: ", format(Sys.time(), tz = "UTC"), " UTC"),
  "Central manuscript figure regeneration: COMPLETE",
  "",
  "Figure 1:",
  " - regenerated from coordinates on a static France outline;",
  " - original satellite/GIS basemap cannot be exactly regenerated because the GIS project/tile source was not supplied;",
  " - submitted raster retained under legacy/submitted_figures/Figure_1_submitted.png.",
  "",
  "Figure 2:",
  " - submitted-style two-panel RSr + density violin/box/jitter figure regenerated from agricultural reference data;",
  " - submitted significance letters retained as explicit visual metadata.",
  "",
  "Figure 3:",
  " - submitted legacy rendering reproduces the historical TIGA abundance < 75000 scoring exception;",
  " - revised candidate uses the full untrimmed TIGA abundance reference.",
  "",
  "Figure 4:",
  " - four-site analysis: QualiAgro, Thil, Yvetot, MetalEurope;",
  " - submitted legacy rendering uses absolute Delta score and historical TIGA density scoring;",
  " - revised candidate retains signed Delta score and uses the full untrimmed reference.",
  "",
  "Submitted manuscript rasters for Figures 1-4 are stored under legacy/submitted_figures/ for visual audit."
)
write_text_lines(manifest, file.path(MANUSCRIPT_OUTPUT_DIR, "FIGURE_REPRODUCTION_MANIFEST.txt"))

capture.output(sessionInfo(), file = file.path(MANUSCRIPT_OUTPUT_DIR, "sessionInfo_manuscript_figures.txt"))

message("[03] Manuscript figures written to: ", MANUSCRIPT_OUTPUT_DIR)
