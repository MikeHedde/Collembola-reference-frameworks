# Reproduction of manuscript Figures 1–4

Version 0.3 adds an explicit manuscript-figure layer to the reproducibility workflow.

Run:

```r
source("analysis/03_regenerate_manuscript_figures.R")
```

or run the full workflow:

```r
source("analysis/run_all.R")
```

Outputs are written to `output/manuscript/`.

## Figure 1 — sampling map

The submitted map used a satellite/GIS basemap. No GIS project, cached basemap or tile-source specification was present in the files supplied by the original analyst. The v0.3 script therefore regenerates the sampling locations without an online dependency, using a static France outline from the R `maps` package. The submitted raster is retained in `legacy/submitted_figures/` for comparison.

For the submitted-figure audit, the plotting population follows the numbers/selection represented in the submitted map: all RMQS-Biodiversity and TIGA locations, all RMQS BioDiv Bretagne locations, agricultural ANDRA locations, and agricultural Bioindicateur2 locations. The exact mapped counts are exported as `Figure_01_map_counts.csv`.

## Figure 2 — reference distributions

The two-panel rarefied-richness/density violin + boxplot + jitter figure is regenerated from the agricultural reference data. The colour palette and the letters shown in the submitted figure are retained explicitly for visual reproduction. Those letters are treated as submitted-figure metadata rather than silently recomputed under a potentially different package version.

## Figure 3 — Bioindicateur2 scores

Two versions are produced:

- `Figure_3_Bioindicateur2_scores_submitted_legacy.*` reproduces the historical scoring sequence, including the ad hoc `TIGA_rural$ab < 75000` filter applied **only for density scoring**.
- `Figure_3_Bioindicateur2_scores_revised_untrimmed.*` uses the reviewer-audited full, untrimmed reference and is the candidate for the revised manuscript.

## Figure 4 — local contrasts

The submitted analysis used **four sites**: QualiAgro, Thil, Yvetot and MetalEurope. The three-site wording in the submitted Methods was an omission. Explicit controls are read from `config/bioindicateur2_controls.csv`.

Two versions are produced:

- `Figure_4_local_contrast_scores_submitted_legacy.*` reconstructs the historical plot using absolute Delta scores and best/co-best emphasis.
- `Figure_4_local_contrast_scores_revised_signed.*` retains the sign of Delta score and uses untrimmed reference distributions, matching the revised interpretation of framework-dependent score translation rather than ecological “discriminant power”.

The ranking summaries and plotting data are also exported as CSV files so manuscript percentages can be recomputed transparently.
