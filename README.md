# Collembola reference-framework revision workflow

Reproducible R workflow rebuilt from the scripts and data supplied by H. Suarez for the *Ecological Indicators* revision on empirical Collembola reference frameworks.

The repository now separates six stages that had previously been mixed together:

1. **provenance audit** of the historical preprocessing;
2. **automatic reconstruction of the merged analysis dataset** from the five separate curated workbooks;
3. **reproduction of the original indicator/scoring core**;
4. **reviewer-requested robustness analyses**;
5. **regeneration of the central manuscript figures**;
6. **regeneration of the historical and reviewer-driven Supplementary Material**.

The historical `Stage_M2_Helio_Suarez.xlsx` is no longer a required input. It is treated only as an optional audit target.

## Quick start

Open `collembola-reference-frameworks.Rproj`, then run:

```r
source("analysis/run_all.R")
```

The workflow will:

- audit the supplied workbooks;
- rebuild the merged reference dataset from the five `Feuille stats` sheets;
- verify the expected reference sample sizes;
- recompute the historical indicators and empirical 0-6 scores;
- expose the historical TIGA abundance filter as a separate scenario;
- run the reviewer robustness analyses;
- export diagnostics for the Bioindicateur2 control mapping;
- regenerate the central manuscript figures;
- regenerate historical Supplementary Materials SM1-SM11 and reviewer-driven candidate supplements;
- save `sessionInfo()` and run parameters.

The first complete run can finish with the LRR analysis intentionally skipped. This is expected until `config/bioindicateur2_controls.csv` has been populated with explicit checked controls.

## Repository structure

```text
.
|-- analysis/
|   |-- 00_audit_input_reconstruction.R
|   |-- 01_build_analysis_dataset.R
|   |-- 02_reproduce_original_results.R
|   |-- reviewer_robustness_analysis.R
|   |-- 03_regenerate_manuscript_figures.R
|   |-- 04_regenerate_supplementary_material.R
|   `-- run_all.R
|-- R/
|   |-- project_paths.R
|   |-- io_helpers.R
|   |-- data_preparation.R
|   |-- indicator_scoring.R
|   `-- reviewer_functions.R
|-- config/
|   |-- bioindicateur2_controls.csv
|   `-- README.md
|-- data/
|   |-- raw/       # supplied workbooks; ignored by Git by default
|   |-- audit/     # optional historical merged file; audit only
|   |-- derived/   # generated RDS / manifests
|   `-- README.md
|-- output/
|   |-- audit/
|   |-- baseline/
|   |-- reviewer/
|   |-- manuscript/
|   `-- supplementary/
|-- legacy/
|   |-- Analyse_statistique_Helio_2025.R
|   |-- GiveScores_Helio_2025.R
|   |-- Mise_en_forme_Helio_2025.R
|   |-- multiCvM_Helio_2025.R
|   `-- statsvar_Helio_2025.R
|-- PROVENANCE_AUDIT.md
|-- REVIEWER_ANALYSIS_MAP.md
|-- CITATION.cff.template
|-- LICENSE
`-- collembola-reference-frameworks.Rproj
```

## Data inputs

The analysis-level workflow uses the current curated `Feuille stats` sheets from:

- `RMQS_2024_COLLEMBOLA.xlsx`
- `RMQS_2021_COLLEMBOLA.xlsx`
- `RMQS_Bretagne.xls`
- `ANDRA_ARTHROPODA.xlsx`
- `TIGA_rural_MESOFAUNA.xlsx`
- `Bioindicateur2_ARTHROPODA.xls`

and uses `TaxRef18_Collembola.csv` for the historical TAXREF quality-control audit.

See `data/README.md` and `PROVENANCE_AUDIT.md` for the distinction between curated analysis inputs and upstream raw-data reconstruction.

## Why not use `Stage_M2_Helio_Suarez.xlsx` directly?

The historical `Mise en forme.R` shows that `Stage_M2_Helio_Suarez.xlsx` is created after the separate datasets have been prepared and harmonised. Treating it as a primary input hides provenance and makes the analysis dependent on an opaque intermediate.

`analysis/01_build_analysis_dataset.R` therefore recreates that merge automatically. If the historical file is present under `data/audit/`, the workflow compares the reconstruction against it but never uses it to generate the revised results.

## Indicator definitions retained from the historical scripts

The baseline reproduction deliberately preserves the historical definitions so reviewer sensitivities are not confounded by a silent methodological change:

- taxon densities are rounded before community calculations;
- taxon columns whose names contain a space define the historical `communities_strict` subset;
- rarefied richness uses `vegan::rarefy(..., sample = 500)`;
- density (`ab`) is the row sum of all taxon columns;
- Shannon is computed on `communities_strict`;
- Pielou is `Shannon / log(RSr)` as in the supplied script;
- 0-6 scores reproduce the quantile logic of `GiveScores.R`.

The submitted manuscript focuses on **rarefied richness (`RSr`) and density (`ab`)**. Shannon and Pielou remain in the baseline audit for provenance but are not propagated through every reviewer sensitivity by default.

## Historical TIGA outlier rule

The supplied analysis script removed `TIGA_rural$ab >= 75000` immediately before density scoring. This was ad hoc and was not part of a general sensitivity framework.

The revised workflow therefore reports it explicitly rather than silently keeping or deleting it:

- `legacy_tiga_ab_lt_75000`;
- `full_untrimmed`;
- `winsorised_01_99`.

Inspect the resulting threshold and score differences before deciding which scenario becomes the revised main analysis.

## Reviewer analyses

The reviewer workflow addresses five questions:

1. **Sample-size imbalance** — 1,000 repeated subsamples of each agricultural framework to `n = 58` without replacement.
2. **Outlier sensitivity** — full untrimmed data versus 1-99% winsorisation, plus transparent >3 IQR counts and the historical TIGA filter comparison.
3. **Cross-framework concordance** — Spearman correlation, exact agreement, agreement within one class, and broad `0-2 / 3-4 / 5-6` agreement for identical Bioindicateur2 observations.
4. **Pooled framework** — raw site-weighted pooling versus balanced pooling with equal `n = 58` contribution from every framework.
5. **Raw ecological contrast versus score translation** — LRR on raw values alongside framework-dependent Delta-score, but only after explicit controls have been mapped.

See `REVIEWER_ANALYSIS_MAP.md` for the reviewer-to-output correspondence.

## Bioindicateur2 control mapping

No automatic control inference is used in the final workflow.

The first run writes:

`output/reviewer/diagnostic_Bioindicateur2_treatment_levels.csv`

Use this to populate:

`config/bioindicateur2_controls.csv`

with exact experiment/control labels. Rerun `analysis/reviewer_robustness_analysis.R` afterwards to obtain the LRR outputs.


## Manuscript Figures 1–4

Version 0.3 adds `analysis/03_regenerate_manuscript_figures.R`, which rebuilds the central manuscript figures into `output/manuscript/`. It creates submitted-style reproduction figures plus revised candidates where reviewer-driven methodological decisions differ from the original analysis. See `MANUSCRIPT_FIGURE_REPRODUCTION.md`.

The submitted Figure 4 is explicitly treated as a **four-site** Bioindicateur2 analysis (QualiAgro, Thil, Yvetot and MetalEurope); the three-site wording in the submitted Methods was an omission.

## Supplementary Material

Version 0.4 adds `analysis/04_regenerate_supplementary_material.R`. It rebuilds the submitted Supplementary Materials SM1-SM11 from the current reproducible workflow and stores reviewer-driven candidate supplementary items separately under `output/supplementary/reviewer_candidates/`. See `SUPPLEMENTARY_MATERIAL_REPRODUCTION.md`.

The revised supplement does **not** retain the historical ad hoc TIGA abundance `< 75000` filter. Shannon and Pielou are retained only to reproduce the historical supplementary scope; the new reviewer robustness analyses remain focused on RSr and density.

## Main output folders

### `output/audit/`

- workbook inventory;
- curated input dimensions;
- reconstructed merge versus historical merged workbook;
- upstream `temporary` -> curated-sheet audit where reconstruction is possible;
- TAXREF diagnostic.

### `output/baseline/`

- agricultural sample sizes;
- legacy indicator summaries;
- full untrimmed thresholds;
- historical TIGA-filter thresholds;
- score changes induced by that historical filter;
- rarefaction support diagnostic.

### `output/reviewer/`

- balanced `n = 58` threshold stability;
- outlier sensitivity;
- cross-framework concordance;
- pooled raw versus balanced results;
- LRR / Delta-score outputs once controls are mapped;
- figures, run metadata and `sessionInfo()`.

### `output/manuscript/`

- reconstructed Figures 1–4 in PDF and PNG;
- submitted-legacy and revised-candidate versions for Figures 3–4;
- plotting data, score counts and four-site ranking summaries;
- figure-reproduction manifest and `sessionInfo()`.

### `output/supplementary/`

- regenerated historical SM1-SM11 figures/tables;
- all-site and agricultural-site four-indicator distributions;
- descriptive statistics and Cramer-von Mises tests;
- seven-class scoring illustration and framework-specific scoring schemes;
- Bioindicateur2 four-indicator score plots plus Friedman/Nemenyi results;
- `reviewer_candidates/` containing equal-n, winsorisation, concordance, pooling and LRR/Delta-score candidates;
- regeneration manifest and `sessionInfo()`.

## R packages

The scripts check for required packages and stop with a clear message if one is absent. Core requirements are:

```r
c(
  "readxl", "dplyr", "tidyr", "purrr", "stringr", "tibble",
  "vegan", "ggplot2", "cowplot", "maps", "scales",
  "e1071", "cramer", "PMCMRplus"
)
```

No `setwd()` is used in the revised workflow, no object is written into `.GlobalEnv`, and no source workbook is modified in place.

## GitHub -> Zenodo

Before public release:

1. run from a clean checkout;
2. resolve any provenance-audit differences;
3. validate `config/bioindicateur2_controls.csv`;
4. decide and document the retained TIGA outlier treatment;
5. check redistribution rights for `data/raw/` and `data/audit/`;
6. commit the final code and selected reproducible outputs;
7. update and rename `CITATION.cff.template` to `CITATION.cff`;
8. create a versioned GitHub release;
9. archive that release in Zenodo and record the DOI.
