# Changelog

## 0.4.0 — Supplementary Material regeneration

- Added `analysis/04_regenerate_supplementary_material.R`.
- Added `output/supplementary/` and a reviewer-candidate subfolder.
- Rebuilds historical Supplementary Materials SM1-SM11 from the current reproducible pipeline.
- Retains Shannon and Pielou only as historical supplementary indicators; reviewer robustness analyses remain focused on RSr and density.
- Recomputes pairwise Cramer-von Mises tests with explicit Bonferroni-adjusted p-values and regenerates compact-letter annotations.
- Recomputes Friedman/Nemenyi analyses for Bioindicateur2 scores.
- Regenerates the score-class method illustration and framework-specific scoring schemes.
- Does not retain the historical ad hoc TIGA abundance `< 75000` filter in revised supplementary results.
- Archives submitted supplementary rasters under `legacy/submitted_supplementary_figures/` for visual audit.
- Archives the recovered historical `delta_score_comp.csv` under `legacy/` as provenance only.
- Adds reviewer-driven candidate supplementary outputs for equal-n resampling, winsorisation, cross-framework concordance, pooling, and LRR versus Delta-score.

## 0.3.0 — manuscript figure reproduction

- Added `analysis/03_regenerate_manuscript_figures.R`.
- Added `output/manuscript/` as a dedicated output layer.
- Regenerates central manuscript Figures 1–4 from the cleaned workflow.
- Stores submitted manuscript figure rasters under `legacy/submitted_figures/` for visual audit only.
- Corrected the Bioindicateur2 contrast design to four sites by adding MetalEurope (`B09MTEC` control).
- Preserved the historical Figure 4 four-site script as `legacy/scores_four_sites_Hedde_2026.R`.
- Corrected the legacy TIGA filter reproduction so `ab < 75000` affects density scoring only, matching the historical script sequence.
- Produces both historical submitted-style and revised untrimmed/signed candidates for Figures 3–4.
- Documents the limitation that the exact Figure 1 satellite basemap cannot be regenerated from the supplied analytical files.


## v0.2.1

- Fixed a dplyr data-mask bug in `get_breaks_from_table()` that mixed thresholds from all frameworks/indicators and caused `findInterval()` to fail with unsorted breaks.
- Added explicit validation that recovered internal score breaks are finite and sorted.
- Added automatic `.xlsx`/`.xls` resolution for all Excel inputs.
- Removed tidyselect `.data` deprecation warnings in the affected `select()`/`pull()` calls.

## v0.2.0-audit — 2026-08-19

- removed `Stage_M2_Helio_Suarez.xlsx` as a required analysis input;
- rebuilds the merged reference dataset from five separate curated workbooks;
- added provenance audit for historical `Mise en forme.R` transformations;
- added input checksums and historical merge comparison;
- separated audit, baseline reproduction and reviewer outputs;
- preserved legacy RSr/Shannon/Pielou/density definitions while restricting reviewer sensitivities to manuscript metrics RSr and density;
- exposed the historical ad hoc `TIGA_rural$ab < 75000` density filter as an explicit scenario;
- replaced automatic Bioindicateur2 control inference with an explicit mapping requirement;
- added one-command `analysis/run_all.R` workflow;
- bundled supplied data in the local test archive while keeping them ignored by Git pending redistribution checks.

## v0.2.2

- Added explicit `control_Site` mapping for Bioindicateur2 LRR analyses.
- Confirmed controls: QualiAgro `B09TEM`, Thil `B10THLT`, Yvetot `B10YVGC`.
- Prevented the second Yvetot row labelled `blé` (`B10YVSIII`) from being silently treated as a control.
- Replaced a few deprecated `.data` tidyselect uses in `pivot_wider()`.
