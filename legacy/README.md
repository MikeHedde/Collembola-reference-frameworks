# Legacy scripts supplied by H. Suarez

These files are preserved **unchanged** for provenance. They are byte-identical to the scripts in the supplied `Stats.zip` archive.

They are not sourced by the revised workflow because they contain local `setwd()` paths, in-place workbook modification, `.GlobalEnv` side effects, duplicated code, and historical/manual preprocessing assumptions.

- `Analyse_statistique_Helio_2025.R` — MD5 `51999e75f90ba441b9f2cc9918926ea5`
- `GiveScores_Helio_2025.R` — MD5 `784f9d0acc29e31c0297e3213ad7bf65`
- `Mise_en_forme_Helio_2025.R` — MD5 `dfc64a6a75b1c00ed6a39416c976f3a8`
- `multiCvM_Helio_2025.R` — MD5 `5d99d7ac94265eb01a03348de578ea52`
- `statsvar_Helio_2025.R` — MD5 `f32041b5792de66fbdcc0d8dd2a660a5`

All new code lives under `analysis/` and `R/`.

## Additional recovered provenance files

- `delta_score_comp_historical.csv` — recovered intermediate used by the later four-site `scores_four_sites_Hedde_2026.R`; retained for audit only and **not** used by the revised workflow.
- `submitted_supplementary_figures/` — raster copies extracted from the submitted manuscript for visual comparison only.
