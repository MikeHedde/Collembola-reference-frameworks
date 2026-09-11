# Supplementary Material regeneration

Version 0.4 adds a complete regeneration layer for the Supplementary Materials that accompanied the submitted *Ecological Indicators* manuscript.

## Principle

The **current reproducible pipeline is the source of truth**. Historical figures/tables are kept only for provenance and visual comparison. In particular, the revised Supplementary Material does not retain the historical ad hoc `TIGA_rural$ab < 75000` density filter; its effect remains documented separately under `output/baseline/` and `output/reviewer/`.

## Historical Supplementary Materials rebuilt

`analysis/04_regenerate_supplementary_material.R` regenerates:

1. **SM1** — site/observation counts by land-cover category and framework;
2. **SM2** — distributions of RSr, Shannon, Pielou and density for all reference sites;
3. **SM3** — descriptive statistics for all sites;
4. **SM4** — pairwise Cramer-von Mises tests for all sites, with Bonferroni-adjusted p-values;
5. **SM5** — distributions of the four indicators for agricultural sites;
6. **SM6** — descriptive statistics for agricultural sites;
7. **SM7** — pairwise Cramer-von Mises tests for agricultural sites;
8. **SM8** — illustration of the seven empirical score classes;
9. **SM9** — framework-specific scoring schemes for the four historical indicators;
10. **SM10** — Bioindicateur2 score distributions for the four historical indicators;
11. **SM11** — Friedman and Nemenyi tests for framework-dependent Bioindicateur2 scores.

The submitted manuscript's supplementary raster figures are archived under `legacy/submitted_supplementary_figures/` for visual audit only.

## Reviewer-driven candidate Supplementary Materials

The script also gathers/creates candidate supplementary items for:

- balanced repeated subsampling to `n = 58`;
- sensitivity of q20-q80 to 1-99% winsorisation;
- cross-framework score concordance;
- pooled raw versus pooled-balanced thresholds;
- raw LRR versus framework-dependent Delta-score.

These are written to `output/supplementary/reviewer_candidates/`. Their final numbering should be decided only after the revised main figures and Results text are frozen.

## Indicator scope

The historical Supplementary Material is rebuilt for all four indicators used in the original analytical scripts (`RSr`, Shannon, Pielou, density), because those plots/tables were part of the submitted supplement.

The **reviewer robustness analyses remain restricted to the two manuscript indicators**, rarefied species richness and density, to avoid silently expanding the scope of the revision.

## Run

From the project root:

```r
source("analysis/04_regenerate_supplementary_material.R")
```

or regenerate the complete project with:

```r
source("analysis/run_all.R")
```

Additional packages needed for the historical statistical tables are `e1071`, `cramer`, and `PMCMRplus`.
