# First run: what to inspect

Run from the opened RStudio project:

```r
source("analysis/run_all.R")
```

Then inspect these files first:

1. `output/audit/audit_merged_verdict.csv`
   - desired result: `PASS`;
   - this verifies that the separate curated reference workbooks recreate the historical merged table without using it as input.

2. `output/audit/audit_upstream_preparation.csv`
   - this is expected to show partial reconstruction and/or differences where undocumented Excel curation occurred;
   - differences here do **not** invalidate the reviewer analyses, but must be documented for provenance.

3. `output/baseline/diagnostic_legacy_tiga_ab_filter.csv`
   - tells us exactly how many TIGA agricultural observations were affected by the historical `ab < 75000` rule.

4. `output/baseline/baseline_legacy_vs_untrimmed_thresholds.csv`
   - quantifies whether that historical rule materially changes density score boundaries.

5. `output/reviewer/diagnostic_Bioindicateur2_treatment_levels.csv`
   - use the exact labels here to fill `config/bioindicateur2_controls.csv` before rerunning the LRR section.


6. `output/manuscript/FIGURE_REPRODUCTION_MANIFEST.txt`
   - confirms regeneration of Figures 1–4;
   - compare the generated PNGs with `legacy/submitted_figures/`;
   - Figure 1 is expected to differ in basemap appearance because the original GIS/satellite source was not supplied.

7. `output/manuscript/Figure_04_control_diagnostics.csv`
   - all four sites should show `used_explicit_mapping`;
   - expected controls: B09TEM, B10THLT, B10YVGC, B09MTEC.

If the workflow stops, send the complete R error plus the name of the script/step shown immediately above it.

8. `output/supplementary/SUPPLEMENTARY_REPRODUCTION_MANIFEST.txt`
   - confirms regeneration of historical SM1-SM11;
   - compare SM02, SM05, SM08, SM09 and SM10 PNGs with `legacy/submitted_supplementary_figures/`.

9. `output/supplementary/SM06_descriptive_statistics_agricultural_sites.csv`
   - RSr and density should reproduce the values already validated in the central Figure 2 workflow;
   - Shannon/Pielou are retained as historical supplementary indicators.

10. `output/supplementary/SM11_Friedman_tests.csv` and `SM11_Nemenyi_pairwise_tests.csv`
    - these replace the historical test table with results recalculated from the current scoring pipeline.

11. `output/supplementary/reviewer_candidates/`
    - candidate material only; final numbering and whether each item remains supplementary or moves to the main text should be decided during manuscript drafting.
