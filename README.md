# Collembola reference frameworks

Reproducible analysis workflow associated with the manuscript:

> **Empirical reference frameworks shape Collembola indicator distributions, scores and contrast sensitivity**

Suarez H., Cortet J., Auclerc A., Bougon N., Brand M., Henon N., Jolivet C.,
Lévêque A., Pouzenc S., Versavel C. & Hedde M.

Target journal: *Ecological Indicators*.

## Overview

Ecological indicators are often interpreted relative to empirical reference
populations. When different monitoring datasets are used to construct these
references, differences in their realised geographic, environmental,
land-use and management composition can affect indicator distributions,
quantile-derived thresholds and the scores subsequently assigned to the same
observations.

This repository reproduces the analyses used to compare four empirical
reference frameworks for Collembola indicators:

- **RMQS-Biodiversité**
- **RMQS BioDiv Bretagne**
- **ANDRA**
- **TIGA rural**

The same independent set of 26 agricultural Bioindicateur2 observations is
then scored against each reference framework.

The principal manuscript analyses focus on:

- rarefied species richness (`RSr`);
- Collembola density (`ab`, individuals m^-2).

Shannon diversity and Pielou evenness are retained in the supplementary
workflow to document the broader historical analysis.

## Main analytical questions

The workflow evaluates three related questions.

1. How do empirical Collembola indicator distributions and quantile-derived
   score boundaries differ among candidate reference frameworks?

2. How do these framework-specific boundaries affect the scores assigned to
   the same independent observations, including their relative ranking and
   classification agreement?

3. How are raw ecological treatment-control contrasts translated through
   framework-specific ordinal score boundaries?

Sensitivity analyses additionally examine:

- unequal reference sample sizes;
- extreme observations and 1–99% winsorisation;
- cross-framework score concordance;
- site-weighted versus equal-weight pooled reference populations.

These analyses are intended as robustness and transferability assessments.
They do not isolate a causal effect of sampling design from geographic extent,
environmental composition, land use or management history.

## Repository structure

```text
.
├── analysis/
│   ├── 01_validate_analysis_data.R
│   ├── 02_main_analysis.R
│   ├── 03_sensitivity_analyses.R
│   ├── 04_manuscript_figures.R
│   ├── 05_supplementary_material.R
│   └── run_all.R
│
├── R/
│   ├── analysis_functions.R
│   ├── indicator_scoring.R
│   ├── io_helpers.R
│   └── project_paths.R
│
├── data/
│   ├── README.md
│   └── analysis/
│       ├── reference_indicator_values.csv
│       ├── bioindicateur2_indicator_values.csv
│       ├── bioindicateur2_contrasts.csv
│       ├── bioindicateur2_controls.csv
│       └── map_points.csv
│
├── output/
│   ├── analysis/
│   ├── sensitivity/
│   ├── manuscript/
│   └── supplementary/
│
├── reproducibility/
│   └── sessionInfo.txt
│
├── CITATION.cff
├── LICENSE
├── VERSION
└── collembola-reference-frameworks.Rproj
```
## Analysis-ready data

The repository starts from frozen site-level indicator values rather than from
the original taxon-by-site abundance matrices.

The four candidate reference datasets contain:

Framework	All observations	Agricultural observations
RMQS-Biodiversité	96	58
RMQS BioDiv Bretagne	98	89
ANDRA	132	90
TIGA rural	430	345
Total	756	582

Agricultural sites are defined from Corine Land Cover level 1.

Three RMQS-Biodiversité observations do not have an assigned CLC level-1
category and are explicitly retained as missing in the supplementary
land-cover summary.

For ANDRA, the 90 agricultural observations correspond to 89 unique mapped
locations because site O17 was sampled in two different years.

See data/README.md for detailed variable definitions and
data-provenance information.

## Scoring system

For each empirical reference framework and indicator, the observed:

minimum (q00);
q20;
q40;
q60;
q80;
maximum (q100);

define seven ordered relative score classes from 0 to 6.

Scores 1–5 partition observations within the empirical reference range.
Score 0 represents values below the observed minimum, whereas score 6
represents values at or above the observed maximum.

The resulting score therefore represents an observation's position relative to
a selected empirical reference population. It should not be interpreted as an
independent measure of complete soil health or as a biological effect-size
scale.

## Reproducing the analyses

From the repository root:

Rscript analysis/run_all.R

run_all.R executes each step in a separate R session so that scripts cannot
depend on objects left in memory by previous steps.

The workflow performs:

analysis-ready data
        ↓
01_validate_analysis_data.R
        ↓
02_main_analysis.R
        ↓
03_sensitivity_analyses.R
        ↓
04_manuscript_figures.R
        ↓
05_supplementary_material.R

A successful complete run ends with:

Complete workflow finished successfully.

## Analysis steps
01 — Data validation

01_validate_analysis_data.R checks:

presence and structure of all analysis-ready files;
unique observation keys;
expected dataset totals;
agricultural reference sample sizes;
Bioindicateur2 sample size;
number of treatment-control contrasts;
indicator validity.

Expected totals are 756 reference observations, 26 Bioindicateur2 observations
and 28 raw treatment-control contrasts.

02 — Main analysis

02_main_analysis.R:

selects agricultural reference populations;
calculates descriptive summaries;
derives framework-specific empirical score boundaries;
scores the same 26 Bioindicateur2 observations against all four frameworks;
records rarefaction-support metadata.
03 — Sensitivity analyses

03_sensitivity_analyses.R implements:

1,000 repeated equal-size subsamples at n = 58 per framework;
sensitivity to extreme observations and 1–99% winsorisation;
pairwise Spearman rank concordance;
exact score agreement;
agreement within ±1 score class;
agreement after grouping scores into broad classes 0–2, 3–4 and 5–6;
site-weighted and equal-weight pooled references;
translation of raw treatment-control log response ratios into signed
framework-specific score differences.

For directional treatment-control comparisons, absolute LRR values below
1e-8 are treated as numerical zero solely to avoid assigning direction to
floating-point noise.

04 — Main manuscript figures

The script generates:

Figure 1 — spatial distribution of datasets;
Figure 2 — agricultural reference distributions;
Figure 3 — scores assigned to Bioindicateur2 observations;
Figure 4 — signed treatment-control score contrasts.

Both PDF and PNG versions are produced.

05 — Supplementary material

The script generates the final supplementary material:

Figures S1–S7;
Tables S1–S11.

Internal consistency checks verify several key published results before the
script terminates.

## Key reproducibility checks

A successful workflow reproduces, among others, the following results:

equal-size sensitivity: 95% resampling intervals exclude zero for
17/24 pairwise internal richness-threshold differences and 15/24
density-threshold differences;
pooled raw versus equal-weight references give exact Bioindicateur2 score
agreement of 96.15% for richness and 88.46% for density;
all pooled-reference scores remain within one class;
26 directional raw contrasts evaluated with four reference frameworks yield
104 framework-by-contrast translations;
among these, 75 (72.1%) yield a non-zero score difference in the same
direction as the raw response, 29 (27.9%) are compressed to
Delta score = 0, and none reverse direction.
Statistical analyses

Distributional differences among empirical reference frameworks are evaluated
using pairwise two-sample Cramér–von Mises tests with 10,000 ordinary
resampling replicates and Bonferroni adjustment.

Differences among scores assigned to the same Bioindicateur2 observations are
evaluated using Friedman tests followed by Nemenyi pairwise comparisons.

Cross-framework Spearman correlations and classification-agreement statistics
are interpreted descriptively.

Treatment-control comparisons are also descriptive because the available
analysis-ready Bioindicateur2 dataset contains one aggregated value per
treatment or control condition rather than replicated treatment-level
observations.

## Software dependencies

The analysis uses R and the following packages:

dplyr
tidyr
tibble
cramer
PMCMRplus
e1071
ggplot2
cowplot
maps
scales

Exact package and R versions used in a run are recorded through sessionInfo()
in the reproducibility and output directories.

## Data provenance and redistribution

The original source workbooks and taxon-by-site abundance matrices are not
redistributed in this repository because dataset-specific redistribution
rights for those files have not been established.

The repository instead provides the frozen analysis-ready site-level indicator
values underlying the statistical analyses presented in the manuscript.

Accordingly, the public workflow reproduces the analyses from the
analysis-ready indicator values onward. It does not independently reconstruct
site-level indicator values from the original taxonomic source files.

The original calculation of rarefied richness used rarefaction implemented in
R with the vegan package; those calculations precede the public workflow.

## Outputs

Generated files are written to:

output/analysis/
output/sensitivity/
output/manuscript/
output/supplementary/

The output directories are regenerated by the workflow and should therefore be
treated as derived products rather than source data.

## Citation

If you use this repository, please cite the associated manuscript:

Suarez H., Cortet J., Auclerc A., Bougon N., Brand M., Henon N., Jolivet C.,
Lévêque A., Pouzenc S., Versavel C. & Hedde M.
Empirical reference frameworks shape Collembola indicator distributions,
scores and contrast sensitivity.
Ecological Indicators.

A version-specific archival citation will be added after deposition of the
public release in Zenodo.

## License

Code in this repository is distributed under the terms of the MIT License.

Reuse of the analysis-ready data remains subject to the provenance and rights
associated with the original contributing datasets.
