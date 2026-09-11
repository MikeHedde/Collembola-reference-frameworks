# Analysis-ready data

This directory contains the analysis-ready datasets required to reproduce the
analyses, figures and supplementary material associated with:

**Empirical reference frameworks shape Collembola indicator distributions,
scores and contrast sensitivity**

The public workflow starts from site-level indicator values rather than from
the original taxon-by-site abundance matrices.

## Directory structure

```text
data/
└── analysis/
    ├── README.md
    ├── reference_indicator_values.csv
    ├── bioindicateur2_indicator_values.csv
    ├── bioindicateur2_contrasts.csv
    ├── bioindicateur2_controls.csv
    └── map_points.csv
Files
reference_indicator_values.csv

Site-level indicator values for the four candidate empirical reference
frameworks:

RMQS-Biodiversité
RMQS BioDiv Bretagne
ANDRA
TIGA rural

The file contains 756 observations in total.

Variables include dataset identity, sampling metadata, coordinates, Corine Land
Cover information and the site-level Collembola metrics used in the analyses:

RSr: rarefied species richness;
Shannon: Shannon diversity;
Pielou: Pielou evenness;
ab: Collembola density (individuals m^-2).

The main manuscript analyses use RSr and ab. Shannon diversity and Pielou
evenness are retained to reproduce the broader supplementary analyses.

The analytical agricultural reference populations contain:

Framework	Agricultural observations
RMQS-Biodiversité	58
RMQS BioDiv Bretagne	89
ANDRA	90
TIGA rural	345

Note that site_id is not globally unique. Some RMQS identifiers occur in both
the national and Brittany datasets, and one ANDRA site (O17) was sampled in
two different years. The combination framework + site_id + year uniquely
identifies reference observations.

bioindicateur2_indicator_values.csv

Site/treatment-level indicator values for the 26 agricultural Bioindicateur2
observations used as an independent evaluation dataset.

These observations are scored against each of the four empirical reference
frameworks.

bioindicateur2_contrasts.csv

The 28 raw treatment-control contrasts used for the four Bioindicateur2 case
studies:

QualiAgro
Thil
Yvetot
MetalEurope

The file contains treatment and control values and their log response ratios
(LRR). One richness LRR is undefined because the treatment value is zero; one
additional richness contrast is treated as numerical zero using the tolerance
defined in analysis/03_sensitivity_analyses.R.

bioindicateur2_controls.csv

Explicit control definitions for the four Bioindicateur2 case studies.

This file documents the treatment-control mapping used when the aggregated
contrast dataset was constructed and is retained for provenance and validation.

map_points.csv

Coordinates used to generate Figure 1.

For ANDRA, 90 agricultural observations correspond to 89 unique mapped
locations because site O17 was sampled in two different years.

Data provenance

The analysis-ready files were derived from the original Collembola datasets
used in the study. They contain aggregated indicator values and analytical
metadata, not the original taxon-by-site abundance matrices.

The original source workbooks are not distributed in this repository because
dataset-specific redistribution rights for those source files have not been
established.

Consequently, this repository reproduces the statistical analyses, sensitivity
analyses, scoring, figures and supplementary material from the frozen
analysis-ready indicator values onward. It does not reproduce the initial
taxonomic-data processing and calculation of site-level indicators from the
original source workbooks.

Reproducibility

From the repository root, run:

Rscript analysis/run_all.R

The workflow validates the supplied datasets before running the analyses.

The expected reference dataset totals are:

Framework	Total observations
RMQS-Biodiversité	96
RMQS BioDiv Bretagne	98
ANDRA	132
TIGA rural	430

The complete reference dataset therefore contains 756 observations.
