# Provenance audit: historical preprocessing versus revision workflow

## Why the workflow was changed

The first repo draft treated `Stage_M2_Helio_Suarez.xlsx` as a required input. Inspection of `legacy/Mise_en_forme_Helio_2025.R` showed that this Excel file is actually a **derived merged product** created at the end of H. Suarez's preprocessing workflow.

The historical script performs more than concatenation:

1. site-level averaging and conversion to individuals m-2 using a 0.002827 m2 core area;
2. dataset-specific site-name cleanup and grouping rules;
3. coordinate joins and CRS conversions;
4. a TAXREF name check;
5. final harmonisation and concatenation of the five reference datasets.

It also explicitly states that some original datasets were manually modified in Excel before the scripted steps.

## Two reproducibility layers

The revised repo therefore separates two questions.

### Analysis-level reproducibility

**Supported from the current curated `Feuille stats` sheets.**

The five reference sheets are read independently, harmonised using the historical long -> bind -> wide logic, and merged automatically. `Stage_M2_Helio_Suarez.xlsx` is not required.

### Upstream raw-data reconstruction

**Partial and dataset-specific.**

`analysis/00_audit_input_reconstruction.R` attempts only transformations explicitly preserved in the historical script. It does not guess taxonomic corrections or undocumented Excel edits.

Known limitations before running the audit:

- RMQS 2024: the historical script starts from a `Feuille stats` sheet and later overwrites that sheet, so the exact pre-script version is no longer recoverable from the workbook/script pair alone.
- RMQS 2021: the historical script contains coordinate transformation but no preserved block reconstructing biological densities from raw counts.
- ANDRA, TIGA rural, RMQS Bretagne, and Bioindicateur2: `temporary` sheets exist and can be compared with the current curated sheets for the explicitly scripted averaging/scaling operations.

## Historical analysis issue made explicit

`legacy/Analyse_statistique_Helio_2025.R` contains an ad hoc removal:

`gamme_TIGA_rural <- gamme_TIGA_rural[gamme_TIGA_rural$ab < 75000, ]`

immediately before density scoring. The revised workflow does not hide this. It produces three explicit scenarios where relevant:

- historical TIGA density filter `< 75000`;
- full untrimmed reference;
- systematic 1-99% winsorisation sensitivity.

The manuscript/rebuttal should state which scenario is retained after inspecting the outputs.

## Bioindicateur2 controls

LRR values are not computed from automatic control heuristics. The revision workflow requires a checked `experiment -> control` mapping in `config/bioindicateur2_controls.csv`. The first run exports all observed labels to help construct that mapping.
