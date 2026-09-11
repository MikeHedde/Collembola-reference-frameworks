# Data layout and provenance

The revision workflow no longer requires the pre-merged `Stage_M2_Helio_Suarez.xlsx` as an input.

## `raw/` — supplied workbooks used by the revision workflow

For local reproduction, place the following files in `data/raw/`:

- `RMQS_2024_COLLEMBOLA.xlsx`
- `RMQS_2021_COLLEMBOLA.xlsx`
- `RMQS_Bretagne.xls`
- `ANDRA_ARTHROPODA.xlsx`
- `TIGA_rural_MESOFAUNA.xlsx`
- `Bioindicateur2_ARTHROPODA.xls`
- `TaxRef18_Collembola.csv`

The primary revision pipeline starts from the current curated `Feuille stats` sheets in these workbooks. This is intentional: H. Suarez's historical `Mise en forme.R` states that some upstream restructuring/correction was performed manually in Excel and the preserved script does not fully reconstruct every curated sheet from raw field/lab exports.

`analysis/00_audit_input_reconstruction.R` separately tests which upstream transformations can be reconstructed from preserved `temporary` sheets without guessing undocumented edits.

## `audit/` — historical derived files, optional

- `Stage_M2_Helio_Suarez.xlsx`

This file is **never used to build the revised analysis dataset**. If present, it is used only as an audit target to verify that the five separate curated reference workbooks recreate the historical merged table.

## `derived/` — generated automatically

`analysis/01_build_analysis_dataset.R` creates:

- `reference_frameworks_merged.rds`
- `bioindicateur2_curated.rds`
- `reference_counts.csv`
- `build_manifest.txt`
- optional `historical_merge_comparison.csv`

These are reproducible intermediate products and should not be edited manually.

## Public GitHub / Zenodo release

The supplied data are bundled in the current working/test archive because they were provided for this analysis. Before a public GitHub/Zenodo release, verify redistribution rights for every workbook. `.gitignore` excludes `data/raw/*` and `data/audit/*` by default to prevent accidental publication.
