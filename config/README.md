# Configuration

## `bioindicateur2_controls.csv`

The LRR analysis deliberately requires an **explicit experiment -> control mapping**.
No automatic control inference is used in the reviewer workflow.

Run the workflow once. It will create:

`output/reviewer/diagnostic_Bioindicateur2_treatment_levels.csv`

Use the exact `Nom_site`, `Traitement`, and (when relevant) `Contamination` labels from that diagnostic to populate one row per experiment to be included.

Columns:

- `Nom_site`: exact experiment/site label.
- `control_Site`: exact site/sample identifier for the local control.
- `control_Traitement`: exact control treatment label, or `*` if treatment is not the discriminating control field.
- `control_Contamination`: exact control contamination label, or `*` if not relevant.
- `include`: `TRUE` only after the mapping has been checked.
- `notes`: provenance / rationale for the control definition.

The manuscript validation uses explicit local references. This file is therefore part of the analytical provenance and should be committed once validated.


## Bioindicateur2 control mapping

`bioindicateur2_controls.csv` identifies the local control explicitly by `control_Site`.
The `control_Traitement` and `control_Contamination` fields are retained as cross-checks.
This is necessary for Yvetot because two observations share the label `blé`; only
`B10YVGC` is the confirmed local control used for the reviewer LRR analysis.

The four manuscript contrast sites are now explicit:

- QualiAgro: `B09TEM` (`TEM`, uncontaminated)
- Thil: `B10THLT` (`labour 30cm`, uncontaminated)
- Yvetot: `B10YVGC` (`blé`, uncontaminated)
- MetalEurope: `B09MTEC` (`T`, uncontaminated)

The submitted Methods accidentally described only three case-study sites; the submitted Figure 4 analysis and the recovered `scores.R` script use all four.
