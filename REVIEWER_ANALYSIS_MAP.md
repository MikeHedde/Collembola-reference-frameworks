# Reviewer-analysis map

| Reviewer concern | Analysis / response | Key outputs | Intended manuscript use |
|---|---|---|---|
| Frameworks have very different sample sizes | Repeated subsampling without replacement of each agricultural framework to `n = 58`, 1,000 repeats | `Table_02_balanced_n58_threshold_stability.csv`, `Table_S1b_balanced_n58_pairwise_threshold_differences.csv`, `Figure_S1_balanced_n58_threshold_stability.pdf` | Main robustness result + detailed Supplement |
| Extreme values may affect quantile score boundaries | Compare full untrimmed thresholds with systematic 1-99% winsorisation; report >3 IQR observations; separately expose historical TIGA `ab < 75000` rule | `Table_S0_historical_TIGA_filter_sensitivity.csv`, `Table_S3_outlier_threshold_sensitivity.csv`, `Table_S3b_outlier_score_sensitivity.csv`, `Table_S4_extreme_observation_diagnostics.csv` | Supplement + explicit Methods clarification |
| Some broad classifications may remain consistent across frameworks | Score the same Bioindicateur2 observations under all four frameworks; quantify rank, exact, +/-1 and broad-class agreement | `Table_04_cross_framework_score_agreement.csv`, `Figure_01_cross_framework_score_agreement.pdf` | Main compact result; full table Supplement |
| A pooled dataset may be informative but can be dominated by large datasets | Compare raw site-weighted pooling with equal-framework pooling based on repeated `n = 58` contribution | `Table_05_pooled_reference_thresholds.csv`, `Table_07_pooled_raw_vs_balanced_agreement.csv`, `Table_S5_pooled_balanced_all_repeats.csv` | Supplement / short Discussion result |
| `|Delta score|` may reflect threshold placement rather than ecological effect magnitude | Report raw-value `ln(treatment/control)` LRR alongside framework-dependent score contrast; require explicit control mapping | `diagnostic_Bioindicateur2_treatment_levels.csv`, then `Table_08_raw_LRR_and_delta_score.csv`, `Figure_02_raw_LRR_vs_delta_score.pdf` | Main validation result after control mapping is checked |
| Sampling design, spatial extent and represented population are confounded | Change interpretation from causal effect of `sampling strategy alone` to dependence on the empirical reference / monitoring framework | wording, not a separate test | Throughout manuscript |
| Soil-health threshold language is too broad | Restrict claims to empirical quantile-derived biodiversity thresholds and framework-dependent scores | wording, not a separate test | Title, Abstract, Introduction, Discussion, Conclusion |
| Historical analysis provenance is unclear | Reconstruct merge from separate curated workbooks; audit upstream scripted preprocessing; do not use historical merged Excel as an input | `output/audit/*`, `PROVENANCE_AUDIT.md` | Reproducibility statement / GitHub-Zenodo documentation |

## Terminology preferred in the revision

Prefer:

- `reference framework`
- `monitoring framework`
- `reference population`
- `empirical quantile-derived thresholds`
- `framework-dependent scores`
- `score contrast sensitivity`
- `framework-dependent score separation`

Avoid or qualify:

- `sampling strategy alone`
- `primary driver`
- `first-order determinant`
- `universal soil-health threshold`
- `reference condition` when only an empirical reference distribution is meant
- `discriminant power` when the metric is separation in ordinal score space
