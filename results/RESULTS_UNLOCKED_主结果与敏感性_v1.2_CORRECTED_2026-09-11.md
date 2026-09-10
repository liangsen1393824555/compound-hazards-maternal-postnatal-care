# RESULTS UNLOCKED：主结果与敏感性 v1.2（更正版，2026-09-11）

This document supersedes the public wording in v1.1 while preserving the v1.1 files and every reported numeric result. See the accompanying erratum for the correction log.

## Mechanical conclusion

**CONFIRMATORY_NULL_OR_INCONCLUSIVE.** Neither pre-specified Nigeria confirmatory estimand passed the Holm familywise threshold of 0.05. The adjusted `RD_11-00` was negative (lower standardized risk in cell 11 than cell 00), but its Holm-adjusted P value was 0.0584 and its 97.5% simultaneous confidence interval included zero. `IC_add` was also negative, imprecise, and compatible with zero. These are associations, not causal effects.

## Integrity check before unlock

- Final sealed inventory state before the original unlock: `READY_FOR_RESULTS_UNLOCK_RESULTS_STILL_SEALED`; its unlock flag was `false`.
- Original runs A/B: all five sealed files matched byte-for-byte, including `sealed_results.json` SHA-256 `9070F0511B88F4FAC2D4E97FA492AC86CD1AB1E48288A9A5459B345632EDE59D`.
- Scheme B runs A/B: all five sealed files matched byte-for-byte, including `scheme_b_sealed_results.json` SHA-256 `56591B4D92DF8FD3BD57C8CEF4DBF33EBEC4F49CBACCA04238A9E6FE11F3BB32`.
- Across the inventory and both run pairs, all 28 checked artifact instances matched their paired/registered SHA-256 values; numeric maximum absolute A/B difference for Scheme B was 0.
- No model was refitted, changed, or selected during unlock or correction. This report only reads and summarizes sealed outputs.

## Nigeria primary analysis

Adjusted survey-weighted logistic model, complete-case `n = 10,643`.

| Joint exposure cell | Adjusted standardized risk | 95% CI |
|---|---:|---:|
| 00 | 58.67% | 57.07% to 60.27% |
| 10 | 60.90% | 55.99% to 65.81% |
| 01 | 53.21% | 50.97% to 55.45% |
| 11 | 52.64% | 47.50% to 57.78% |

Confirmatory estimands:

| Estimand | Estimate | 95% CI | Raw P | Holm P | 97.5% simultaneous CI | Holm pass? |
|---|---:|---:|---:|---:|---:|---:|
| `RD_11-00` | -6.03 percentage points | -11.45 to -0.61 pp | 0.0292 | 0.0584 | -12.23 to 0.17 pp | No |
| `IC_add` | -2.80 percentage points | -10.04 to 4.44 pp | 0.4484 | 0.4484 | -11.08 to 5.48 pp | No |

The ordinary 95% interval and raw P value for `RD_11-00` exclude the null, but the pre-specified multiplicity-adjusted decision does not. The simultaneous interval also includes zero. Therefore this is not a confirmatory-positive result.

Frozen secondary estimands (descriptive/supportive, not substitutes for the confirmatory tests):

| Estimand | Estimate | 95% CI |
|---|---:|---:|
| `RD_10-00` | 2.23 pp | -2.92 to 7.38 pp |
| `RD_01-00` | -5.46 pp | -8.26 to -2.66 pp |
| `RR_10/00` | 1.038 | 0.954 to 1.130 |
| `RR_01/00` | 0.907 | 0.862 to 0.954 |
| `RR_11/00` | 0.897 | 0.810 to 0.994 |
| Risk-scale ratio interaction | 0.953 | 0.838 to 1.084 |

## Sensitivity analyses, in the pre-specified execution order

These analyses are supportive. None replaces or redefines the primary estimands.

| Order | Sensitivity | `RD_11-00` (95% CI), pp | Direction/precision versus primary | `IC_add` (95% CI), pp | Interaction consistency |
|---:|---|---:|---|---:|---|
| 1 | 25 km radius | -7.02 (-13.41, -0.63) | Same direction; slightly larger magnitude; ordinary 95% CI excludes 0, but Holm P=0.0628 and simultaneous CI includes 0 | -1.46 (-9.29, 6.36) | Same negative direction; imprecise |
| 2 | 100 km radius | -4.79 (-9.61, 0.03) | Same direction; smaller magnitude; CI narrowly includes 0 | -0.71 (-9.60, 8.18) | Same negative direction; imprecise |
| 3 | Scheme B | -4.59 (-8.42, -0.77) | Same direction; smaller magnitude; supportive adjusted estimate, P=0.0185 unadjusted for multiplicity | 1.79 (-2.40, 5.97) | Direction reverses; CI includes 0; P=0.4025 unadjusted for multiplicity |
| 4 | Livebirth-only | -5.95 (-11.38, -0.51) | Same direction and magnitude; ordinary 95% CI excludes 0, but Holm P=0.0639 and simultaneous CI includes 0 | -2.84 (-10.11, 4.44) | Same direction; imprecise |
| 5 | Strict fixed effects | -3.00 (-8.21, 2.21) | Same direction; attenuated; CI includes 0 | 1.51 (-5.63, 8.64) | Direction reverses; CI includes 0 |
| 6 | `v190` specification | -4.97 (-10.42, 0.47) | Same direction; modestly attenuated; CI includes 0 | -2.83 (-10.02, 4.36) | Same direction; imprecise |
| 7 | PNC outcome unknown/DK assigned `Y=1` | -5.59 (-10.96, -0.22) | Same direction and similar magnitude; ordinary 95% CI excludes 0, but Holm P=0.0828 and simultaneous CI includes 0 | -2.81 (-10.02, 4.39) | Same direction; imprecise |
| 8 | PNC outcome unknown/DK assigned `Y=0` | -6.15 (-11.55, -0.75) | Same direction and similar magnitude; ordinary 95% CI excludes 0, but Holm P=0.0512 and simultaneous CI includes 0 | -2.66 (-9.89, 4.56) | Same direction; imprecise |

Across all executed sensitivity analyses, the risk-difference estimate remained negative, ranging from -7.02 to -3.00 percentage points. Interval-level support was not uniform: several ordinary 95% intervals excluded zero, but all available multiplicity-adjusted/simultaneous evaluations remained non-confirmatory. The additive-interaction estimate was small and imprecise throughout, and its sign reversed under Scheme B and strict fixed effects. Sensitivity robustness is therefore directional for `RD_11-00`, but not confirmatory or uniformly precise; it is weak for `IC_add`.

The **200-replicate displacement sensitivity** was withdrawn without replacement because the frozen inputs did not mechanically support recomputing the required exposure measures at 200 perturbed locations. It is not a 200 km radius analysis. The executed radius analyses were 25 km and 100 km around the frozen primary 50 km definition.

## Burkina Faso and pooled descriptive supplements

Burkina Faso is **direct-only** for the four joint cells because the adjusted C×K support/overlap gate failed. Its direct survey-weighted descriptive risks were: 00, 17.38% (95% CI 15.56%–19.21%); 10, 10.05% (6.06%–14.04%); 01, 21.52% (17.75%–25.29%); and 11, 26.07% (11.76%–40.37%). The direct survey-weighted descriptive `RD_11-00` was 8.68 pp (95% CI -5.52 to 22.88 pp). These four-cell values are not adjusted or model-standardized joint risks.

The reduced, pre-specified Burkina Faso C+K supplement (`n = 4,539`) produced: C contrast RD -5.41 pp (95% CI -9.73 to -1.09; supportive P=0.0143 unadjusted for multiplicity; RR 0.705, 95% CI 0.508–0.977), and K contrast RD -1.51 pp (95% CI -5.47 to 2.45; supportive P=0.4534 unadjusted for multiplicity; RR 0.917, 95% CI 0.728–1.156). These reduced-model associations do not estimate the joint C×K interaction.

For Burkina Faso Scheme B, the direct survey-weighted descriptive cell risks were 17.86%, 15.58%, 18.16%, and 20.65% for cells 00, 10, 01, and 11, respectively. The reduced C+K model passed its QC gates, but the sealed aggregate output did not supply a standardized contrast; it is not used to replace the primary Burkina Faso supplement.

The pooled analysis is **equal-country, direct-only** and descriptive. The risks were 41.38%, 35.72%, 30.83%, and 29.34% for cells 00, 10, 01, and 11. Its direct `RD_11-00` was -12.04 pp (95% CI -19.67 to -4.40 pp). It is not a pooled adjusted causal effect or a substitute for Nigeria's confirmatory analysis.

## Diagnostics and limitations

- Nigeria's primary model converged in four iterations, was full rank (35/35), had full-rank covariance, no boundary flag, and no warnings. Design degrees of freedom were 1,288. Standardized predictions ranged from 0.526 to 0.609.
- Pearson residuals ranged from -3.04 to 4.62; the 1st and 99th percentiles were -2.49 and 2.30. Anonymous PSU influence norms had median 0.0119, 99th percentile 0.109, and maximum 0.613. No high-influence PSU was deleted.
- Nigeria had 10,694 eligible records, 10,645 with known PNC outcome, 49 with PNC outcome unknown/DK, and 10,643 complete cases. Weighted effective sample size was 7,591.8. Burkina Faso had 4,771 eligible records, 4,547 with known PNC outcome, 224 with PNC outcome unknown/DK, and 4,539 reduced-model complete cases; weighted effective sample size was 3,641.9.
- Burkina Faso's adjusted joint model was not fitted because regional and calendar-month overlap failed the pre-specified gate; direct survey-weighted descriptive cells and the reduced C+K supplement are the appropriate reporting level.
- Cross-sectional/observational measurement, exposure classification, residual confounding, spatial/temporal linkage uncertainty, PNC outcome unknown/DK, and finite support limit interpretation. Associations must not be described as causal.
- Scheme B source-event counts were 1,723 for Nigeria and 545 for Burkina Faso. These are **UCDP source-event records**, not outcome events, not exposed women, and not participant counts.

## Disclosure control

This public summary contains only safe aggregates. It omits identifiable coefficient tables, region identifiers, PSU identifiers, coordinates, microdata, and point-level climate or conflict data.

