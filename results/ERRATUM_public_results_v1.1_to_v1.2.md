# Erratum: public results v1.1 → v1.2

Date: 2026-09-11

This erratum corrects labels and wording in the public v1.1 materials. The sealed source results, estimates, confidence intervals, P values, model specifications, mechanical conclusion, and v1.1 files are unchanged.

## Corrections

1. **Outcome unknown/DK:** Nigeria's 49 and Burkina Faso's 224 records are PNC outcome unknown/don't know, not unknown exposure. The sample-flow labels, limitations, and sensitivity labels now state outcome unknown/DK. `unknown_as_1` means the PNC outcome unknown/DK records were assigned `Y=1`; `unknown_as_0` means they were assigned `Y=0`.
2. **Scheme B Nigeria:** the sealed analysis is an adjusted survey-weighted logistic model with standardized risks. Its reported sensitivity estimates are supportive adjusted estimates. The Scheme B P values are unadjusted **for multiplicity**; the model itself is not unadjusted.
3. **Burkina Faso four-cell results:** these are direct survey-weighted descriptive risks, not adjusted or model-standardized joint risks.
4. **Displacement sensitivity:** the withdrawn item is the 200-replicate displacement sensitivity, not a 200 km radius sensitivity. The 25 km, primary 50 km, and 100 km exposure radii are distinct from the 200-replicate procedure.
5. **Machine-readable labels:** corresponding JSON keys and limitation text were corrected from exposure-unknown wording to outcome unknown/DK; the withdrawn analysis identifier is now `200_replicate_displacement`.

## Unchanged conclusions and numbers

- Mechanical conclusion remains **CONFIRMATORY_NULL_OR_INCONCLUSIVE**.
- Neither confirmatory estimand passes Holm α=0.05.
- All numeric leaves in `public_results_v1.2.json` match `public_results_v1.1.json` item-for-item after semantic key alignment and version-label exclusions.
- No model was refitted, reselected, or modified.

