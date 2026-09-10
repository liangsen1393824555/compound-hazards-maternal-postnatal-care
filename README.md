# Compound hazards and maternal postnatal care

This repository accompanies the manuscript **“Compound heat–drought and nearby organized violence in relation to nonreceipt of maternal postnatal care within two days: a prospectively locked survey-weighted analysis in Nigeria with descriptive Burkina Faso supplements.”**

## Study status and interpretation

The analysis plan and estimands were locked before results were unsealed. Neither of the two Nigeria confirmatory estimands passed the prespecified Holm-adjusted 0.05 familywise threshold. The negative joint-versus-neither risk-difference estimate is therefore descriptive rather than confirmatory evidence of a protective or causal effect. Burkina Faso and equal-country results are descriptive or supportive only.

## Repository contents

- `code/`: data-construction, survey-analysis, sensitivity-analysis, and dual-run comparison programs.
- `protocol/`: the locked statistical analysis plan, pre-unlock amendment, and decision log.
- `results/`: disclosure-controlled aggregate results, the corrected results narrative, and the erratum documenting the v1.1-to-v1.2 correction.
- `figures/`: aggregate manuscript figures.
- `DATA_ACCESS.md`: data-access and privacy restrictions.
- `SHA256SUMS.txt`: checksums for the public release.

## Reproducibility boundary

The programs operate on inputs that the researcher must obtain lawfully from their original providers. This repository intentionally excludes DHS survey microdata, DHS geographic data, cluster and household identifiers, coordinates, row-level linked data, point-level climate values, and event-level conflict records. A reader without separately authorized inputs can inspect the locked methods and reproduce all reported manuscript numbers from `results/public_results_v1.2.json`, but cannot reconstruct the restricted row-level analysis file from this repository alone.

Python dependencies are listed in `requirements.txt`; R dependencies are listed in `requirements-r.txt`. The original execution environment used Python 3.12 and R with the `survey` and `jsonlite` packages. Paths are supplied at runtime; no credentials or user-specific paths are embedded in the release.

## Main aggregate result

For Nigeria, the adjusted joint-versus-neither risk difference was −6.03 percentage points (ordinary 95% CI −11.45 to −0.61; Holm-adjusted P=0.0584; Bonferroni 97.5% simultaneous CI −12.23 to 0.17). The additive interaction contrast was −2.80 percentage points (ordinary 95% CI −10.04 to 4.44; Holm-adjusted P=0.4484; simultaneous CI −11.08 to 5.48). These findings are confirmatory-null or inconclusive under the prespecified decision rule.

## AI-use disclosure

During preparation of this work, the author used OpenAI Codex to assist with programming, quality-control documentation, language editing, and manuscript preparation. The author reviewed and verified the outputs and takes full responsibility for the content.

## Author and citation

Sen Liang (ORCID: 0009-0009-2911-1795). See `CITATION.cff` for machine-readable citation metadata.

## License status

No software or content license has yet been assigned. Public visibility does not by itself grant reuse rights; third-party datasets remain governed by their providers’ terms.
