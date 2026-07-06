# Data & code availability

All primary data are public and **fully sourced below** (release, exact filename, DOI/URL). Nothing
about the inputs is unspecified — the only step not captured in git is the local download of the large
files themselves (redistributed per DepMap terms and kept out of version control; see
`data/raw/README.md`). This section mirrors, and sharpens, the manuscript's Data Availability statement.

## Code
- **Method module:** https://github.com/open-biosciences/sprime
- **This analysis:** https://github.com/open-biosciences/sprime-lung-sl  (mint a Zenodo DOI on tagged release)

## Primary data (DepMap / Broad Institute)

| Input | Release | File | Source / DOI |
|---|---|---|---|
| Drug-response 4PL parameters | **PRISM Repurposing 19Q4** (secondary screen) | `secondary-screen-dose-response-curve-parameters.csv` | DepMap portal <https://depmap.org/portal/download/all>; method Corsello et al. 2020, *Nat Cancer*, PMID 32613204, doi:10.1038/s43018-019-0018-6; DepMap 19Q4 Public figshare doi:10.6084/m9.figshare.11384241.v2 |
| Damaging-mutation matrix | **DepMap Public 24Q2** *(see discrepancy note)* | `OmicsSomaticMutationsMatrixDamaging.csv` | DepMap portal <https://depmap.org/portal/download/all> |
| Model / cell-line metadata | DepMap (matching release) | model file / `ccle_tissue`, lineage | DepMap portal |
| RNAi gene-dependency (validation) | **DEMETER2** (combined Broad/Novartis/Marcotte) | combined RNAi dependency | DepMap portal / Broad DEMETER2 |
| Screen selection & QC | HTS002 (base) + MTS010 (overlay of corrected values); rows filtered on `passed_str_profiling == TRUE` | — | documented in Methods / preprocessing flow (Concept Figure A/B) |
| Cell-line data-source notes | — | — | https://github.com/mocomakers/nf_streamlit/ |

## Curated reference (committed in this repo)
- `data/reference/validated_sl_reference.csv` — known genotype-specific vulnerabilities used for the
  S′-concordance check, each row grounded (BioGRID ORCS essentiality, ChEMBL MoA, PubMed SL evidence),
  seeded 2026-07-05. See `GROUNDING_REPORT_R1.7_validation.md` in the paper workspace.

## Methodological antecedent (peer-reviewed)
- Zamora et al., *Cancers* 2023, 15(24):5811, doi:10.3390/cancers15245811 (log-based predecessor of S′).

## ⚠ Discrepancy to reconcile before submission
The submitted manuscript's Data Availability statement cites the damaging-mutation matrix from the
**"24Q4 Releases,"** while the R1.8 preprocessing flow (Concept Figure A/B) and the pipeline config
(`config/definitions.yaml`) cite **DepMap Public 24Q2**. Confirm which release was actually used and
make the manuscript, this file, `data/raw/README.md`, and `definitions.yaml` agree. (PRISM 19Q4 is
consistent across all sources.)

---
This repository is the processed-data-and-code deposit requested by Reviewer 1 (comment 11) and supports
the metric-transparency point raised by Reviewer 2 (comment 1).
