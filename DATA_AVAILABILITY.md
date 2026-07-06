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

## ⚠ Mutation-matrix release: 24Q2 vs 24Q4 — resolve by checksum
The submitted manuscript's Data Availability statement cites the damaging-mutation matrix from the
**"24Q4 Releases,"** while the R1.8 preprocessing flow (Concept Figure A/B) and `config/definitions.yaml`
cite **DepMap Public 24Q2**. The file `OmicsSomaticMutationsMatrixDamaging.csv` exists in **both**
releases, but they are **different files** (24Q4 added cell lines), so the actual downloaded copy is
self-identifying — no need to rely on memory:

| Release | Figshare | Size (bytes) | MD5 |
|---|---|---|---|
| DepMap Public **24Q2** (2024-05-23) | doi:10.25452/figshare.plus.25880521.v1 | 135,585,447 | `02f3568b71af0ca3e8d10e681eefac86` |
| DepMap Public **24Q4** (2024-12-10) | doi:10.25452/figshare.plus.27993248.v1 | 147,655,356 | `cb20fdbe1cf3b9b0d8ed4f53e1f399b6` |

**To settle it**, on the machine that ran the analysis:
```
md5sum OmicsSomaticMutationsMatrixDamaging.csv    # or: certutil -hashfile <file> MD5   (Windows)
ls -l  OmicsSomaticMutationsMatrixDamaging.csv    # size alone also distinguishes the two
```
**Complementary value-based check** (tells you whether the release choice even *matters*):
`python src/check_release.py` reproduces the paper's own RB1 cohort from a matrix on disk — the release
whose calls yield the 7 named RB1-mutant lung lines (all damaging value == 2) plus the reported counts
is the one used. With both matrices present, `--compare` diffs the analyzed-gene calls: identical ⇒ the
discrepancy is cosmetic (fix the manuscript string only); different ⇒ regenerate the counts on the
confirmed release. Expectations are pinned in `config/definitions.yaml → expected_cohorts`.

Then set that release consistently across: the manuscript Data Availability statement, this file,
`data/raw/README.md`, and `config/definitions.yaml`. (Because the two versions differ in cell-line
coverage, if it turns out 24Q4 was used, the per-genotype n and counts should be regenerated on 24Q4.)
*Verified against the DepMap Figshare+ manifests, 2026-07-05.* (PRISM Repurposing 19Q4 is consistent
across all sources; no conflict there.)

---
This repository is the processed-data-and-code deposit requested by Reviewer 1 (comment 11) and supports
the metric-transparency point raised by Reviewer 2 (comment 1).
