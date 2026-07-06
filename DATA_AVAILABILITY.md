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
| Damaging-mutation matrix | **DepMap Public 24Q2** *(resolved — see below)* | `OmicsSomaticMutationsMatrixDamaging.csv` | DepMap portal <https://depmap.org/portal/download/all> |
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

## Mutation-matrix release: RESOLVED — DepMap Public 24Q2
**The release used is 24Q2** (resolved 2026-07-05). Evidence, without needing the matrix itself:
- The analysis working files on disk (`prism_mts010_filtered_data.csv`, `Model.csv`, CRISPR exports) date
  **Oct–Nov 2024** — *before* DepMap 24Q4 was released (**2024-12-10**), so 24Q4 could not have been used.
- The local DepMap `Model.csv` has **1,959 models**, consistent with **24Q2** (24Q3/24Q4 have more; 23Q4 fewer).
- The R1.8 preprocessing flow (Concept Figure A/B) and `config/definitions.yaml` already say 24Q2.

The **manuscript must be corrected to 24Q2**: its reference **[24] cites "DepMap 23Q4 Public"** and the
Data Availability statement names no quarter (some materials also said "24Q4") — all should read 24Q2.

If the matrix ever needs a *checksum-level* confirmation, the two releases are different files and
self-identifying (`certutil -hashfile <file> MD5`, or size alone):

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

This repo (this file, `data/raw/README.md`, `config/definitions.yaml`) is now consistent on **24Q2**;
the manuscript Data Availability statement and ref [24] remain to be corrected to 24Q2 by the authors.
*Checksums verified against the DepMap Figshare+ manifests, 2026-07-05.* (PRISM Repurposing 19Q4 is
consistent across all sources; no conflict there.)

---
This repository is the processed-data-and-code deposit requested by Reviewer 1 (comment 11) and supports
the metric-transparency point raised by Reviewer 2 (comment 1).
