# Reproducibility confirmations

Specific values the manuscript references, confirmed against the local DepMap/PRISM files by
`src/confirm_claims.py`. The queries use the paper's **named entities** (its 7 RB1-mutant lines) as
keys, so they need no re-derivation of genotype (the damaging-mutations matrix is not required for
confirmation). Machine-readable results: `data/derived/claim_confirmations.csv`.

**Reproduce:**
```
python src/confirm_claims.py \
  --prism  <PRISM working file, e.g. prism_mts010_filtered_data.csv> \
  --model  <DepMap Model.csv> \
  --crispr <CRISPRGeneEffect.csv> \
  --common-essential <AchillesCommonEssentialControls.csv>
```
Run 2026-07-05 against DepMap files dated Oct–Nov 2024 (PRISM MTS010 working file; DepMap Model.csv
= 1,959 models, consistent with 24Q2; CRISPRGeneEffect.csv Chronos).

## Q1 — RB1-mutant lung cohort identity  ✅ CONFIRMED
All 7 named lines resolve to ModelIDs and are Lung:

| Line | ModelID | Subtype |
|---|---|---|
| NCIH2228 | ACH-000447 | NSCLC |
| HCC44 | ACH-000667 | NSCLC |
| DMS273 | ACH-000749 | Lung neuroendocrine (SCLC) |
| NCIH446 | ACH-000800 | Lung neuroendocrine (SCLC) |
| T3M10 | ACH-000813 | NSCLC |
| NCIH1048 | ACH-000866 | Lung neuroendocrine (SCLC) |
| HCC15 | ACH-000878 | NSCLC |

The cohort is **4 NSCLC + 3 SCLC/neuroendocrine** — a mix, not purely SCLC.

## Q2 — Common-essential status  ⚠ SOURCES DISAGREE on AURKB
DepMap's Achilles common-essential control list vs the R1.7 grounding (which used BioGRID ORCS):

| Gene | DepMap Achilles | BioGRID ORCS (grounding) |
|---|---|---|
| PLK1 | **common-essential** | common-essential (929) — agree |
| CHEK1 | **common-essential** | common-essential (928) — agree |
| **AURKB** | **selective (not on list)** | common-essential (912) — **disagree** |
| AURKA | selective | — |
| PARP1 | selective | context-selective (67) — agree |

DepMap does **not** classify AURKB as common-essential. This is *favorable* to the paper (AURKB reads
as a selective vulnerability), but the letter's current "PLK1/AURKB/CHEK1 are common-essential" phrasing
overstates it — per DepMap, only PLK1 and CHEK1 are.

## Q3 — RB1-mut vs WT-lung CRISPR (Chronos) differential  ⚠ SMALL, AURKB-led
Real gene-effect, RB1-mut (n=5 present in CRISPR data) vs WT-lung (n=65); more negative = more dependent:

| Gene | mut mean | WT mean | differential | direction |
|---|---|---|---|---|
| **AURKB** | −2.376 | −2.228 | **−0.148** | mutant more dependent (strongest) |
| PLK1 | −2.861 | −2.760 | −0.101 | mutant more dependent (small) |
| PARP1 | −0.174 | −0.160 | −0.015 | negligible |
| CHEK1 | −1.784 | −1.790 | +0.007 | none |
| AURKA | −1.256 | −1.317 | +0.061 | wrong direction |

The differentials are **modest and not uniform**. **AURKB is the clearest RB1-selective dependency**
(and is not common-essential, Q2) — consistent with Oser 2019 (RB1→AURKB). PARP1 shows ~no CRISPR
dependency differential, expected because PARP inhibitors act by **trapping**, not by creating a
knockout dependency — so CRISPR/RNAi loss-of-function *understates* PARP1's synthetic lethality.

**Caveats:** this is CRISPR (Chronos), whereas the manuscript's R1.7 validation used **RNAi (DEMETER2)** —
a different modality that can disagree; n=5 RB1-mut lines are in the CRISPR data (2 of 7 absent); "WT-lung"
is lung∩PRISM minus the 7 named lines (approximate, since the damaging-mutations matrix is not on disk).

## Q4 — Worked Example 3 (doxorubicin pS′ = 8.15)  ◻ UNAVAILABLE
Doxorubicin is not present for these lines in the filtered MTS010 working file, so Example 3's
pS′ = 8.15 cannot be reproduced from the data currently on disk (needs the compound's rows).

## Bottom line for R1.7
The real DepMap data supports the *direction* of the RB1→Aurora story but reorders it: **lead with
AURKB** (selective per DepMap + mutant-more-dependent in CRISPR + Oser SL literature), treat **PLK1/CHEK1**
as common-essential with small differentials, and cite **PARP1** on its **trapping** mechanism rather than
a dependency differential. The current letter follows the BioGRID-ORCS grounding, which led with PARP1 and
called AURKB common-essential — worth reconciling against these on-disk DepMap results.
