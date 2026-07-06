# data/raw/ — canonical inputs (hand-managed, not committed)

The **sources are fully specified** (release + exact filename + download location below and in
`DATA_AVAILABILITY.md`). These large public files are gitignored — download them from the DepMap portal
into this folder locally, then `python src/run_pipeline.py` computes everything into `data/derived/`.
Downloading them is the only remaining step; nothing about the data is unknown or unspecified.

| local file | release | portal file to download | where | key columns |
|---|---|---|---|---|
| `prism_ssdrc.csv` | PRISM Repurposing **19Q4** (secondary) | `secondary-screen-dose-response-curve-parameters.csv` | depmap.org/portal/download/all (19Q4) | depmap_id, compound (broad_id/name), upper_asymptote, lower_asymptote, hill_slope, ec50, auc, ic50, passed_str_profiling, screen_id (HTS002/MTS010) |
| `omics_damaging_mutations.csv` | DepMap Public **24Q2** (resolved 2026-07-05; manuscript to be corrected) | `OmicsSomaticMutationsMatrixDamaging.csv` | depmap.org/portal/download/all (24Q2) | ModelID (ACH-######), gene, damaging value {0,1,2} |
| `model_metadata.csv` | DepMap (matching release) | model / sample-info table | depmap.org/portal | ModelID (ACH), ccle_tissue, lineage |
| `demeter2_rnai.csv` | DEMETER2 (combined RNAi) | combined RNAi gene-dependency | depmap.org/portal (DEMETER2) | ModelID, gene, dependency_score |

Curated reference tables that ARE committed live in `../reference/` (e.g. `validated_sl_reference.csv`) —
do not put those here.

**Provenance mirrors the manuscript's Data Availability statement** (portal, PRISM 19Q4 SSDRC, damaging
mutations matrix) and Corsello et al. 2020 (PMID 32613204) for PRISM methodology. See
`DATA_AVAILABILITY.md` for DOIs and the release resolution (24Q2).
