# data/raw/ — canonical inputs (hand-managed)

Large public files are **not committed**; place them here locally (see `.gitignore`) and record
provenance in DATA_AVAILABILITY.md. Expected files & key columns:

| file | source | key columns |
|---|---|---|
| `prism_ssdrc.csv` | PRISM Repurposing 19Q4 | depmap_id, compound, upper_asymptote, lower_asymptote, hill, ec50, auc, ic50, passed_str_profiling, screen_id |
| `omics_damaging_mutations.csv` | DepMap Public 24Q2 | ModelID (ACH), gene, damaging_value(0/1/2) |
| `model_metadata.csv` | DepMap | ModelID (ACH), ccle_tissue, lineage |
| `demeter2_rnai.csv` | DEMETER2 | ModelID, gene, dependency_score |
| `validated_sl_reference.csv` | curated (this study) | gene, compound, moa, interaction, source_pmid |
