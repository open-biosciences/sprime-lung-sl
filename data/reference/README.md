# data/reference/ — curated, version-controlled reference tables

Small hand-curated tables that ARE committed (unlike the large external inputs in `data/raw/`,
which are gitignored). The pipeline reads these as ground truth for validation/concordance.

## validated_sl_reference.csv
Known genotype-specific vulnerabilities used for the S′-concordance check (R1.7 / C9). Each row is
grounded against a public resource and dated in `grounded`:
- `orcs_hit_screens` / `essentiality_class` — BioGRID ORCS CRISPR essentiality (common-essential vs
  context-selective). **Common-essential targets require reporting the mutant-vs-WT *differential*,
  not raw dependency.**
- `compound_chembl_id` / `compound_moa` — ChEMBL mechanism-of-action.
- `sl_evidence_pmid` / `sl_evidence_doi` — published synthetic-lethality evidence.

Seeded 2026-07-05 from the bio-research grounding pass (see `GROUNDING_REPORT_R1.7_validation.md` in
the paper workspace). Extend as more hits are validated; re-verify identifiers before manuscript use.
