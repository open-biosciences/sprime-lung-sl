# CLAUDE.md — sprime-lung-sl

Orientation for a Claude Code session in this repo.

## What this is
The **paper companion** to Zamora et al. (S′ synthetic-lethality in lung): the reproducible
analysis + data deposit that Reviewer 1 (comment 11) asked for. It imports the `sprime` method
module (sibling repo) and turns raw DepMap/PRISM inputs into every number the manuscript reports.

## The single-source-of-truth contract
`config/definitions.yaml` pins **all** conventions and thresholds. `src/run_pipeline.py` reads
ONLY `config/definitions.yaml` + `data/raw/` and writes `data/derived/*.csv`. No number in the
paper is hand-typed — it comes out of the pipeline. `docs/concerns.csv` (C1–C19) is the project
spine: each reviewer concern points at the artifact that resolves it.

## Layout
- `config/definitions.yaml` — conventions: E_max=percent, C_ref=1µM, SL rule (WT pS′>0 & mut pS′>0
  & ΔpS′≤−2), high-stringency ΔpS′≤−4, genotype rule (damaging-AF-sum>0.95 ⇒ mutant), FDR run-not-gated.
- `src/run_pipeline.py` — **STUB**. Needs raw CSVs in `data/raw/`, then wire in `sprime` + emit derived.
- `src/figures.py` — figure generation from `data/derived/` (Fig 1/4/5 fixes = concern C12).
- `data/raw/` — DepMap/PRISM inputs (empty until you drop them in; see `data/raw/README.md`).
- `data/derived/` — the "numbers ledger": psprime_by_genotype.csv, genotype_counts.csv,
  metric_correlations.csv, rnai_validation.csv, concordance_summary.csv, ...
- `docs/METHODS_CANONICAL.md`, `DATA_AVAILABILITY.md`, `docs/concerns.csv`.

## Current state / next work
The pipeline is a stub and `data/raw/` is empty — the blocker is dropping in the DepMap 24Q2
`OmicsSomaticMutationsMatrixDamaging.csv` and PRISM 19Q4 SSDRC parameters, then implementing
`run_pipeline.main()` per the numbered contract in the file. Until then, counts cited in the
response letter (e.g. RB1 7/75) are placeholders to be replaced by real derived values (C4, C15).

## E_max convention (concerns.csv C1)
Pipeline uses the **percent** E_max scale (matches `sprime`). OPEN: confirm against the reference
SPrime implementation before release. Do not change the scale without updating `sprime/tests` too.

## Run
    pip install -r requirements.txt
    pip install -e ../sprime         # method module
    python src/run_pipeline.py       # once data/raw/ is populated and main() is implemented
    python src/figures.py

## License
MIT (code) + CC-BY-4.0 (processed data). `CITATION.cff` lists the manuscript author team.
