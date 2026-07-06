# sprime-lung-sl

Reproducible analysis for **"Use of S′ Indices to Evaluate Synthetic Lethality as a Systems
Phenomenon in Lung Cancer Cells"** (Zamora et al., submitted to *Precision Oncology*).

This is the paper-companion repository. It turns canonical DepMap/PRISM inputs into the S′-based
synthetic-lethality results, using the [`sprime`](https://github.com/open-biosciences/sprime)
method module. **Everything downstream — tables, figures, and the numbers in the manuscript and
response letter — regenerates from `data/raw/` + `config/definitions.yaml`, so values cannot drift.**

## Principle: one source of truth
1. `config/definitions.yaml` pins every convention (E_max, C_ref, SL thresholds, genotype rule).
2. `src/run_pipeline.py` reads `data/raw/` + definitions → writes `data/derived/*.csv`.
3. `src/figures.py` renders figures from `data/derived/` only.
4. The manuscript/response cite cells in `data/derived/` (a "numbers ledger"), not hand-typed values.

## Layout
    config/definitions.yaml      # pinned conventions — the single source of truth
    data/raw/                    # canonical inputs (see data/raw/README.md; not committed if large)
    data/derived/                # generated outputs (see data/derived/README.md)
    src/run_pipeline.py          # raw + definitions -> derived
    src/figures.py               # derived -> figures/
    docs/METHODS_CANONICAL.md    # frozen method text (mirrors definitions.yaml)
    docs/concerns.csv            # reviewer-concern tracker (project spine)
    DATA_AVAILABILITY.md         # what to cite; DOI once minted

## Reproduce
    pip install -r requirements.txt
    python src/run_pipeline.py            # writes data/derived/*.csv
    python src/figures.py                 # writes figures/*
    pytest                                # (via sprime) reconciliation tests

## ⚠ Blocking open item
The **E_max convention** (percent vs fraction) is not yet locked in `definitions.yaml`; it sets the
S′ scale and the ±2/±4 thresholds. Resolve against the `sprime` reference before generating
final numbers. See `docs/concerns.csv` (C1).

## License
Code: MIT. Processed data: CC-BY-4.0 (confirm with Open Biosciences). Raw DepMap/PRISM data are
redistributed per their original terms — see DATA_AVAILABILITY.md.
