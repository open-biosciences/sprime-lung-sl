"""End-to-end smoke test for run_pipeline against tiny synthetic DepMap/PRISM inputs.

Locks the pipeline's behavior without needing the real (gitignored) matrices: genotype counting,
S' computation, pS'/dpS' classification (a planted Class-A hit must be recovered), the STR-profiling
filter, and that every derived CSV is written. Requires the analysis deps + sprime:
    pip install -r requirements.txt && pip install -e ../sprime
Run with:  PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 pytest tests/test_pipeline_smoke.py
"""
import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

SRC = Path(__file__).resolve().parents[1] / "src"


def _load_run_pipeline():
    sys.path.insert(0, str(SRC))
    spec = importlib.util.spec_from_file_location("run_pipeline", SRC / "run_pipeline.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _synth_inputs(raw: Path):
    """5 lung lines (RB1 mut/wt/excluded mix) + non-lung; 3 compounds incl. a planted RB1 Class-A hit."""
    raw.mkdir(parents=True, exist_ok=True)
    lines = [f"ACH-{i:04d}" for i in range(10)]
    rb1 = [2, 2, 2, 0, 0, 1, 0, 0, 0, 0]                # lung 0..5: mut=3, wt=2, excluded=1
    pten = [0, 0, 2, 2, 0, 0, 0, 0, 0, 0]
    model = pd.DataFrame({
        "ModelID": lines,
        "StrippedCellLineName": [f"CELL{i}" for i in range(10)],
        "OncotreeLineage": ["Lung"] * 6 + ["Skin"] * 4,
    })
    matrix = pd.DataFrame({"RB1 (5925)": rb1, "PTEN (5728)": pten,
                           "CDKN2A (1029)": [0] * 10, "TP53 (7157)": [0] * 10}, index=lines)
    matrix.index.name = "ModelID"

    rows = []
    for i in range(6):
        if rb1[i] == 2:   # RB1 mutant: high S' (potent, efficacious)
            rows.append(dict(depmap_id=lines[i], name="DRUGX", upper_limit=1.0, lower_limit=0.01,
                             ec50=0.1, auc=0.3, passed_str_profiling="TRUE"))
        else:             # non-mutant: positive but lower S'
            rows.append(dict(depmap_id=lines[i], name="DRUGX", upper_limit=0.7, lower_limit=0.1,
                             ec50=2.0, auc=0.7, passed_str_profiling="TRUE"))
        rows.append(dict(depmap_id=lines[i], name="DRUGY", upper_limit=0.5, lower_limit=0.1,
                         ec50=1.0, auc=0.6, passed_str_profiling="TRUE"))
        rows.append(dict(depmap_id=lines[i], name="DRUGZ", upper_limit=0.9, lower_limit=0.2,
                         ec50=0.5, auc=0.5, passed_str_profiling="FALSE"))   # must be dropped by STR filter

    matrix.to_csv(raw / "omics_damaging_mutations.csv")
    model.to_csv(raw / "model_metadata.csv", index=False)
    pd.DataFrame(rows).to_csv(raw / "prism_ssdrc.csv", index=False)


def test_pipeline_smoke(tmp_path):
    rp = _load_run_pipeline()
    raw, derived = tmp_path / "raw", tmp_path / "derived"
    _synth_inputs(raw)
    rp.cr.RAW = raw           # so check_release._resolve finds the synthetic inputs
    rp.DERIVED = derived

    assert rp.main([]) == 0

    counts = pd.read_csv(derived / "genotype_counts.csv").set_index("gene")
    assert counts.loc["RB1", "n_mutant"] == 3
    assert counts.loc["RB1", "n_wild_type"] == 2
    assert counts.loc["RB1", "n_excluded"] == 1

    ps = pd.read_csv(derived / "psprime_by_genotype.csv")
    hit = ps[(ps.gene == "RB1") & (ps.compound == "DRUGX")].iloc[0]
    assert bool(hit["class_A"]) is True                 # planted Class-A recovered
    assert hit["delta_ps"] <= -2

    assert "DRUGZ" not in set(ps.compound)              # STR-failed compound excluded
    assert (derived / "metric_correlations.csv").exists()
