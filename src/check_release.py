"""Reconcile the DepMap mutations release (24Q2 vs 24Q4) from the paper's own reported cohort.

The manuscript's Data Availability says the damaging-mutation matrix is "24Q4"; the R1.8 flow
diagram and config/definitions.yaml say "24Q2" (see the RECONCILE note in definitions.yaml).
Because DepMap re-runs mutation calling each quarter, the same line can carry a different
RB1/PTEN/TP53 damaging value across releases — so the release actually used is the one whose
matrix reproduces the paper's reported RB1 cohort.

Two paper-derived tests (expectations pinned in config/definitions.yaml -> expected_cohorts):
  1. NAMED COHORT — the 7 named RB1-mutant lung lines must all carry RB1 damaging value == mutant_value.
  2. COUNT        — # RB1-mutant lines (lung ∩ PRISM-19Q4) must equal the reported n_mutant (7).

Nothing here fetches from the web: the matrices are 135–147 MB local files (redistribution-restricted,
gitignored). Drop them into data/raw/ and run:

    python src/check_release.py                                  # checks the matrix in data/raw/
    python src/check_release.py --matrix path/to/24Q2.csv --label 24Q2
    python src/check_release.py --matrix 24Q2.csv --label 24Q2 --compare 24Q4.csv --label2 24Q4

If both releases give identical calls for the analyzed genes/lines, the discrepancy is cosmetic
(just fix the manuscript's release string); if they differ, re-run the pipeline on the confirmed release.
"""
from __future__ import annotations
import argparse
import re
import sys
from pathlib import Path

import yaml

# pandas is only needed once real files are on disk; import lazily so the "where to put the files"
# help still prints in a broken/absent-pandas environment.
try:
    import pandas as pd
except Exception:  # pragma: no cover - environment-dependent
    pd = None

ROOT = Path(__file__).resolve().parents[1]
DEFS = yaml.safe_load((ROOT / "config" / "definitions.yaml").read_text())
RAW = ROOT / "data" / "raw"

# Accept either the short local name (data/raw/README.md) or the canonical DepMap download name.
MATRIX_CANDIDATES = ["omics_damaging_mutations.csv", "OmicsSomaticMutationsMatrixDamaging.csv"]
MODEL_CANDIDATES = ["model_metadata.csv", "Model.csv", "sample_info.csv"]
PRISM_CANDIDATES = ["prism_ssdrc.csv", "secondary-screen-dose-response-curve-parameters.csv"]


def _norm(s) -> str:
    """Normalize a cell-line name for matching: uppercase, strip non-alphanumerics."""
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def _resolve(explicit, candidates, what):
    if explicit:
        p = Path(explicit)
        return p if p.exists() else None
    for name in candidates:
        p = RAW / name
        if p.exists():
            return p
    return None


def _find_gene_column(columns, symbol, entrez):
    """DepMap matrix columns look like 'RB1 (5925)'. Match by symbol first, then Entrez."""
    want = _norm(symbol)
    for c in columns:
        head = re.split(r"[\s(]", str(c), 1)[0]
        if _norm(head) == want:
            return c
    if entrez is not None:
        for c in columns:
            if str(entrez) in str(c):
                return c
    return None


def _name_column(model: pd.DataFrame):
    for c in ["StrippedCellLineName", "CellLineName", "CCLEName", "stripped_cell_line_name", "CCLE_Name"]:
        if c in model.columns:
            return c
    return None


def _model_id_column(model: pd.DataFrame):
    for c in ["ModelID", "DepMap_ID", "depmap_id", "depMapID"]:
        if c in model.columns:
            return c
    return None


def _lineage_mask_lung(model: pd.DataFrame):
    for c in ["OncotreeLineage", "lineage", "Lineage", "ccle_tissue", "primary_tissue"]:
        if c in model.columns:
            return model[c].astype(str).str.contains("lung", case=False, na=False)
    return None


def _map_names(model: pd.DataFrame, names):
    """Return {input_name: ModelID} plus a list of unmatched names."""
    idc, namec = _model_id_column(model), _name_column(model)
    if idc is None or namec is None:
        raise SystemExit(f"! model table missing ID/name columns (have: {list(model.columns)[:8]}...)")
    lut = {_norm(n): mid for n, mid in zip(model[namec], model[idc])}
    out, missing = {}, []
    for n in names:
        mid = lut.get(_norm(n))
        (out.__setitem__(n, mid) if mid is not None else missing.append(n))
    return out, missing


def _load_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0, low_memory=False)


def _rb1_cohort_config():
    ec = (DEFS.get("expected_cohorts") or {}).get("RB1")
    if not ec:
        raise SystemExit("! definitions.yaml has no expected_cohorts.RB1 block")
    g = DEFS["genotype"]
    return ec, g["mutant_value"], g["wild_type_value"], g["excluded_value"]


def test_named_cohort(matrix, model):
    """Test 1 — every named RB1-mutant line carries value == mutant_value. Prints per-line evidence."""
    ec, mut, wt, exc = _rb1_cohort_config()
    col = _find_gene_column(matrix.columns, ec["gene_symbol"], ec.get("gene_entrez"))
    if col is None:
        print(f"  [NAMED]  UNKNOWN — no {ec['gene_symbol']} column in matrix"); return None
    ids, missing = _map_names(model, ec["mutant_lines"])
    print(f"  [NAMED]  gene column = {col!r}; mutant_value = {mut}")
    ok = True
    for name, mid in ids.items():
        if mid is None or mid not in matrix.index:
            print(f"           {name:10s} {mid or '(no ModelID)'}: not in matrix"); ok = False; continue
        v = matrix.at[mid, col]
        flag = "OK" if v == mut else ("WT" if v == wt else ("ambiguous" if v == exc else "OTHER"))
        print(f"           {name:10s} {mid}: value={v}  [{flag}]")
        ok = ok and (v == mut)
    if missing:
        print(f"           unmapped names (fix via model table): {missing}"); ok = False
    print(f"  [NAMED]  {'PASS — all 7 are mutant in this release' if ok else 'FAIL — cohort not reproduced'}")
    return ok


def test_count(matrix, model, prism):
    """Test 2 — RB1-mutant count over (lung ∩ PRISM-19Q4) equals reported n_mutant."""
    ec, mut, wt, exc = _rb1_cohort_config()
    col = _find_gene_column(matrix.columns, ec["gene_symbol"], ec.get("gene_entrez"))
    if col is None:
        print("  [COUNT]  UNKNOWN — no RB1 column"); return None
    lung = _lineage_mask_lung(model)
    if lung is None:
        print("  [COUNT]  UNKNOWN — no lineage column in model table"); return None
    idc = _model_id_column(model)
    lung_ids = set(model.loc[lung, idc])

    prism_ids = None
    if prism is not None:
        pcol = next((c for c in ["depmap_id", "DepMap_ID", "ModelID"] if c in prism.columns), None)
        if pcol:
            prism_ids = set(prism[pcol].dropna().unique())
    universe = lung_ids & prism_ids if prism_ids else lung_ids
    scope = "lung ∩ PRISM-19Q4" if prism_ids else "lung (PRISM file absent — not intersected)"

    col_vals = matrix[col]
    in_scope = [i for i in universe if i in matrix.index]
    n_mut = sum(col_vals[i] == mut for i in in_scope)
    n_wt = sum(col_vals[i] == wt for i in in_scope)
    exp_mut, exp_wt = ec["n_mutant"], ec.get("n_wild_type")
    print(f"  [COUNT]  scope = {scope}; lines in matrix = {len(in_scope)}")
    print(f"  [COUNT]  RB1-mutant = {n_mut} (expect {exp_mut});  WT = {n_wt} (expect {exp_wt})")
    ok = (n_mut == exp_mut)
    print(f"  [COUNT]  {'PASS' if ok else 'MISMATCH'} on mutant count")
    return ok


def run_one(label, matrix_path, model, prism):
    print(f"\n=== Release under test: {label or matrix_path.name} ({matrix_path.name}) ===")
    matrix = _load_matrix(matrix_path)
    print(f"  matrix: {matrix.shape[0]} lines × {matrix.shape[1]} genes")
    r1 = test_named_cohort(matrix, model)
    r2 = test_count(matrix, model, prism)
    verdict = ("MATCHES the paper cohort — likely the release used" if (r1 and r2)
               else "does NOT reproduce the paper cohort" if (r1 is False or r2 is False)
               else "inconclusive (missing inputs)")
    print(f"  VERDICT [{label or matrix_path.name}]: {verdict}")
    return matrix


def diff_releases(m_a, label_a, m_b, label_b):
    """If both matrices are present, diff the analyzed-gene columns across shared lines."""
    print(f"\n=== Diff {label_a} vs {label_b} (analyzed genes) ===")
    genes = DEFS["genotype"]["genes"]
    shared = m_a.index.intersection(m_b.index)
    any_diff = False
    for g in genes:
        ca = _find_gene_column(m_a.columns, g, None)
        cb = _find_gene_column(m_b.columns, g, None)
        if ca is None or cb is None:
            print(f"  {g:7s}: column missing in one release — skipped"); continue
        d = (m_a.loc[shared, ca] != m_b.loc[shared, cb]).sum()
        any_diff = any_diff or d > 0
        print(f"  {g:7s}: {d} of {len(shared)} shared lines differ")
    print("  => " + ("calls DIFFER — release choice moves the numbers; re-run on the confirmed release."
                     if any_diff else
                     "calls IDENTICAL for analyzed genes — discrepancy is cosmetic; just fix the manuscript string."))


def main(argv=None):
    ap = argparse.ArgumentParser(description="Reconcile DepMap mutations release from the paper's RB1 cohort.")
    ap.add_argument("--matrix", help="damaging-mutations matrix CSV (default: search data/raw/)")
    ap.add_argument("--label", help="release label for the primary matrix, e.g. 24Q2")
    ap.add_argument("--compare", help="second matrix CSV to diff against (e.g. the other release)")
    ap.add_argument("--label2", help="release label for the second matrix, e.g. 24Q4")
    ap.add_argument("--model", help="model/sample-info table (default: search data/raw/)")
    ap.add_argument("--prism", help="PRISM SSDRC parameters CSV (default: search data/raw/)")
    args = ap.parse_args(argv)

    try:  # emit UTF-8 so ∩/×/— render on Windows consoles
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

    matrix_path = _resolve(args.matrix, MATRIX_CANDIDATES, "mutations matrix")
    model_path = _resolve(args.model, MODEL_CANDIDATES, "model table")
    prism_path = _resolve(args.prism, PRISM_CANDIDATES, "PRISM SSDRC")

    if matrix_path is None or model_path is None:
        print("Nothing to check yet — required local files not found.\n")
        print(f"  Put the DepMap files in {RAW} (see data/raw/README.md and DATA_AVAILABILITY.md):")
        print(f"    - mutations matrix: one of {MATRIX_CANDIDATES}  {'FOUND' if matrix_path else 'MISSING'}")
        print(f"    - model table:      one of {MODEL_CANDIDATES}  {'FOUND' if model_path else 'MISSING'}")
        print(f"    - PRISM SSDRC:      one of {PRISM_CANDIDATES}  {'(optional) FOUND' if prism_path else '(optional) MISSING'}")
        print("\nThen re-run. This is the same trigger as wiring up run_pipeline.py.")
        return 2

    if pd is None:
        print("Files are present, but pandas failed to import in this environment.")
        print("  (Seen here: a numpy/pandas binary mismatch — pandas built for numpy 1.x vs installed 2.x.)")
        print("  Fix with a matched pair, e.g.:  pip install 'numpy<2' 'pandas>=1.4'   (or upgrade pandas>=2.1).")
        return 3

    model = pd.read_csv(model_path, low_memory=False)
    prism = pd.read_csv(prism_path, low_memory=False) if prism_path else None
    if prism is None:
        print("(note: PRISM file absent — COUNT test runs on lung lines without the PRISM intersection.)")

    m_a = run_one(args.label, matrix_path, model, prism)

    if args.compare:
        cmp_path = Path(args.compare)
        if not cmp_path.exists():
            print(f"\n! --compare file not found: {cmp_path}"); return 1
        m_b = run_one(args.label2, cmp_path, model, prism)
        diff_releases(m_a, args.label or matrix_path.name, m_b, args.label2 or cmp_path.name)

    print("\nRule of thumb: the release whose matrix PASSES both tests is the one used. If the diff shows"
          "\nidentical analyzed-gene calls, the 24Q2/24Q4 discrepancy is cosmetic — fix the manuscript string only.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
