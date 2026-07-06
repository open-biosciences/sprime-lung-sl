"""Canonical pipeline: data/raw/ + config/definitions.yaml -> data/derived/*.csv

Contract: every number the paper reports comes out of this script into data/derived/, computed under
the pinned definitions - no hand-typed values. Fires the moment the DepMap/PRISM files land in data/raw/.

Order of operations:
  0. preflight - resolve local files; reconcile the 24Q2/24Q4 release from the paper's RB1 cohort
     (reuses src/check_release.py) so the run is stamped with the release it actually used.
  1. genotype  - assign WT/mutant/excluded per gene from the damaging matrix (definitions.genotype),
     restricted to lung & PRISM-19Q4 lines            -> genotype_counts.csv
  2. S'        - per (compound, line) from the PRISM 4PL fit via sprime.s_prime (percent E_max)
  3. pS'/dpS'  - pool by genotype, bootstrap CI, classify Class-A SL (definitions.thresholds)
                                                        -> psprime_by_genotype.csv
  4. metrics   - S' vs AUC vs pEC50 rank correlations   -> metric_correlations.csv
  5. RNAi      - (optional, if DEMETER2 present) RB1-mut-vs-WT dRNAi differentials -> rnai_validation.csv

NOTE: validated end-to-end on synthetic inputs (tests/test_pipeline_smoke.py) but not yet run against the
real matrices (they are gitignored). Column lookups are defensive (real PRISM 19Q4 uses
upper_limit/lower_limit/slope, not the README's idealized names); the first real run may still surface a
column the resolvers don't know - add it to the *_COLS lists below.
"""
from __future__ import annotations
import sys
from pathlib import Path

import yaml

# check_release lives beside this file; reuse its tested resolvers rather than duplicating them.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_release as cr  # noqa: E402

try:
    import numpy as np
    import pandas as pd
    from scipy.stats import spearmanr
    from sprime import s_prime, pooled_ps, delta_ps, bootstrap_delta_ps_ci
    _DEPS_OK = True
except Exception as _e:  # pragma: no cover - environment-dependent
    _DEPS_OK = False
    _IMPORT_ERR = _e

ROOT = Path(__file__).resolve().parents[1]
DEFS = yaml.safe_load((ROOT / "config" / "definitions.yaml").read_text())
RAW = ROOT / "data" / "raw"
DERIVED = ROOT / "data" / "derived"
DEMETER_CANDIDATES = ["demeter2_rnai.csv", "D2_combined_gene_dep_scores.csv"]

# PRISM 19Q4 secondary-screen column aliases (idealized README name -> real-file variants).
PRISM_ID_COLS = ["depmap_id", "DepMap_ID", "ModelID"]
PRISM_COMPOUND_COLS = ["broad_id", "name", "compound", "Compound"]
PRISM_UPPER_COLS = ["upper_limit", "upper_asymptote", "upper", "UpperAsymptote"]
PRISM_LOWER_COLS = ["lower_limit", "lower_asymptote", "lower", "LowerAsymptote"]
PRISM_EC50_COLS = ["ec50", "EC50", "ec50_uM"]
PRISM_AUC_COLS = ["auc", "AUC"]
PRISM_STR_COLS = ["passed_str_profiling", "passed_STR_profiling"]


def _col(df, names, required=True, what=""):
    for n in names:
        if n in df.columns:
            return n
    if required:
        raise SystemExit(f"! column not found for {what or names[0]} (looked for {names}); "
                         f"have: {list(df.columns)[:12]}...  - add the real name to the *_COLS list.")
    return None


def _genotypes(matrix, model, prism):
    """Per analyzed gene: sets of mutant / WT / excluded ModelIDs over lung & PRISM-19Q4."""
    g = DEFS["genotype"]
    mut, wt, exc = g["mutant_value"], g["wild_type_value"], g["excluded_value"]
    lung_mask = cr._lineage_mask_lung(model)
    idc = cr._model_id_column(model)
    lung_ids = set(model.loc[lung_mask, idc]) if lung_mask is not None else set(model[idc])
    pcol = _col(prism, PRISM_ID_COLS, what="PRISM depmap_id")
    universe = (lung_ids & set(prism[pcol].dropna().unique())) & set(matrix.index)

    out = {}
    for gene in g["genes"]:
        col = cr._find_gene_column(matrix.columns, gene, None)
        if col is None:
            print(f"  ! {gene}: no column in matrix - skipped"); continue
        vals = matrix.loc[list(universe), col]
        out[gene] = {
            "mutant": set(vals.index[vals == mut]),
            "wt": set(vals.index[vals == wt]),
            "excluded": set(vals.index[vals == exc]),
        }
    return out, universe


def _write_counts(genos):
    rows = [{"gene": gene, "n_mutant": len(d["mutant"]), "n_wild_type": len(d["wt"]),
             "n_excluded": len(d["excluded"])} for gene, d in genos.items()]
    df = pd.DataFrame(rows)
    # cross-check against the paper's pinned expectation (drift-catcher)
    for gene, ec in (DEFS.get("expected_cohorts") or {}).items():
        row = df[df.gene == gene]
        if not row.empty and "n_mutant" in ec:
            got, exp = int(row.n_mutant.iloc[0]), ec["n_mutant"]
            print(f"  [{gene}] n_mutant = {got} (paper expects {exp}) "
                  f"{'OK' if got == exp else '! MISMATCH - check release/genotype rule'}")
    _emit(df, "genotype_counts.csv")
    return df


def _s_prime_table(prism):
    """Long table of S' per (compound, line), computed on the percent E_max convention."""
    conv = DEFS["s_prime"]
    assert conv["emax_convention"] == "percent", "pipeline computes S' on the percent scale"
    idc = _col(prism, PRISM_ID_COLS, what="PRISM depmap_id")
    ccol = _col(prism, PRISM_COMPOUND_COLS, what="PRISM compound id")
    up, lo = _col(prism, PRISM_UPPER_COLS, what="upper asymptote"), _col(prism, PRISM_LOWER_COLS, what="lower asymptote")
    ec = _col(prism, PRISM_EC50_COLS, what="ec50")
    auc = _col(prism, PRISM_AUC_COLS, required=False)
    strc = _col(prism, PRISM_STR_COLS, required=False)

    df = prism.copy()
    if strc:  # inherit PRISM QC - keep only STR-verified rows (definitions.data_sources.screens)
        df = df[df[strc].astype(str).str.upper().isin(["TRUE", "1", "YES"])]
    df = df.dropna(subset=[up, lo, ec])
    emax_pct = (df[up].astype(float) - df[lo].astype(float)) * 100.0   # fraction asymptotes -> percent, sign preserved
    df = df.assign(
        emax_pct=emax_pct,
        s_prime=s_prime(emax_pct.values, df[ec].astype(float).values,
                        c_ref_uM=conv["c_ref_uM"], emax_as_percent=True),
        pec50=6.0 - np.log10(df[ec].astype(float)),   # ec50 in uM -> pEC50 (M)
    )
    keep = {idc: "line", ccol: "compound", ec: "ec50_uM", "emax_pct": "emax_pct", "s_prime": "s_prime", "pec50": "pec50"}
    if auc:
        keep[auc] = "auc"
    tbl = df[list(keep)].rename(columns=keep)
    # collapse any duplicate (compound, line) pairs (e.g. HTS002/MTS010 overlap) to the mean S'
    return tbl.groupby(["compound", "line"], as_index=False).mean(numeric_only=True)


def _psprime_by_genotype(stbl, genos):
    """For every (gene, compound): pooled pS' per genotype, dpS', bootstrap CI, Class-A flag."""
    th = DEFS["thresholds"]["synthetic_lethal"]
    hs = DEFS["thresholds"]["high_stringency_delta_ps_le"]
    nboot = DEFS["statistics"]["bootstrap_iterations"]
    alpha = 1 - DEFS["statistics"]["ci"]
    s_by_line = {c: g.set_index("line")["s_prime"] for c, g in stbl.groupby("compound")}

    rows = []
    for gene, d in genos.items():
        wt_ids, mut_ids = d["wt"], d["mutant"]
        for compound, s in s_by_line.items():
            s_wt = s.reindex(list(wt_ids)).dropna().values
            s_mut = s.reindex(list(mut_ids)).dropna().values
            if len(s_wt) < 2 or len(s_mut) < 2:      # need both cohorts to contrast
                continue
            ps_wt, ps_mut = pooled_ps(s_wt), pooled_ps(s_mut)
            dps, lo, hi = bootstrap_delta_ps_ci(s_wt, s_mut, n_boot=nboot, alpha=alpha,
                                                seed=DEFS["statistics"].get("seed", 0))
            class_a = (ps_wt > th["wt_ps_gt"]) and (ps_mut > th["mut_ps_gt"]) and (dps <= th["delta_ps_le"])
            rows.append({
                "gene": gene, "compound": compound, "n_wt": len(s_wt), "n_mut": len(s_mut),
                "ps_wt": round(ps_wt, 4), "ps_mut": round(ps_mut, 4), "delta_ps": round(dps, 4),
                "ci_lo": round(lo, 4), "ci_hi": round(hi, 4),
                "class_A": class_a, "high_stringency": bool(dps <= hs),
            })
    df = pd.DataFrame(rows).sort_values(["gene", "delta_ps"]) if rows else pd.DataFrame(
        columns=["gene", "compound", "n_wt", "n_mut", "ps_wt", "ps_mut", "delta_ps", "ci_lo", "ci_hi", "class_A", "high_stringency"])
    _emit(df, "psprime_by_genotype.csv")
    print(f"  Class-A SL candidates (dpS'<={th['delta_ps_le']}, both cohorts inhibited): {int(df.get('class_A', pd.Series(dtype=bool)).sum())}")
    return df


def _metric_correlations(stbl):
    """Rank correlations across all (compound, line) pairs - the C8 benchmark backing."""
    rows = []
    if "auc" in stbl.columns:
        pairs = [("s_prime", "auc"), ("s_prime", "pec50"), ("auc", "pec50")]
    else:
        pairs = [("s_prime", "pec50")]
    for a, b in pairs:
        sub = stbl[[a, b]].dropna()
        if len(sub) >= 3:
            r, p = spearmanr(sub[a], sub[b])
            rows.append({"metric_a": a, "metric_b": b, "spearman_r": round(float(r), 4),
                         "p_value": float(p), "n_pairs": len(sub)})
    _emit(pd.DataFrame(rows), "metric_correlations.csv")


def _rnai_validation(genos):
    """Optional: RB1-mut-vs-WT dRNAi differential per analyzed gene (DEMETER2). Best-effort on schema."""
    path = cr._resolve(None, DEMETER_CANDIDATES, "DEMETER2")
    if path is None:
        print("  (DEMETER2 absent - rnai_validation.csv skipped)"); return
    try:
        d2 = pd.read_csv(path, index_col=0, low_memory=False)  # expect lines x genes (or transpose)
        rows = []
        for gene, d in genos.items():
            col = cr._find_gene_column(d2.columns, gene, None)
            if col is None:
                continue
            wt = d2.loc[d2.index.intersection(d["wt"]), col].dropna()
            mut = d2.loc[d2.index.intersection(d["mutant"]), col].dropna()
            if len(wt) and len(mut):
                rows.append({"gene": gene, "n_wt": len(wt), "n_mut": len(mut),
                             "mean_dep_wt": round(float(wt.mean()), 4),
                             "mean_dep_mut": round(float(mut.mean()), 4),
                             "delta_rnai": round(float(mut.mean() - wt.mean()), 4)})
        _emit(pd.DataFrame(rows), "rnai_validation.csv")
        print("  (dRNAi < 0 => mutant cells more dependent; interpret on a common-essential baseline - see reference table)")
    except Exception as e:
        print(f"  ! DEMETER2 schema not recognized ({e}); rnai_validation.csv skipped - inspect columns and extend.")


def _emit(df, name):
    DERIVED.mkdir(parents=True, exist_ok=True)
    df.to_csv(DERIVED / name, index=False)
    print(f"  wrote data/derived/{name}  ({len(df)} rows)")


def main(argv=None):
    if not _DEPS_OK:
        print("Dependencies not importable in this environment:")
        print(f"  {_IMPORT_ERR}")
        print("  Fix: pip install -r requirements.txt ; pip install -e ../sprime")
        print("  (If it's a numpy/pandas binary mismatch, pin a matched pair, e.g. 'numpy<2' 'pandas>=1.4'.)")
        return 3

    matrix_path = cr._resolve(None, cr.MATRIX_CANDIDATES, "mutations matrix")
    model_path = cr._resolve(None, cr.MODEL_CANDIDATES, "model table")
    prism_path = cr._resolve(None, cr.PRISM_CANDIDATES, "PRISM SSDRC")
    if not (matrix_path and model_path and prism_path):
        print("Pipeline not ready - required inputs missing from data/raw/:")
        print(f"  mutations matrix: {'FOUND' if matrix_path else 'MISSING'} ({cr.MATRIX_CANDIDATES})")
        print(f"  model table:      {'FOUND' if model_path else 'MISSING'} ({cr.MODEL_CANDIDATES})")
        print(f"  PRISM SSDRC:      {'FOUND' if prism_path else 'MISSING'} ({cr.PRISM_CANDIDATES})")
        print("\nSee data/raw/README.md + DATA_AVAILABILITY.md. Downloading the files is the only remaining step.")
        return 2

    print("Loading inputs...")
    matrix = pd.read_csv(matrix_path, index_col=0, low_memory=False)
    model = pd.read_csv(model_path, low_memory=False)
    prism = pd.read_csv(prism_path, low_memory=False)

    print("\n[0] Release reconciliation (from the paper's RB1 cohort):")
    cr.test_named_cohort(matrix, model)
    cr.test_count(matrix, model, prism)

    print("\n[1] Genotype assignment (lung & PRISM-19Q4):")
    genos, universe = _genotypes(matrix, model, prism)
    print(f"  analysis universe: {len(universe)} lung & PRISM lines")
    _write_counts(genos)

    print("\n[2-3] S' / pS' / dpS' by genotype:")
    stbl = _s_prime_table(prism)
    print(f"  computed S' for {len(stbl)} (compound, line) pairs")
    _psprime_by_genotype(stbl, genos)

    print("\n[4] Metric correlations:")
    _metric_correlations(stbl)

    print("\n[5] RNAi validation (optional):")
    _rnai_validation(genos)

    print(f"\nDone. Derived tables in {DERIVED}. These replace the placeholder counts in the response letter (C4/C15).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
