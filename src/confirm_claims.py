"""Succinct reproducibility queries — confirm the specific values the manuscript references
against local DepMap/PRISM files, using the paper's NAMED entities as query keys (so it needs
no re-derivation of genotype). Reads only the columns/rows required, never whole matrices.

Checks:
  Q1  the 7 named RB1-mutant lung lines resolve to ModelIDs and are Lung           (Model.csv)
  Q2  common-essential status of the RB1 target genes                              (Achilles control list)
  Q3  RB1-mut vs WT-lung CRISPR (Chronos) gene-effect differential per target      (CRISPRGeneEffect.csv)
  Q4  doxorubicin pS' across the 7 RB1-mut lines vs the paper's Example 3 (8.15)    (PRISM working file)

Point it at the files wherever they live (they are large / redistribution-restricted, not in git):
    python src/confirm_claims.py --prism PRISM.csv --model Model.csv \
        --crispr CRISPRGeneEffect.csv --common-essential AchillesCommonEssentialControls.csv
Writes data/derived/claim_confirmations.csv. Expectations are read from config/definitions.yaml.
"""
from __future__ import annotations
import argparse
import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
DEFS = yaml.safe_load((ROOT / "config" / "definitions.yaml").read_text())
DERIVED = ROOT / "data" / "derived"
TARGETS = ["AURKA", "AURKB", "PLK1", "CHEK1", "PARP1"]   # RB1 dependency targets discussed in R1.7


def _norm(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def _gene_col(cols, sym):
    for c in cols:
        if _norm(re.split(r"[\s(]", str(c), 1)[0]) == _norm(sym):
            return c
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--prism", required=True, help="PRISM working/SSDRC csv")
    ap.add_argument("--model", required=True, help="DepMap Model.csv")
    ap.add_argument("--crispr", help="CRISPRGeneEffect.csv (Chronos); optional")
    ap.add_argument("--common-essential", dest="ce", help="AchillesCommonEssentialControls.csv; optional")
    args = ap.parse_args(argv)

    import numpy as np
    import pandas as pd
    sys.path.insert(0, str(ROOT.parent / "sprime" / "src"))
    from sprime import s_prime, pooled_ps

    ec = (DEFS.get("expected_cohorts") or {}).get("RB1", {})
    rb1_lines = ec.get("mutant_lines", [])
    results = []

    # Q1 ---------------------------------------------------------------------
    print("Q1  RB1-mutant lung cohort identity (Model.csv)")
    model = pd.read_csv(args.model, usecols=lambda c: c in
                        {"ModelID", "StrippedCellLineName", "OncotreeLineage", "OncotreePrimaryDisease"})
    line2id = {}
    for n in rb1_lines:
        row = model[model.StrippedCellLineName.map(_norm) == _norm(n)]
        if len(row):
            r = row.iloc[0]; line2id[n] = r.ModelID
            print(f"    {n:9s} -> {r.ModelID}  {r.OncotreeLineage} / {r.OncotreePrimaryDisease}")
    all_lung = all(model.loc[model.ModelID == m, "OncotreeLineage"].iloc[0] == "Lung" for m in line2id.values())
    status = "CONFIRMED" if (len(line2id) == len(rb1_lines) and all_lung) else "PARTIAL"
    print(f"    => {len(line2id)}/{len(rb1_lines)} resolved, all Lung={all_lung}  [{status}]")
    results.append(dict(query="Q1_cohort", item="7 named RB1-mut lung lines",
                        observed=f"{len(line2id)}/{len(rb1_lines)} resolved, all_lung={all_lung}",
                        reference="paper names 7 RB1-mut lung lines", status=status))

    # Q2 ---------------------------------------------------------------------
    if args.ce:
        print("Q2  Common-essential status (Achilles control list)")
        cel = pd.read_csv(args.ce)
        ce_genes = {_norm(re.split(r"[\s(]", str(x), 1)[0]) for x in cel.iloc[:, 0]}
        for g in TARGETS:
            hit = _norm(g) in ce_genes
            print(f"    {g:6s}: {'common-essential' if hit else 'selective (not on list)'}")
            results.append(dict(query="Q2_common_essential", item=g,
                                observed="common-essential" if hit else "selective",
                                reference="R1.7 grounding (BioGRID ORCS)", status="OBSERVED"))

    # Q3 ---------------------------------------------------------------------
    if args.crispr:
        print("Q3  RB1-mut vs WT-lung CRISPR (Chronos) gene-effect differential")
        hdr = pd.read_csv(args.crispr, nrows=0)
        idcol = hdr.columns[0]
        cols = {g: _gene_col(hdr.columns, g) for g in TARGETS}
        use = [idcol] + [c for c in cols.values() if c]
        ge = pd.read_csv(args.crispr, usecols=use).set_index(idcol)
        lung_ids = set(model.loc[model.OncotreeLineage == "Lung", "ModelID"])
        prism_ids = set(pd.read_csv(args.prism, usecols=["depmap_id"]).depmap_id.unique())
        lung_prism = lung_ids & prism_ids
        mut_in = set(line2id.values()) & set(ge.index)
        wt_ids = (lung_prism - set(line2id.values())) & set(ge.index)
        print(f"    cohorts: RB1-mut n={len(mut_in)}, WT-lung n={len(wt_ids)} (more negative = more dependent)")
        for g in TARGETS:
            c = cols[g]
            if not c:
                continue
            mm = float(ge.loc[list(mut_in), c].mean()); wm = float(ge.loc[list(wt_ids), c].mean())
            d = mm - wm
            direction = "mut more dependent" if d < 0 else "no SL-direction differential"
            print(f"    {g:6s} mut={mm:7.3f} wt={wm:7.3f} diff={d:+.3f}  {direction}")
            results.append(dict(query="Q3_crispr_differential", item=g,
                                observed=f"mut={mm:.3f} wt={wm:.3f} diff={d:+.3f}",
                                reference="R1.7 RB1-mut-more-dependent claim", status=direction))

    # Q4 ---------------------------------------------------------------------
    print("Q4  Doxorubicin pS' across the 7 RB1-mut lines vs paper Example 3 (8.15)")
    prism = pd.read_csv(args.prism, usecols=lambda c: c in
                        {"depmap_id", "name", "upper_limit", "lower_limit", "ec50"})
    dox = prism[prism.name.astype(str).str.lower().str.contains("doxorubicin", na=False)]
    dox = dox[dox.depmap_id.isin(set(line2id.values()))].dropna(subset=["upper_limit", "lower_limit", "ec50"])
    if len(dox):
        emax = (dox.upper_limit.astype(float) - dox.lower_limit.astype(float)) * 100
        sp = s_prime(emax.values, dox.ec50.astype(float).values, emax_as_percent=True)
        ps = pooled_ps(sp)
        print(f"    pS'(doxorubicin, RB1-mut, n={len(sp)}) = {ps:.2f}  vs paper 8.15")
        results.append(dict(query="Q4_worked_example", item="doxorubicin pS' RB1-mut",
                            observed=f"{ps:.2f} (n={len(sp)})", reference="paper Example 3 = 8.15",
                            status="CONFIRMED" if abs(ps - 8.15) < 0.3 else "DIFFERS"))
    else:
        print("    doxorubicin not present for these lines in this file")
        results.append(dict(query="Q4_worked_example", item="doxorubicin pS' RB1-mut",
                            observed="not in file", reference="paper Example 3 = 8.15", status="UNAVAILABLE"))

    DERIVED.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(results).to_csv(DERIVED / "claim_confirmations.csv", index=False)
    print(f"\nwrote data/derived/claim_confirmations.csv ({len(results)} checks)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
