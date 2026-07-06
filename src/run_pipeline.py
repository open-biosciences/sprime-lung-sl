"""Canonical pipeline: data/raw/ + config/definitions.yaml -> data/derived/*.csv

STUB. Drop the real analysis in; the contract is that every number the paper reports comes
out of this script into data/derived/, computed under the pinned definitions (no hand-typed values).
"""
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[1]
DEFS = yaml.safe_load((ROOT / "config" / "definitions.yaml").read_text())


def main():
    assert DEFS["s_prime"]["emax_convention"] in {"percent", "fraction"}, \
        "Lock the E_max convention in definitions.yaml before running (see concerns.csv C1)."
    # 1. load data/raw/*  2. assign genotype (per definitions.genotype)
    # 3. compute S′ (sprime.s_prime under definitions.s_prime)  4. pool pS′/ΔpS′ + bootstrap CIs
    # 5. classify (definitions.thresholds)  6. correlations, OAT sensitivity, RNAi, concordance
    # 7. write data/derived/*.csv
    raise NotImplementedError("Wire in raw data + sprime module, then emit data/derived/*.csv")


if __name__ == "__main__":
    main()
