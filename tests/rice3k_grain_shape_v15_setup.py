"""v0.1.5 setup: build λ_GC-deflated .ma files for SBayesRC.

For each trait, multiplies SE by sqrt(λ_GC) so that the implied z = BETA/SE'
becomes z' = z / sqrt(λ_GC). This is equivalent to dividing chi-square by
λ_GC (the standard genomic-control correction), and is applied to all six
fine-mapping methods in v0.1.5.

Inputs:
  data/rice_3k/results/grain_gwas/grain_lambda_gc.tsv
  data/rice_3k/sbayesrc_rice/sumstats/grain_<TRAIT>.ma  (v0.1.4)

Outputs:
  data/rice_3k/sbayesrc_rice/sumstats_v15/grain_<TRAIT>.ma  (deflated)
"""
from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

DATA = Path("/mnt/data/GraphGWAS/data/rice_3k")
SBAYES = DATA / "sbayesrc_rice"
SUMS_OLD = SBAYES / "sumstats"
SUMS_NEW = SBAYES / "sumstats_v15"
SUMS_NEW.mkdir(parents=True, exist_ok=True)

LAMBDA_TSV = DATA / "results" / "grain_gwas" / "grain_lambda_gc.tsv"

lam = pd.read_csv(LAMBDA_TSV, sep="\t")
lam_map = dict(zip(lam["trait"], lam["lambda_gc"]))
print("λ_GC by trait:")
for t, v in lam_map.items():
    print(f"  {t}: {v:.4f}  →  sqrt={math.sqrt(v):.4f}")

for trait in ["TGW", "GL", "GW", "RLW"]:
    inp = SUMS_OLD / f"grain_{trait}.ma"
    out = SUMS_NEW / f"grain_{trait}.ma"
    df = pd.read_csv(inp, sep="\t")
    factor = math.sqrt(lam_map[trait])
    df["se"] = df["se"] * factor
    # Recompute p-value from deflated z using χ²₁ tail
    from scipy import stats as sst
    z_def = df["b"] / df["se"]
    df["p"] = sst.norm.sf(z_def.abs()) * 2
    df.to_csv(out, sep="\t", index=False)
    print(f"  wrote {out.name}  N_snps={len(df):,}  factor={factor:.4f}")
