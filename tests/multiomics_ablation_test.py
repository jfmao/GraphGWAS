"""Multi-omics subsampling test — does more annotation coverage help?

For a single GWAS locus with a known causal variant, we run HBP at
5 levels of annotation coverage (10%, 25%, 50%, 75%, 100% of the
multi-omics edges randomly kept) × 20 bootstrap seeds. We quantify:
  - credible-set (CS) size
  - posterior inclusion probability (PIP) on the true causal
  - rank of the true causal by PIP

Expected pattern (per the paper's relational-prior claim):
  - CS size should decrease with increasing edge coverage
  - Causal PIP should rise with increasing edge coverage
  - Baseline (flat-prior L1 without graph) is flat across coverage

The script runs species-by-species on any config that provides:
  - a GWAS sumstats TSV (+ per-variant LD)
  - a causal_variant (variant_id known ground truth)
  - a full multi-omics graph cache
  - a locus window

This is a reusable scaffold — specific species configs are passed in.
"""
from __future__ import annotations

import json
import random
import subprocess
import sys
from copy import deepcopy
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

sys.path.insert(0, "/mnt/data/GraphGWAS/src/python")

from graphgwas.finemapping_v2 import (
    hbp_finemap_from_sumstats,
    l1_finemap_from_sumstats,
)


# --------------------------------------------------------------------
# Graph-cache subsampling helpers
# --------------------------------------------------------------------

def subsample_cache(full_cache: dict, fraction: float, seed: int) -> dict:
    """Retain a random `fraction` of (variant, edge) entries per variant.

    Strategy:
      - With probability (1-fraction) drop each gene/pathway/ppi entry
      - Variants with no surviving edges are pruned from the cache
    """
    rng = random.Random(seed)
    out: dict[str, dict] = {}
    for vid, info in full_cache.items():
        new_info = {}
        for kind in ("genes", "pathways", "ppi"):
            items = info.get(kind, [])
            kept = [x for x in items if rng.random() < fraction]
            if kept:
                new_info[kind] = kept
        if new_info:
            out[vid] = new_info
    return out


# --------------------------------------------------------------------
# Per-species runner
# --------------------------------------------------------------------

def run_ablation(
    species: str,
    locus_name: str,
    variants: list[dict],
    z: np.ndarray,
    R_sq: np.ndarray,
    full_cache: dict,
    causal_variant_id: str,
    fractions: list[float] = (0.1, 0.25, 0.5, 0.75, 1.0),
    n_seeds: int = 20,
    output_dir: Path | None = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """Run HBP at each coverage fraction with bootstrap seeds.
    Returns long-format DataFrame with columns:
      species, locus, fraction, seed, method, cs_size, causal_pip, causal_rank.
    """
    rows = []
    for frac in fractions:
        if frac == 1.0:
            seeds = [0]   # full cache is deterministic
        else:
            seeds = list(range(n_seeds))
        for seed in seeds:
            sub = subsample_cache(full_cache, frac, seed) if frac < 1.0 else full_cache
            # Run HBP
            cand = hbp_finemap_from_sumstats(
                variants, z, R_sq, graph_cache=sub,
                r2_smooth=0.3, credible_set_coverage=0.95,
                chr_name=variants[0]["chr"] if variants else "",
            )
            # Compute metrics
            cs = [c for c in cand if c.in_credible_set]
            cs_size = len(cs)
            causal_match = [c for c in cand if c.variant_id == causal_variant_id]
            causal_pip = float(causal_match[0].pip) if causal_match else 0.0
            ranked = sorted(cand, key=lambda c: -c.pip)
            causal_rank = next((i+1 for i, c in enumerate(ranked)
                                if c.variant_id == causal_variant_id), None)
            rows.append({
                "species": species, "locus": locus_name,
                "fraction": frac, "seed": seed,
                "method": "HBP",
                "cs_size": cs_size, "causal_pip": causal_pip,
                "causal_rank": causal_rank,
                "cache_size": len(sub),
            })
            if verbose and seed == 0:
                print(f"    [{frac*100:>5.0f}%] cache={len(sub):>7,}  "
                      f"CS={cs_size:>5}  causal_pip={causal_pip:.4f}  "
                      f"rank={causal_rank}", flush=True)

        # Also run L1 BASELINE (no graph prior) once per fraction for
        # comparison — it should be flat across fractions.
        l1_cand = l1_finemap_from_sumstats(
            variants, z, R_sq, z_func=None, alpha=1.0,
            r2_smooth=0.3, credible_set_coverage=0.95,
            chr_name=variants[0]["chr"] if variants else "",
        )
        cs_l1 = [c for c in l1_cand if c.in_credible_set]
        cm_l1 = [c for c in l1_cand if c.variant_id == causal_variant_id]
        rank_l1 = next((i+1 for i, c in enumerate(sorted(l1_cand, key=lambda c: -c.pip))
                       if c.variant_id == causal_variant_id), None)
        rows.append({
            "species": species, "locus": locus_name,
            "fraction": frac, "seed": -1,
            "method": "L1_baseline",
            "cs_size": len(cs_l1),
            "causal_pip": float(cm_l1[0].pip) if cm_l1 else 0.0,
            "causal_rank": rank_l1, "cache_size": 0,
        })

    df = pd.DataFrame(rows)
    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_dir / f"ablation_{species}_{locus_name}.tsv",
                  sep="\t", index=False)
    return df
