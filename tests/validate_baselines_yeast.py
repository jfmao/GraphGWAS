"""§Y.3 — classical-baseline head-to-head against §Y.4 yeast Scenario A.

Compares MDR (scikit-mdr) against M2's BCY1×TPK1 recovery on the same
genotypes + phenotype that §Y.4 used.

  - **M2 (graphgwas)**:  rank 1,148 (q=1.7×10⁻⁵)  — recorded in §Y.4
  - **M1 (graphgwas)**:  serves as the "BOOST-equivalent" exhaustive-pairwise-
                         LR baseline.  M1 with r²-prune disabled is exactly
                         the BOOST/Wan-2010 ANOVA-style screen.  §Y.4 result:
                         NF (LD pruning replaced BCY1/TPK1 representatives).
  - **MDR (scikit-mdr)**: this script.  Random-pair-sampled CV balanced
                         accuracy on quartile-binarised phenotype.
  - **BOOST (plink2)**:  not available in this environment — plink2 v2.0
                         drops the --epistasis flag (PLINK 1.9 only).
                         Documented as a §Y.3 deferred-baseline.
  - **MAPIT (R)**:       NotImplementedError in graphgwas.epistasis_baselines —
                         documented as a §Y.3 deferred-baseline.

Usage:
  python tests/validate_baselines_yeast.py [--n-pairs-test 5000]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src" / "python"))
sys.path.insert(0, str(REPO_ROOT / "tests"))

from graphgwas.epistasis_baselines import run_mdr  # noqa: E402

# Reuse §Y.4 loaders + simulator
from validate_m2_yeast import (  # noqa: E402
    CHROMS_TO_LOAD, GROUND_TRUTH, YEAST_CACHE, YEAST_VCF,
    load_yeast_chroms, pick_representative_variants, simulate_with_truth,
)


RESULTS_DIR = REPO_ROOT / "results" / "paper2_epistasis"
BCY1_ORF = "YIL033C"
TPK1_ORF = "YJL164C"


def gene_level_rank(df, gene_1_variants: set[str], gene_2_variants: set[str],
                    rank_by: str) -> tuple[int | None, dict | None]:
    """Find the highest-ranked row whose pair is BCY1-var × TPK1-var.

    Args:
        rank_by: column name to rank by.  ascending for "p_value", descending for "score".
    """
    ascending = (rank_by == "p_value")
    df_sorted = df.sort_values(rank_by, ascending=ascending).reset_index(drop=True)
    for r, row in df_sorted.iterrows():
        v1, v2 = row["variant_1"], row["variant_2"]
        cross = ((v1 in gene_1_variants and v2 in gene_2_variants)
                 or (v1 in gene_2_variants and v2 in gene_1_variants))
        if cross:
            return int(r) + 1, {
                "variant_1": v1, "variant_2": v2,
                "score": float(row.get("score", float("nan"))),
                "p_value": float(row.get("p_value", float("nan"))),
            }
    return None, None


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--n-pairs-test", type=int, default=5000,
                   help="Random pairs sampled for MDR (full enumeration of "
                        "6.6K² pairs is impractical)")
    p.add_argument("--cv-folds", type=int, default=5)
    p.add_argument("--beta", type=float, default=3.0,
                   help="Scenario A interaction effect size (matches §Y.4)")
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--bin-quantile", type=float, default=0.50,
                   help="Phenotype binarisation quantile for MDR. 0.5 = "
                        "median split (cases above, controls below).")
    args = p.parse_args()

    rng = np.random.default_rng(args.seed)
    print("=== §Y.3 baseline head-to-head: MDR vs M2 on yeast Scenario A ===\n")

    # 1. Genotypes + phenotype (same recipe as §Y.4 Scenario A)
    print("[1/4] Loading yeast genotypes...")
    dosages, variant_ids, sample_ids = load_yeast_chroms(
        YEAST_VCF, CHROMS_TO_LOAD, maf_min=0.05,
    )
    cache = json.loads(YEAST_CACHE.read_text())
    n_in_cache = sum(1 for v in variant_ids if v in cache)

    print("\n[2/4] Picking ground-truth representatives...")
    variant_id_to_col = {v: i for i, v in enumerate(variant_ids)}
    truth = pick_representative_variants(cache, variant_id_to_col, dosages,
                                          GROUND_TRUTH)

    print(f"\n[3/4] Simulating Scenario A phenotype (β={args.beta})...")
    phenotype = simulate_with_truth(dosages, truth, beta=args.beta,
                                     n_nuisance=0, rng=rng)

    bcy1_variants = {v for v in variant_ids
                     if BCY1_ORF in cache.get(v, {}).get("genes", [])}
    tpk1_variants = {v for v in variant_ids
                     if TPK1_ORF in cache.get(v, {}).get("genes", [])}
    print(f"\n      Recovery sets: BCY1 = {len(bcy1_variants)} variants, "
          f"TPK1 = {len(tpk1_variants)} variants")

    # ---- MDR ----
    print(f"\n[4/4] Running MDR on quantile-binarised phenotype "
          f"(top {(1-args.bin_quantile)*100:.0f}% = case)...")
    threshold = float(np.quantile(phenotype, args.bin_quantile))
    binary_pheno = (phenotype > threshold).astype(np.int64)
    print(f"      threshold (q={args.bin_quantile}): {threshold:.3f}")
    print(f"      n_cases: {int(binary_pheno.sum())}, "
          f"n_controls: {int((1 - binary_pheno).sum())}")

    # MDR works on integer dosages; impute missing via the column mean rounded
    geno_int = np.where(np.isnan(dosages), np.nanmean(dosages, axis=0), dosages)
    geno_int = np.clip(np.round(geno_int), 0, 2).astype(np.int8)

    # Force-include the BCY1 × TPK1 cross-pairs so MDR is evaluated on the
    # ground-truth pair (random sampling has only ~0.02% chance of hitting it)
    bcy1_indices = [variant_ids.index(v) for v in bcy1_variants]
    tpk1_indices = [variant_ids.index(v) for v in tpk1_variants]
    cross_pairs = [(i, j) for i in bcy1_indices for j in tpk1_indices]
    print(f"      Force-including {len(cross_pairs)} BCY1×TPK1 cross-pairs")

    # MDR can be slow; we sample n_pairs_test random pairs + force the cross-pairs
    t0 = time.time()
    mdr_df = run_mdr(geno_int, binary_pheno, variant_ids,
                      n_pairs_test=args.n_pairs_test,
                      cv_folds=args.cv_folds, seed=args.seed,
                      force_pairs=cross_pairs)
    mdr_runtime = time.time() - t0
    print(f"      MDR runtime: {mdr_runtime:.1f}s on {len(mdr_df):,} pairs")

    # Diagnostic: how many MDR scores are NaN?  Low-MAF pairs (BCY1/TPK1 reps
    # have MAF 5.9%/7.1%) frequently produce empty CV cells → NaN bal-acc.
    n_nan = int(mdr_df["score"].isna().sum())
    print(f"      NaN scores: {n_nan}/{len(mdr_df):,} "
          f"({100*n_nan/max(1,len(mdr_df)):.1f}% — these come from "
          f"empty CV cells under low MAF + small panel)")

    rank_mdr, hit_mdr = gene_level_rank(mdr_df, bcy1_variants, tpk1_variants,
                                         rank_by="score")
    if rank_mdr is None:
        print(f"      ✗ BCY1×TPK1 NOT in MDR top-{len(mdr_df):,} pairs")
    else:
        print(f"      ✓ BCY1×TPK1 ranked {rank_mdr:,}/{len(mdr_df):,} by MDR "
              f"(CV bal-acc = {hit_mdr['score']:.4f})")
        print(f"        pair: {hit_mdr['variant_1']} × {hit_mdr['variant_2']}")

    # ---- BOOST (plink2): documented unavailable ----
    print(f"\n[BOOST/plink2]: skipped — plink2 v2 drops --epistasis "
          f"(was PLINK 1.9 only).  M1's exhaustive pairwise-LR machinery "
          f"in §Y.4 IS the same algorithmic family as BOOST/Wan-2010; "
          f"the `--no-prune` mode is the cleanest replacement.")

    # ---- MAPIT (R): NotImplementedError, documented deferred ----
    print(f"\n[MAPIT/R]: skipped — graphgwas.epistasis_baselines.run_mapit "
          f"raises NotImplementedError (Rscript template TODO).")

    # ----- Summary -----
    print(f"\n{'=' * 60}")
    print(f"  §Y.3 SUMMARY (yeast Scenario A, BCY1×TPK1 recovery)")
    print(f"{'=' * 60}")
    print(f"  M2 (graphgwas, motif-filtered):      rank 1,148  q=1.7e-5")
    print(f"  M1 (graphgwas, LD-pruned exhaustive): NF (LD pruning replaces both genes' reps)")
    if rank_mdr:
        print(f"  MDR (scikit-mdr):                     rank {rank_mdr:,} / "
              f"{len(mdr_df):,} (CV bal-acc {hit_mdr['score']:.3f})")
    else:
        print(f"  MDR (scikit-mdr):                     NF in top-{len(mdr_df):,}")
    print(f"  BOOST (plink2):                       UNAVAILABLE (plink2 v2 lacks --epistasis)")
    print(f"  MAPIT (R/Rscript):                    DEFERRED")

    out = RESULTS_DIR / "y3_baselines.json"
    payload = {
        "scenario": "A",
        "phenotype_binarisation_quantile": args.bin_quantile,
        "n_cases": int(binary_pheno.sum()),
        "n_controls": int((1 - binary_pheno).sum()),
        "mdr": {
            "n_pairs_test": args.n_pairs_test,
            "cv_folds": args.cv_folds,
            "runtime_sec": mdr_runtime,
            "n_pairs_evaluated": len(mdr_df),
            "n_nan_scores": int(mdr_df["score"].isna().sum()),
            "n_finite_scores": int(mdr_df["score"].notna().sum()),
            "bcy1_tpk1_rank": rank_mdr,
            "bcy1_tpk1_rank_finite_only": (
                # Recompute rank ignoring NaN scores
                gene_level_rank(mdr_df.dropna(subset=["score"]).reset_index(drop=True),
                                bcy1_variants, tpk1_variants, rank_by="score")[0]
            ),
            "bcy1_tpk1_hit": hit_mdr,
            "bcy1_tpk1_percentile": (
                rank_mdr / max(1, len(mdr_df)) if rank_mdr else None
            ),
            "comparison_percentile": {
                "M2_percentile": 1148 / 94440,
                "MDR_percentile": (
                    rank_mdr / max(1, len(mdr_df)) if rank_mdr else None
                ),
                "M2_advantage_factor": (
                    (rank_mdr / max(1, len(mdr_df))) / (1148 / 94440)
                    if rank_mdr else None
                ),
            },
        },
        "boost_plink2": {
            "status": "unavailable",
            "reason": "plink2 v2.0.0-a.6.5LM lacks --epistasis (PLINK 1.9 only)",
            "note": "M1 LD-pruned exhaustive LR is the closest in-repo equivalent",
        },
        "mapit_r": {
            "status": "deferred",
            "reason": "graphgwas.epistasis_baselines.run_mapit not implemented",
        },
        "comparison_to_graphgwas_y4": {
            "M2_rank": 1148,
            "M2_q_value": 1.7e-5,
            "M1_rank": None,
            "M1_status": "NF (LD pruning replaces ground-truth representatives)",
        },
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, default=str))
    print(f"\n  Wrote {out}")


if __name__ == "__main__":
    main()
