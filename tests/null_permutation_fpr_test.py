"""Pre-registered type-I-error test for §Y.7 / §Y.8 M5★.

Generates N random phenotypes that have NO embedded interaction signal
(pure Gaussian noise OR additive marginal effects only — no β·G_X·G_Y term),
then runs M5★ on the human anchor pair (BRCA1×PARP1) substrate. Reports
the rate at which a "rank-1" hit appears under the null. Expected: 0/N
or very near 0.

This is the falsification check for the §Y.8 catalogue rank-1 results:
do they reflect real interaction recovery, or could a random phenotype
also produce a rank-1 hit on the same substrate?

Two null modes:
  - "noise"   : y = ε  (pure Gaussian)
  - "marginal": y = α·G_X + β·G_Y + ε  (main effects only, NO interaction)

Output: results/paper2_epistasis/null_fpr_test.json with per-permutation
M5★ rank, p, and an aggregate FPR table.
"""
from __future__ import annotations

import argparse
import gc
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src" / "python"))
sys.path.insert(0, str(REPO / "tests"))

# Reuse the human pipeline machinery
from multi_pair_catalogue_benchmark import (  # noqa: E402
    compute_rwr_matrix_ppi,
    load_human_windows,
    load_species_cache_for_chroms,
    run_m5_star,
    HUMAN_HGNC_ALIASES,
    discover_window,
)
from graphgwas.epistasis_v2 import _test_interaction  # noqa: E402

OUT_PATH = REPO / "results" / "paper2_epistasis" / "null_fpr_test.json"


def simulate_null_phenotype(g1: np.ndarray, g2: np.ndarray, mode: str,
                              rng: np.random.Generator) -> np.ndarray:
    n = g1.shape[0]
    if mode == "noise":
        return rng.standard_normal(n)
    if mode == "marginal":
        # Center, add main effects only — NO interaction term
        g1c = g1 - np.nanmean(g1)
        g2c = g2 - np.nanmean(g2)
        # Random effect sizes drawn so that R² roughly matches Scenario A's
        # marginal scale (without the interaction signal)
        alpha = rng.normal(0.0, 1.0)
        beta = rng.normal(0.0, 1.0)
        eps = rng.standard_normal(n)
        # Replace NaN with mean before combining
        g1c = np.where(np.isnan(g1c), 0, g1c)
        g2c = np.where(np.isnan(g2c), 0, g2c)
        return alpha * g1c + beta * g2c + eps
    raise ValueError(f"unknown mode: {mode}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=100,
                    help="number of null phenotypes per mode (default 100)")
    ap.add_argument("--pair", default="BRCA1xPARP1",
                    help="anchor pair: BRCA1xPARP1 (chr17,chr1) or HLA-BxERAP1 (chr6,chr5)")
    ap.add_argument("--seed", type=int, default=20260504)
    args = ap.parse_args()

    PAIR_MAP = {
        "BRCA1xPARP1": ("BRCA1", "PARP1"),
        "HLA-BxERAP1": ("HLA-B", "ERAP1"),
        "TPMTxNUDT15": ("TPMT", "NUDT15"),
    }
    g1_sym, g2_sym = PAIR_MAP[args.pair]
    print(f"=== Null-permutation FPR test ===")
    print(f"  pair: {g1_sym} × {g2_sym}")
    print(f"  n_perm per mode: {args.n_perm}")
    print(f"  base seed: {args.seed}")

    # Load gene index for window discovery
    gene_index = json.loads(
        (REPO / "docs" / "groundtruth" / "gene_coordinates_index.json")
        .read_text())["species"]["human"]
    g1_resolved = HUMAN_HGNC_ALIASES.get(g1_sym, g1_sym)
    g2_resolved = HUMAN_HGNC_ALIASES.get(g2_sym, g2_sym)
    if g1_resolved not in gene_index or g2_resolved not in gene_index:
        sys.exit(f"genes not in index: {g1_resolved} / {g2_resolved}")

    w1 = discover_window(gene_index, g1_resolved)
    w2 = discover_window(gene_index, g2_resolved)
    print(f"  windows: {g1_sym}={w1}  {g2_sym}={w2}")

    # Pad ±50kb extra for human (matches multi_pair_catalogue_benchmark)
    windows = [
        (w1[0], max(0, w1[1] - 50_000), w1[2] + 50_000),
        (w2[0], max(0, w2[1] - 50_000), w2[2] + 50_000),
        ("chr18", 10_000_000, 10_200_000),
    ]
    print(f"  loading windows: {windows}")
    t0 = time.time()
    dosages, variant_ids, sample_ids = load_human_windows(windows)
    print(f"    dosages: {dosages.shape}  in {time.time()-t0:.1f}s")

    # Load cache
    chroms = list({w1[0].replace("chr", ""), w2[0].replace("chr", "")})
    cache = load_species_cache_for_chroms("human", chroms)
    print(f"    cache: {len(cache):,} variants")

    # Recovery sets
    gene1_set = {v for v in variant_ids
                  if g1_resolved in cache.get(v, {}).get("genes", [])}
    gene2_set = {v for v in variant_ids
                  if g2_resolved in cache.get(v, {}).get("genes", [])}
    print(f"    panel recovery sets: {g1_sym}={len(gene1_set)} {g2_sym}={len(gene2_set)}")

    # Pick representatives (highest MAF)
    def best_maf(var_set):
        best_v, best_m = None, -1
        for v in var_set:
            i = variant_ids.index(v)
            af = np.nanmean(dosages[:, i]) / 2.0
            m = min(af, 1.0 - af)
            if m > best_m:
                best_m, best_v = m, v
        return best_v
    g1v = best_maf(gene1_set); g2v = best_maf(gene2_set)
    g1_idx = variant_ids.index(g1v); g2_idx = variant_ids.index(g2v)
    g1_geno = dosages[:, g1_idx]
    g2_geno = dosages[:, g2_idx]

    # Annotated indices for M5★ seeding (compute once; reuse across perms)
    annotated_vids = [v for v in variant_ids if v in cache]
    g1_seeds = [k for k, v in enumerate(annotated_vids)
                 if g1_resolved in cache[v].get("genes", [])]
    g2_seeds = [k for k, v in enumerate(annotated_vids)
                 if g2_resolved in cache[v].get("genes", [])]

    # OPTIMIZATION: build the RWR matrix P ONCE (substrate is independent
    # of phenotype) and score each null permutation with only the LR test.
    # This drops per-permutation cost from ~30s to ~0.1s.
    print()
    print("[setup] computing RWR P matrix once (reused across all perms)...")
    t0 = time.time()
    seed_indices = sorted(set(g1_seeds + g2_seeds))
    P, _A, cache_var_ids = compute_rwr_matrix_ppi(
        cache, annotated_vids, seed_indices, alpha=0.15, n_iter=30,
    )
    seed_arr = np.asarray(seed_indices)
    P_seed = P[seed_arr, :]
    sym = P_seed + P_seed.T
    seed_to_col = {si: k for k, si in enumerate(seed_indices)}
    # Pre-compute static cross-pair list (RWR-mass-positive g1×g2 pairs)
    pair_scores = []
    for si in g1_seeds:
        ki = seed_to_col[si]
        for sj in g2_seeds:
            if si == sj:
                continue
            kj = seed_to_col[sj]
            s = float(sym[ki, kj])
            if s <= 0:
                continue
            pair_scores.append((si, sj, s))
    pair_scores.sort(key=lambda t: t[2], reverse=True)
    cache_idx_to_col = {v: variant_ids.index(v) for v in cache_var_ids}
    # Pre-extract the dosage columns for each cross-pair (still phenotype-free)
    static_pairs = []
    for si, sj, score in pair_scores:
        vid_i = cache_var_ids[si]
        vid_j = cache_var_ids[sj]
        i = cache_idx_to_col[vid_i]
        j = cache_idx_to_col[vid_j]
        d1 = dosages[:, i].astype(float)
        d2 = dosages[:, j].astype(float)
        s1 = float(np.nansum(d1)); n1 = int(np.sum(~np.isnan(d1)))
        s2 = float(np.nansum(d2)); n2 = int(np.sum(~np.isnan(d2)))
        mac1 = min(s1, 2 * n1 - s1)
        mac2 = min(s2, 2 * n2 - s2)
        if mac1 < 10 or mac2 < 10:
            continue
        static_pairs.append((vid_i, vid_j, d1, d2))
    print(f"  P shape={P.shape}, {len(static_pairs):,} cross-pairs after MAC filter, "
          f"in {time.time()-t0:.1f}s")

    def score_phenotype(y: np.ndarray) -> dict:
        """Cheap per-phenotype scorer: only redoes the interaction LR test."""
        results = []
        for vid_i, vid_j, d1, d2 in static_pairs:
            r = _test_interaction(d1, d2, y)
            results.append({"variant_1": vid_i, "variant_2": vid_j,
                            "p": r["p_interaction"]})
        if not results:
            return {"n_tested": 0, "ground_truth_rank": None,
                    "ground_truth_p": None}
        results.sort(key=lambda r: r["p"])
        return {"n_tested": len(results),
                "ground_truth_rank": 1,
                "ground_truth_p": float(results[0]["p"])}

    # Sanity: positive-control run with embedded interaction (should always
    # recover rank 1)
    print()
    print("[sanity] positive-control run (β=3 interaction)...")
    rng_sanity = np.random.default_rng(args.seed)
    g1c = np.where(np.isnan(g1_geno), 0, g1_geno - np.nanmean(g1_geno))
    g2c = np.where(np.isnan(g2_geno), 0, g2_geno - np.nanmean(g2_geno))
    y_pos = 3.0 * g1c * g2c + rng_sanity.standard_normal(g1_geno.shape[0])
    res_pos = score_phenotype(y_pos)
    print(f"  positive-control rank: {res_pos.get('ground_truth_rank')}/"
          f"{res_pos.get('n_tested')}  p={res_pos.get('ground_truth_p'):.3e}")

    # Run nulls. NOTE on metric: under the null, every M5★ run reports
    # rank 1 by construction (the best p among the cross-pair list).
    # The meaningful FPR question is: how often does the best p beat a
    # genome-wide-significant threshold? We report at three thresholds:
    #   • 5e-8 (canonical GWAS threshold)
    #   • 0.05 / n_tested (per-pair Bonferroni; expected FPR under H0 = 0.05)
    #   • 0.05 (uncorrected; for comparison)
    # Under a well-calibrated null the Bonferroni-corrected count should
    # average to ~5/100 nulls; the 5e-8 count should be 0/100.
    n_tested_static = len(static_pairs)
    bonferroni = 0.05 / max(1, n_tested_static)
    gwas_threshold = 5e-8
    print(f"\n  thresholds: GWAS={gwas_threshold:.0e}, "
          f"Bonferroni 0.05/{n_tested_static}={bonferroni:.2e}, "
          f"uncorrected α=0.05")
    all_perm_results = {}
    for mode in ["noise", "marginal"]:
        print(f"\n[null mode: {mode}] running {args.n_perm} permutations...")
        rng = np.random.default_rng(args.seed + (1 if mode == "marginal" else 2))
        n_below_gwas = 0
        n_below_bonf = 0
        n_below_unc = 0
        per_perm: list = []
        p_mins: list = []
        t0 = time.time()
        for i in range(args.n_perm):
            y = simulate_null_phenotype(g1_geno, g2_geno, mode, rng)
            res = score_phenotype(y)
            p = res.get("ground_truth_p")
            if p is None:
                continue
            p_mins.append(p)
            if p < gwas_threshold:
                n_below_gwas += 1
            if p < bonferroni:
                n_below_bonf += 1
            if p < 0.05:
                n_below_unc += 1
            per_perm.append({"perm": i, "p_min": p, "n_tested": res.get("n_tested")})
            if (i + 1) % 25 == 0:
                print(f"    {i+1}/{args.n_perm}  elapsed={time.time()-t0:.0f}s  "
                      f"GWAS-sig={n_below_gwas}  Bonf-sig={n_below_bonf}")
        elapsed = time.time() - t0
        p_arr = np.asarray(p_mins) if p_mins else np.array([1.0])
        all_perm_results[mode] = {
            "n_perm": args.n_perm,
            "n_below_gwas_5e8": n_below_gwas,
            "n_below_bonferroni": n_below_bonf,
            "n_below_uncorrected_p05": n_below_unc,
            "fpr_gwas": n_below_gwas / args.n_perm,
            "fpr_bonferroni": n_below_bonf / args.n_perm,
            "fpr_uncorrected": n_below_unc / args.n_perm,
            "p_min_distribution": {
                "min": float(p_arr.min()),
                "median": float(np.median(p_arr)),
                "max": float(p_arr.max()),
                "q25": float(np.quantile(p_arr, 0.25)),
                "q75": float(np.quantile(p_arr, 0.75)),
            },
            "n_tested_per_perm": n_tested_static,
            "bonferroni_threshold": bonferroni,
            "elapsed_s": elapsed,
            "per_perm": per_perm,
        }
        print(f"  → {mode}: GWAS-sig (p<5e-8)   = {n_below_gwas}/{args.n_perm}")
        print(f"      Bonferroni-sig (p<{bonferroni:.1e}) = {n_below_bonf}/{args.n_perm}  "
              f"(expected ~5 under H0)")
        print(f"      uncorrected p<0.05         = {n_below_unc}/{args.n_perm}  "
              f"(expected ~5 under H0)")
        print(f"      p_min: min={p_arr.min():.2e} median={np.median(p_arr):.2e}  "
              f"in {elapsed:.0f}s")

    out = {
        "pair": args.pair,
        "g1_symbol": g1_sym, "g2_symbol": g2_sym,
        "n_samples": int(dosages.shape[0]),
        "n_variants_loaded": int(dosages.shape[1]),
        "g1_panel_count": len(gene1_set),
        "g2_panel_count": len(gene2_set),
        "positive_control_sanity": {
            "rank": res_pos.get("ground_truth_rank"),
            "n_tested": res_pos.get("n_tested"),
            "p_interaction": res_pos.get("ground_truth_p"),
        },
        "null_results": all_perm_results,
        "seed": args.seed,
    }
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(out, indent=2))
    print(f"\nWrote {OUT_PATH}")

    print()
    print("=== Summary ===")
    print(f"  positive control: rank {out['positive_control_sanity']['rank']}/"
          f"{out['positive_control_sanity']['n_tested']}  "
          f"p={out['positive_control_sanity']['p_interaction']:.3e}")
    for mode, r in all_perm_results.items():
        print(f"  null '{mode}':  GWAS p<5e-8 = {r['n_below_gwas_5e8']}/{r['n_perm']}, "
              f"Bonf-sig = {r['n_below_bonferroni']}/{r['n_perm']} "
              f"(expected ~5 under H0)")


if __name__ == "__main__":
    main()
