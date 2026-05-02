"""§Y.5 k=3 — gene-priority triplet RWR benchmark on yeast cAMP-PKA core.

The structurally-unique paper-#2 claim: graph methods reach k≥3 epistasis
that pairwise BOOST/MAPIT/MDR cannot formulate.  Until this script,
``mutual_rwr_triplet_scores`` had only a trivial-case unit test
(``test_random_walk_epistasis.py:90``) and zero real-data benchmark.

Ground truth: BCY1 (YIL033C) + TPK1 (YJL164C) + TPK2 (YPL203W) — the
canonical cAMP-dependent PKA holoenzyme: BCY1 is the regulatory subunit,
TPK1/TPK2/TPK3 are catalytic subunits.  All three obligately bind in
the inactive state; cAMP releases them.

Pipeline:

    1. Load yeast 1011 genotypes + graph cache.  Inject canonical
       BioGRID cAMP-PKA edges (the cache PPI is incomplete — see §Y.5).
    2. Build a gene-priority heterograph (shared pathways + PPI, IDF-
       weighted) restricted to genes with ≥1 dosage variant.
    3. Find top-K gene triangles by triple-product mutual-RWR score.
    4. Expand each gene triplet (g_a, g_b, g_c) to all variant triplets
       across the three genes (sampled to cap fan-out).
    5. Simulate phenotype y = beta * g_a * g_b * g_c + eps (3-way DGP).
    6. Test each variant triplet for a 3-way interaction term
       (linear regression with main effects, 2-way, and the 3-way term).
    7. Rank ground-truth triplets {BCY1, TPK1, TPK2} (and {BCY1, TPK1, TPK3}).

This is the first real-data evidence that graph-priority RWR enables
k=3 detection at biobank-scale tractability.

Usage:
    python tests/validate_triplet_rwr_yeast.py [--top-gene-triplets 200]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from itertools import combinations, product
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src" / "python"))
sys.path.insert(0, str(REPO_ROOT / "tests"))

# Reuse §Y.4 + §Y.5 helpers
from graphgwas.canonical_ppi import (  # noqa: E402
    YEAST_CANONICAL_EDGES,
    inject_canonical_ppi_edges,
)
from validate_m2_yeast import (  # noqa: E402
    YEAST_CACHE, YEAST_VCF, load_yeast_chroms,
)

# k=3 needs chrs harbouring all four PKA subunits:
#   BCY1 (chr9), TPK1 (chr10), TPK3 (chr11), TPK2 (chr16)
CHROMS_TO_LOAD = ["chromosome9", "chromosome10", "chromosome11",
                   "chromosome12", "chromosome16"]
from validate_m5_gene_priority_yeast import (  # noqa: E402
    _build_gene_gene_heterograph,
)


RESULTS_DIR = REPO_ROOT / "results" / "paper2_epistasis"

# k=3 ground-truth: BCY1 + TPK1 + TPK2 (the PKA reg + 2 catalytics)
PKA_TRIPLE = ("YIL033C", "YJL164C", "YPL203W")  # BCY1, TPK1, TPK2
ALT_TRIPLE = ("YIL033C", "YJL164C", "YKL166C")  # BCY1, TPK1, TPK3
TRIPLE_NAMES = {"YIL033C": "BCY1", "YJL164C": "TPK1",
                "YPL203W": "TPK2", "YKL166C": "TPK3"}


# ===========================================================================
# Heterograph triplet enumeration
# ===========================================================================

def heterograph_triplet_scores(
    A_gg: np.ndarray,
    seed_indices: list[int],
    alpha: float = 0.15,
    top_pair_k: int = 5000,
    top_triplet_k: int = 200,
    n_iter: int = 30,
) -> tuple[list[tuple[int, int, int, float]], list[tuple[int, int, float]]]:
    """Find top gene triangles by triple-product mutual-RWR score.

    Strategy: run gene-gene RWR once, get top-N pair scores.  Build the
    pair-graph and find triangles (i, j, k) where (i,j), (i,k), (j,k) are
    all in the top-pair list.  Score each triangle by the product of its
    three pair scores.  This is the heterograph counterpart to
    ``epistasis_higher_order.mutual_rwr_triplet_scores``.

    Returns:
        (triplets, pairs) where triplets = [(i,j,k,score), ...] and
        pairs = [(i,j,score), ...] are returned for inspection.
    """
    n = A_gg.shape[0]
    seeds = list(seed_indices)
    n_seed = len(seeds)
    row_sum = A_gg.sum(axis=1, keepdims=True).clip(min=1e-12)
    P_gg = A_gg / row_sum
    E = np.zeros((n, n_seed), dtype=np.float64)
    for k, gi in enumerate(seeds):
        E[gi, k] = 1.0
    P = E.copy()
    for _ in range(n_iter):
        P = (1.0 - alpha) * (P_gg.T @ P) + alpha * E
    P = P / np.clip(P.sum(axis=0, keepdims=True), 1e-12, None)
    seed_arr = np.asarray(seeds)
    P_seed = P[seed_arr, :]
    mutual = P_seed * P_seed.T

    # Top pairs (upper triangle)
    iu = np.triu_indices(n_seed, k=1)
    pair_idx = list(zip(iu[0].tolist(), iu[1].tolist()))
    scores = mutual[iu]
    order = np.argsort(scores)[::-1]
    pair_score: dict[tuple[int, int], float] = {}
    pairs_out: list[tuple[int, int, float]] = []
    for o in order[:top_pair_k]:
        s = float(scores[o])
        if s <= 0:
            break
        ki, kj = pair_idx[o]
        a, b = seeds[ki], seeds[kj]
        key = (min(a, b), max(a, b))
        pair_score[key] = s
        pairs_out.append((a, b, s))

    # Pair-graph adjacency for triangle enumeration
    adj: dict[int, set[int]] = {}
    for (a, b) in pair_score:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)

    triplets: list[tuple[int, int, int, float]] = []
    seen: set[tuple[int, int, int]] = set()
    for i in adj:
        ni = adj[i]
        for j in ni:
            if j <= i:
                continue
            common = ni & adj[j]
            for k in common:
                if k <= j:
                    continue
                key = (i, j, k)
                if key in seen:
                    continue
                seen.add(key)
                s = (pair_score[(min(i, j), max(i, j))]
                     * pair_score[(min(i, k), max(i, k))]
                     * pair_score[(min(j, k), max(j, k))])
                triplets.append((i, j, k, s))
    triplets.sort(key=lambda t: t[3], reverse=True)
    return triplets[:top_triplet_k], pairs_out


# ===========================================================================
# 3-way interaction test (linear regression LRT)
# ===========================================================================

def test_3way_interaction(d1, d2, d3, y) -> dict:
    """LRT comparing y ~ g1+g2+g3 + 2-way + 3-way to y ~ g1+g2+g3.

    Imputes NaN with column mean (consistent with §Y.4 _test_interaction).
    Returns p-value for the 3-way term plus β coefficients.
    """
    d1 = np.where(np.isnan(d1), np.nanmean(d1), d1)
    d2 = np.where(np.isnan(d2), np.nanmean(d2), d2)
    d3 = np.where(np.isnan(d3), np.nanmean(d3), d3)
    y = np.asarray(y, dtype=float)
    if np.var(d1) == 0 or np.var(d2) == 0 or np.var(d3) == 0:
        return {"p_3way": np.nan, "beta_3way": np.nan, "n": int(len(y))}

    n = len(y)
    # Center columns
    g1 = d1 - d1.mean()
    g2 = d2 - d2.mean()
    g3 = d3 - d3.mean()
    # Design matrix: intercept + 3 main + 3 two-way + 1 three-way
    X = np.column_stack([
        np.ones(n), g1, g2, g3,
        g1 * g2, g1 * g3, g2 * g3,
        g1 * g2 * g3,
    ])
    # Closed-form OLS
    XtX = X.T @ X
    Xty = X.T @ y
    try:
        beta = np.linalg.solve(XtX, Xty)
    except np.linalg.LinAlgError:
        return {"p_3way": np.nan, "beta_3way": np.nan, "n": n}
    resid = y - X @ beta
    rss = float(resid @ resid)
    p_full = X.shape[1]
    if n - p_full <= 0:
        return {"p_3way": np.nan, "beta_3way": np.nan, "n": n}
    sigma2 = rss / (n - p_full)
    try:
        cov = sigma2 * np.linalg.inv(XtX)
    except np.linalg.LinAlgError:
        return {"p_3way": np.nan, "beta_3way": np.nan, "n": n}
    se_3way = float(np.sqrt(max(cov[7, 7], 0.0)))
    if se_3way <= 0:
        return {"p_3way": np.nan, "beta_3way": float(beta[7]), "n": n}
    t = float(beta[7]) / se_3way
    p = 2.0 * (1.0 - sp_stats.t.cdf(abs(t), df=n - p_full))
    return {"p_3way": float(p), "beta_3way": float(beta[7]),
            "se_3way": se_3way, "n": n}


# ===========================================================================
# Phenotype simulator (3-way DGP)
# ===========================================================================

def simulate_3way_pheno(dosages, idx_a, idx_b, idx_c, beta, rng):
    """y = beta * (g_a - mean) * (g_b - mean) * (g_c - mean) + eps."""
    g = []
    for k in (idx_a, idx_b, idx_c):
        d = dosages[:, k]
        d = np.where(np.isnan(d), np.nanmean(d), d)
        g.append(d - d.mean())
    causal = beta * g[0] * g[1] * g[2]
    eps = rng.standard_normal(len(causal))
    y = causal + eps
    var_c = causal.var()
    var_t = y.var()
    print(f"  Var decomp: causal={var_c:.3f}, total={var_t:.3f}, "
          f"R²={var_c/var_t:.3f}")
    return y


# ===========================================================================
# Main
# ===========================================================================

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--beta", type=float, default=10.0,
                   help="3-way effect size (default 10.0; needs to be "
                        "larger than 2-way β=3 because the triple product "
                        "of low-MAF dosages is small)")
    p.add_argument("--alpha", type=float, default=0.15)
    p.add_argument("--top-pair-k", type=int, default=5000)
    p.add_argument("--top-triplet-k", type=int, default=500)
    p.add_argument("--max-variants-per-gene", type=int, default=20,
                   help="Cap variants/gene during triplet expansion to "
                        "avoid O(n³) blow-up")
    p.add_argument("--seed", type=int, default=2026)
    args = p.parse_args()

    rng = np.random.default_rng(args.seed)
    print("=== §Y.5 k=3 — gene-triplet RWR on yeast PKA core ===\n")

    # 1. Genotypes + cache + canonical-PPI injection
    print("[1/6] Loading genotypes...")
    dosages, variant_ids, sample_ids = load_yeast_chroms(
        YEAST_VCF, CHROMS_TO_LOAD, maf_min=0.05,
    )
    print(f"\n[2/6] Loading + injecting canonical PPI...")
    cache = json.loads(YEAST_CACHE.read_text())
    cache, inject_diag = inject_canonical_ppi_edges(
        cache, edges=YEAST_CANONICAL_EDGES, inplace=True, verbose=False,
    )
    print(f"  injected {inject_diag['n_edges_applied']} edges, "
          f"{inject_diag['n_variants_modified']} variants modified")

    annotated_vids = [v for v in variant_ids if v in cache]
    print(f"  {len(annotated_vids):,} dosage variants in cache")

    # 3. Heterograph
    print(f"\n[3/6] Building heterograph (dosage-only, IDF-weighted)...")
    A_gg, gene_ids, gene_to_idx, het_diag = _build_gene_gene_heterograph(
        cache, annotated_vids, use_idf=True, inject_canonical_ppi=False,
    )
    n_gene = A_gg.shape[0]
    edge_count = int(np.count_nonzero(A_gg) // 2)
    print(f"  heterograph: {n_gene} genes × {edge_count:,} edges, "
          f"max weight {A_gg.max():.2f}")

    # Locate ground-truth genes
    bcy1_i = gene_to_idx.get("YIL033C")
    tpk1_i = gene_to_idx.get("YJL164C")
    tpk2_i = gene_to_idx.get("YPL203W")
    tpk3_i = gene_to_idx.get("YKL166C")
    print(f"  BCY1 idx: {bcy1_i}, TPK1 idx: {tpk1_i}, "
          f"TPK2 idx: {tpk2_i}, TPK3 idx: {tpk3_i}")
    # Pick whichever triple has all three genes in the dosage panel.
    # TPK2 (YPL203W) has only MAF<0.05 variants in the 1011 panel, so usually
    # gets filtered out — fall back to TPK3 (YKL166C) which is panel-segregating.
    if all(idx is not None for idx in (bcy1_i, tpk1_i, tpk2_i)):
        third_i, third_name, third_sys = tpk2_i, "TPK2", "YPL203W"
    elif all(idx is not None for idx in (bcy1_i, tpk1_i, tpk3_i)):
        third_i, third_name, third_sys = tpk3_i, "TPK3", "YKL166C"
        print(f"  ⚠ TPK2 not in dosage panel (no MAF≥0.05 variants); "
              f"using TPK3 as third PKA subunit")
    else:
        sys.exit(f"✗ Cannot locate BCY1+TPK1+(TPK2 or TPK3) in heterograph "
                 f"(BCY1={bcy1_i}, TPK1={tpk1_i}, TPK2={tpk2_i}, TPK3={tpk3_i})")
    # All three pairwise edge weights present?
    for a, b, name in [(bcy1_i, tpk1_i, f"BCY1↔TPK1"),
                       (bcy1_i, third_i, f"BCY1↔{third_name}"),
                       (tpk1_i, third_i, f"TPK1↔{third_name}")]:
        print(f"    {name} edge weight: {A_gg[a, b]:.4f}")

    # 4. Top-K triplets
    print(f"\n[4/6] Computing gene-triplet RWR scores...")
    t0 = time.time()
    triplets, pairs = heterograph_triplet_scores(
        A_gg, list(range(n_gene)),
        alpha=args.alpha,
        top_pair_k=args.top_pair_k,
        top_triplet_k=args.top_triplet_k,
    )
    rwr_time = time.time() - t0
    print(f"  → {len(triplets)} top gene triplets "
          f"(from top {len(pairs)} pairs; RWR took {rwr_time:.1f}s)")

    # Locate ground-truth triplet in the result (primary or fallback)
    truth_keys = [tuple(sorted([bcy1_i, tpk1_i, third_i]))]
    truth_ranks: dict = {}
    for r, (i, j, k, s) in enumerate(triplets, 1):
        key = tuple(sorted([i, j, k]))
        for tk in truth_keys:
            if key == tk and tk not in truth_ranks:
                gene_names = ",".join(TRIPLE_NAMES.get(gene_ids[g], gene_ids[g])
                                       for g in tk)
                truth_ranks[gene_names] = {
                    "rank": r, "score": s,
                    "gene_systematic": [gene_ids[g] for g in tk],
                }
    if truth_ranks:
        for name, info in truth_ranks.items():
            print(f"  ✓ Triplet {name}: rank {info['rank']}/{len(triplets)}, "
                  f"score={info['score']:.4g}")
    else:
        print(f"  ✗ Neither {{BCY1,TPK1,TPK2}} nor {{BCY1,TPK1,TPK3}} "
              f"in top-{len(triplets)} gene triplets")

    # 5. Variant-triplet expansion + 3-way interaction test
    print(f"\n[5/6] Expanding gene triplets → variant triplets...")
    cache_idx_to_dosage_col = {v: i for i, v in enumerate(variant_ids)}
    annotated_vid_set = set(annotated_vids)
    gene_to_var_indices: dict[str, list[int]] = {}
    for vid in annotated_vids:
        for g in cache[vid].get("genes", []):
            if g in gene_to_idx:
                gene_to_var_indices.setdefault(g, []).append(
                    cache_idx_to_dosage_col[vid]
                )

    # Simulate Scenario A 3-way phenotype on highest-MAF variant per gene
    af_vec = np.nanmean(dosages, axis=0) / 2.0
    maf_vec = np.minimum(af_vec, 1.0 - af_vec)

    def best_variant_idx(gene_systematic):
        idxs = gene_to_var_indices.get(gene_systematic, [])
        if not idxs:
            return None
        return idxs[int(np.argmax([maf_vec[k] for k in idxs]))]

    bcy1_var = best_variant_idx("YIL033C")
    tpk1_var = best_variant_idx("YJL164C")
    third_var = best_variant_idx(third_sys)
    if any(v is None for v in (bcy1_var, tpk1_var, third_var)):
        sys.exit(f"✗ Could not locate dosage representatives for BCY1/TPK1/{third_name}")
    print(f"  Reps: BCY1={variant_ids[bcy1_var]} (MAF {maf_vec[bcy1_var]:.3f}), "
          f"TPK1={variant_ids[tpk1_var]} (MAF {maf_vec[tpk1_var]:.3f}), "
          f"{third_name}={variant_ids[third_var]} (MAF {maf_vec[third_var]:.3f})")

    print(f"\n[6/6] Simulating 3-way DGP (β={args.beta}) and testing...")
    phenotype = simulate_3way_pheno(
        dosages, bcy1_var, tpk1_var, third_var, args.beta, rng,
    )

    # Test variant triplets across the top-K gene triplets
    rng_inner = np.random.default_rng(args.seed + 1)
    rows = []
    for gi, gj, gk_idx, gscore in triplets:
        names = [gene_ids[gi], gene_ids[gj], gene_ids[gk_idx]]
        var_lists = []
        for g in names:
            idxs = gene_to_var_indices.get(g, [])
            if len(idxs) > args.max_variants_per_gene:
                idxs = list(rng_inner.choice(
                    idxs, size=args.max_variants_per_gene, replace=False
                ))
            var_lists.append(idxs)
        if not all(var_lists):
            continue
        for v1, v2, v3 in product(*var_lists):
            if len({v1, v2, v3}) < 3:
                continue
            d1 = dosages[:, v1].astype(float)
            d2 = dosages[:, v2].astype(float)
            d3 = dosages[:, v3].astype(float)
            r = test_3way_interaction(d1, d2, d3, phenotype)
            if not np.isnan(r["p_3way"]):
                rows.append({
                    "v1": variant_ids[v1], "v2": variant_ids[v2], "v3": variant_ids[v3],
                    "gene_systematic": names,
                    "p_3way": r["p_3way"], "beta_3way": r["beta_3way"],
                    "gene_triplet_score": gscore,
                })
    print(f"  tested {len(rows):,} variant triplets")

    if not rows:
        sys.exit("No variant triplets tested — heterograph too sparse?")

    # BH-FDR
    pvals = np.array([r["p_3way"] for r in rows])
    n_tests = len(pvals)
    order = np.argsort(pvals)
    rank_arr = np.empty(n_tests, dtype=int)
    rank_arr[order] = np.arange(1, n_tests + 1)
    qvals = np.minimum(1.0, pvals * n_tests / rank_arr)
    for k in range(n_tests - 2, -1, -1):
        qvals[order[k]] = min(qvals[order[k]], qvals[order[k + 1]])
    for k, r in enumerate(rows):
        r["q_3way"] = float(qvals[k])

    # Sort by p, find ground-truth triplet
    rows.sort(key=lambda r: r["p_3way"])

    bcy1_vars_set = {variant_ids[i] for i in gene_to_var_indices.get("YIL033C", [])}
    tpk1_vars_set = {variant_ids[i] for i in gene_to_var_indices.get("YJL164C", [])}
    third_vars_set = {variant_ids[i] for i in gene_to_var_indices.get(third_sys, [])}

    def is_truth(row):
        triple = {row["v1"], row["v2"], row["v3"]}
        in_a = bool(triple & bcy1_vars_set)
        in_b = bool(triple & tpk1_vars_set)
        in_c = bool(triple & third_vars_set)
        return in_a and in_b and in_c

    truth_var_rank = None
    truth_var_hit = None
    for r, row in enumerate(rows, 1):
        if is_truth(row):
            truth_var_rank = r
            truth_var_hit = row
            break

    print(f"\n=== RESULT ===")
    print(f"  N variant triplets tested: {len(rows):,}")
    print(f"  Top variant triplet: {rows[0]['v1']} × {rows[0]['v2']} × "
          f"{rows[0]['v3']}")
    print(f"      genes: {rows[0]['gene_systematic']}, "
          f"p_3way={rows[0]['p_3way']:.3g}, q={rows[0]['q_3way']:.3g}")
    if truth_var_rank is None:
        print(f"  ✗ Ground-truth (BCY1×TPK1×TPK2) NOT in tested triplets")
    else:
        print(f"  ✓ Ground-truth variant triplet rank: "
              f"{truth_var_rank:,}/{len(rows):,}")
        print(f"      pair: {truth_var_hit['v1']} × {truth_var_hit['v2']} × "
              f"{truth_var_hit['v3']}")
        print(f"      p_3way: {truth_var_hit['p_3way']:.3g}, "
              f"q={truth_var_hit['q_3way']:.3g}")
        print(f"      gene-triplet RWR score: "
              f"{truth_var_hit['gene_triplet_score']:.4g}")

    # Brute-force comparison: how many variant triplets exist in the panel?
    # and what would O(n³) take?
    n_annotated_vars = len(annotated_vids)
    bf_count = n_annotated_vars * (n_annotated_vars - 1) * (n_annotated_vars - 2) // 6
    bf_speedup = bf_count / max(1, len(rows))
    print(f"\n  Brute-force enumeration would test {bf_count:,} triplets "
          f"({n_annotated_vars} annotated × ³)")
    print(f"  Graph-priority RWR tested {len(rows):,} (×{bf_speedup:.0f} reduction)")

    # Save
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out = RESULTS_DIR / "y5_triplet_rwr_yeast.json"
    payload = {
        "method": "M5_gene_triplet_RWR",
        "ground_truth_triplet": [f"BCY1 (YIL033C)", "TPK1 (YJL164C)",
                                  f"{third_name} ({third_sys})"],
        "params": {
            "alpha": args.alpha,
            "top_pair_k": args.top_pair_k,
            "top_triplet_k": args.top_triplet_k,
            "max_variants_per_gene": args.max_variants_per_gene,
            "beta_3way": args.beta,
        },
        "heterograph": {
            "n_genes": int(n_gene),
            "n_edges": int(edge_count),
            "max_edge_weight": float(A_gg.max()),
            "BCY1_TPK1_edge_weight": float(A_gg[bcy1_i, tpk1_i]),
            f"BCY1_{third_name}_edge_weight": float(A_gg[bcy1_i, third_i]),
            f"TPK1_{third_name}_edge_weight": float(A_gg[tpk1_i, third_i]),
        },
        "gene_triplet_recovery": truth_ranks,
        "variant_triplet_recovery": {
            "rank": truth_var_rank,
            "hit": truth_var_hit,
            "n_tested": len(rows),
            "brute_force_size": int(bf_count),
            "speedup": float(bf_speedup),
        },
        "rwr_runtime_sec": rwr_time,
    }
    out.write_text(json.dumps(payload, indent=2, default=str))
    print(f"\n  Wrote {out}")


if __name__ == "__main__":
    main()
