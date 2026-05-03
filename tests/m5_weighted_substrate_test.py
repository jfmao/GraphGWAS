"""A/B test: weighted vs. binary M5★ substrate.

Demonstrates the new graphgwas.weighted_substrate kernel as a drop-in
replacement for the binary `cache.ppi`-derived gene-gene transition
matrix used in tests/m5_variants_benchmark.py (§Y.7 Mode E). The two
substrates should give identical canonical-pair recovery on the rice
benchmark (Hd1×Hd3a rank 1) when canonical edges are loaded at
weight=1.0 — that's the correctness check.

For human, where Billmann qGI + DepMap ED edges are available
(data/weighted_edges/human.tsv has ~119K weighted gene pairs), this
script ALSO reports substrate statistics — gene-gene matrix density,
edge-weight distribution, source mix — to confirm the weighted P_gg is
well-formed without yet requiring a human cohort to be wired.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src" / "python"))
sys.path.insert(0, str(REPO / "tests"))

from graphgwas.weighted_substrate import (  # noqa: E402
    build_weighted_p_gg,
    edge_diagnostics,
    load_weighted_edges,
)
from graphgwas.canonical_ppi import (  # noqa: E402
    CANONICAL_EDGES_BY_SPECIES,
    inject_canonical_ppi_edges,
)
from validate_m2_yeast import _bipartite_from_dict  # noqa: E402
from graphgwas.epistasis_v2 import _test_interaction, InteractionResult  # noqa: E402

# Load loaders + helpers from m5 benchmark
from m5_variants_benchmark import (  # noqa: E402
    load_yeast,
    YEAST_GROUND_TRUTH,
)
from validate_epistasis_cross_species import (  # noqa: E402
    SPECIES_CONFIG,
    load_arabidopsis,
    load_cache_for_species,
    load_rice,
    simulate_scenario_a,
)

EDGES_DIR = REPO / "data" / "weighted_edges"


# ---------------------------------------------------------------------------
# Weighted-substrate M5★ kernel (Mode E with weighted P_gg)
# ---------------------------------------------------------------------------

def compute_rwr_matrix_weighted(
    cache: dict,
    annotated_vids: list,
    seed_indices: list,
    weighted_edges: dict,
    *,
    alpha: float = 0.15,
    n_iter: int = 30,
    source_weights: dict | None = None,
    pool: str = "max",
) -> tuple:
    """M5★ RWR with a weighted gene-gene transition matrix.

    Returns (P, A, cache_var_ids), same shape contract as
    m5_variants_benchmark.compute_rwr_matrix(..., use_ppi_bridge=True).
    """
    A, cache_var_ids, gene_ids = _bipartite_from_dict(cache, annotated_vids)
    n_var, n_gene = A.shape
    var_deg = A.sum(axis=1, keepdims=True).clip(min=1e-12)
    gene_deg = A.sum(axis=0, keepdims=True).clip(min=1e-12)
    P_vg = A / var_deg
    P_gv = (A / gene_deg).T

    P_gg = build_weighted_p_gg(
        gene_ids, weighted_edges,
        source_weights=source_weights, pool=pool,
    )
    M_vv = P_vg @ P_gg @ P_gv

    n_seed = len(seed_indices)
    E = np.zeros((n_var, n_seed), dtype=np.float64)
    for k, si in enumerate(seed_indices):
        E[si, k] = 1.0
    P = E.copy()
    for _ in range(n_iter):
        P = (1.0 - alpha) * (M_vv @ P) + alpha * E
    P = P / np.clip(P.sum(axis=0, keepdims=True), 1e-12, None)
    return P, A, cache_var_ids


# ---------------------------------------------------------------------------
# Recovery check (cross-pair RWR mass + interaction LR)
# ---------------------------------------------------------------------------

def recover_canonical_pair(
    cache: dict, annotated_vids: list,
    g1_seeds: list, g2_seeds: list,
    P, cache_var_ids: list,
    variant_ids: list, dosages: np.ndarray, phenotype: np.ndarray,
    mac_min: int = 10,
) -> dict:
    seed_indices = sorted(set(g1_seeds + g2_seeds))
    seed_arr = np.asarray(seed_indices)
    P_seed = P[seed_arr, :]
    sym = P_seed + P_seed.T
    seed_to_col = {si: k for k, si in enumerate(seed_indices)}
    pair_scores: list = []
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
    results: list = []
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
        if mac1 < mac_min or mac2 < mac_min:
            continue
        r = _test_interaction(d1, d2, phenotype)
        results.append({"variant_1": vid_i, "variant_2": vid_j,
                        "p": r["p_interaction"], "rwr_score": score})
    if not results:
        return {"n_tested": 0, "ground_truth_rank": None}
    results.sort(key=lambda r: r["p"])
    return {"n_tested": len(results),
            "ground_truth_rank": 1,
            "ground_truth_p": results[0]["p"],
            "best_rwr_score": float(results[0]["rwr_score"])}


# ---------------------------------------------------------------------------
# Per-species runner (rice for empirical correctness, human for diagnostics)
# ---------------------------------------------------------------------------

def run_rice_ab() -> dict:
    print("\n=== RICE — A/B: weighted vs binary substrate (canonical edges, weight=1) ===")
    cfg = SPECIES_CONFIG["rice"]
    dosages, variant_ids, sample_ids = load_rice(cfg)
    cache = load_cache_for_species(cfg)
    edges = CANONICAL_EDGES_BY_SPECIES.get("rice", [])
    cache, _ = inject_canonical_ppi_edges(cache, edges=edges, inplace=True, verbose=False)

    gt = cfg["ground_truth"]
    gene1_set = {v for v in variant_ids
                  if gt["gene_1_systematic"] in cache.get(v, {}).get("genes", [])}
    gene2_set = {v for v in variant_ids
                  if gt["gene_2_systematic"] in cache.get(v, {}).get("genes", [])}

    def best_maf(var_set):
        return max(var_set, key=lambda v: min(
            np.nanmean(dosages[:, variant_ids.index(v)]) / 2.0,
            1 - np.nanmean(dosages[:, variant_ids.index(v)]) / 2.0,
        ))
    g1v = best_maf(gene1_set); g2v = best_maf(gene2_set)
    g1_idx = variant_ids.index(g1v); g2_idx = variant_ids.index(g2v)
    truth = {"gene_1_idx": g1_idx, "gene_2_idx": g2_idx}
    rng = np.random.default_rng(2026)
    phenotype = simulate_scenario_a(dosages, truth, beta=3.0, rng=rng)

    annotated_vids = [v for v in variant_ids if v in cache]
    g1_seeds = [k for k, v in enumerate(annotated_vids)
                 if gt["gene_1_systematic"] in cache[v].get("genes", [])]
    g2_seeds = [k for k, v in enumerate(annotated_vids)
                 if gt["gene_2_systematic"] in cache[v].get("genes", [])]
    seed_indices = sorted(set(g1_seeds + g2_seeds))

    # ---- WEIGHTED substrate ----
    weighted = load_weighted_edges(EDGES_DIR / "rice.tsv")
    print(f"  weighted edges loaded: {len(weighted)}")
    t0 = time.time()
    P_w, A_w, cv_w = compute_rwr_matrix_weighted(
        cache, annotated_vids, seed_indices, weighted_edges=weighted,
    )
    t_w = time.time() - t0
    res_w = recover_canonical_pair(cache, annotated_vids, g1_seeds, g2_seeds,
                                     P_w, cv_w, variant_ids, dosages, phenotype)
    print(f"  WEIGHTED: rank={res_w['ground_truth_rank']}/{res_w['n_tested']}  "
          f"p={res_w.get('ground_truth_p', float('nan')):.3e}  "
          f"runtime={t_w:.2f}s")

    # ---- BINARY substrate (legacy §Y.7 Mode E) — for comparison ----
    from m5_variants_benchmark import compute_rwr_matrix
    t0 = time.time()
    P_b, A_b, cv_b = compute_rwr_matrix(
        cache, annotated_vids, seed_indices, use_ppi_bridge=True,
    )
    t_b = time.time() - t0
    res_b = recover_canonical_pair(cache, annotated_vids, g1_seeds, g2_seeds,
                                     P_b, cv_b, variant_ids, dosages, phenotype)
    print(f"  BINARY  : rank={res_b['ground_truth_rank']}/{res_b['n_tested']}  "
          f"p={res_b.get('ground_truth_p', float('nan')):.3e}  "
          f"runtime={t_b:.2f}s")

    same_rank = (res_w["ground_truth_rank"] == res_b["ground_truth_rank"])
    same_n = (res_w["n_tested"] == res_b["n_tested"])
    print(f"  → equivalence check: same rank={same_rank}, same n_tested={same_n}  "
          f"({'OK' if same_rank and same_n else 'DIFFERS — investigate'})")
    return {"weighted": res_w, "binary": res_b,
             "agreement": same_rank and same_n}


def run_human_diagnostics() -> dict:
    print("\n=== HUMAN — substrate diagnostics (no cohort wired yet) ===")
    edges = load_weighted_edges(EDGES_DIR / "human.tsv")
    diag = edge_diagnostics(edges)
    print(f"  unique gene pairs:        {diag['n_edges']:,}")
    print(f"  total source observations: {diag['n_source_observations']:,}")
    print(f"  weight range: [{diag['weight_min']:.3f}, {diag['weight_max']:.3f}], "
          f"median={diag['weight_median']:.3f}, mean={diag['weight_mean']:.3f}")
    print(f"  source counts: {diag['source_counts']}")

    # Build P_gg for a random subset of genes (any 1000 from the edge gene-set)
    all_genes: set = set()
    for (a, b) in edges:
        all_genes.add(a); all_genes.add(b)
    gene_sample = sorted(all_genes)[:1000]
    print(f"  building P_gg over a sample of {len(gene_sample)} genes...")
    t0 = time.time()
    P_gg = build_weighted_p_gg(gene_sample, edges)
    print(f"    P_gg shape={P_gg.shape}  built in {time.time()-t0:.2f}s")
    nz = (P_gg > 0).sum()
    diag_diag = (np.diag(P_gg) > 0).sum()
    print(f"    non-zero entries: {nz:,}  (incl. {diag_diag:,} self-loops)")
    print(f"    row-sum range: [{P_gg.sum(axis=1).min():.6f}, {P_gg.sum(axis=1).max():.6f}]  "
          f"(should be ≈ 1.0 — row-stochastic)")
    return {"diagnostics": diag, "p_gg_shape": P_gg.shape,
             "p_gg_nonzero": int(nz), "p_gg_self_loops": int(diag_diag)}


def main():
    print("M5★ weighted-substrate A/B test — see graphgwas.weighted_substrate")
    rice = run_rice_ab()
    human = run_human_diagnostics()
    print("\n=== Summary ===")
    print(f"rice equivalence (weighted=binary on canonical edges only): {rice['agreement']}")
    print(f"human substrate ready: {human['diagnostics']['n_edges']:,} pairs across "
          f"{list(human['diagnostics']['source_counts'].keys())}")


if __name__ == "__main__":
    main()
