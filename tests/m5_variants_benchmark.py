"""Paper #2 §Y.5/§Y.6 follow-up — M5 RWR variant benchmark.

Investigates three remediations for the M5 RWR failure that left Hd1 × Hd3a
(rice), FT × FLC (Arabidopsis), and BCY1 × TPK1 (yeast) all NF in the
top-2,000 mutual-RWR pair set:

  (A) **High-cap**       — keep mutual-RWR scoring but raise top_k from
                            2,000 to a much larger value (default 50,000)
                            so the canonical pair has more room to surface.
  (B) **Hub-deflated**   — replace mutual = P_a[b] * P_b[a] with
                            score = (P_a[b] + P_b[a]) / (deg_a * deg_b)
                            where deg_x is the variant's annotated-gene
                            count in the bipartite graph (resource-
                            allocation-style normalisation).
  (C) **Focused seeds**  — seed RWR only from variants annotated to the
                            canonical pair's two genes (no other
                            annotated variants in the seed set), then
                            score every (g1_var × g2_var) cross-pair.
                            Models a two-stage workflow: a coarse
                            gene-pair selector (from §Y.5 gene-priority
                            RWR or M2 motif-filtering) hands a candidate
                            gene-pair to M5, which finds the best
                            variant-pair representative.

Runs each variant on all three species and reports the canonical-pair rank.
A variant is considered to recover the pair if the rank lands in the
top-1% of pairs scored. Output:

    results/paper2_epistasis/m5_variants_benchmark.json

This is a pure-method-tuning script: phenotype DGP, MAF filter, and
canonical-PPI injection are identical to validate_epistasis_cross_species.py.
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

from graphgwas.epistasis_v2 import (  # noqa: E402
    InteractionResult,
    _test_interaction,
)
from validate_m2_yeast import _bipartite_from_dict  # noqa: E402

# Reuse the loaders + species config from the cross-species runner.
from validate_epistasis_cross_species import (  # noqa: E402
    SPECIES_CONFIG,
    load_arabidopsis,
    load_cache_for_species,
    load_rice,
    pick_representatives_cross_species,
    simulate_scenario_a,
)

RESULTS_DIR = REPO_ROOT / "results" / "paper2_epistasis"

# ---------------------------------------------------------------------------
# Yeast loader (script-local — borrowed from validate_m2_yeast)
# ---------------------------------------------------------------------------

YEAST_GROUND_TRUTH = {
    "gene_1_systematic": "YIL033C",
    "gene_1_symbol":     "BCY1",
    "gene_2_systematic": "YJL164C",
    "gene_2_symbol":     "TPK1",
}


def load_yeast() -> tuple[np.ndarray, list, list, dict]:
    """Load yeast chr9/10/12 dosages + cache by reusing validate_m2_yeast."""
    from validate_m2_yeast import load_yeast_chroms, YEAST_VCF, YEAST_CACHE
    dosages, vids, sample_ids = load_yeast_chroms(
        YEAST_VCF, chroms=["chromosome9", "chromosome10", "chromosome12"],
        maf_min=0.05,
    )
    cache = json.loads(Path(YEAST_CACHE).read_text())
    return dosages, vids, sample_ids, cache


# ---------------------------------------------------------------------------
# Core RWR (shared by all 4 modes)
# ---------------------------------------------------------------------------

def compute_rwr_matrix(cache: dict,
                        annotated_vids: list,
                        seed_indices: list,
                        alpha: float = 0.15,
                        n_iter: int = 30,
                        use_ppi_bridge: bool = False) -> tuple:
    """Build bipartite, run batched RWR seeded from `seed_indices`.

    If ``use_ppi_bridge`` is True, the variant→variant transition is
    augmented with a gene-gene PPI hop: each (g_a, g_b) PPI annotation
    in the cache becomes a gene-gene edge, and the walk becomes
    variant → gene_a → gene_b → variant. This bridges variant pairs
    whose genes are in different bipartite components but share a
    canonical PPI edge (e.g. Hd1 → Hd3a via the rice flowering edge).

    Returns (P, A, cache_var_ids) where:
      P: shape (n_var, n_seed); P[i, k] = stationary mass at variant i
         when the walk was seeded from variant seed_indices[k].
      A: bipartite adjacency, shape (n_var, n_gene).
      cache_var_ids: list aligning variant index k to its variantId.
    """
    A, cache_var_ids, gene_ids = _bipartite_from_dict(cache, annotated_vids)
    n_var, n_gene = A.shape
    var_deg = A.sum(axis=1, keepdims=True).clip(min=1e-12)
    gene_deg = A.sum(axis=0, keepdims=True).clip(min=1e-12)
    P_vg = A / var_deg
    P_gv = (A / gene_deg).T

    if use_ppi_bridge:
        # Build gene-gene PPI adjacency from the cache's `ppi` field.
        # Symmetrise; include self-loops so the walk can stay at the gene.
        gene_to_idx = {g: i for i, g in enumerate(gene_ids)}
        G = np.zeros((n_gene, n_gene), dtype=np.float64)
        np.fill_diagonal(G, 1.0)
        for vid in cache_var_ids:
            entry = cache[vid]
            entry_genes = [g for g in entry.get("genes", []) if g in gene_to_idx]
            entry_ppis = [g for g in entry.get("ppi", []) if g in gene_to_idx]
            for ga in entry_genes:
                ia = gene_to_idx[ga]
                for gb in entry_ppis:
                    if ga == gb:
                        continue
                    ib = gene_to_idx[gb]
                    G[ia, ib] = 1.0
                    G[ib, ia] = 1.0
        # Row-normalise G to a transition matrix
        G_row_sum = G.sum(axis=1, keepdims=True).clip(min=1e-12)
        P_gg = G / G_row_sum
        # 4-step transition: variant → gene → gene' → variant
        M_vv = P_vg @ P_gg @ P_gv
    else:
        # Original 2-step transition: variant → gene → variant
        M_vv = P_vg @ P_gv

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
# Test interaction on a list of (cache_idx_i, cache_idx_j, score) triples
# ---------------------------------------------------------------------------

def test_pairs_and_rank(pair_scores: list,
                          cache_var_ids: list,
                          variant_ids: list,
                          dosages: np.ndarray,
                          phenotype: np.ndarray,
                          gene1_variants: set,
                          gene2_variants: set,
                          mac_min: int) -> dict:
    """Run interaction LR on each (i_cache, j_cache, score) and rank by p.

    Returns dict with: ground_truth_rank (None if no g1×g2 pair surfaced),
    ground_truth_p, ground_truth_q, n_tested, sig_at_q05, best_p.
    """
    cache_idx_to_col = {v: variant_ids.index(v) for v in cache_var_ids}

    results: list = []
    for i_cache, j_cache, score in pair_scores:
        vid_i = cache_var_ids[i_cache]
        vid_j = cache_var_ids[j_cache]
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
        results.append({
            "variant_1": vid_i, "variant_2": vid_j,
            "rwr_score": float(score),
            "p_interaction": r["p_interaction"],
        })
    if not results:
        return {"n_tested": 0, "ground_truth_rank": None,
                "ground_truth_p": None, "ground_truth_q": None,
                "sig_at_q05": 0, "best_p": None}

    # Sort by p_interaction ascending (smallest = best)
    results.sort(key=lambda r: r["p_interaction"])
    n = len(results)

    # BH-FDR
    pvals = np.array([r["p_interaction"] for r in results])
    order = np.argsort(pvals)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)
    qvals = np.minimum(1.0, pvals * n / ranks)
    for k in range(n - 2, -1, -1):
        qvals[order[k]] = min(qvals[order[k]], qvals[order[k + 1]])
    sig_at_q05 = int(np.sum(qvals < 0.05))

    # Find canonical pair: any (g1, g2) cross-pair in results
    gt_rank = None
    gt_p = None
    gt_q = None
    for rk, r in enumerate(results, start=1):
        v1, v2 = r["variant_1"], r["variant_2"]
        is_gt = ((v1 in gene1_variants and v2 in gene2_variants) or
                  (v1 in gene2_variants and v2 in gene1_variants))
        if is_gt:
            gt_rank = rk
            gt_p = r["p_interaction"]
            # find q for this index
            gt_q = float(qvals[rk - 1])
            break

    return {
        "n_tested": n,
        "ground_truth_rank": gt_rank,
        "ground_truth_p": gt_p,
        "ground_truth_q": gt_q,
        "sig_at_q05": sig_at_q05,
        "best_p": float(results[0]["p_interaction"]),
        "top1pct_threshold": int(np.ceil(n * 0.01)),
    }


# ---------------------------------------------------------------------------
# Mode runners
# ---------------------------------------------------------------------------

def run_baseline(P, seed_indices, top_k: int) -> list:
    """Mode 0: original mutual-RWR scoring with given top_k cap."""
    seed_arr = np.asarray(seed_indices)
    P_seed = P[seed_arr, :]
    mutual = P_seed * P_seed.T
    iu = np.triu_indices(len(seed_indices), k=1)
    scores = mutual[iu]
    pair_idx_pairs = list(zip(iu[0], iu[1]))
    order = np.argsort(scores)[::-1]
    pair_scores: list = []
    for o in order:
        s = float(scores[o])
        if s <= 0:
            break
        ki, kj = pair_idx_pairs[o]
        pair_scores.append((seed_indices[ki], seed_indices[kj], s))
        if len(pair_scores) >= top_k:
            break
    return pair_scores


def run_hub_deflated(P, A, seed_indices, top_k: int) -> list:
    """Mode B: (P[i,j_col] + P[j,i_col]) / (deg_i * deg_j) where deg is
    the variant's bipartite degree (number of annotated genes)."""
    seed_arr = np.asarray(seed_indices)
    P_seed = P[seed_arr, :]                        # (n_seed, n_seed)
    sym = (P_seed + P_seed.T)
    deg = np.asarray(A.sum(axis=1)).flatten()      # variant bipartite degree
    deg_seed = np.maximum(deg[seed_arr], 1.0)      # avoid div-by-zero
    norm = np.outer(deg_seed, deg_seed)
    score_mat = sym / norm
    iu = np.triu_indices(len(seed_indices), k=1)
    scores = score_mat[iu]
    pair_idx_pairs = list(zip(iu[0], iu[1]))
    order = np.argsort(scores)[::-1]
    pair_scores: list = []
    for o in order:
        s = float(scores[o])
        if s <= 0:
            break
        ki, kj = pair_idx_pairs[o]
        pair_scores.append((seed_indices[ki], seed_indices[kj], s))
        if len(pair_scores) >= top_k:
            break
    return pair_scores


def run_focused_seeds(cache, annotated_vids, g1_seeds, g2_seeds,
                       alpha: float = 0.15, n_iter: int = 30,
                       use_ppi_bridge: bool = False) -> tuple:
    """Mode C/E: seed RWR ONLY from g1+g2 variants; score all g1×g2 cross-pairs.

    With ``use_ppi_bridge=True`` (Mode E), the walk is variant→gene→gene'→variant,
    which allows the canonical PPI edge (e.g. Hd1↔Hd3a) to bridge cross-pairs.
    """
    seed_indices = sorted(set(g1_seeds + g2_seeds))
    P, A, cache_var_ids = compute_rwr_matrix(
        cache, annotated_vids, seed_indices, alpha=alpha, n_iter=n_iter,
        use_ppi_bridge=use_ppi_bridge,
    )
    seed_arr = np.asarray(seed_indices)
    P_seed = P[seed_arr, :]
    sym = (P_seed + P_seed.T)
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
    return pair_scores, P, A, cache_var_ids


# ---------------------------------------------------------------------------
# Per-species runner
# ---------------------------------------------------------------------------

def run_species(species: str, beta: float, seed: int,
                 top_k_high: int, mac_min: int = 10) -> dict:
    print(f"\n{'=' * 70}\n  SPECIES: {species.upper()}\n{'=' * 70}")
    rng = np.random.default_rng(seed)

    # Load
    if species == "yeast":
        dosages, variant_ids, sample_ids, cache = load_yeast()
        gt = YEAST_GROUND_TRUTH
    elif species == "arabidopsis":
        config = SPECIES_CONFIG["arabidopsis"]
        dosages, variant_ids, sample_ids = load_arabidopsis(config)
        cache = load_cache_for_species(config)
        gt = config["ground_truth"]
    elif species == "rice":
        config = SPECIES_CONFIG["rice"]
        dosages, variant_ids, sample_ids = load_rice(config)
        cache = load_cache_for_species(config)
        gt = config["ground_truth"]
    else:
        sys.exit(f"unknown species: {species}")

    # Inject canonical PPI for ALL species. The yeast cache also has the
    # BCY1↔TPK1 curation gap (per §Y.5), so injection is required for the
    # PPI-bridged substrate (Mode D / E) to recover the canonical pair.
    from graphgwas.canonical_ppi import (
        CANONICAL_EDGES_BY_SPECIES,
        inject_canonical_ppi_edges,
    )
    edges = CANONICAL_EDGES_BY_SPECIES.get(species, [])
    cache, inj_diag = inject_canonical_ppi_edges(
        cache, edges=edges, inplace=True, verbose=False,
    )
    print(f"      canonical-PPI injection: "
          f"{inj_diag['n_edges_applied']}/{inj_diag['n_edges_supplied']} "
          f"edges, {inj_diag['n_variants_modified']} variants modified")

    # Recovery sets
    gene1_variants = {
        vid for vid in variant_ids
        if gt["gene_1_systematic"] in cache.get(vid, {}).get("genes", [])
    }
    gene2_variants = {
        vid for vid in variant_ids
        if gt["gene_2_systematic"] in cache.get(vid, {}).get("genes", [])
    }
    print(f"      Recovery sets: {gt['gene_1_symbol']}={len(gene1_variants)} "
          f"{gt['gene_2_symbol']}={len(gene2_variants)}")

    # Phenotype simulation (Scenario A)
    g1v_pick = next(iter(gene1_variants))
    g2v_pick = next(iter(gene2_variants))
    g1_idx = variant_ids.index(g1v_pick)
    g2_idx = variant_ids.index(g2v_pick)
    # Pick highest-MAF representatives instead of arbitrary first
    def best_maf(var_set):
        af = lambda v: np.nanmean(dosages[:, variant_ids.index(v)]) / 2.0
        maf = lambda v: min(af(v), 1.0 - af(v))
        return max(var_set, key=maf)
    g1v_pick = best_maf(gene1_variants)
    g2v_pick = best_maf(gene2_variants)
    g1_idx = variant_ids.index(g1v_pick)
    g2_idx = variant_ids.index(g2v_pick)
    truth = {"gene_1_idx": g1_idx, "gene_2_idx": g2_idx}
    phenotype = simulate_scenario_a(dosages, truth, beta, rng)

    # Annotated variant set + g1/g2 seed indices in cache_var_ids
    annotated_vids = [v for v in variant_ids if v in cache]
    print(f"      annotated: {len(annotated_vids):,}")

    # ================ MODE 0: BASELINE (current; top_k=2000) ================
    print(f"\n  [Mode 0] Baseline mutual-RWR, top_k=2,000")
    t0 = time.time()
    g1_seeds_baseline = [k for k, v in enumerate(annotated_vids)
                          if gt["gene_1_systematic"] in cache[v].get("genes", [])]
    g2_seeds_baseline = [k for k, v in enumerate(annotated_vids)
                          if gt["gene_2_systematic"] in cache[v].get("genes", [])]
    seed_indices_all = sorted(set(g1_seeds_baseline + g2_seeds_baseline +
                                    list(range(len(annotated_vids)))))
    P_full, A_full, cache_var_ids_full = compute_rwr_matrix(
        cache, annotated_vids, seed_indices_all, alpha=0.15, n_iter=30,
    )
    # Re-derive g1/g2 indices in the seed_indices_all space
    g1_in_cache_var = {cache_var_ids_full[k] for k in g1_seeds_baseline}
    g2_in_cache_var = {cache_var_ids_full[k] for k in g2_seeds_baseline}

    baseline_pairs = run_baseline(P_full, seed_indices_all, top_k=2000)
    print(f"      pairs to test: {len(baseline_pairs):,}  (in {time.time()-t0:.1f}s)")
    t0 = time.time()
    baseline_res = test_pairs_and_rank(
        baseline_pairs, cache_var_ids_full, variant_ids, dosages, phenotype,
        gene1_variants, gene2_variants, mac_min,
    )
    print(f"      tested in {time.time()-t0:.1f}s; "
          f"ground_truth_rank={baseline_res['ground_truth_rank']} "
          f"of {baseline_res['n_tested']}")

    # ================ MODE A: HIGH-CAP top_k ================
    print(f"\n  [Mode A] Mutual-RWR, top_k={top_k_high:,}")
    t0 = time.time()
    high_cap_pairs = run_baseline(P_full, seed_indices_all, top_k=top_k_high)
    print(f"      pairs to test: {len(high_cap_pairs):,}  (RWR cached)")
    t0 = time.time()
    high_cap_res = test_pairs_and_rank(
        high_cap_pairs, cache_var_ids_full, variant_ids, dosages, phenotype,
        gene1_variants, gene2_variants, mac_min,
    )
    print(f"      tested in {time.time()-t0:.1f}s; "
          f"ground_truth_rank={high_cap_res['ground_truth_rank']} "
          f"of {high_cap_res['n_tested']}")

    # ================ MODE B: HUB-DEFLATED ================
    print(f"\n  [Mode B] Hub-deflated (P_a+P_b)/(deg_a*deg_b), top_k=2,000")
    t0 = time.time()
    hub_pairs = run_hub_deflated(P_full, A_full, seed_indices_all, top_k=2000)
    print(f"      pairs to test: {len(hub_pairs):,}")
    t0 = time.time()
    hub_res = test_pairs_and_rank(
        hub_pairs, cache_var_ids_full, variant_ids, dosages, phenotype,
        gene1_variants, gene2_variants, mac_min,
    )
    print(f"      tested in {time.time()-t0:.1f}s; "
          f"ground_truth_rank={hub_res['ground_truth_rank']} "
          f"of {hub_res['n_tested']}")

    # ================ MODE C: FOCUSED SEEDS (g1+g2 only) ================
    print(f"\n  [Mode C] Focused seeds (g1+g2 only); score g1×g2 cross-pairs")
    t0 = time.time()
    focused_pairs, _, _, cache_var_ids_focused = run_focused_seeds(
        cache, annotated_vids, g1_seeds_baseline, g2_seeds_baseline,
        alpha=0.15, n_iter=30,
    )
    print(f"      pairs to test: {len(focused_pairs):,}  (in {time.time()-t0:.1f}s)")
    t0 = time.time()
    focused_res = test_pairs_and_rank(
        focused_pairs, cache_var_ids_focused, variant_ids, dosages, phenotype,
        gene1_variants, gene2_variants, mac_min,
    )
    print(f"      tested in {time.time()-t0:.1f}s; "
          f"ground_truth_rank={focused_res['ground_truth_rank']} "
          f"of {focused_res['n_tested']}")

    # ================ MODE D: PPI-AUGMENTED (full-seeded) ================
    # Diagnostic for the ROOT cause: the bipartite walk can't reach across
    # gene components without a PPI bridge. Re-run baseline with the
    # ``use_ppi_bridge`` substrate that adds gene-gene PPI hops.
    print(f"\n  [Mode D] PPI-augmented substrate (var→gene→gene'→var), "
          f"all annotated as seeds, mutual scoring, top_k=2,000")
    t0 = time.time()
    P_ppi, A_ppi, cache_var_ids_ppi = compute_rwr_matrix(
        cache, annotated_vids, seed_indices_all, alpha=0.15, n_iter=30,
        use_ppi_bridge=True,
    )
    print(f"      RWR (PPI-bridged) computed in {time.time()-t0:.1f}s")
    t0 = time.time()
    ppi_pairs = run_baseline(P_ppi, seed_indices_all, top_k=2000)
    print(f"      pairs to test: {len(ppi_pairs):,}")
    t0 = time.time()
    ppi_res = test_pairs_and_rank(
        ppi_pairs, cache_var_ids_ppi, variant_ids, dosages, phenotype,
        gene1_variants, gene2_variants, mac_min,
    )
    print(f"      tested in {time.time()-t0:.1f}s; "
          f"ground_truth_rank={ppi_res['ground_truth_rank']} "
          f"of {ppi_res['n_tested']}")

    # ================ MODE E: PPI-AUGMENTED + FOCUSED SEEDS ================
    # Combination of Mode C + Mode D: seed only g1+g2 variants AND use the
    # PPI-bridged substrate. Should be the strongest mode if the canonical
    # PPI edge exists in the cache (after canonical-PPI injection).
    print(f"\n  [Mode E] PPI-augmented + focused seeds; score g1×g2 cross-pairs")
    t0 = time.time()
    e_pairs, _, _, cache_var_ids_e = run_focused_seeds(
        cache, annotated_vids, g1_seeds_baseline, g2_seeds_baseline,
        alpha=0.15, n_iter=30, use_ppi_bridge=True,
    )
    print(f"      pairs to test: {len(e_pairs):,}  (in {time.time()-t0:.1f}s)")
    t0 = time.time()
    e_res = test_pairs_and_rank(
        e_pairs, cache_var_ids_e, variant_ids, dosages, phenotype,
        gene1_variants, gene2_variants, mac_min,
    )
    print(f"      tested in {time.time()-t0:.1f}s; "
          f"ground_truth_rank={e_res['ground_truth_rank']} "
          f"of {e_res['n_tested']}")

    return {
        "species": species,
        "n_samples": int(dosages.shape[0]),
        "n_variants": int(dosages.shape[1]),
        "n_annotated": len(annotated_vids),
        "g1_count": len(gene1_variants),
        "g2_count": len(gene2_variants),
        "modes": {
            "0_baseline":            baseline_res,
            "A_high_cap":            high_cap_res,
            "B_hub_deflated":        hub_res,
            "C_focused_seeds":       focused_res,
            "D_ppi_bridge":          ppi_res,
            "E_ppi_bridge_focused":  e_res,
        },
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--species", choices=["yeast", "arabidopsis", "rice", "all"],
                   default="all")
    p.add_argument("--beta", type=float, default=3.0)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--top-k-high", type=int, default=50_000)
    args = p.parse_args()

    if args.species == "all":
        species_list = ["yeast", "arabidopsis", "rice"]
    else:
        species_list = [args.species]

    all_results: dict = {}
    for sp in species_list:
        all_results[sp] = run_species(sp, args.beta, args.seed,
                                        top_k_high=args.top_k_high)

    # Final summary
    print(f"\n{'=' * 70}\n  M5 VARIANT BENCHMARK SUMMARY\n{'=' * 70}\n")
    modes = ["0_baseline", "A_high_cap", "B_hub_deflated", "C_focused_seeds",
             "D_ppi_bridge", "E_ppi_bridge_focused"]
    print(f"  {'Species':<14}  " + "  ".join(f"{m:>17}" for m in modes))
    print(f"  {'-' * 14}  " + "  ".join("-" * 17 for _ in modes))
    for sp, sp_data in all_results.items():
        cells = []
        for m in modes:
            r = sp_data["modes"][m]
            rk = r.get("ground_truth_rank")
            n = r.get("n_tested", 0)
            top1pct = r.get("top1pct_threshold", 1)
            in_top1 = rk is not None and rk <= top1pct
            mark = "✓" if in_top1 else (" " if rk is None else "·")
            cell = f"{mark}{rk}/{n}" if rk is not None else f" NF/{n}"
            cells.append(cell)
        print(f"  {sp:<14}  " + "  ".join(f"{c:>17}" for c in cells))
    print(f"\n  ✓ = rank in top 1% of tested pairs (winning configuration)\n")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out = RESULTS_DIR / "m5_variants_benchmark.json"
    out.write_text(json.dumps(all_results, indent=2))
    print(f"  Wrote {out}")


if __name__ == "__main__":
    main()
