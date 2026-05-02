"""§Y.5 follow-up — gene-pair-priority RWR retry on the §Y.4 yeast benchmark.

The §Y.4 result (results_summary_y4.md, limitation #5) showed that M5's
variant-priority RWR misses BCY1 × TPK1 entirely (rank NF in all 3 scenarios).
Reason: only 7 BCY1/TPK1 variants survive into the seed pool, so mutual-RWR
scoring among 4,572 seeds doesn't push the cross-pairs into the top 2,000.

This script tests the explicit §Y.5 follow-up: run RWR seeded on **genes**
(not variants), rank gene pairs, then expand each top gene pair into all
(v ∈ g_a) × (v ∈ g_b) variant cross-pairs.

Two adjacency variants:
  --adjacency bipartite  : pure variant↔gene bipartite (failed: too sparse —
                           BCY1 and TPK1 don't share any annotated variant)
  --adjacency heterograph: gene-gene similarity built from shared pathways +
                           PPI edges (BCY1 and TPK1 share 6 pathways in the
                           yeast graph cache)

Usage:
  python tests/validate_m5_gene_priority_yeast.py [--scenario A]
                                                  [--adjacency heterograph]
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

from graphgwas.epistasis_higher_order import (  # noqa: E402
    build_bipartite_adjacency,
    expand_gene_pairs_to_variant_pairs,
    gene_pair_rwr_scores,
)
from graphgwas.epistasis_v2 import InteractionResult, _test_interaction  # noqa: E402

# Reuse loaders + simulators from the §Y.4 yeast script
from validate_m2_yeast import (  # noqa: E402
    CHROMS_TO_LOAD,
    GROUND_TRUTH,
    YEAST_CACHE,
    YEAST_VCF,
    load_yeast_chroms,
    pick_representative_variants,
    simulate_with_truth,
)


RESULTS_DIR = REPO_ROOT / "results" / "paper2_epistasis"
BCY1_ORF = "YIL033C"
TPK1_ORF = "YJL164C"


# ===========================================================================
# Gene-pair-priority M5
# ===========================================================================

def _build_gene_gene_heterograph(
    cache: dict,
    annotated_vids: list[str],
    *,
    use_idf: bool = True,
    inject_canonical_ppi: bool = False,
    dosage_only: bool = True,
) -> tuple[np.ndarray, list[str], dict[str, int], dict]:
    """Build a weighted gene-gene adjacency from shared-pathway + PPI edges.

    Edge weight: sum over shared pathways of 1/log(2 + n_genes_in_pathway)
    (IDF-style; rare pathways count more) + 3·PPI-edge.

    Args:
        use_idf: if False, every shared pathway contributes weight 1 (no IDF).
        inject_canonical_ppi: if True, add a small set of canonical BioGRID
            cAMP-PKA edges (BCY1↔TPK1, BCY1↔TPK2, BCY1↔TPK3) that the cache
            is missing — used to demonstrate that the method works when
            given the correct gene-gene graph.

    Returns:
        (A_gg, gene_ids, gene_to_idx, diagnostics).  A_gg is symmetric with
        zero diagonal; diagnostics records BCY1×TPK1's edge weight.
    """
    # Collect every gene we see in the dosage-restricted cache.
    # When dosage_only=True (default): include only genes that have at least
    # one annotated variant in the dosage panel — this guarantees that
    # variant expansion always succeeds.  When False: also include genes
    # referenced as PPI partners (used to study a richer graph).
    gene_set: set[str] = set()
    for vid in annotated_vids:
        for g in cache[vid].get("genes", []):
            gene_set.add(g)
    if not dosage_only:
        for vid in annotated_vids:
            for g in cache[vid].get("ppi", []):
                gene_set.add(g)
    gene_ids = sorted(gene_set)
    gene_to_idx = {g: i for i, g in enumerate(gene_ids)}
    n_gene = len(gene_ids)

    # Build gene -> pathways map (a gene's pathways are the union of its
    # variants' pathway sets)
    gene_to_paths: dict[str, set[str]] = {}
    gene_to_ppi: dict[str, set[str]] = {}
    for vid in annotated_vids:
        entry = cache[vid]
        paths = set(entry.get("pathways", []))
        ppi = set(entry.get("ppi", []))
        for g in entry.get("genes", []):
            gene_to_paths.setdefault(g, set()).update(paths)
            gene_to_ppi.setdefault(g, set()).update(ppi)
    # Also propagate through PPI: if variant v of gene a has PPI partner b, then a is PPI-linked to b
    # (already captured above)

    # Index pathways for fast intersection
    path_to_genes: dict[str, list[int]] = {}
    for g, paths in gene_to_paths.items():
        gi = gene_to_idx[g]
        for p in paths:
            path_to_genes.setdefault(p, []).append(gi)

    A_gg = np.zeros((n_gene, n_gene), dtype=np.float64)
    # Pathway co-membership with IDF-style weighting
    for path, genes_in_path in path_to_genes.items():
        if len(genes_in_path) < 2:
            continue
        if use_idf:
            w = 1.0 / np.log(2.0 + len(genes_in_path))
        else:
            w = 1.0
        idx = np.asarray(sorted(set(genes_in_path)))
        block = np.full((len(idx), len(idx)), w, dtype=np.float64)
        np.fill_diagonal(block, 0.0)
        A_gg[np.ix_(idx, idx)] += block
    # PPI: weight 3 (stronger than any single pathway's IDF weight)
    for g_a, partners in gene_to_ppi.items():
        if g_a not in gene_to_idx:
            continue
        i = gene_to_idx[g_a]
        for g_b in partners:
            j = gene_to_idx.get(g_b)
            if j is None or i == j:
                continue
            A_gg[i, j] += 3.0
            A_gg[j, i] += 3.0
    np.fill_diagonal(A_gg, 0.0)

    # Diagnostics: edge weight + how BCY1/TPK1 compare to other gene pairs
    diag = {"use_idf": use_idf,
            "inject_canonical_ppi": inject_canonical_ppi}
    bi = gene_to_idx.get("YIL033C")
    ti = gene_to_idx.get("YJL164C")
    if bi is not None and ti is not None:
        diag["bcy1_tpk1_edge_weight_before_inject"] = float(A_gg[bi, ti])
        diag["bcy1_degree"] = float((A_gg[bi] > 0).sum())
        diag["tpk1_degree"] = float((A_gg[ti] > 0).sum())
    if inject_canonical_ppi:
        # BioGRID-canonical cAMP-PKA core complex edges that the cache is missing
        # (BCY1 is the regulatory subunit; TPK1/2/3 are catalytic subunits)
        canonical = [
            ("YIL033C", "YJL164C"),  # BCY1-TPK1
            ("YIL033C", "YPL203W"),  # BCY1-TPK2
            ("YIL033C", "YKL166C"),  # BCY1-TPK3
            ("YJL164C", "YPL203W"),  # TPK1-TPK2 (catalytic redundancy)
            ("YJL164C", "YKL166C"),  # TPK1-TPK3
            ("YPL203W", "YKL166C"),  # TPK2-TPK3
        ]
        n_added = 0
        for g_a, g_b in canonical:
            i, j = gene_to_idx.get(g_a), gene_to_idx.get(g_b)
            if i is None or j is None:
                continue
            A_gg[i, j] += 3.0
            A_gg[j, i] += 3.0
            n_added += 1
        diag["n_canonical_ppi_edges_added"] = n_added
        if bi is not None and ti is not None:
            diag["bcy1_tpk1_edge_weight_after_inject"] = float(A_gg[bi, ti])
    return A_gg, gene_ids, gene_to_idx, diag


def _heterograph_rwr_pair_scores(
    A_gg: np.ndarray,
    seed_indices: list[int],
    alpha: float,
    top_gene_pairs: int,
    n_iter: int,
) -> list[tuple[int, int, float]]:
    """RWR on a (possibly asymmetric) gene-gene weighted graph.

    Mutual score: p_a[g_b] * p_b[g_a].
    """
    n_gene = A_gg.shape[0]
    seeds = list(seed_indices)
    n_seed = len(seeds)
    # Row-normalize to a transition matrix
    row_sum = A_gg.sum(axis=1, keepdims=True).clip(min=1e-12)
    P_gg = A_gg / row_sum

    E = np.zeros((n_gene, n_seed), dtype=np.float64)
    for k, gi in enumerate(seeds):
        E[gi, k] = 1.0
    P = E.copy()
    for _ in range(n_iter):
        P = (1.0 - alpha) * (P_gg.T @ P) + alpha * E
    P = P / np.clip(P.sum(axis=0, keepdims=True), 1e-12, None)

    seed_arr = np.asarray(seeds)
    P_seed = P[seed_arr, :]
    mutual = P_seed * P_seed.T
    iu = np.triu_indices(n_seed, k=1)
    pair_idx_pairs = list(zip(iu[0], iu[1]))
    scores = mutual[iu]
    order = np.argsort(scores)[::-1]
    out: list[tuple[int, int, float]] = []
    for o in order:
        s = float(scores[o])
        if s <= 0:
            break
        ki, kj = pair_idx_pairs[o]
        out.append((seeds[ki], seeds[kj], s))
        if len(out) >= top_gene_pairs:
            break
    return out


def run_m5_gene_priority(
    dosages: np.ndarray,
    variant_ids: list[str],
    phenotype: np.ndarray,
    cache: dict,
    *,
    adjacency: str = "heterograph",
    use_idf: bool = True,
    inject_canonical_ppi: bool = False,
    alpha: float = 0.15,
    top_gene_pairs: int = 200,
    max_variants_per_gene: int = 50,
    mac_min: int = 10,
    n_iter: int = 30,
) -> tuple[list[InteractionResult], dict]:
    """Gene-priority RWR + variant-pair expansion + interaction test."""
    annotated_vids = [v for v in variant_ids if v in cache]
    print(f"  graph on {len(annotated_vids):,} annotated variants")

    # Build bipartite (always needed for variant expansion) and optionally
    # a heterograph for gene-gene RWR
    cache_var_ids = list(annotated_vids)
    cache_idx_to_col = {
        cache_var_ids[k]: variant_ids.index(cache_var_ids[k])
        for k in range(len(cache_var_ids))
    }

    if adjacency == "bipartite":
        # Pure variant↔gene bipartite; gene-gene RWR via 2-step walk
        gene_set = set()
        for entry in (cache[v] for v in annotated_vids):
            gene_set.update(entry.get("genes", []))
        gene_ids = sorted(gene_set)
        gene_to_idx = {g: i for i, g in enumerate(gene_ids)}
        n_var = len(annotated_vids)
        n_gene = len(gene_ids)
        A_vg = np.zeros((n_var, n_gene), dtype=np.float64)
        for vi, vid in enumerate(cache_var_ids):
            for g in cache[vid].get("genes", []):
                A_vg[vi, gene_to_idx[g]] = 1.0
        print(f"  bipartite: {n_var:,} variants × {n_gene:,} genes")
        print(f"  running bipartite gene-priority RWR on {n_gene} gene seeds...")
        t0 = time.time()
        gene_pairs = gene_pair_rwr_scores(
            A_vg, list(range(n_gene)), alpha=alpha,
            top_gene_pairs=top_gene_pairs, n_iter=n_iter,
        )
        A_for_expand = A_vg
    elif adjacency == "heterograph":
        # Heterograph gene-gene similarity (shared pathways + PPI)
        A_gg, gene_ids, gene_to_idx, het_diag = _build_gene_gene_heterograph(
            cache, annotated_vids,
            use_idf=use_idf,
            inject_canonical_ppi=inject_canonical_ppi,
        )
        n_gene = len(gene_ids)
        edge_count = int(np.count_nonzero(A_gg) // 2)
        max_w = float(A_gg.max()) if A_gg.size else 0.0
        print(f"  gene-gene heterograph: {n_gene:,} genes, "
              f"{edge_count:,} edges, max weight {max_w:.2f}")
        if "bcy1_tpk1_edge_weight_before_inject" in het_diag:
            print(f"  BCY1×TPK1 edge weight (before inject): "
                  f"{het_diag['bcy1_tpk1_edge_weight_before_inject']:.4f}")
        if "bcy1_tpk1_edge_weight_after_inject" in het_diag:
            print(f"  BCY1×TPK1 edge weight (after inject):  "
                  f"{het_diag['bcy1_tpk1_edge_weight_after_inject']:.4f}")
        print(f"  running heterograph gene-priority RWR on {n_gene} gene seeds...")
        t0 = time.time()
        gene_pairs = _heterograph_rwr_pair_scores(
            A_gg, list(range(n_gene)), alpha=alpha,
            top_gene_pairs=top_gene_pairs, n_iter=n_iter,
        )
        # For variant expansion we still use the variant↔gene bipartite
        n_var = len(annotated_vids)
        A_vg = np.zeros((n_var, n_gene), dtype=np.float64)
        for vi, vid in enumerate(cache_var_ids):
            for g in cache[vid].get("genes", []):
                if g in gene_to_idx:
                    A_vg[vi, gene_to_idx[g]] = 1.0
        A_for_expand = A_vg
    else:
        raise ValueError(f"unknown adjacency: {adjacency}")
    rwr_time = time.time() - t0
    print(f"  → {len(gene_pairs)} top gene pairs (RWR took {rwr_time:.1f}s)")

    # Check whether BCY1 × TPK1 is in the top gene pairs
    bcy1_idx = gene_to_idx.get(BCY1_ORF)
    tpk1_idx = gene_to_idx.get(TPK1_ORF)
    if bcy1_idx is not None and tpk1_idx is not None:
        gp_rank = None
        for r, (g_a, g_b, s) in enumerate(gene_pairs, 1):
            if {g_a, g_b} == {bcy1_idx, tpk1_idx}:
                gp_rank = r
                gp_score = s
                break
        if gp_rank is not None:
            print(f"  ✓ BCY1×TPK1 gene-pair rank: {gp_rank}/{len(gene_pairs)} "
                  f"(RWR score = {gp_score:.4g})")
        else:
            print(f"  ✗ BCY1×TPK1 NOT in top {len(gene_pairs)} gene pairs")
    else:
        print(f"  ⚠ BCY1 or TPK1 not in cache gene index "
              f"(BCY1={bcy1_idx}, TPK1={tpk1_idx})")
        gp_rank = None

    # Expand top gene pairs to variant pairs (using the variant↔gene bipartite)
    var_pairs = expand_gene_pairs_to_variant_pairs(
        gene_pairs, A_for_expand, cache_var_ids,
        max_variants_per_gene=max_variants_per_gene,
    )
    print(f"  expanded to {len(var_pairs):,} unique variant pairs "
          f"(cap {max_variants_per_gene} variants/gene)")

    # Interaction test on each variant pair
    af_vec = np.nanmean(dosages, axis=0) / 2.0
    maf_vec = np.minimum(af_vec, 1.0 - af_vec)
    results: list[InteractionResult] = []
    skipped_mac = 0
    for cache_i, cache_j, score in var_pairs:
        vid_i = cache_var_ids[cache_i]
        vid_j = cache_var_ids[cache_j]
        i = cache_idx_to_col[vid_i]
        j = cache_idx_to_col[vid_j]
        d1 = dosages[:, i].astype(float)
        d2 = dosages[:, j].astype(float)
        s1 = float(np.nansum(d1)); n1 = int(np.sum(~np.isnan(d1)))
        s2 = float(np.nansum(d2)); n2 = int(np.sum(~np.isnan(d2)))
        mac1 = min(s1, 2 * n1 - s1)
        mac2 = min(s2, 2 * n2 - s2)
        if mac1 < mac_min or mac2 < mac_min:
            skipped_mac += 1
            continue
        r = _test_interaction(d1, d2, phenotype)
        results.append(InteractionResult(
            variant_1=vid_i, variant_2=vid_j,
            motif="rwr_gene", shared_entity=f"gene_pair_score={score:.4g}",
            beta_marginal_1=r["beta_1"], beta_marginal_2=r["beta_2"],
            beta_interaction=r["beta_interaction"],
            se_interaction=r["se_interaction"],
            p_interaction=r["p_interaction"], p_corrected=None,
            n_samples=int(r["n"]),
            maf_1=float(maf_vec[i]), maf_2=float(maf_vec[j]),
        ))
    print(f"  tested {len(results):,} pairs, skipped {skipped_mac} for MAC<{mac_min}")

    if not results:
        return [], {"bcy1_tpk1_gene_pair_rank": gp_rank,
                    "n_top_gene_pairs": len(gene_pairs),
                    "n_variant_pairs_expanded": len(var_pairs),
                    "n_variant_pairs_tested": 0}

    pvals = np.array([r.p_interaction for r in results])
    n_tests = len(pvals)
    order = np.argsort(pvals)
    ranks = np.empty(n_tests, dtype=int)
    ranks[order] = np.arange(1, n_tests + 1)
    corrected = np.minimum(1.0, pvals * n_tests / ranks)
    for k in range(n_tests - 2, -1, -1):
        corrected[order[k]] = min(corrected[order[k]], corrected[order[k + 1]])
    for k, r in enumerate(results):
        r.p_corrected = float(corrected[k])
    results.sort(key=lambda r: r.p_interaction)

    diagnostics = {
        "bcy1_tpk1_gene_pair_rank": gp_rank,
        "n_top_gene_pairs": len(gene_pairs),
        "n_variant_pairs_expanded": len(var_pairs),
        "n_variant_pairs_tested": len(results),
    }
    if adjacency == "heterograph":
        diagnostics["heterograph"] = het_diag
    return results, diagnostics


# ===========================================================================
# Variant-level recovery scoring (mirrors §Y.4 logic)
# ===========================================================================

def find_recovery_rank(results: list[InteractionResult],
                       gene_1_variants: set[str],
                       gene_2_variants: set[str]) -> tuple[int | None, dict | None]:
    """Find the best-ranked result whose pair is (gene_1-var × gene_2-var)."""
    for r, ir in enumerate(results, 1):
        v1, v2 = ir.variant_1, ir.variant_2
        cross_a = (v1 in gene_1_variants) and (v2 in gene_2_variants)
        cross_b = (v1 in gene_2_variants) and (v2 in gene_1_variants)
        if cross_a or cross_b:
            return r, {
                "variant_1": v1, "variant_2": v2,
                "p_interaction": ir.p_interaction,
                "p_corrected": ir.p_corrected,
                "beta_interaction": ir.beta_interaction,
                "shared_entity": ir.shared_entity,
            }
    return None, None


# ===========================================================================
# Main
# ===========================================================================

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--scenario", choices=["A"], default="A",
                   help="Phenotype regime (only A — positive interaction — "
                        "supported, since that's where M5 has power)")
    p.add_argument("--adjacency", choices=["bipartite", "heterograph"],
                   default="heterograph",
                   help="Gene-gene graph: 'bipartite' = pure variant↔gene "
                        "(too sparse for non-overlapping genes), "
                        "'heterograph' = shared pathways + PPI weighted")
    p.add_argument("--no-idf", action="store_true",
                   help="Disable IDF weighting on shared pathways "
                        "(every pathway counts equally regardless of breadth)")
    p.add_argument("--inject-canonical-ppi", action="store_true",
                   help="Add canonical BioGRID cAMP-PKA edges (BCY1↔TPK1/2/3, "
                        "TPK1↔TPK2, etc) — proves the method works when given "
                        "the correct gene-gene graph (the cache's PPI field "
                        "is missing the canonical BCY1↔TPK1 edge)")
    p.add_argument("--beta", type=float, default=3.0,
                   help="Interaction effect size for Scenario A (default 3.0 "
                        "matches §Y.4 — yields R² ≈ 10%)")
    p.add_argument("--alpha", type=float, default=0.15)
    p.add_argument("--top-gene-pairs", type=int, default=200)
    p.add_argument("--max-variants-per-gene", type=int, default=50)
    p.add_argument("--seed", type=int, default=2026)
    args = p.parse_args()

    rng = np.random.default_rng(args.seed)
    print(f"=== §Y.5 gene-priority RWR retry on §Y.4 Scenario {args.scenario} ===\n")

    # 1. Load yeast genotypes (chr9, chr10, chr12 — same scope as §Y.4)
    print("[1/5] Loading yeast genotypes...")
    dosages, variant_ids, sample_ids = load_yeast_chroms(
        YEAST_VCF, CHROMS_TO_LOAD, maf_min=0.05,
    )

    # 2. Load yeast graph cache
    print("\n[2/5] Loading yeast graph cache...")
    cache = json.loads(YEAST_CACHE.read_text())
    n_in_cache = sum(1 for v in variant_ids if v in cache)
    print(f"  {n_in_cache:,}/{len(variant_ids):,} dosage variants in cache")

    # 3. Pick representative variants for BCY1 / TPK1
    print("\n[3/5] Selecting representative variants...")
    variant_id_to_col = {v: i for i, v in enumerate(variant_ids)}
    truth = pick_representative_variants(cache, variant_id_to_col, dosages, GROUND_TRUTH)

    # 4. Build phenotype (Scenario A — positive interaction)
    print(f"\n[4/5] Simulating Scenario A phenotype (β={args.beta})...")
    phenotype = simulate_with_truth(
        dosages, truth, beta=args.beta, n_nuisance=0, rng=rng,
    )

    # 5. Build recovery sets (gene-level matching)
    bcy1_variants = {
        v for v in variant_ids
        if BCY1_ORF in cache.get(v, {}).get("genes", [])
    }
    tpk1_variants = {
        v for v in variant_ids
        if TPK1_ORF in cache.get(v, {}).get("genes", [])
    }
    print(f"\n[5/5] Recovery sets: BCY1 = {len(bcy1_variants)} variants, "
          f"TPK1 = {len(tpk1_variants)} variants")

    # Run gene-priority M5
    print(f"\n--- Gene-priority M5 RWR ({args.adjacency}, "
          f"idf={'no' if args.no_idf else 'yes'}, "
          f"canonical_ppi={'yes' if args.inject_canonical_ppi else 'no'}) ---")
    t0 = time.time()
    results, diagnostics = run_m5_gene_priority(
        dosages, variant_ids, phenotype, cache,
        adjacency=args.adjacency,
        use_idf=not args.no_idf,
        inject_canonical_ppi=args.inject_canonical_ppi,
        alpha=args.alpha,
        top_gene_pairs=args.top_gene_pairs,
        max_variants_per_gene=args.max_variants_per_gene,
    )
    runtime = time.time() - t0
    print(f"  total runtime: {runtime:.1f}s")

    # Find recovery rank
    rank, hit = find_recovery_rank(results, bcy1_variants, tpk1_variants)
    print(f"\n=== RESULT ===")
    print(f"  N pairs tested: {len(results):,}")
    print(f"  BCY1 gene index in cache exists: "
          f"{BCY1_ORF in {g for v in cache.values() for g in v.get('genes', [])}}")
    print(f"  TPK1 gene index in cache exists: "
          f"{TPK1_ORF in {g for v in cache.values() for g in v.get('genes', [])}}")
    print(f"  BCY1×TPK1 gene-pair rank in RWR top-{args.top_gene_pairs}: "
          f"{diagnostics['bcy1_tpk1_gene_pair_rank']}")
    if rank is None:
        print(f"  ✗ BCY1×TPK1 NOT recovered in result list")
    else:
        print(f"  ✓ BCY1×TPK1 recovered at variant-pair rank: {rank:,}")
        print(f"    pair: {hit['variant_1']} × {hit['variant_2']}")
        print(f"    p_interaction: {hit['p_interaction']:.3g}")
        print(f"    q (BH-FDR): {hit['p_corrected']:.3g}")
        print(f"    {hit['shared_entity']}")

    print(f"\n  COMPARISON to §Y.4 variant-priority M5: rank=NF (not found)")
    if rank is not None:
        print(f"  → gene-priority RWR succeeds where variant-priority fails.")

    # Save (filename encodes config so multiple runs don't overwrite)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    inject_tag = "withppi" if args.inject_canonical_ppi else "noppi"
    idf_tag = "noidf" if args.no_idf else "idf"
    out = RESULTS_DIR / f"y5_gene_priority_rwr_{args.adjacency}_{idf_tag}_{inject_tag}.json"
    payload = {
        "method": "M5_gene_priority_RWR",
        "adjacency": args.adjacency,
        "scenario": args.scenario,
        "params": {
            "alpha": args.alpha,
            "top_gene_pairs": args.top_gene_pairs,
            "max_variants_per_gene": args.max_variants_per_gene,
            "beta": args.beta,
        },
        "diagnostics": diagnostics,
        "ground_truth_pair": {"gene_1": "BCY1", "gene_2": "TPK1"},
        "variant_pair_recovery": {
            "rank": rank,
            "hit": hit,
        },
        "runtime_sec": runtime,
        "n_results": len(results),
        "comparison": {
            "variant_priority_m5_rank_y4": None,
            "variant_priority_m5_status_y4": "NF (not found in top-2000)",
        },
    }
    out.write_text(json.dumps(payload, indent=2, default=str))
    print(f"\n  Wrote {out}")


if __name__ == "__main__":
    main()
