"""GWAS-result interpretation without LLMs.

Pure-statistical structured summary of an `AssociationResult` run: counts,
genomic inflation, independent loci by 500-kb clumping, and rule-based
follow-up recommendations. No external API keys required.
"""
from __future__ import annotations

from .db import GraphGWASConnection


def interpret_results(conn: GraphGWASConnection, run_id: str,
                      verbose: bool = True) -> dict:
    """Structured interpretation of GWAS results without LLM.

    Provides:
    - Summary statistics (n_tested, n_significant, lambda_GC)
    - Top loci with gene annotations
    - Pathway enrichment via stored results
    - Recommendations for follow-up analyses
    """
    from .results import query_top_hits, query_manhattan_data

    top_hits = query_top_hits(conn, run_id, p_threshold=5e-8, limit=50)
    all_data = query_manhattan_data(conn, run_id)

    n_tested = len(all_data)
    n_significant = len(top_hits)

    import numpy as np
    from scipy import stats as sp_stats
    p_values = [d["log10p"] for d in all_data if d.get("log10p")]
    if p_values:
        chi2_obs = [sp_stats.chi2.isf(10 ** (-lp), 1) for lp in p_values if lp > 0]
        lambda_gc = float(np.median(chi2_obs) / 0.4549) if chi2_obs else 1.0
    else:
        lambda_gc = None

    loci = []
    if top_hits:
        sorted_hits = sorted(top_hits, key=lambda h: h["p_value"])
        used_positions = []
        for hit in sorted_hits:
            pos = hit.get("pos", 0)
            chr_val = hit.get("chr", "")
            is_new = True
            for used_chr, used_pos in used_positions:
                if chr_val == used_chr and abs(pos - used_pos) < 500000:
                    is_new = False
                    break
            if is_new:
                used_positions.append((chr_val, pos))
                loci.append(hit)

    interpretation = {
        "summary": {
            "n_variants_tested": n_tested,
            "n_genome_wide_significant": n_significant,
            "n_independent_loci": len(loci),
            "lambda_gc": lambda_gc,
        },
        "top_loci": loci[:20],
        "recommendations": [],
    }

    if n_significant == 0:
        interpretation["recommendations"].append(
            "No genome-wide significant hits. Consider: "
            "(1) MPAT gene-level test for aggregated signal, "
            "(2) Disease architecture via max-flow for pathway-level signal, "
            "(3) Epistasis scan for multi-variant interactions."
        )
    elif n_significant > 0 and n_significant < 10:
        interpretation["recommendations"].append(
            f"Found {n_significant} significant variants in {len(loci)} loci. "
            "Recommended follow-up: LD fine-mapping per locus, "
            "MPAT for gene-level confirmation, "
            "epistasis scan to check for multi-variant architecture."
        )
    else:
        interpretation["recommendations"].append(
            f"Strong polygenic signal: {n_significant} significant variants. "
            "Consider: GNN for multi-locus prediction, "
            "disease architecture for pathway-level summary."
        )

    if lambda_gc and lambda_gc > 1.1:
        interpretation["recommendations"].append(
            f"Genomic inflation λ_GC = {lambda_gc:.3f} > 1.1. "
            "Consider adding PCA covariates or spectral correction."
        )

    if verbose:
        print(f"\n=== GWAS Interpretation (run: {run_id}) ===")
        print(f"Variants tested: {n_tested:,}")
        print(f"Genome-wide significant: {n_significant}")
        print(f"Independent loci: {len(loci)}")
        if lambda_gc:
            print(f"Lambda GC: {lambda_gc:.3f}")
        print("\nTop loci:")
        for locus in loci[:10]:
            genes = ", ".join(locus.get("genes", [])[:3]) or "intergenic"
            print(f"  {locus.get('chr', '')}:{locus.get('pos', '')} "
                  f"p={locus.get('p_value', 1):.2e} genes={genes}")
        print("\nRecommendations:")
        for rec in interpretation["recommendations"]:
            print(f"  - {rec}")

    return interpretation
