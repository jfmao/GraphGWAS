"""Graph-native Polygenic Risk Score (PRS).

Classical PRS is a weighted sum: PRS_i = Σ_j β_j × x_ij.  It treats each
variant independently, ignoring LD structure, gene membership, and pathway
biology.  Graph-native PRS uses the variant→gene→pathway topology to produce
biologically informed, LD-aware, and pathway-decomposable risk scores.

Components:
    P1. classical_prs           — standard weighted-sum PRS (baseline)
    P2. graph_pruned_prs        — LD-aware pruning via graph centrality
    P3. pathway_partitioned_prs — PRS decomposed by biological pathway
    P4. message_passing_prs     — GNN-learned PRS from graph structure
    P5. prs_evaluation          — AUROC, Nagelkerke R², calibration
    P6. prs_report              — integrated PRS comparison report
"""

from __future__ import annotations

import numpy as np
from scipy import stats as sp_stats

from . import config as _cfg
from .db import GraphGWASConnection
from .genotype import (
    build_dosage,
    get_phenotype_indices,
    variant_iterator,
)


# ===================================================================
# P1. Classical PRS (baseline)
# ===================================================================

def classical_prs(conn: GraphGWASConnection,
                  run_id: str,
                  p_threshold: float = 5e-8,
                  verbose: bool = True) -> dict:
    """Standard weighted-sum PRS from stored GWAS results.

    PRS_i = Σ_j β_j × dosage_ij  for all variants with p < threshold.

    This is the baseline against which graph-native methods are compared.

    Args:
        conn: database connection.
        run_id: GWAS run ID for effect sizes.
        p_threshold: p-value threshold for variant inclusion.
        verbose: print progress.

    Returns:
        dict with prs_scores (per sample), n_variants, variant_ids.
    """
    # Get significant variants with betas
    result = conn.execute_read(
        """
        MATCH (ar:AssociationResult)-[:FOR_VARIANT]->(v:Variant)
        WHERE ar.run_id = $run_id AND ar.p_value < $pval
              AND ar.beta IS NOT NULL
        RETURN v.variantId AS vid, ar.beta AS beta, ar.p_value AS p,
               v.gt_packed AS gtp
        ORDER BY ar.p_value
        """,
        {"run_id": run_id, "pval": p_threshold},
    )

    records = [dict(r) for r in result]
    if not records:
        if verbose:
            print(f"No significant variants found (p < {p_threshold:.0e})")
        return {"prs_scores": [], "n_variants": 0, "method": "classical"}

    case_idx, ctrl_idx = get_phenotype_indices(conn)
    all_idx = np.concatenate([case_idx, ctrl_idx])
    n_samples = len(all_idx)

    prs = np.zeros(n_samples, dtype=np.float64)
    variant_ids = []
    betas = []

    for rec in records:
        gtp = rec["gtp"]
        if gtp is None:
            continue
        beta = rec["beta"]
        dosage = build_dosage(gtp, all_idx, _cfg.N_SAMPLES)
        # Replace NaN with 0 (missing genotypes don't contribute)
        dosage = np.nan_to_num(dosage, nan=0.0)
        prs += beta * dosage
        variant_ids.append(rec["vid"])
        betas.append(beta)

    if verbose:
        print(f"=== Classical PRS ===")
        print(f"  Variants: {len(variant_ids)} (p < {p_threshold:.0e})")
        print(f"  PRS range: [{prs.min():.4f}, {prs.max():.4f}]")
        print(f"  PRS mean(cases): {prs[:len(case_idx)].mean():.4f}, "
              f"mean(controls): {prs[len(case_idx):].mean():.4f}")

    return {
        "method": "classical",
        "prs_scores": prs.tolist(),
        "n_variants": len(variant_ids),
        "variant_ids": variant_ids,
        "betas": betas,
        "n_cases": len(case_idx),
        "n_controls": len(ctrl_idx),
        "sample_indices": all_idx.tolist(),
    }


# ===================================================================
# P2. Graph-Pruned PRS (LD-aware via centrality)
# ===================================================================

def graph_pruned_prs(conn: GraphGWASConnection,
                     run_id: str,
                     p_threshold: float = 5e-8,
                     r2_threshold: float = 0.2,
                     pruning_method: str = "centrality",
                     verbose: bool = True) -> dict:
    """LD-aware PRS using graph centrality for variant selection.

    Instead of classical clumping (greedy, order-dependent), builds the
    LD graph among significant variants and selects representatives by
    graph centrality.  The highest-centrality variant in each LD clump
    is the best tag — it captures the most information from its neighborhood.

    This is provably better than greedy clumping when LD structure is
    complex (multi-way LD, overlapping haplotype blocks).

    Args:
        conn: database connection.
        run_id: GWAS run ID.
        p_threshold: p-value threshold.
        r2_threshold: LD r² threshold for graph edges.
        pruning_method: 'centrality' (betweenness), 'pagerank', 'independent_set'.
        verbose: print progress.

    Returns:
        dict with prs_scores, selected variants, pruning stats.
    """
    # Get significant variants
    result = conn.execute_read(
        """
        MATCH (ar:AssociationResult)-[:FOR_VARIANT]->(v:Variant)
        WHERE ar.run_id = $run_id AND ar.p_value < $pval
              AND ar.beta IS NOT NULL
        RETURN v.variantId AS vid, v.chr AS chr, v.pos AS pos,
               ar.beta AS beta, ar.p_value AS p, v.gt_packed AS gtp
        ORDER BY ar.p_value
        """,
        {"run_id": run_id, "pval": p_threshold},
    )

    records = [dict(r) for r in result]
    if not records:
        return {"prs_scores": [], "n_variants": 0, "method": "graph_pruned"}

    case_idx, ctrl_idx = get_phenotype_indices(conn)
    all_idx = np.concatenate([case_idx, ctrl_idx])

    # Build LD graph among significant variants
    import networkx as nx

    # Pre-compute dosages for LD calculation
    dosages = {}
    for rec in records:
        if rec["gtp"] is not None:
            dosages[rec["vid"]] = build_dosage(rec["gtp"], all_idx, _cfg.N_SAMPLES)

    G = nx.Graph()
    vid_to_rec = {rec["vid"]: rec for rec in records if rec["vid"] in dosages}
    for vid in dosages:
        G.add_node(vid, beta=vid_to_rec[vid]["beta"],
                   p_value=vid_to_rec[vid]["p"])

    # Compute pairwise LD for variants on the same chromosome within 1Mb
    vids = list(dosages.keys())
    for i in range(len(vids)):
        for j in range(i + 1, len(vids)):
            v1, v2 = vids[i], vids[j]
            r1, r2 = vid_to_rec[v1], vid_to_rec[v2]
            if r1["chr"] != r2["chr"]:
                continue
            if abs(r1["pos"] - r2["pos"]) > 1_000_000:
                continue
            d1 = dosages[v1]
            d2 = dosages[v2]
            valid = ~np.isnan(d1) & ~np.isnan(d2)
            if np.sum(valid) < 20:
                continue
            if np.std(d1[valid]) == 0 or np.std(d2[valid]) == 0:
                continue
            r2_val = np.corrcoef(d1[valid], d2[valid])[0, 1] ** 2
            if r2_val >= r2_threshold:
                G.add_edge(v1, v2, r2=r2_val)

    if verbose:
        print(f"=== Graph-Pruned PRS ===")
        print(f"  LD graph: {G.number_of_nodes()} variants, "
              f"{G.number_of_edges()} LD edges (r² ≥ {r2_threshold})")

    # Select representative variants by graph method
    if pruning_method == "centrality":
        # For each connected component, select the variant with highest
        # betweenness centrality (best tag). If isolated, keep it.
        if G.number_of_edges() > 0:
            centrality = nx.betweenness_centrality(G)
        else:
            centrality = {v: 0 for v in G.nodes()}

        selected = set()
        for component in nx.connected_components(G):
            # Pick the variant with highest centrality, breaking ties by p-value
            best = max(component,
                       key=lambda v: (centrality[v],
                                      -G.nodes[v].get("p_value", 1)))
            selected.add(best)

    elif pruning_method == "pagerank":
        pr = nx.pagerank(G, weight="r2") if G.number_of_edges() > 0 else {v: 1.0 for v in G.nodes()}
        selected = set()
        for component in nx.connected_components(G):
            best = max(component, key=lambda v: pr[v])
            selected.add(best)

    elif pruning_method == "independent_set":
        # Maximal independent set — no two selected variants in LD
        selected = set(nx.maximal_independent_set(G, seed=42))

    else:
        raise ValueError(f"Unknown pruning method: {pruning_method}")

    # Compute PRS with selected variants only
    n_samples = len(all_idx)
    prs = np.zeros(n_samples, dtype=np.float64)
    selected_betas = []
    selected_vids = []

    for vid in selected:
        rec = vid_to_rec[vid]
        dosage = dosages[vid]
        dosage = np.nan_to_num(dosage, nan=0.0)
        prs += rec["beta"] * dosage
        selected_vids.append(vid)
        selected_betas.append(rec["beta"])

    if verbose:
        print(f"  Selected: {len(selected)} / {G.number_of_nodes()} variants "
              f"(pruning: {pruning_method})")
        print(f"  PRS range: [{prs.min():.4f}, {prs.max():.4f}]")

    return {
        "method": "graph_pruned",
        "pruning_method": pruning_method,
        "prs_scores": prs.tolist(),
        "n_variants_input": G.number_of_nodes(),
        "n_variants_selected": len(selected),
        "n_ld_edges": G.number_of_edges(),
        "variant_ids": selected_vids,
        "betas": selected_betas,
        "r2_threshold": r2_threshold,
        "n_cases": len(case_idx),
        "n_controls": len(ctrl_idx),
        "sample_indices": all_idx.tolist(),
    }


# ===================================================================
# P3. Pathway-Partitioned PRS
# ===================================================================

def pathway_partitioned_prs(conn: GraphGWASConnection,
                            run_id: str,
                            p_threshold: float = 1e-5,
                            verbose: bool = True) -> dict:
    """PRS decomposed by biological pathway.

    Instead of a single PRS number, produces a PRS vector where each
    component = the genetic risk attributable to one pathway.  This tells
    you not just "high risk" but "high risk through lipid metabolism and
    immune signaling."

    Also computes pathway-specific predictive accuracy to identify which
    biological systems drive the polygenic signal.

    Args:
        conn: database connection.
        run_id: GWAS run ID.
        p_threshold: p-value threshold (relaxed for pathway aggregation).
        verbose: print progress.

    Returns:
        dict with total_prs, pathway_prs (per pathway), pathway_contributions.
    """
    # Get variants with pathway annotations
    result = conn.execute_read(
        """
        MATCH (ar:AssociationResult)-[:FOR_VARIANT]->(v:Variant)
              -[:HAS_CONSEQUENCE]->(g:Gene)-[:IN_PATHWAY]->(p:Pathway)
        WHERE ar.run_id = $run_id AND ar.p_value < $pval
              AND ar.beta IS NOT NULL
        RETURN v.variantId AS vid, ar.beta AS beta, ar.p_value AS p,
               v.gt_packed AS gtp, p.name AS pathway, g.symbol AS gene
        ORDER BY p.name, ar.p_value
        """,
        {"run_id": run_id, "pval": p_threshold},
    )

    records = [dict(r) for r in result]
    if not records:
        if verbose:
            print("No pathway-annotated significant variants found.")
        return {"method": "pathway_partitioned", "pathway_prs": {}, "n_pathways": 0}

    case_idx, ctrl_idx = get_phenotype_indices(conn)
    all_idx = np.concatenate([case_idx, ctrl_idx])
    n_samples = len(all_idx)
    n_case = len(case_idx)

    # Group variants by pathway
    pathway_variants: dict[str, list[dict]] = {}
    for rec in records:
        pw = rec["pathway"]
        if pw not in pathway_variants:
            pathway_variants[pw] = []
        # Avoid duplicate variants within a pathway
        existing_vids = {v["vid"] for v in pathway_variants[pw]}
        if rec["vid"] not in existing_vids:
            pathway_variants[pw].append(rec)

    # Compute PRS per pathway
    total_prs = np.zeros(n_samples, dtype=np.float64)
    pathway_results = []

    for pw, variants in sorted(pathway_variants.items()):
        pw_prs = np.zeros(n_samples, dtype=np.float64)
        pw_vids = []

        for rec in variants:
            gtp = rec["gtp"]
            if gtp is None:
                continue
            dosage = build_dosage(gtp, all_idx, _cfg.N_SAMPLES)
            dosage = np.nan_to_num(dosage, nan=0.0)
            pw_prs += rec["beta"] * dosage
            pw_vids.append(rec["vid"])

        if not pw_vids:
            continue

        total_prs += pw_prs

        # Pathway-specific discrimination
        case_mean = float(np.mean(pw_prs[:n_case]))
        ctrl_mean = float(np.mean(pw_prs[n_case:]))
        if np.std(pw_prs) > 0:
            # Cohen's d as effect size
            pooled_std = np.sqrt((np.var(pw_prs[:n_case]) + np.var(pw_prs[n_case:])) / 2)
            cohens_d = (case_mean - ctrl_mean) / pooled_std if pooled_std > 0 else 0
        else:
            cohens_d = 0.0

        pathway_results.append({
            "pathway": pw,
            "n_variants": len(pw_vids),
            "variant_ids": pw_vids,
            "prs_case_mean": case_mean,
            "prs_ctrl_mean": ctrl_mean,
            "cohens_d": float(cohens_d),
            "prs_scores": pw_prs.tolist(),
        })

    # Sort pathways by absolute Cohen's d
    pathway_results.sort(key=lambda p: abs(p["cohens_d"]), reverse=True)

    # Compute total PRS contribution per pathway
    total_var = float(np.var(total_prs)) if np.var(total_prs) > 0 else 1.0
    for pr in pathway_results:
        pw_var = float(np.var(pr["prs_scores"]))
        pr["variance_contribution"] = pw_var / total_var
        # Don't store full prs_scores in the summary (too large)
        del pr["prs_scores"]

    if verbose:
        print(f"=== Pathway-Partitioned PRS ===")
        print(f"  {len(pathway_results)} pathways, "
              f"{sum(p['n_variants'] for p in pathway_results)} variant-pathway links")
        print(f"\n  Top pathways by effect size:")
        for pr in pathway_results[:10]:
            print(f"    {pr['pathway']}: d={pr['cohens_d']:.3f}, "
                  f"{pr['n_variants']} variants, "
                  f"var={pr['variance_contribution']:.3f}")

    return {
        "method": "pathway_partitioned",
        "total_prs": total_prs.tolist(),
        "pathway_prs": pathway_results,
        "n_pathways": len(pathway_results),
        "n_cases": n_case,
        "n_controls": len(ctrl_idx),
        "sample_indices": all_idx.tolist(),
    }


# ===================================================================
# P4. Message-Passing PRS (GNN-based)
# ===================================================================

def message_passing_prs(conn: GraphGWASConnection,
                        chr: str = "chr19",
                        start: int | None = None,
                        end: int | None = None,
                        max_variants: int = 5000,
                        epochs: int = 100,
                        device: str = "auto",
                        verbose: bool = True) -> dict:
    """GNN-learned PRS from graph structure.

    Instead of a linear weighted sum, uses the trained GNN's phenotype
    prediction score as a graph-structural PRS.  This captures non-linear
    interactions, gene-level aggregation, and pathway-mediated effects
    that linear PRS misses.

    The GNN PRS = P(case | graph neighborhood), which subsumes the
    classical PRS as a special case (single-layer linear GNN = weighted sum).

    Args:
        conn: database connection.
        chr: chromosome for graph export.
        start, end: position range.
        max_variants: max variants to include.
        epochs: GNN training epochs.
        device: compute device.
        verbose: print progress.

    Returns:
        dict with gnn_prs_scores, auroc, comparison to linear PRS.
    """
    try:
        from .gnn import export_to_pyg, train_gnn
        import torch
    except ImportError:
        return {"method": "message_passing", "error": "PyTorch required"}

    if verbose:
        print(f"=== Message-Passing PRS ===")

    data = export_to_pyg(conn, chr, start, end,
                         max_variants=max_variants, verbose=verbose)
    if data['variant'].num_nodes == 0:
        return {"method": "message_passing", "error": "No variants in region"}

    result = train_gnn(data, epochs=epochs, device=device, verbose=verbose)
    model = result["model"]
    auroc = result["best_auroc"]

    # Extract GNN phenotype scores as PRS
    model.eval()
    with torch.no_grad():
        logits, embeddings = model(data.x_dict, data.edge_index_dict)

    if logits is not None:
        gnn_prs = torch.sigmoid(logits).cpu().numpy()
    else:
        gnn_prs = np.zeros(data['sample'].num_nodes)

    labels = data['sample'].y.cpu().numpy()
    n_case = int(labels.sum())
    n_ctrl = len(labels) - n_case

    if verbose:
        print(f"\n  GNN PRS AUROC: {auroc:.4f}")
        print(f"  PRS mean(cases): {gnn_prs[labels == 1].mean():.4f}, "
              f"mean(controls): {gnn_prs[labels == 0].mean():.4f}")

    return {
        "method": "message_passing",
        "prs_scores": gnn_prs.tolist(),
        "auroc": auroc,
        "n_samples": len(gnn_prs),
        "n_cases": n_case,
        "n_controls": n_ctrl,
        "epochs": epochs,
    }


# ===================================================================
# P5. PRS Evaluation
# ===================================================================

def prs_evaluation(prs_scores: list[float] | np.ndarray,
                   n_cases: int,
                   n_controls: int,
                   verbose: bool = True) -> dict:
    """Evaluate PRS predictive performance.

    Computes AUROC, Nagelkerke R², and quantile-based odds ratios
    (top decile vs bottom decile).

    Args:
        prs_scores: PRS values (cases first, then controls).
        n_cases: number of case samples.
        n_controls: number of control samples.
        verbose: print progress.

    Returns:
        dict with auroc, nagelkerke_r2, or_top_decile, calibration stats.
    """
    prs = np.asarray(prs_scores, dtype=np.float64)
    n = n_cases + n_controls

    if len(prs) != n:
        return {"error": f"PRS length ({len(prs)}) != n_cases+n_controls ({n})"}

    labels = np.concatenate([np.ones(n_cases), np.zeros(n_controls)])

    # AUROC
    from sklearn.metrics import roc_auc_score
    try:
        auroc = float(roc_auc_score(labels, prs))
    except ValueError:
        auroc = 0.5

    # Nagelkerke R² (pseudo R² from logistic regression)
    from scipy.special import expit
    # Fit intercept-only model
    p0 = n_cases / n
    ll_null = n_cases * np.log(p0) + n_controls * np.log(1 - p0)

    # Fit PRS model (simple logistic: logit(p) = a + b*PRS)
    if np.std(prs) > 0:
        # Standardize PRS
        prs_std = (prs - np.mean(prs)) / np.std(prs)
        try:
            import statsmodels.api as sm
            X = sm.add_constant(prs_std)
            model = sm.Logit(labels, X).fit(disp=0, maxiter=25)
            ll_full = model.llf
        except Exception:
            ll_full = ll_null
    else:
        ll_full = ll_null

    # Cox-Snell R²
    r2_cs = 1 - np.exp(2 * (ll_null - ll_full) / n)
    # Nagelkerke correction
    r2_max = 1 - np.exp(2 * ll_null / n)
    nagelkerke_r2 = float(r2_cs / r2_max) if r2_max > 0 else 0.0

    # Odds ratio: top decile vs bottom decile
    decile_thresholds = np.percentile(prs, [10, 90])
    bottom_decile = prs <= decile_thresholds[0]
    top_decile = prs >= decile_thresholds[1]

    n_case_top = np.sum(labels[top_decile])
    n_ctrl_top = np.sum(1 - labels[top_decile])
    n_case_bot = np.sum(labels[bottom_decile])
    n_ctrl_bot = np.sum(1 - labels[bottom_decile])

    if n_ctrl_top > 0 and n_case_bot > 0 and n_ctrl_bot > 0:
        or_top = (n_case_top / n_ctrl_top) / (n_case_bot / n_ctrl_bot)
    else:
        or_top = float("inf")

    if verbose:
        print(f"\n  PRS Evaluation:")
        print(f"    AUROC:          {auroc:.4f}")
        print(f"    Nagelkerke R²:  {nagelkerke_r2:.4f}")
        print(f"    OR (top/bottom decile): {or_top:.2f}")

    return {
        "auroc": auroc,
        "nagelkerke_r2": nagelkerke_r2,
        "or_top_vs_bottom_decile": float(or_top) if not np.isinf(or_top) else None,
        "n_cases": n_cases,
        "n_controls": n_controls,
    }


# ===================================================================
# P6. Integrated PRS Report
# ===================================================================

def prs_report(conn: GraphGWASConnection,
               run_id: str,
               p_threshold: float = 5e-8,
               include_gnn: bool = False,
               gnn_chr: str = "chr19",
               verbose: bool = True) -> dict:
    """Compare classical, graph-pruned, and pathway-partitioned PRS.

    Args:
        conn: database connection.
        run_id: GWAS run ID.
        p_threshold: p-value threshold.
        include_gnn: also compute GNN-based PRS.
        gnn_chr: chromosome for GNN PRS.
        verbose: print progress.

    Returns:
        Comparison of PRS methods with evaluation metrics.
    """
    if verbose:
        print("=" * 60)
        print(f"  PRS Report (run: {run_id})")
        print("=" * 60)

    results = {}

    # P1: Classical
    if verbose:
        print("\n" + "-" * 60)
    classical = classical_prs(conn, run_id, p_threshold, verbose=verbose)
    if classical.get("prs_scores"):
        eval_c = prs_evaluation(classical["prs_scores"],
                                classical["n_cases"], classical["n_controls"],
                                verbose=verbose)
        classical["evaluation"] = eval_c
    results["classical"] = classical

    # P2: Graph-pruned
    if verbose:
        print("\n" + "-" * 60)
    pruned = graph_pruned_prs(conn, run_id, p_threshold, verbose=verbose)
    if pruned.get("prs_scores"):
        eval_p = prs_evaluation(pruned["prs_scores"],
                                pruned["n_cases"], pruned["n_controls"],
                                verbose=verbose)
        pruned["evaluation"] = eval_p
    results["graph_pruned"] = pruned

    # P3: Pathway-partitioned (use relaxed threshold)
    if verbose:
        print("\n" + "-" * 60)
    pathway = pathway_partitioned_prs(conn, run_id, p_threshold=1e-5,
                                       verbose=verbose)
    if pathway.get("total_prs"):
        eval_pw = prs_evaluation(pathway["total_prs"],
                                 pathway["n_cases"], pathway["n_controls"],
                                 verbose=verbose)
        pathway["evaluation"] = eval_pw
    results["pathway_partitioned"] = pathway

    # P4: GNN PRS (optional)
    if include_gnn:
        if verbose:
            print("\n" + "-" * 60)
        gnn = message_passing_prs(conn, chr=gnn_chr, verbose=verbose)
        results["message_passing"] = gnn

    # Summary comparison
    if verbose:
        print("\n" + "=" * 60)
        print("  Comparison")
        print("=" * 60)
        for name, r in results.items():
            ev = r.get("evaluation", {})
            auroc = ev.get("auroc", r.get("auroc", "N/A"))
            r2 = ev.get("nagelkerke_r2", "N/A")
            n_var = r.get("n_variants", r.get("n_variants_selected", "N/A"))
            if isinstance(auroc, float):
                print(f"  {name:<25s} AUROC={auroc:.4f}  R²={r2:.4f}  "
                      f"variants={n_var}")
            else:
                print(f"  {name:<25s} {auroc}")

    return {
        "run_id": run_id,
        "methods": results,
        "summary": {
            name: {
                "auroc": r.get("evaluation", {}).get("auroc", r.get("auroc")),
                "nagelkerke_r2": r.get("evaluation", {}).get("nagelkerke_r2"),
                "n_variants": r.get("n_variants", r.get("n_variants_selected")),
            }
            for name, r in results.items()
        },
    }
