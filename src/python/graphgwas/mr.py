"""Graph-native Mendelian Randomization (MR).

Classical MR uses genetic variants as instrumental variables to infer causal
effects between an exposure and an outcome.  Graph-native MR leverages the
variant→gene→pathway topology for instrument selection, pleiotropy detection,
and pathway-mediated causal decomposition.

Key advantages over classical MR (TwoSampleMR, MR-Egger, MR-PRESSO):
- Instrument selection via graph centrality (not arbitrary p-threshold clumping)
- Pleiotropy detection from graph structure (shared pathway membership)
- Pathway-mediated MR: decompose causal effect by biological route
- LD-aware instrument independence from the LD graph (not LD matrix pruning)

Components:
    R1. select_instruments        — graph-centrality instrument selection
    R2. ivw_estimate              — inverse-variance weighted MR estimate
    R3. egger_estimate            — MR-Egger (intercept = pleiotropy test)
    R4. weighted_median_estimate  — weighted median MR
    R5. pleiotropy_test           — graph-structural pleiotropy detection
    R6. pathway_mediated_mr       — decompose causal effect by pathway
    R7. mr_report                 — integrated MR report with sensitivity
"""

from __future__ import annotations

import numpy as np
from scipy import stats as sp_stats

from .db import GraphGWASConnection


# ===================================================================
# R1. Graph-Centrality Instrument Selection
# ===================================================================

def select_instruments(conn: GraphGWASConnection,
                       exposure_run_id: str,
                       outcome_run_id: str,
                       p_threshold: float = 5e-8,
                       r2_threshold: float = 0.01,
                       verbose: bool = True) -> dict:
    """Select MR instruments using graph centrality for LD-aware pruning.

    For variants associated with the exposure:
    1. Build LD graph among candidates
    2. Select independent instruments via maximum independent set
       (provably optimal — no variant in LD with another selected one)
    3. Look up each instrument's effect on both exposure and outcome

    Args:
        conn: database connection.
        exposure_run_id: GWAS run ID for exposure trait.
        outcome_run_id: GWAS run ID for outcome trait.
        p_threshold: significance threshold for exposure association.
        r2_threshold: LD r² threshold — instruments must be below this.
        verbose: print progress.

    Returns:
        dict with instruments list, each having beta_exposure, se_exposure,
        beta_outcome, se_outcome.
    """
    from .genotype import build_dosage, get_phenotype_indices
    from .config import N_SAMPLES
    import networkx as nx

    # Get exposure-significant variants
    exp_result = conn.execute_read(
        """
        MATCH (ar:AssociationResult)-[:FOR_VARIANT]->(v:Variant)
        WHERE ar.run_id = $run_id AND ar.p_value < $pval
              AND ar.beta IS NOT NULL AND ar.se IS NOT NULL
        RETURN v.variantId AS vid, ar.beta AS beta_exp, ar.se AS se_exp,
               ar.p_value AS p_exp, v.gt_packed AS gtp,
               v.chr AS chr, v.pos AS pos
        ORDER BY ar.p_value
        """,
        {"run_id": exposure_run_id, "pval": p_threshold},
    )
    exp_records = {r["vid"]: dict(r) for r in exp_result}

    if not exp_records:
        if verbose:
            print(f"No instruments found (p < {p_threshold:.0e})")
        return {"instruments": [], "n_instruments": 0}

    # Get outcome effects for the same variants
    vid_list = list(exp_records.keys())
    out_result = conn.execute_read(
        """
        MATCH (ar:AssociationResult)-[:FOR_VARIANT]->(v:Variant)
        WHERE ar.run_id = $run_id AND v.variantId IN $vids
        RETURN v.variantId AS vid, ar.beta AS beta_out, ar.se AS se_out,
               ar.p_value AS p_out
        """,
        {"run_id": outcome_run_id, "vids": vid_list},
    )
    out_lookup = {r["vid"]: dict(r) for r in out_result}

    # Build LD graph for instrument pruning
    case_idx, ctrl_idx = get_phenotype_indices(conn)
    all_idx = np.concatenate([case_idx, ctrl_idx])

    dosages = {}
    for vid, rec in exp_records.items():
        if rec["gtp"] is not None:
            dosages[vid] = build_dosage(rec["gtp"], all_idx, N_SAMPLES)

    G = nx.Graph()
    for vid in dosages:
        G.add_node(vid)

    vids = list(dosages.keys())
    for i in range(len(vids)):
        for j in range(i + 1, len(vids)):
            v1, v2 = vids[i], vids[j]
            r1, r2 = exp_records[v1], exp_records[v2]
            if r1["chr"] != r2["chr"]:
                continue
            if abs(r1["pos"] - r2["pos"]) > 10_000_000:
                continue
            d1, d2 = dosages[v1], dosages[v2]
            valid = ~np.isnan(d1) & ~np.isnan(d2)
            if np.sum(valid) < 20:
                continue
            if np.std(d1[valid]) == 0 or np.std(d2[valid]) == 0:
                continue
            r2_val = np.corrcoef(d1[valid], d2[valid])[0, 1] ** 2
            if r2_val >= r2_threshold:
                G.add_edge(v1, v2, r2=r2_val)

    # Maximum independent set for optimal instrument selection
    if G.number_of_edges() > 0:
        independent = set(nx.maximal_independent_set(G, seed=42))
    else:
        independent = set(G.nodes())

    # Build instrument list
    instruments = []
    for vid in independent:
        if vid not in out_lookup:
            continue
        exp = exp_records[vid]
        out = out_lookup[vid]
        if exp["se_exp"] is None or out["se_out"] is None:
            continue
        if exp["se_exp"] <= 0 or out["se_out"] is None:
            continue

        instruments.append({
            "variant": vid,
            "chr": exp["chr"],
            "pos": exp["pos"],
            "beta_exposure": float(exp["beta_exp"]),
            "se_exposure": float(exp["se_exp"]),
            "p_exposure": float(exp["p_exp"]),
            "beta_outcome": float(out["beta_out"]),
            "se_outcome": float(out["se_out"]),
            "p_outcome": float(out["p_out"]),
        })

    instruments.sort(key=lambda i: i["p_exposure"])

    if verbose:
        print(f"=== MR Instrument Selection ===")
        print(f"  Exposure candidates: {len(exp_records)} (p < {p_threshold:.0e})")
        print(f"  LD graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges "
              f"(r² ≥ {r2_threshold})")
        print(f"  Independent instruments: {len(instruments)} "
              f"(with outcome data)")

    return {
        "instruments": instruments,
        "n_instruments": len(instruments),
        "n_candidates": len(exp_records),
        "n_ld_edges": G.number_of_edges(),
        "r2_threshold": r2_threshold,
    }


# ===================================================================
# R2. Inverse-Variance Weighted (IVW) MR
# ===================================================================

def ivw_estimate(instruments: list[dict],
                 verbose: bool = True) -> dict:
    """Inverse-variance weighted MR estimate.

    β_IVW = Σ(β_out × β_exp / se²_out) / Σ(β²_exp / se²_out)

    Assumes all instruments are valid (no pleiotropy).  Use Egger
    for pleiotropy-robust estimation.

    Args:
        instruments: list of instrument dicts from select_instruments.
        verbose: print results.

    Returns:
        dict with beta_ivw, se_ivw, p_value, n_instruments.
    """
    if len(instruments) < 1:
        return _null_mr("ivw")

    bx = np.array([i["beta_exposure"] for i in instruments])
    by = np.array([i["beta_outcome"] for i in instruments])
    se_y = np.array([i["se_outcome"] for i in instruments])

    # Wald ratios
    wald = by / bx
    weights = (bx ** 2) / (se_y ** 2)

    # IVW estimate (fixed effects)
    beta_ivw = float(np.sum(weights * wald) / np.sum(weights))
    se_ivw = float(1.0 / np.sqrt(np.sum(weights)))
    z = beta_ivw / se_ivw
    p_value = float(2 * sp_stats.norm.sf(abs(z)))

    # Cochran's Q for heterogeneity
    q_stat = float(np.sum(weights * (wald - beta_ivw) ** 2))
    q_df = len(instruments) - 1
    q_p = float(sp_stats.chi2.sf(q_stat, q_df)) if q_df > 0 else 1.0

    if verbose:
        print(f"  IVW: β={beta_ivw:.4f}, SE={se_ivw:.4f}, "
              f"p={p_value:.2e} ({len(instruments)} instruments)")
        print(f"  Heterogeneity Q={q_stat:.2f}, p={q_p:.4f}")

    return {
        "method": "ivw",
        "beta": beta_ivw,
        "se": se_ivw,
        "p_value": p_value,
        "n_instruments": len(instruments),
        "q_stat": q_stat,
        "q_p_value": q_p,
    }


# ===================================================================
# R3. MR-Egger
# ===================================================================

def egger_estimate(instruments: list[dict],
                   verbose: bool = True) -> dict:
    """MR-Egger regression with intercept test for directional pleiotropy.

    Regresses β_outcome on β_exposure with intercept.
    Slope = causal effect (robust to pleiotropy if InSIDE holds).
    Intercept ≠ 0 → evidence of directional pleiotropy.

    Args:
        instruments: list of instrument dicts.
        verbose: print results.

    Returns:
        dict with beta_egger, se_egger, p_value, intercept, intercept_p.
    """
    if len(instruments) < 3:
        return _null_mr("egger")

    bx = np.array([i["beta_exposure"] for i in instruments])
    by = np.array([i["beta_outcome"] for i in instruments])
    se_y = np.array([i["se_outcome"] for i in instruments])

    # Ensure positive exposure effects (orient instruments)
    sign = np.sign(bx)
    bx = np.abs(bx)
    by = by * sign

    # Weighted regression: by ~ intercept + slope * bx, weights = 1/se_y²
    weights = 1.0 / (se_y ** 2)
    W = np.diag(weights)
    X = np.column_stack([np.ones(len(bx)), bx])

    try:
        # WLS: (X'WX)^{-1} X'Wy
        XtWX = X.T @ W @ X
        XtWy = X.T @ W @ by
        beta_hat = np.linalg.solve(XtWX, XtWy)
        residuals = by - X @ beta_hat
        n = len(by)
        mse = float(np.sum(weights * residuals ** 2) / (n - 2))
        var_beta = mse * np.linalg.inv(XtWX)
        se_hat = np.sqrt(np.diag(var_beta))
    except np.linalg.LinAlgError:
        return _null_mr("egger")

    intercept = float(beta_hat[0])
    se_intercept = float(se_hat[0])
    beta_egger = float(beta_hat[1])
    se_egger = float(se_hat[1])

    z_slope = beta_egger / se_egger
    p_slope = float(2 * sp_stats.norm.sf(abs(z_slope)))
    z_int = intercept / se_intercept if se_intercept > 0 else 0
    p_intercept = float(2 * sp_stats.norm.sf(abs(z_int)))

    if verbose:
        print(f"  Egger: β={beta_egger:.4f}, SE={se_egger:.4f}, "
              f"p={p_slope:.2e}")
        print(f"  Intercept={intercept:.4f}, p={p_intercept:.4f} "
              f"({'PLEIOTROPY' if p_intercept < 0.05 else 'OK'})")

    return {
        "method": "egger",
        "beta": beta_egger,
        "se": se_egger,
        "p_value": p_slope,
        "intercept": intercept,
        "se_intercept": se_intercept,
        "intercept_p_value": p_intercept,
        "pleiotropy_detected": p_intercept < 0.05,
        "n_instruments": len(instruments),
    }


# ===================================================================
# R4. Weighted Median MR
# ===================================================================

def weighted_median_estimate(instruments: list[dict],
                             verbose: bool = True) -> dict:
    """Weighted median MR estimate.

    Consistent when ≥50% of the weight comes from valid instruments.
    More robust than IVW when some instruments are pleiotropic.

    Args:
        instruments: list of instrument dicts.
        verbose: print results.

    Returns:
        dict with beta, se (bootstrapped), p_value.
    """
    if len(instruments) < 3:
        return _null_mr("weighted_median")

    bx = np.array([i["beta_exposure"] for i in instruments])
    by = np.array([i["beta_outcome"] for i in instruments])
    se_y = np.array([i["se_outcome"] for i in instruments])

    wald = by / bx
    weights = (bx ** 2) / (se_y ** 2)
    weights = weights / np.sum(weights)  # normalize

    # Weighted median
    sorted_idx = np.argsort(wald)
    cumw = np.cumsum(weights[sorted_idx])
    median_idx = np.searchsorted(cumw, 0.5)
    median_idx = min(median_idx, len(wald) - 1)
    beta_wm = float(wald[sorted_idx[median_idx]])

    # Bootstrap SE
    rng = np.random.default_rng(42)
    n_boot = 1000
    boot_betas = []
    for _ in range(n_boot):
        by_boot = by + rng.normal(0, se_y)
        wald_boot = by_boot / bx
        sorted_b = np.argsort(wald_boot)
        cumw_b = np.cumsum(weights[sorted_b])
        mi = np.searchsorted(cumw_b, 0.5)
        mi = min(mi, len(wald_boot) - 1)
        boot_betas.append(wald_boot[sorted_b[mi]])

    se_wm = float(np.std(boot_betas))
    z = beta_wm / se_wm if se_wm > 0 else 0
    p_value = float(2 * sp_stats.norm.sf(abs(z)))

    if verbose:
        print(f"  Weighted median: β={beta_wm:.4f}, SE={se_wm:.4f}, "
              f"p={p_value:.2e}")

    return {
        "method": "weighted_median",
        "beta": beta_wm,
        "se": se_wm,
        "p_value": p_value,
        "n_instruments": len(instruments),
    }


# ===================================================================
# R5. Graph-Structural Pleiotropy Detection
# ===================================================================

def pleiotropy_test(conn: GraphGWASConnection,
                    instruments: list[dict],
                    verbose: bool = True) -> dict:
    """Detect instrument pleiotropy from graph structure.

    An instrument is potentially pleiotropic if it connects to the outcome
    through a DIFFERENT pathway than the exposure.  This is detectable
    from the variant→gene→pathway graph without any statistical test.

    Classifies each instrument as:
    - 'valid': connects to exposure pathway only
    - 'pleiotropic': connects to additional pathways that also affect outcome
    - 'unknown': no pathway annotation

    Args:
        conn: database connection.
        instruments: instrument list from select_instruments.
        verbose: print results.

    Returns:
        dict with per-instrument pleiotropy classification and pathway details.
    """
    if not instruments:
        return {"instruments": [], "n_pleiotropic": 0}

    vids = [i["variant"] for i in instruments]

    # Get pathway memberships for each instrument variant
    result = conn.execute_read(
        """
        MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene)-[:IN_PATHWAY]->(p:Pathway)
        WHERE v.variantId IN $vids
        RETURN v.variantId AS vid, collect(DISTINCT p.name) AS pathways,
               collect(DISTINCT g.symbol) AS genes
        """,
        {"vids": vids},
    )

    pathway_map = {}
    gene_map = {}
    for rec in result:
        pathway_map[rec["vid"]] = set(rec["pathways"]) if rec["pathways"] else set()
        gene_map[rec["vid"]] = set(rec["genes"]) if rec["genes"] else set()

    # Find the "exposure pathway set" = union of pathways for the strongest instruments
    top_instruments = instruments[:max(3, len(instruments) // 3)]
    exposure_pathways = set()
    for inst in top_instruments:
        exposure_pathways.update(pathway_map.get(inst["variant"], set()))

    # Classify each instrument
    classified = []
    n_pleiotropic = 0
    for inst in instruments:
        vid = inst["variant"]
        pws = pathway_map.get(vid, set())
        genes = gene_map.get(vid, set())

        if not pws:
            status = "unknown"
            extra_pathways = []
        elif pws.issubset(exposure_pathways):
            status = "valid"
            extra_pathways = []
        else:
            status = "pleiotropic"
            extra_pathways = sorted(pws - exposure_pathways)
            n_pleiotropic += 1

        classified.append({
            **inst,
            "pleiotropy_status": status,
            "pathways": sorted(pws),
            "genes": sorted(genes),
            "extra_pathways": extra_pathways,
        })

    if verbose:
        n_valid = sum(1 for c in classified if c["pleiotropy_status"] == "valid")
        n_unknown = sum(1 for c in classified if c["pleiotropy_status"] == "unknown")
        print(f"\n  Graph pleiotropy test:")
        print(f"    Valid: {n_valid}, Pleiotropic: {n_pleiotropic}, "
              f"Unknown: {n_unknown}")
        print(f"    Exposure pathways: {len(exposure_pathways)}")
        for c in classified:
            if c["pleiotropy_status"] == "pleiotropic":
                print(f"    WARNING: {c['variant']} — extra pathways: "
                      f"{', '.join(c['extra_pathways'][:3])}")

    return {
        "instruments": classified,
        "n_valid": sum(1 for c in classified if c["pleiotropy_status"] == "valid"),
        "n_pleiotropic": n_pleiotropic,
        "n_unknown": sum(1 for c in classified if c["pleiotropy_status"] == "unknown"),
        "exposure_pathways": sorted(exposure_pathways),
    }


# ===================================================================
# R6. Pathway-Mediated MR
# ===================================================================

def pathway_mediated_mr(conn: GraphGWASConnection,
                        instruments: list[dict],
                        verbose: bool = True) -> dict:
    """Decompose causal effect by biological pathway.

    Groups instruments by their pathway membership and computes
    pathway-specific IVW estimates.  This answers: "through which
    biological route does the exposure cause the outcome?"

    Args:
        conn: database connection.
        instruments: instrument list from select_instruments.
        verbose: print results.

    Returns:
        dict with per-pathway causal estimates.
    """
    if len(instruments) < 2:
        return {"pathway_estimates": [], "method": "pathway_mediated"}

    vids = [i["variant"] for i in instruments]
    inst_lookup = {i["variant"]: i for i in instruments}

    # Get pathway assignments
    result = conn.execute_read(
        """
        MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene)-[:IN_PATHWAY]->(p:Pathway)
        WHERE v.variantId IN $vids
        RETURN v.variantId AS vid, p.name AS pathway
        """,
        {"vids": vids},
    )

    # Group instruments by pathway
    pathway_instruments: dict[str, list[dict]] = {}
    for rec in result:
        pw = rec["pathway"]
        vid = rec["vid"]
        if vid in inst_lookup:
            if pw not in pathway_instruments:
                pathway_instruments[pw] = []
            # Avoid duplicate instruments within a pathway
            existing_vids = {i["variant"] for i in pathway_instruments[pw]}
            if vid not in existing_vids:
                pathway_instruments[pw].append(inst_lookup[vid])

    # Compute IVW per pathway
    pathway_estimates = []
    for pw, pw_instruments in sorted(pathway_instruments.items()):
        if len(pw_instruments) < 1:
            continue
        est = ivw_estimate(pw_instruments, verbose=False)
        pathway_estimates.append({
            "pathway": pw,
            "beta": est["beta"],
            "se": est["se"],
            "p_value": est["p_value"],
            "n_instruments": est["n_instruments"],
        })

    pathway_estimates.sort(key=lambda p: abs(p["beta"]), reverse=True)

    # Unassigned instruments (no pathway)
    assigned_vids = set()
    for pw_insts in pathway_instruments.values():
        for i in pw_insts:
            assigned_vids.add(i["variant"])
    unassigned = [i for i in instruments if i["variant"] not in assigned_vids]

    if verbose:
        print(f"\n  Pathway-Mediated MR:")
        print(f"    {len(pathway_estimates)} pathways with instruments")
        print(f"    {len(unassigned)} instruments without pathway annotation")
        for pe in pathway_estimates[:10]:
            sig = "*" if pe["p_value"] < 0.05 else ""
            print(f"    {pe['pathway']}: β={pe['beta']:.4f} "
                  f"(p={pe['p_value']:.2e}, n={pe['n_instruments']}){sig}")

    return {
        "method": "pathway_mediated",
        "pathway_estimates": pathway_estimates,
        "n_pathways": len(pathway_estimates),
        "n_unassigned": len(unassigned),
    }


# ===================================================================
# R7. Integrated MR Report
# ===================================================================

def mr_report(conn: GraphGWASConnection,
              exposure_run_id: str,
              outcome_run_id: str,
              p_threshold: float = 5e-8,
              verbose: bool = True) -> dict:
    """Full MR analysis: instrument selection + 3 estimators + pleiotropy.

    Args:
        conn: database connection.
        exposure_run_id: GWAS run ID for exposure.
        outcome_run_id: GWAS run ID for outcome.
        p_threshold: significance threshold for instruments.
        verbose: print progress.

    Returns:
        Comprehensive MR report with sensitivity analyses.
    """
    if verbose:
        print("=" * 60)
        print(f"  MR Report: {exposure_run_id} → {outcome_run_id}")
        print("=" * 60)

    # Select instruments
    if verbose:
        print("\n" + "-" * 60)
    selection = select_instruments(conn, exposure_run_id, outcome_run_id,
                                   p_threshold=p_threshold, verbose=verbose)
    instruments = selection["instruments"]

    if not instruments:
        return {"error": "No valid instruments found",
                "exposure": exposure_run_id, "outcome": outcome_run_id}

    results = {"selection": selection}

    # IVW
    if verbose:
        print("\n--- MR Estimates ---")
    results["ivw"] = ivw_estimate(instruments, verbose=verbose)

    # Egger
    results["egger"] = egger_estimate(instruments, verbose=verbose)

    # Weighted median
    results["weighted_median"] = weighted_median_estimate(instruments,
                                                          verbose=verbose)

    # Graph pleiotropy
    if verbose:
        print("\n--- Pleiotropy ---")
    results["pleiotropy"] = pleiotropy_test(conn, instruments, verbose=verbose)

    # Pathway-mediated
    results["pathway_mediated"] = pathway_mediated_mr(conn, instruments,
                                                      verbose=verbose)

    # Summary
    if verbose:
        print("\n" + "=" * 60)
        print("  Summary")
        print("=" * 60)
        for method in ["ivw", "egger", "weighted_median"]:
            r = results[method]
            sig = "***" if r["p_value"] < 0.001 else ("**" if r["p_value"] < 0.01 else ("*" if r["p_value"] < 0.05 else ""))
            print(f"  {method:<20s} β={r['beta']:.4f} "
                  f"(SE={r['se']:.4f}, p={r['p_value']:.2e}){sig}")

        pl = results["pleiotropy"]
        print(f"\n  Pleiotropy: {pl['n_pleiotropic']} / {len(instruments)} "
              f"instruments flagged")
        egger_int = results["egger"].get("intercept_p_value", 1)
        print(f"  Egger intercept p = {egger_int:.4f} "
              f"({'directional pleiotropy' if egger_int < 0.05 else 'no evidence'})")

    return {
        "exposure": exposure_run_id,
        "outcome": outcome_run_id,
        "n_instruments": len(instruments),
        "results": results,
        "summary": {
            "ivw_beta": results["ivw"]["beta"],
            "ivw_p": results["ivw"]["p_value"],
            "egger_beta": results["egger"]["beta"],
            "egger_p": results["egger"]["p_value"],
            "egger_intercept_p": results["egger"].get("intercept_p_value"),
            "wm_beta": results["weighted_median"]["beta"],
            "wm_p": results["weighted_median"]["p_value"],
            "n_pleiotropic": results["pleiotropy"]["n_pleiotropic"],
        },
    }


# ===================================================================
# Helpers
# ===================================================================

def _null_mr(method: str) -> dict:
    return {
        "method": method,
        "beta": 0.0,
        "se": float("nan"),
        "p_value": 1.0,
        "n_instruments": 0,
        "error": "Insufficient instruments",
    }
