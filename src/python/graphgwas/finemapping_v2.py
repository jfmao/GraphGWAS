"""Graph-Native LD Fine-Mapping v2.

Four methods for identifying causal variants within LD blocks:
  L1: Dual-Graph Fine-Mapping (LD + functional annotations)
  L2: LD-Regularized Association (spectral decomposition)
  L3: Haplotype Graph Association
  L4: Recombination-Aware Embedding

This module implements L1 and L4 first.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import numpy as np
from scipy import stats as sp_stats
from scipy.linalg import cho_factor, cho_solve

from . import config as _cfg
from .db import GraphGWASConnection
from .genotype import (
    get_all_indices,
    get_phenotype_values,
    build_dosage,
    variant_iterator,
)


# ===================================================================
# Data classes
# ===================================================================

@dataclass
class FinemapCandidate:
    """A variant with fine-mapping scores."""
    variant_id: str
    chr: str
    pos: int
    ref: str
    alt: str
    af: float
    z_stat: float           # statistical signal (-log10 p or z-score)
    z_functional: float     # functional annotation score
    unique_stat: float      # signal after LD deconvolution
    combined_score: float   # final combined score
    pip: float              # posterior inclusion probability
    in_credible_set: bool   # in 95% credible set
    annotations: list[str]  # gene/pathway annotations that contributed
    n_ld_neighbors: int     # number of LD neighbors (r² > threshold)


# ===================================================================
# Shared utilities
# ===================================================================

def _load_locus_variants(
    conn: GraphGWASConnection,
    chr: str,
    center: int,
    window: int,
    all_idx: np.ndarray,
) -> list[dict]:
    """Load variants in a locus window with dosages and p-values."""
    start = center - window // 2
    end = center + window // 2

    variants = []
    for v in variant_iterator(conn, chr, start, end):
        af = v.get("af_total", 0)
        if af < 0.001 or af > 0.999:
            continue
        gtp = v.get("gt_packed")
        if gtp is None:
            continue
        dosage = build_dosage(gtp, all_idx, _cfg.N_SAMPLES)
        v["dosage"] = dosage
        variants.append(v)

    return variants


def _load_locus_variants_bgen(
    reader,
    chr: str,
    center: int,
    window: int,
    min_af: float = 0.001,
) -> list[dict]:
    """BGEN-backed locus loader. Returns dicts matching the Neo4j loader shape."""
    start = center - window // 2
    end = center + window // 2
    variant_df, dosage_matrix = reader.load_locus(chr, start, end, format="dosage")
    variants: list[dict] = []
    for i, row in variant_df.iterrows():
        d = dosage_matrix[:, i]
        af = float(d.mean() / 2.0)
        if af < min_af or af > (1 - min_af):
            continue
        variants.append(
            {
                "variantId": f"{row['chr']}:{row['pos']}:{row['a1']}:{row['a2']}",
                "chr": row["chr"],
                "pos": int(row["pos"]),
                "ref": row["a1"],
                "alt": row["a2"],
                "af_total": af,
                "dosage": d.astype(np.float64),
            }
        )
    return variants


def load_locus_variants(
    chr: str,
    center: int,
    window: int,
    source,
    all_idx: np.ndarray | None = None,
) -> list[dict]:
    """Unified locus loader dispatching on source type.

    If `source` is a GraphGWASConnection, uses Neo4j-backed _load_locus_variants
    (requires `all_idx`). If `source` is a BgenReader, uses BGEN-backed loader.
    Both return the same dict shape for downstream methods.
    """
    from .bgen_reader import BgenReader  # local import to avoid heavy import cost
    if isinstance(source, BgenReader):
        return _load_locus_variants_bgen(source, chr, center, window)
    if all_idx is None:
        raise ValueError("all_idx required when source is a Neo4j connection")
    return _load_locus_variants(source, chr, center, window, all_idx)


def _compute_ld_matrix(variants: list[dict]) -> np.ndarray:
    """Compute pairwise r² matrix for locus variants.

    Vectorized: single matrix multiplication via BLAS instead of O(n²) Python loop.
    For 700 variants × 3202 samples: ~0.5s (vs ~160s with loop).
    """
    n = len(variants)
    # Stack dosages into matrix (n_variants × n_samples)
    D = np.array([v["dosage"] for v in variants], dtype=np.float64)
    # Mean-impute NaN
    row_means = np.nanmean(D, axis=1, keepdims=True)
    nan_mask = np.isnan(D)
    D = np.where(nan_mask, row_means, D)
    # Standardize each variant (mean=0, std=1)
    stds = D.std(axis=1, keepdims=True)
    stds = np.where(stds < 1e-10, 1.0, stds)  # avoid division by zero
    D = (D - D.mean(axis=1, keepdims=True)) / stds
    # Correlation matrix = (D @ D.T) / n_samples
    R = (D @ D.T) / D.shape[1]
    # r² = correlation²
    R = R ** 2
    np.fill_diagonal(R, 1.0)
    # Clamp to [0, 1]
    R = np.clip(R, 0, 1)
    return R


def _compute_association_stats(
    variants: list[dict],
    phenotype: np.ndarray,
) -> np.ndarray:
    """Compute -log10(p) for each variant via vectorized regression.

    Vectorized: computes all variant associations in one matrix operation
    instead of per-variant loop. For 700 variants: ~0.01s vs ~1s.
    """
    n_var = len(variants)
    # Stack dosages (n_variants × n_samples)
    D = np.array([v["dosage"] for v in variants], dtype=np.float64)
    y = phenotype.copy()

    # Handle NaN: mean-impute dosages, use valid phenotype samples
    valid_pheno = ~np.isnan(y)
    y_clean = np.where(np.isnan(y), 0, y)
    D_clean = np.where(np.isnan(D), 0, D)

    # Use only valid phenotype samples
    if valid_pheno.sum() < 20:
        return np.zeros(n_var)

    D_v = D_clean[:, valid_pheno]
    y_v = y_clean[valid_pheno]
    n = int(valid_pheno.sum())

    # Vectorized simple linear regression: beta = cov(G, Y) / var(G)
    y_centered = y_v - y_v.mean()
    D_centered = D_v - D_v.mean(axis=1, keepdims=True)

    var_g = np.var(D_v, axis=1)  # (n_var,)
    cov_gy = (D_centered @ y_centered) / n  # (n_var,)

    # Avoid division by zero
    safe_var = np.where(var_g > 1e-10, var_g, 1.0)
    beta = cov_gy / safe_var  # (n_var,)

    # Residuals and standard error (vectorized)
    predicted = D_centered * beta[:, np.newaxis]  # (n_var, n)
    residuals = y_centered[np.newaxis, :] - predicted  # (n_var, n)
    mse = np.sum(residuals**2, axis=1) / (n - 2)  # (n_var,)
    se = np.sqrt(mse / (n * safe_var))  # (n_var,)

    # t-statistics and p-values
    safe_se = np.where(se > 1e-10, se, 1.0)
    t_stats = beta / safe_se  # (n_var,)
    p_values = 2 * sp_stats.t.sf(np.abs(t_stats), n - 2)

    # -log10(p), capped at 300
    z_stats = -np.log10(np.maximum(p_values, 1e-300))
    z_stats = np.where(var_g > 1e-10, z_stats, 0)  # zero for monomorphic

    return z_stats


# ===================================================================
# L1: Dual-Graph Fine-Mapping (LD + Functional)
# ===================================================================

def dual_graph_finemap(
    conn: GraphGWASConnection,
    chr: str,
    lead_pos: int,
    window: int = 500_000,
    r2_smooth: float = 0.3,
    alpha: float = 0.5,
    credible_set_coverage: float = 0.95,
    verbose: bool = True,
) -> list[FinemapCandidate]:
    """L1: Dual-Graph Fine-Mapping combining LD + functional annotations.

    For each variant, computes:
    - z_stat: statistical signal from GWAS
    - z_functional: traversal-based functional annotation score
    - unique_stat: signal after LD deconvolution (subtracting LD neighbors)
    - combined: α * unique_stat + (1-α) * z_functional

    Variants with both statistical signal AND functional context rank higher.
    Proxy variants (LD-only signal, no functional annotations) get suppressed.

    Args:
        conn: database connection.
        chr: chromosome of lead variant.
        lead_pos: position of lead GWAS variant.
        window: window size around lead (bp).
        r2_smooth: LD threshold for smoothing/deconvolution.
        alpha: balance statistical vs functional (0=pure functional, 1=pure statistical).
        credible_set_coverage: coverage for credible set (default 95%).
        verbose: print progress.

    Returns:
        list of FinemapCandidate sorted by combined_score.
    """
    all_idx = get_all_indices(conn)
    pheno = get_phenotype_values(conn, all_idx)

    if verbose:
        print(f"L1 Dual-Graph Fine-Mapping: {chr}:{lead_pos} ±{window//1000}kb")

    # Step 1: Load locus variants
    variants = _load_locus_variants(conn, chr, lead_pos, window, all_idx)
    n_var = len(variants)

    if verbose:
        print(f"  {n_var} variants in locus")

    if n_var < 3:
        return []

    # Step 2: Compute association statistics
    z_stats = _compute_association_stats(variants, pheno)

    # Step 3: Compute LD matrix
    R = _compute_ld_matrix(variants)

    # Step 4: LD deconvolution — unique contribution after removing neighbors
    unique_stats = np.zeros(n_var)
    n_neighbors = np.zeros(n_var, dtype=int)
    for i in range(n_var):
        neighbors = np.where(R[i] > r2_smooth)[0]
        neighbors = neighbors[neighbors != i]
        n_neighbors[i] = len(neighbors)
        if len(neighbors) > 0:
            ld_contribution = np.sum(R[i, neighbors] * z_stats[neighbors])
            unique_stats[i] = z_stats[i] - ld_contribution / len(neighbors)
        else:
            unique_stats[i] = z_stats[i]

    # Normalize unique_stats to [0, max(z_stats)]
    if unique_stats.max() > 0:
        unique_stats = unique_stats / unique_stats.max() * z_stats.max()
    unique_stats = np.maximum(unique_stats, 0)

    # Step 5: Multi-omics functional annotation score via graph traversal
    z_func = np.zeros(n_var)

    # Batch query: get all annotation layers for locus variants at once
    vids = [v["variantId"] for v in variants]

    # Layer 1: Gene + Pathway (graph traversal)
    gene_scores = {}
    result = conn.execute_read("""
        UNWIND $vids AS vid
        MATCH (v:Variant {variantId: vid})-[:HAS_CONSEQUENCE]->(g:Gene)
        OPTIONAL MATCH (g)-[:IN_PATHWAY]->(p:Pathway)
        OPTIONAL MATCH (g)-[:INTERACTS_WITH]-(g2:Gene)
        RETURN vid, g.symbol AS gene, collect(DISTINCT p.name) AS pathways,
               count(DISTINCT g2) AS ppi_partners
    """, {"vids": vids})
    for rec in result:
        vid = rec["vid"]
        gene_scores.setdefault(vid, {"genes": [], "pathways": [], "ppi": 0})
        if rec["gene"]:
            gene_scores[vid]["genes"].append(rec["gene"])
            gene_scores[vid]["ppi"] = max(gene_scores[vid]["ppi"], rec["ppi_partners"])
        for pw in rec.get("pathways", []):
            if pw:
                gene_scores[vid]["pathways"].append(pw)

    # Layer 2: eQTL score (from variant properties)
    eqtl_scores = {}
    result = conn.execute_read("""
        UNWIND $vids AS vid
        MATCH (v:Variant {variantId: vid})
        WHERE v.eqtl_score IS NOT NULL
        RETURN vid, v.eqtl_score AS eqtl, v.eqtl_gene AS eqtl_gene
    """, {"vids": vids})
    for rec in result:
        eqtl_scores[rec["vid"]] = {"score": rec["eqtl"], "gene": rec.get("eqtl_gene", "")}

    # Layer 3: Conservation score (from variant properties)
    cons_scores = {}
    result = conn.execute_read("""
        UNWIND $vids AS vid
        MATCH (v:Variant {variantId: vid})
        WHERE v.conservation_score IS NOT NULL
        RETURN vid, v.conservation_score AS cons
    """, {"vids": vids})
    for rec in result:
        cons_scores[rec["vid"]] = rec["cons"]

    # Compute combined functional score per variant
    for i, v in enumerate(variants):
        vid = v["variantId"]
        annotations = []
        score = 0.0

        # Gene/Pathway layer (weight: 1.0 per gene, 0.5 per pathway)
        gs = gene_scores.get(vid, {"genes": [], "pathways": [], "ppi": 0})
        if gs["genes"]:
            score += 1.0
            annotations.extend(f"gene:{g}" for g in gs["genes"])
        if gs["pathways"]:
            score += len(set(gs["pathways"])) * 0.5
            annotations.extend(f"pathway:{p}" for p in set(gs["pathways"]))

        # PPI layer (weight: 0.3 per interaction partner, capped at 3.0)
        if gs["ppi"] > 0:
            score += min(gs["ppi"] * 0.3, 3.0)
            annotations.append(f"ppi_partners:{gs['ppi']}")

        # eQTL layer (weight: 2.0 × eqtl_score — strong signal)
        eq = eqtl_scores.get(vid)
        if eq:
            score += 2.0 * eq["score"]
            annotations.append(f"eqtl:{eq['gene']}({eq['score']:.1f})")

        # Conservation layer (weight: 1.5 × conservation_score)
        cons = cons_scores.get(vid, 0)
        if cons > 0.5:  # only count if moderately conserved
            score += 1.5 * cons
            annotations.append(f"conservation:{cons:.2f}")

        z_func[i] = np.log1p(score)
        v["_annotations"] = annotations

    # Normalize functional scores
    if z_func.max() > 0:
        z_func = z_func / z_func.max() * z_stats.max()

    if verbose:
        n_annotated = np.sum(z_func > 0)
        print(f"  {n_annotated}/{n_var} variants with functional annotations")

    # Step 6: Combined score
    combined = alpha * unique_stats + (1 - alpha) * z_func

    # Step 7: PIP (softmax of combined score)
    exp_scores = np.exp(combined - combined.max())  # numerically stable
    pip = exp_scores / exp_scores.sum()

    # Step 8: Credible set (cumulative PIP to coverage)
    order = np.argsort(-pip)
    cumsum = np.cumsum(pip[order])
    in_cs = np.zeros(n_var, dtype=bool)
    for k, idx in enumerate(order):
        in_cs[idx] = True
        if cumsum[k] >= credible_set_coverage:
            break

    # Build results
    candidates = []
    for i in range(n_var):
        candidates.append(FinemapCandidate(
            variant_id=variants[i]["variantId"],
            chr=chr,
            pos=variants[i]["pos"],
            ref=variants[i].get("ref", ""),
            alt=variants[i].get("alt", ""),
            af=variants[i].get("af_total", 0),
            z_stat=float(z_stats[i]),
            z_functional=float(z_func[i]),
            unique_stat=float(unique_stats[i]),
            combined_score=float(combined[i]),
            pip=float(pip[i]),
            in_credible_set=bool(in_cs[i]),
            annotations=variants[i].get("_annotations", []),
            n_ld_neighbors=int(n_neighbors[i]),
        ))

    candidates.sort(key=lambda c: -c.combined_score)

    if verbose:
        cs_size = sum(1 for c in candidates if c.in_credible_set)
        print(f"  Credible set ({credible_set_coverage*100:.0f}%): {cs_size} variants")
        print("\n  Top 5:")
        print(f"  {'Variant':<50s} {'z_stat':>7s} {'z_func':>7s} {'unique':>7s} {'PIP':>6s} {'CS':>3s}")
        print(f"  {'-'*82}")
        for c in candidates[:5]:
            cs_mark = "***" if c.in_credible_set else ""
            print(f"  {c.variant_id:<50s} {c.z_stat:>7.2f} {c.z_functional:>7.2f} "
                  f"{c.unique_stat:>7.2f} {c.pip:>6.4f} {cs_mark:>3s}")

    return candidates


# ===================================================================
# L4: Recombination-Aware Embedding
# ===================================================================

def recombination_embedding_finemap(
    conn: GraphGWASConnection,
    chr: str,
    lead_pos: int,
    window: int = 500_000,
    r2_threshold: float = 0.1,
    n_components: int = 2,
    cluster_method: str = "dbscan",
    eps_cluster: float = 0.5,
    credible_set_coverage: float = 0.95,
    verbose: bool = True,
) -> list[FinemapCandidate]:
    """L4: Recombination-Aware Embedding for multi-signal detection.

    Embeds variants in a metric space where distance reflects recombination
    probability (not physical distance). Then clusters to identify independent
    signals and selects the top candidate per cluster.

    Key insight: variants in LD cluster together. Recombination hotspots
    create gaps. Number of clusters = number of independent causal signals.

    Args:
        conn: database connection.
        chr: chromosome.
        lead_pos: position of lead variant.
        window: window size (bp).
        r2_threshold: LD threshold for distance computation.
        n_components: embedding dimensions (2 for visualization).
        cluster_method: "dbscan" or "hierarchical".
        eps_cluster: DBSCAN epsilon parameter.
        credible_set_coverage: coverage for credible set.
        verbose: print progress.

    Returns:
        list of FinemapCandidate with cluster assignments.
    """
    from sklearn.manifold import MDS
    from sklearn.cluster import DBSCAN

    all_idx = get_all_indices(conn)
    pheno = get_phenotype_values(conn, all_idx)

    if verbose:
        print(f"L4 Recombination Embedding: {chr}:{lead_pos} ±{window//1000}kb")

    # Load variants
    variants = _load_locus_variants(conn, chr, lead_pos, window, all_idx)
    n_var = len(variants)

    if verbose:
        print(f"  {n_var} variants in locus")

    if n_var < 5:
        return []

    # Compute association stats
    z_stats = _compute_association_stats(variants, pheno)

    # Compute LD matrix
    R = _compute_ld_matrix(variants)

    # Step 1: Distance metric = -log(r² + ε)
    # Perfect LD → distance 0, independent → large distance
    eps = 1e-6
    distance = -np.log(R + eps)
    np.fill_diagonal(distance, 0)
    # Cap extreme distances
    max_dist = np.percentile(distance[distance > 0], 95) if (distance > 0).any() else 10
    distance = np.minimum(distance, max_dist)

    # Step 2: MDS embedding
    if verbose:
        print(f"  Computing MDS embedding ({n_components}D)...")
    mds = MDS(n_components=n_components, dissimilarity="precomputed",
              random_state=42, n_init=2, max_iter=300, normalized_stress="auto")
    coords = mds.fit_transform(distance)

    # Step 3: Clustering
    if verbose:
        print(f"  Clustering ({cluster_method}, eps={eps_cluster})...")

    if cluster_method == "dbscan":
        clustering = DBSCAN(eps=eps_cluster, min_samples=3)
        labels = clustering.fit_predict(coords)
    else:
        from sklearn.cluster import AgglomerativeClustering
        clustering = AgglomerativeClustering(
            n_clusters=None, distance_threshold=eps_cluster,
            metric="euclidean", linkage="average")
        labels = clustering.fit_predict(coords)

    n_clusters = len(set(labels) - {-1})
    noise_count = (labels == -1).sum()

    if verbose:
        print(f"  {n_clusters} clusters, {noise_count} noise variants")

    # Step 4: For each cluster, find best candidate
    # Primary: highest z-score
    # Secondary: most central (closest to cluster centroid)
    cluster_info = {}
    for c in set(labels):
        if c == -1:
            continue
        members = np.where(labels == c)[0]
        centroid = coords[members].mean(axis=0)
        distances_to_center = np.sqrt(np.sum((coords[members] - centroid) ** 2, axis=1))

        best_z_idx = members[np.argmax(z_stats[members])]
        best_center_idx = members[np.argmin(distances_to_center)]

        cluster_info[c] = {
            "members": members,
            "n_members": len(members),
            "best_z": best_z_idx,
            "best_center": best_center_idx,
            "centroid": centroid,
        }

    # Step 5: PIP based on z-scores (per-cluster softmax)
    pip = np.zeros(n_var)
    for c, info in cluster_info.items():
        members = info["members"]
        z_m = z_stats[members]
        if z_m.max() > 0:
            exp_z = np.exp(z_m - z_m.max())
            pip[members] = exp_z / exp_z.sum() / n_clusters  # share PIP across clusters
        else:
            pip[members] = 1.0 / len(members) / n_clusters

    # Noise variants get tiny PIP
    noise = np.where(labels == -1)[0]
    if len(noise) > 0:
        pip[noise] = 1e-6

    # Normalize
    pip = pip / pip.sum()

    # Credible set
    order = np.argsort(-pip)
    cumsum = np.cumsum(pip[order])
    in_cs = np.zeros(n_var, dtype=bool)
    for k, idx in enumerate(order):
        in_cs[idx] = True
        if cumsum[k] >= credible_set_coverage:
            break

    # Build results
    candidates = []
    for i in range(n_var):
        cluster_label = int(labels[i])
        is_primary = any(info["best_z"] == i for info in cluster_info.values())

        candidates.append(FinemapCandidate(
            variant_id=variants[i]["variantId"],
            chr=chr,
            pos=variants[i]["pos"],
            ref=variants[i].get("ref", ""),
            alt=variants[i].get("alt", ""),
            af=variants[i].get("af_total", 0),
            z_stat=float(z_stats[i]),
            z_functional=float(cluster_label),  # repurpose as cluster ID
            unique_stat=float(coords[i, 0]),  # embedding x
            combined_score=float(pip[i]),
            pip=float(pip[i]),
            in_credible_set=bool(in_cs[i]),
            annotations=[f"cluster:{cluster_label}",
                          "primary" if is_primary else "secondary"],
            n_ld_neighbors=int(np.sum(R[i] > r2_threshold)) - 1,
        ))

    candidates.sort(key=lambda c: -c.pip)

    if verbose:
        cs_size = sum(1 for c in candidates if c.in_credible_set)
        print(f"  Credible set ({credible_set_coverage*100:.0f}%): {cs_size} variants")
        print("\n  Clusters:")
        for c, info in sorted(cluster_info.items()):
            best = variants[info["best_z"]]
            print(f"    Cluster {c}: {info['n_members']} variants, "
                  f"top={best['variantId'][:40]} z={z_stats[info['best_z']]:.2f}")
        print("\n  Top 5 candidates:")
        print(f"  {'Variant':<50s} {'z_stat':>7s} {'PIP':>6s} {'Cluster':>7s} {'CS':>3s}")
        print(f"  {'-'*75}")
        for c in candidates[:5]:
            cluster = [a for a in c.annotations if a.startswith("cluster:")][0]
            cs_mark = "***" if c.in_credible_set else ""
            print(f"  {c.variant_id:<50s} {c.z_stat:>7.2f} {c.pip:>6.4f} {cluster:>7s} {cs_mark:>3s}")

    return candidates


# ===================================================================
# Shared helpers (v3 methods)
# ===================================================================

def _compute_ld_correlation(variants: list[dict]) -> np.ndarray:
    """Compute signed LD correlation matrix (Pearson r, not r²)."""
    n = len(variants)
    D = np.array([v["dosage"] for v in variants], dtype=np.float64)
    row_means = np.nanmean(D, axis=1, keepdims=True)
    D = np.where(np.isnan(D), row_means, D)
    stds = D.std(axis=1, keepdims=True)
    stds = np.where(stds < 1e-10, 1.0, stds)
    D = (D - D.mean(axis=1, keepdims=True)) / stds
    R = (D @ D.T) / D.shape[1]
    np.fill_diagonal(R, 1.0)
    return R


def _compute_z_scores(
    variants: list[dict], phenotype: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute z-scores (t-statistics = beta/se) per variant.

    Returns (z_scores, betas, se) arrays of length n_var.
    """
    n_var = len(variants)
    D = np.array([v["dosage"] for v in variants], dtype=np.float64)
    y = phenotype.copy()
    valid = ~np.isnan(y)
    y_clean = np.where(np.isnan(y), 0, y)
    D_clean = np.where(np.isnan(D), 0, D)
    if valid.sum() < 20:
        return np.zeros(n_var), np.zeros(n_var), np.ones(n_var)
    D_v = D_clean[:, valid]
    y_v = y_clean[valid]
    n = int(valid.sum())
    y_c = y_v - y_v.mean()
    D_c = D_v - D_v.mean(axis=1, keepdims=True)
    var_g = np.var(D_v, axis=1)
    safe_var = np.where(var_g > 1e-10, var_g, 1.0)
    beta = (D_c @ y_c) / n / safe_var
    predicted = D_c * beta[:, np.newaxis]
    residuals = y_c[np.newaxis, :] - predicted
    mse = np.sum(residuals**2, axis=1) / (n - 2)
    se = np.sqrt(mse / (n * safe_var))
    safe_se = np.where(se > 1e-10, se, 1.0)
    t_stats = np.where(var_g > 1e-10, beta / safe_se, 0.0)
    return t_stats, beta, se


def _ld_deconvolve(
    z_stats: np.ndarray, R_sq: np.ndarray, r2_threshold: float = 0.3,
) -> tuple[np.ndarray, np.ndarray]:
    """LD deconvolution: unique contribution after removing neighbor signal.

    Returns (unique_stats, n_neighbors).
    """
    n = len(z_stats)
    unique = np.zeros(n)
    n_nb = np.zeros(n, dtype=int)
    for i in range(n):
        nb = np.where(R_sq[i] > r2_threshold)[0]
        nb = nb[nb != i]
        n_nb[i] = len(nb)
        if len(nb) > 0:
            unique[i] = z_stats[i] - np.sum(R_sq[i, nb] * z_stats[nb]) / len(nb)
        else:
            unique[i] = z_stats[i]
    if unique.max() > 0:
        unique = unique / unique.max() * z_stats.max()
    unique = np.maximum(unique, 0)
    return unique, n_nb


def _softmax(scores: np.ndarray) -> np.ndarray:
    """Numerically stable softmax → PIP."""
    s = scores - scores.max()
    e = np.exp(s)
    return e / e.sum()


def _build_credible_set(
    pip: np.ndarray, coverage: float = 0.95,
) -> np.ndarray:
    """Build credible set from PIPs. Returns boolean mask."""
    order = np.argsort(-pip)
    cumsum = np.cumsum(pip[order])
    in_cs = np.zeros(len(pip), dtype=bool)
    for k, idx in enumerate(order):
        in_cs[idx] = True
        if cumsum[k] >= coverage:
            break
    return in_cs


def _load_eqtl_cache(
    path: str = "/mnt/data/GraphGWAS/data/annotations/gtex_chr22_enhanced_cache.json",
) -> dict:
    """Load eQTL annotation cache from file."""
    with open(path) as f:
        return json.load(f)


def _soft_threshold(v: np.ndarray, kappa: float) -> np.ndarray:
    """Soft-thresholding (proximal operator for L1 norm)."""
    return np.sign(v) * np.maximum(np.abs(v) - kappa, 0)


def _build_candidates(
    variants: list[dict],
    pip: np.ndarray,
    z_stats: np.ndarray,
    z_func: np.ndarray,
    unique_stats: np.ndarray,
    combined: np.ndarray,
    n_neighbors: np.ndarray,
    chr_name: str,
    credible_set_coverage: float = 0.95,
    annotations_map: dict | None = None,
) -> list[FinemapCandidate]:
    """Build sorted FinemapCandidate list from arrays."""
    in_cs = _build_credible_set(pip, credible_set_coverage)
    candidates = []
    for i in range(len(variants)):
        candidates.append(FinemapCandidate(
            variant_id=variants[i]["variantId"],
            chr=chr_name,
            pos=variants[i]["pos"],
            ref=variants[i].get("ref", ""),
            alt=variants[i].get("alt", ""),
            af=variants[i].get("af_total", 0),
            z_stat=float(z_stats[i]),
            z_functional=float(z_func[i]),
            unique_stat=float(unique_stats[i]),
            combined_score=float(combined[i]),
            pip=float(pip[i]),
            in_credible_set=bool(in_cs[i]),
            annotations=(annotations_map or {}).get(variants[i]["variantId"], []),
            n_ld_neighbors=int(n_neighbors[i]),
        ))
    candidates.sort(key=lambda c: -c.combined_score)
    return candidates


# ===================================================================
# Design C: Graph-Regularized Sparse Deconvolution (GRSD)
# ===================================================================

def graph_regularized_sparse_finemap(
    conn: "GraphGWASConnection",
    chr: str,
    lead_pos: int,
    window: int = 500_000,
    lambda_ridge: float = 0.1,
    lambda_graph: float = 0.5,
    r2_edge_threshold: float = 0.1,
    prior_var: float = 0.04,
    credible_set_coverage: float = 0.95,
    verbose: bool = True,
) -> list[FinemapCandidate]:
    """GRSD: Graph-Regularized Sparse Deconvolution.

    Solves: β* = argmin ||z - R·β||² + λ_r·||β||² + λ_g·β^T·L·β
    where L is the LD graph Laplacian (encourages different betas for
    LD-connected variants = anti-smoothing / sharpening).

    This has a closed-form solution: β* = (R^T R + λ_r I + λ_g L)^{-1} R^T z
    PIPs computed via Wakefield approximate Bayes factors on deconvolved betas.
    """
    all_idx = get_all_indices(conn)
    pheno = get_phenotype_values(conn, all_idx)

    if verbose:
        print(f"GRSD Fine-Mapping: {chr}:{lead_pos} ±{window // 1000}kb")

    variants = _load_locus_variants(conn, chr, lead_pos, window, all_idx)
    n_var = len(variants)
    if n_var < 3:
        return []

    R = _compute_ld_correlation(variants)
    z, beta_marginal, se_marginal = _compute_z_scores(variants, pheno)
    R_sq = R ** 2
    n_samples = int((~np.isnan(pheno)).sum())

    if verbose:
        print(f"  {n_var} variants, max |z|={np.max(np.abs(z)):.1f}")

    # Build LD graph Laplacian
    # Adjacency: A[i,j] = r²(i,j) if > threshold, else 0
    A = np.where(R_sq > r2_edge_threshold, R_sq, 0)
    np.fill_diagonal(A, 0)
    D_diag = A.sum(axis=1)
    L = np.diag(D_diag) - A  # graph Laplacian

    n_edges = int((A > 0).sum() // 2)
    if verbose:
        print(f"  {n_edges} LD edges, λ_r={lambda_ridge}, λ_g={lambda_graph}")

    # Closed-form solution: β* = (R^T R + λ_r I + λ_g L)^{-1} R^T z
    RtR = R.T @ R
    Q = RtR + lambda_ridge * np.eye(n_var) + lambda_graph * L
    Rtz = R.T @ z

    try:
        Q_factor = cho_factor(Q)
        beta_deconv = cho_solve(Q_factor, Rtz)
    except np.linalg.LinAlgError:
        Q = Q + 0.1 * np.eye(n_var)
        Q_factor = cho_factor(Q)
        beta_deconv = cho_solve(Q_factor, Rtz)

    # Wakefield approximate Bayes factors for PIPs
    safe_se = np.where(se_marginal > 1e-10, se_marginal, 1.0)
    se2 = safe_se ** 2
    W = prior_var
    z_deconv = beta_deconv / safe_se
    ratio = W / (se2 + W)
    log_bf = 0.5 * np.log(se2 / (se2 + W)) + 0.5 * z_deconv**2 * ratio
    pip = _softmax(log_bf)

    # -log10(p) for z_stat field
    p_vals = 2 * sp_stats.t.sf(np.abs(z), n_samples - 2)
    z_stat_log = -np.log10(np.maximum(p_vals, 1e-300))

    n_nb = np.array([(R_sq[i] > r2_edge_threshold).sum() - 1 for i in range(n_var)])

    return _build_candidates(
        variants, pip, z_stat_log, np.zeros(n_var), np.abs(beta_deconv),
        log_bf, n_nb, chr, credible_set_coverage)


# ===================================================================
# Design B: Hierarchical Belief Propagation (HBP)
# ===================================================================

def hierarchical_bp_finemap(
    conn: "GraphGWASConnection",
    chr: str,
    lead_pos: int,
    window: int = 500_000,
    n_rounds: int = 5,
    alpha: float = 0.6,
    damping: float = 0.5,
    r2_smooth: float = 0.3,
    credible_set_coverage: float = 0.95,
    eqtl_cache: dict | None = None,
    verbose: bool = True,
) -> list[FinemapCandidate]:
    """HBP: Hierarchical Belief Propagation fine-mapping.

    Models variant–gene–pathway as a 3-layer factor graph.
    Messages propagate upward (evidence aggregation) then downward
    (prior refinement) through the biological hierarchy.
    """
    all_idx = get_all_indices(conn)
    pheno = get_phenotype_values(conn, all_idx)

    if verbose:
        print(f"HBP Fine-Mapping: {chr}:{lead_pos} ±{window // 1000}kb")

    variants = _load_locus_variants(conn, chr, lead_pos, window, all_idx)
    n_var = len(variants)
    if n_var < 3:
        return []

    z_stats = _compute_association_stats(variants, pheno)
    R_sq = _compute_ld_matrix(variants)
    unique_stats, n_nb = _ld_deconvolve(z_stats, R_sq, r2_smooth)

    # --- Build bipartite matrices from Neo4j ---
    vids = [v["variantId"] for v in variants]
    vid_idx = {vid: i for i, vid in enumerate(vids)}

    result = conn.execute_read("""
        UNWIND $vids AS vid
        MATCH (v:Variant {variantId: vid})-[:HAS_CONSEQUENCE]->(g:Gene)
        OPTIONAL MATCH (g)-[:IN_PATHWAY]->(p:Pathway)
        OPTIONAL MATCH (g)-[:INTERACTS_WITH]-(g2:Gene)
        RETURN vid, g.symbol AS gene,
               collect(DISTINCT p.name) AS pathways,
               collect(DISTINCT g2.symbol) AS ppi_neighbors
    """, {"vids": vids})

    # Collect genes and pathways
    gene_set = set()
    pathway_set = set()
    vid_genes = {}  # vid -> [(gene, [pathways], [ppi])]
    for rec in result:
        vid = rec["vid"]
        gene = rec["gene"]
        if not gene:
            continue
        gene_set.add(gene)
        pathways = [p for p in rec.get("pathways", []) if p]
        pathway_set.update(pathways)
        ppi = [g for g in rec.get("ppi_neighbors", []) if g]
        vid_genes.setdefault(vid, []).append((gene, pathways, ppi))

    genes = sorted(gene_set)
    pathways = sorted(pathway_set)
    n_genes = len(genes)
    n_pathways = len(pathways)
    gene_idx = {g: i for i, g in enumerate(genes)}
    pw_idx = {p: i for i, p in enumerate(pathways)}

    if verbose:
        print(f"  {n_var} variants, {n_genes} genes, {n_pathways} pathways")

    if n_genes == 0:
        # No gene annotations — fall back to L1-style
        pip = _softmax(unique_stats)
        return _build_candidates(
            variants, pip, z_stats, np.zeros(n_var), unique_stats,
            unique_stats, n_nb, chr, credible_set_coverage)

    # B_vg: variant-gene matrix (weighted by eQTL if available)
    B_vg = np.zeros((n_var, n_genes))
    for vid, entries in vid_genes.items():
        i = vid_idx.get(vid)
        if i is None:
            continue
        eqtl_weight = 1.0
        if eqtl_cache:
            eq = eqtl_cache.get(vid)
            if eq:
                eqtl_weight = 1.0 + np.log1p(eq.get("composite", 0))
        for gene, _, _ in entries:
            j = gene_idx[gene]
            B_vg[i, j] = max(B_vg[i, j], eqtl_weight)

    # B_gp: gene-pathway matrix
    B_gp = np.zeros((n_genes, n_pathways))
    for vid, entries in vid_genes.items():
        for gene, pws, _ in entries:
            gi = gene_idx[gene]
            for pw in pws:
                pi = pw_idx[pw]
                B_gp[gi, pi] = 1.0

    # W_gg: gene-gene PPI adjacency (within locus genes only)
    W_gg = np.zeros((n_genes, n_genes))
    for vid, entries in vid_genes.items():
        for gene, _, ppis in entries:
            gi = gene_idx[gene]
            for ppi_gene in ppis:
                if ppi_gene in gene_idx:
                    gj = gene_idx[ppi_gene]
                    W_gg[gi, gj] = 1.0
                    W_gg[gj, gi] = 1.0

    # --- Message passing ---
    variant_belief = _softmax(unique_stats)

    for rnd in range(n_rounds):
        # UPWARD: variant → gene → pathway
        gene_score = B_vg.T @ variant_belief
        # PPI diffusion (one step)
        if W_gg.any():
            row_sum = W_gg.sum(axis=1)
            row_sum = np.where(row_sum > 0, row_sum, 1.0)
            W_norm = W_gg / row_sum[:, np.newaxis]
            gene_score = 0.7 * gene_score + 0.3 * (W_norm @ gene_score)
        pathway_score = B_gp.T @ gene_score

        # DOWNWARD: pathway → gene → variant
        gene_prior = B_gp @ pathway_score
        variant_prior = B_vg @ gene_prior

        # Normalize
        vp_sum = variant_prior.sum()
        if vp_sum > 0:
            variant_prior = variant_prior / vp_sum

        # Combine with statistical evidence
        new_belief = alpha * _softmax(unique_stats) + (1 - alpha) * variant_prior
        new_belief = new_belief / new_belief.sum()

        # Damping
        variant_belief = damping * variant_belief + (1 - damping) * new_belief
        variant_belief = variant_belief / variant_belief.sum()

    pip = variant_belief

    # z_functional = variant_prior contribution (for reporting)
    z_func = variant_prior * z_stats.max() if variant_prior.max() > 0 else np.zeros(n_var)

    return _build_candidates(
        variants, pip, z_stats, z_func, unique_stats,
        alpha * unique_stats + (1 - alpha) * z_func,
        n_nb, chr, credible_set_coverage)


def _cache_chr_graph_structure(conn, chr_name: str) -> dict:
    """Pre-cache variant→gene→pathway graph for an entire chromosome.

    Returns dict: vid -> {"genes": [str], "pathways": [str], "ppi": [str]}
    Query runs once (~30s on 70M-node DB), then all HBP calls use the cache.
    """
    result = conn.execute_read("""
        MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene)
        WHERE v.chr = $chr
        OPTIONAL MATCH (g)-[:IN_PATHWAY]->(p:Pathway)
        OPTIONAL MATCH (g)-[:INTERACTS_WITH]-(g2:Gene)
        RETURN v.variantId AS vid, g.symbol AS gene,
               collect(DISTINCT p.name) AS pathways,
               collect(DISTINCT g2.symbol) AS ppi
    """, {"chr": chr_name})

    cache = {}
    for rec in result:
        vid = rec["vid"]
        gene = rec["gene"]
        if not gene:
            continue
        entry = cache.setdefault(vid, {"genes": [], "pathways": [], "ppi": []})
        if gene not in entry["genes"]:
            entry["genes"].append(gene)
        for p in rec.get("pathways", []):
            if p and p not in entry["pathways"]:
                entry["pathways"].append(p)
        for g in rec.get("ppi", []):
            if g and g not in entry["ppi"]:
                entry["ppi"].append(g)
    return cache


def fast_hbp_finemap(
    variants: list[dict],
    phenotype: np.ndarray,
    graph_cache: dict,
    n_rounds: int = 5,
    alpha: float = 0.6,
    damping: float = 0.5,
    r2_smooth: float = 0.3,
    credible_set_coverage: float = 0.95,
    eqtl_cache: dict | None = None,
    chr_name: str = "chr22",
) -> list[FinemapCandidate]:
    """Fast HBP using pre-loaded variants and pre-cached graph structure.

    Avoids Neo4j queries entirely. Runs in ~0.1s per locus.
    """
    n_var = len(variants)
    if n_var < 3:
        return []

    z_stats = _compute_association_stats(variants, phenotype)
    R_sq = _compute_ld_matrix(variants)
    unique_stats, n_nb = _ld_deconvolve(z_stats, R_sq, r2_smooth)

    vids = [v["variantId"] for v in variants]
    vid_idx = {vid: i for i, vid in enumerate(vids)}

    # Build bipartite matrices from graph_cache
    gene_set = set()
    pathway_set = set()
    vid_genes = {}
    for vid in vids:
        info = graph_cache.get(vid)
        if not info:
            continue
        for gene in info.get("genes", []):
            gene_set.add(gene)
            vid_genes.setdefault(vid, []).append(
                (gene, info.get("pathways", []), info.get("ppi", [])))
        pathway_set.update(info.get("pathways", []))

    genes = sorted(gene_set)
    pathways = sorted(pathway_set)
    n_genes = len(genes)
    n_pathways = len(pathways)

    if n_genes == 0:
        pip = _softmax(unique_stats)
        return _build_candidates(
            variants, pip, z_stats, np.zeros(n_var), unique_stats,
            unique_stats, n_nb, chr_name, credible_set_coverage)

    gene_idx = {g: i for i, g in enumerate(genes)}
    pw_idx = {p: i for i, p in enumerate(pathways)}

    B_vg = np.zeros((n_var, n_genes))
    for vid, entries in vid_genes.items():
        i = vid_idx.get(vid)
        if i is None:
            continue
        eqtl_weight = 1.0
        if eqtl_cache:
            eq = eqtl_cache.get(vid)
            if eq:
                eqtl_weight = 1.0 + np.log1p(eq.get("composite", 0))
        for gene, _, _ in entries:
            j = gene_idx[gene]
            B_vg[i, j] = max(B_vg[i, j], eqtl_weight)

    B_gp = np.zeros((n_genes, n_pathways))
    for vid, entries in vid_genes.items():
        for gene, pws, _ in entries:
            gi = gene_idx[gene]
            for pw in pws:
                if pw in pw_idx:
                    B_gp[gi, pw_idx[pw]] = 1.0

    W_gg = np.zeros((n_genes, n_genes))
    for vid, entries in vid_genes.items():
        for gene, _, ppis in entries:
            gi = gene_idx[gene]
            for ppi_gene in ppis:
                if ppi_gene in gene_idx:
                    gj = gene_idx[ppi_gene]
                    W_gg[gi, gj] = 1.0
                    W_gg[gj, gi] = 1.0

    # Message passing
    variant_belief = _softmax(unique_stats)
    variant_prior = np.zeros(n_var)

    for rnd in range(n_rounds):
        gene_score = B_vg.T @ variant_belief
        if W_gg.any():
            row_sum = W_gg.sum(axis=1)
            row_sum = np.where(row_sum > 0, row_sum, 1.0)
            W_norm = W_gg / row_sum[:, np.newaxis]
            gene_score = 0.7 * gene_score + 0.3 * (W_norm @ gene_score)
        pathway_score = B_gp.T @ gene_score

        gene_prior = B_gp @ pathway_score
        variant_prior = B_vg @ gene_prior
        vp_sum = variant_prior.sum()
        if vp_sum > 0:
            variant_prior = variant_prior / vp_sum

        new_belief = alpha * _softmax(unique_stats) + (1 - alpha) * variant_prior
        new_belief = new_belief / new_belief.sum()
        variant_belief = damping * variant_belief + (1 - damping) * new_belief
        variant_belief = variant_belief / variant_belief.sum()

    pip = variant_belief
    z_func = variant_prior * z_stats.max() if variant_prior.max() > 0 else np.zeros(n_var)

    return _build_candidates(
        variants, pip, z_stats, z_func, unique_stats,
        alpha * unique_stats + (1 - alpha) * z_func,
        n_nb, chr_name, credible_set_coverage)


# ===================================================================
# Sumstats-only entry paths (for Pan-UKB and other summary-stat sources)
# ===================================================================

def hbp_finemap_from_sumstats(
    variants: list[dict],
    z_stats: np.ndarray,
    R_sq: np.ndarray,
    graph_cache: dict,
    n_rounds: int = 5,
    alpha: float = 0.6,
    damping: float = 0.5,
    r2_smooth: float = 0.3,
    credible_set_coverage: float = 0.95,
    eqtl_cache: dict | None = None,
    chr_name: str = "chr22",
) -> list[FinemapCandidate]:
    """HBP fine-mapping from pre-computed z-scores and LD matrix.

    Use when individual dosages are not available (e.g. Pan-UKB or any
    other summary-statistics release). Downstream logic is identical to
    `fast_hbp_finemap`; only the statistical inputs are externally supplied.

    Args:
        variants: list of dicts with at least `variantId`, `chr`, `pos`;
            optionally `ref`, `alt`, `af_total`.
        z_stats: length-n array of per-variant z-scores (BETA / SE).
        R_sq: n × n squared-correlation (r²) matrix aligned to `variants`.
            If your LD source gives signed correlation, pass `r**2`.
        graph_cache: {variantId: {genes, pathways, ppi}} — the HBP
            message-passing adjacency source, same format used by
            `fast_hbp_finemap`.
        Other args: see `fast_hbp_finemap`.

    Returns:
        list of FinemapCandidate sorted by combined score.
    """
    n_var = len(variants)
    if n_var < 3:
        return []
    z_stats = np.asarray(z_stats, dtype=np.float64)
    R_sq = np.asarray(R_sq, dtype=np.float64)
    assert z_stats.shape == (n_var,), "z_stats length must match variants"
    assert R_sq.shape == (n_var, n_var), "R_sq must be n × n"

    unique_stats, n_nb = _ld_deconvolve(z_stats, R_sq, r2_smooth)

    vids = [v["variantId"] for v in variants]
    vid_idx = {vid: i for i, vid in enumerate(vids)}

    gene_set: set = set()
    pathway_set: set = set()
    vid_genes: dict[str, list] = {}
    for vid in vids:
        info = graph_cache.get(vid)
        if not info:
            continue
        for gene in info.get("genes", []):
            gene_set.add(gene)
            vid_genes.setdefault(vid, []).append(
                (gene, info.get("pathways", []), info.get("ppi", []))
            )
        pathway_set.update(info.get("pathways", []))

    genes = sorted(gene_set)
    pathways = sorted(pathway_set)
    n_genes = len(genes)
    n_pathways = len(pathways)

    if n_genes == 0:
        pip = _softmax(unique_stats)
        return _build_candidates(
            variants, pip, z_stats, np.zeros(n_var), unique_stats,
            unique_stats, n_nb, chr_name, credible_set_coverage,
        )

    gene_idx = {g: i for i, g in enumerate(genes)}
    pw_idx = {p: i for i, p in enumerate(pathways)}

    B_vg = np.zeros((n_var, n_genes))
    for vid, entries in vid_genes.items():
        i = vid_idx.get(vid)
        if i is None:
            continue
        eqtl_weight = 1.0
        if eqtl_cache:
            eq = eqtl_cache.get(vid)
            if eq:
                eqtl_weight = 1.0 + np.log1p(eq.get("composite", 0))
        for gene, _, _ in entries:
            j = gene_idx[gene]
            B_vg[i, j] = max(B_vg[i, j], eqtl_weight)

    B_gp = np.zeros((n_genes, n_pathways))
    for vid, entries in vid_genes.items():
        for gene, pws, _ in entries:
            gi = gene_idx[gene]
            for pw in pws:
                if pw in pw_idx:
                    B_gp[gi, pw_idx[pw]] = 1.0

    W_gg = np.zeros((n_genes, n_genes))
    for vid, entries in vid_genes.items():
        for gene, _, ppis in entries:
            gi = gene_idx[gene]
            for ppi_gene in ppis:
                if ppi_gene in gene_idx:
                    gj = gene_idx[ppi_gene]
                    W_gg[gi, gj] = 1.0
                    W_gg[gj, gi] = 1.0

    variant_belief = _softmax(unique_stats)
    variant_prior = np.zeros(n_var)

    for _ in range(n_rounds):
        gene_score = B_vg.T @ variant_belief
        if W_gg.any():
            row_sum = W_gg.sum(axis=1)
            row_sum = np.where(row_sum > 0, row_sum, 1.0)
            W_norm = W_gg / row_sum[:, np.newaxis]
            gene_score = 0.7 * gene_score + 0.3 * (W_norm @ gene_score)
        pathway_score = B_gp.T @ gene_score

        gene_prior = B_gp @ pathway_score
        variant_prior = B_vg @ gene_prior
        vp_sum = variant_prior.sum()
        if vp_sum > 0:
            variant_prior = variant_prior / vp_sum

        new_belief = alpha * _softmax(unique_stats) + (1 - alpha) * variant_prior
        new_belief = new_belief / new_belief.sum()
        variant_belief = damping * variant_belief + (1 - damping) * new_belief
        variant_belief = variant_belief / variant_belief.sum()

    pip = variant_belief
    z_func = variant_prior * z_stats.max() if variant_prior.max() > 0 else np.zeros(n_var)

    # Phase 0 schema extension: per-variant continuous prior_score
    # (CADD, MotifBreakR Δ-PWM, DeepSEA effect, ABC contact strength, ...)
    # If any variant has prior_score, re-mix into the final belief and z_func.
    has_prior = False
    prior_extra = np.zeros(n_var)
    for i, vid in enumerate(vids):
        info = graph_cache.get(vid)
        if info is not None:
            ps = info.get("prior_score")
            if ps is not None:
                prior_extra[i] = float(ps)
                has_prior = True
    if has_prior and prior_extra.max() > 0:
        prior_extra = prior_extra / prior_extra.max() * z_stats.max()
        z_func = 0.7 * z_func + 0.3 * prior_extra
        # Also re-mix into PIP so the prior actually changes the credible set
        boosted = alpha * unique_stats + (1 - alpha) * (z_func + 0.5 * prior_extra)
        new_pip = _softmax(boosted)
        # Blend with BP-derived belief so we don't fully overwrite
        pip = 0.5 * pip + 0.5 * new_pip
        pip = pip / pip.sum()

    return _build_candidates(
        variants, pip, z_stats, z_func, unique_stats,
        alpha * unique_stats + (1 - alpha) * z_func,
        n_nb, chr_name, credible_set_coverage,
    )


def l1_finemap_from_sumstats(
    variants: list[dict],
    z_stats: np.ndarray,
    R_sq: np.ndarray,
    z_func: np.ndarray | None = None,
    alpha: float = 0.5,
    r2_smooth: float = 0.3,
    credible_set_coverage: float = 0.95,
    chr_name: str = "chr22",
    annotations_map: dict | None = None,
    graph_cache: dict | None = None,
) -> list[FinemapCandidate]:
    """GAFM (Graph-Augmented Fine-Mapping) from pre-computed z-scores and LD matrix.

    Paper-facing name: GAFM. Python prefix: l1 (historical, preserved for
    backward compatibility in JSON keys and benchmark scripts).

    Summary-stats entry point mirroring `dual_graph_finemap` without the
    Neo4j/phenotype dependencies. The functional annotation score
    `z_func` is supplied by the caller (compute once from the graph or
    from any other annotation source).

    Args:
        variants: list of dicts with at least `variantId`, `chr`, `pos`;
            optionally `ref`, `alt`, `af_total`.
        z_stats: length-n array of per-variant z-scores.
        R_sq: n × n squared-correlation (r²) matrix.
        z_func: length-n array of functional annotation scores, already
            rescaled to match z_stats's dynamic range. Pass None or zeros
            to fall back to a pure-statistical softmax (α = 1 effectively).
        alpha: blend of statistical (α) vs functional (1−α). Default 0.5.
        r2_smooth: LD deconvolution threshold (default 0.3).
        credible_set_coverage: credible-set PIP sum target (default 0.95).
        chr_name: chromosome label for output.
        annotations_map: {variantId: [annotation_strings]} (optional).

    Returns:
        list of FinemapCandidate sorted by combined score.
    """
    n_var = len(variants)
    if n_var < 3:
        return []
    z_stats = np.asarray(z_stats, dtype=np.float64)
    R_sq = np.asarray(R_sq, dtype=np.float64)
    assert z_stats.shape == (n_var,), "z_stats length must match variants"
    assert R_sq.shape == (n_var, n_var), "R_sq must be n × n"

    if z_func is None:
        z_func = np.zeros(n_var)
    else:
        z_func = np.asarray(z_func, dtype=np.float64)
        assert z_func.shape == (n_var,), "z_func length must match variants"

    # Phase 0 schema extension: per-variant prior_score from graph_cache
    # (CADD, MotifBreakR ΔPWM, DeepSEA effect, ABC contact strength, ...)
    if graph_cache:
        prior_extra = np.zeros(n_var)
        has = False
        for i, v in enumerate(variants):
            info = graph_cache.get(v.get("variantId"))
            if info is not None:
                ps = info.get("prior_score")
                if ps is not None:
                    prior_extra[i] = float(ps)
                    has = True
        if has and prior_extra.max() > 0:
            prior_extra = prior_extra / prior_extra.max() * z_stats.max()
            z_func = z_func + 0.5 * prior_extra

    unique_stats, n_nb = _ld_deconvolve(z_stats, R_sq, r2_smooth)
    combined = alpha * unique_stats + (1 - alpha) * z_func
    pip = _softmax(combined)

    return _build_candidates(
        variants, pip, z_stats, z_func, unique_stats, combined,
        n_nb, chr_name, credible_set_coverage,
        annotations_map=annotations_map,
    )


# ===================================================================
# Design A: Cross-Locus Graph Fine-Mapping (CLGF)
# ===================================================================

def cross_locus_graph_finemap(
    conn: "GraphGWASConnection",
    loci: list[dict],
    phenotype: np.ndarray,
    all_idx: np.ndarray,
    n_iterations: int = 3,
    pathway_weight: float = 0.3,
    alpha: float = 0.7,
    r2_smooth: float = 0.3,
    credible_set_coverage: float = 0.95,
    eqtl_cache: dict | None = None,
    verbose: bool = True,
) -> dict[str, list[FinemapCandidate]]:
    """CLGF: Cross-Locus Graph Fine-Mapping.

    Shares evidence across multiple GWAS loci via the biological pathway
    graph. Loci connected through shared pathways reinforce each other's
    causal variant identification — information that single-locus methods
    like SuSiE and FINEMAP structurally cannot use.

    Args:
        loci: list of dicts with keys 'chr', 'lead_pos', 'window'.
        phenotype: pre-loaded phenotype vector aligned to all_idx.
        all_idx: sample packed_index array.
        n_iterations: EM rounds for cross-locus refinement.
        pathway_weight: weight for cross-locus pathway prior vs local annotation.
        alpha: balance statistical vs functional (higher = more statistical).

    Returns:
        dict mapping locus key ("chr:pos") to sorted FinemapCandidate list.
    """
    if verbose:
        print(f"CLGF: Cross-Locus Graph Fine-Mapping ({len(loci)} loci)")

    # --- Phase 1: Initialize each locus ---
    locus_data = {}
    for loc in loci:
        key = f"{loc['chr']}:{loc['lead_pos']}"
        variants = _load_locus_variants(
            conn, loc["chr"], loc["lead_pos"],
            loc.get("window", 500_000), all_idx)
        if len(variants) < 3:
            continue

        z_stats = _compute_association_stats(variants, phenotype)
        R_sq = _compute_ld_matrix(variants)
        unique, n_nb = _ld_deconvolve(z_stats, R_sq, r2_smooth)

        # eQTL annotation scores
        z_func = np.zeros(len(variants))
        if eqtl_cache:
            for i, v in enumerate(variants):
                eq = eqtl_cache.get(v["variantId"])
                if eq:
                    score = eq.get("composite", 0)
                    nt = eq.get("n_tissues", 0)
                    if 1 <= nt <= 3:
                        score *= 1.5
                    elif nt >= 8:
                        score *= 0.3
                    z_func[i] = np.log1p(score)
            if z_func.max() > 0:
                z_func = z_func / z_func.max() * z_stats.max()

        pip = _softmax(unique)

        # Query variant→pathway mapping
        vids = [v["variantId"] for v in variants]
        result = conn.execute_read("""
            UNWIND $vids AS vid
            MATCH (v:Variant {variantId: vid})-[:HAS_CONSEQUENCE]->(g:Gene)
            OPTIONAL MATCH (g)-[:IN_PATHWAY]->(p:Pathway)
            RETURN vid, collect(DISTINCT p.name) AS pathways
        """, {"vids": vids})

        vid_pathways = {}
        for rec in result:
            pws = [p for p in rec.get("pathways", []) if p]
            if pws:
                vid_pathways[rec["vid"]] = pws

        locus_data[key] = {
            "variants": variants, "z_stats": z_stats, "unique": unique,
            "z_func": z_func, "n_nb": n_nb, "pip": pip,
            "vid_pathways": vid_pathways,
        }

    if verbose:
        n_with_pw = sum(1 for ld in locus_data.values() if ld["vid_pathways"])
        print(f"  {len(locus_data)} loci loaded, {n_with_pw} with pathway info")

    # --- Phase 2: Iterative cross-locus refinement ---
    for iteration in range(n_iterations):
        # E-step: Accumulate pathway evidence across ALL loci
        pathway_scores = {}
        for key, ld in locus_data.items():
            for i, v in enumerate(ld["variants"]):
                vid = v["variantId"]
                for pw in ld["vid_pathways"].get(vid, []):
                    pathway_scores.setdefault(pw, 0.0)
                    pathway_scores[pw] += ld["pip"][i]

        # M-step: Update PIPs with cross-locus prior
        for key, ld in locus_data.items():
            n_var = len(ld["variants"])
            cross_prior = np.zeros(n_var)

            for i, v in enumerate(ld["variants"]):
                vid = v["variantId"]
                for pw in ld["vid_pathways"].get(vid, []):
                    total = pathway_scores.get(pw, 0)
                    # Subtract own-locus contribution to prevent self-reinforcement
                    own = sum(
                        ld["pip"][j]
                        for j in range(n_var)
                        if pw in ld["vid_pathways"].get(
                            ld["variants"][j]["variantId"], [])
                    )
                    cross_prior[i] += max(total - own, 0)

            if cross_prior.max() > 0:
                cross_prior = cross_prior / cross_prior.max()

            # Combined score
            func_score = (
                pathway_weight * cross_prior
                + (1 - pathway_weight) * (ld["z_func"] / max(ld["z_func"].max(), 1e-10))
            )
            combined = alpha * ld["unique"] + (1 - alpha) * func_score * ld["z_stats"].max()

            ld["pip"] = _softmax(combined)
            ld["combined"] = combined
            ld["cross_prior"] = cross_prior

        if verbose:
            top_pw = sorted(pathway_scores.items(), key=lambda x: -x[1])[:3]
            print(f"  Iteration {iteration + 1}: top pathways = "
                  + ", ".join(f"{p}({s:.2f})" for p, s in top_pw))

    # --- Phase 3: Build results ---
    results = {}
    for key, ld in locus_data.items():
        chr_name = key.split(":")[0]
        z_func_final = ld.get("cross_prior", np.zeros(len(ld["variants"])))
        combined = ld.get("combined", ld["unique"])
        results[key] = _build_candidates(
            ld["variants"], ld["pip"], ld["z_stats"], z_func_final,
            ld["unique"], combined, ld["n_nb"], chr_name,
            credible_set_coverage)

    return results


# ===================================================================
# TSV output
# ===================================================================

def finemap_to_tsv(results: list[FinemapCandidate], path: str):
    """Write fine-mapping results to TSV."""
    import csv
    cols = ["variant_id", "chr", "pos", "ref", "alt", "af",
            "z_stat", "z_functional", "unique_stat", "combined_score",
            "pip", "in_credible_set", "annotations", "n_ld_neighbors"]

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        writer.writeheader()
        for r in results:
            writer.writerow({
                "variant_id": r.variant_id,
                "chr": r.chr,
                "pos": r.pos,
                "ref": r.ref,
                "alt": r.alt,
                "af": f"{r.af:.4f}",
                "z_stat": f"{r.z_stat:.4f}",
                "z_functional": f"{r.z_functional:.4f}",
                "unique_stat": f"{r.unique_stat:.4f}",
                "combined_score": f"{r.combined_score:.6f}",
                "pip": f"{r.pip:.6f}",
                "in_credible_set": r.in_credible_set,
                "annotations": ";".join(r.annotations),
                "n_ld_neighbors": r.n_ld_neighbors,
            })
