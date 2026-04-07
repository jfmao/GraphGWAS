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
import numpy as np
from scipy import stats as sp_stats
from scipy.linalg import eigh

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


def _compute_ld_matrix(variants: list[dict]) -> np.ndarray:
    """Compute pairwise r² matrix for locus variants."""
    n = len(variants)
    R = np.zeros((n, n))
    for i in range(n):
        d_i = variants[i]["dosage"]
        valid_i = ~np.isnan(d_i)
        for j in range(i, n):
            d_j = variants[j]["dosage"]
            valid = valid_i & ~np.isnan(d_j)
            if valid.sum() < 20:
                continue
            r = np.corrcoef(d_i[valid], d_j[valid])[0, 1]
            R[i, j] = R[j, i] = r ** 2
    np.fill_diagonal(R, 1.0)
    return R


def _compute_association_stats(
    variants: list[dict],
    phenotype: np.ndarray,
) -> np.ndarray:
    """Compute -log10(p) for each variant via simple regression."""
    n_var = len(variants)
    z_stats = np.zeros(n_var)

    for i, v in enumerate(variants):
        d = v["dosage"]
        valid = ~np.isnan(d) & ~np.isnan(phenotype)
        n = valid.sum()
        if n < 20 or np.std(d[valid]) == 0:
            z_stats[i] = 0
            continue

        g = d[valid]
        y = phenotype[valid]
        # Simple linear regression
        X = np.column_stack([np.ones(n), g])
        try:
            beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
            resid = y - X @ beta
            mse = np.sum(resid**2) / (n - 2)
            se = np.sqrt(mse * np.linalg.inv(X.T @ X)[1, 1])
            t = beta[1] / se
            p = 2 * sp_stats.t.sf(abs(t), n - 2)
            z_stats[i] = -np.log10(max(p, 1e-300))
        except Exception:
            z_stats[i] = 0

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

    # Step 5: Functional annotation score via graph traversal
    z_func = np.zeros(n_var)

    # Query: which variants fall in annotated genes? which genes are in pathways?
    for i, v in enumerate(variants):
        vid = v["variantId"]
        result = conn.execute_read("""
            MATCH (v:Variant {variantId: $vid})-[:HAS_CONSEQUENCE]->(g:Gene)
            OPTIONAL MATCH (g)-[:IN_PATHWAY]->(p:Pathway)
            RETURN g.symbol AS gene, collect(DISTINCT p.name) AS pathways
        """, {"vid": vid})

        annotations = []
        score = 0.0
        for rec in result:
            gene = rec.get("gene")
            pathways = rec.get("pathways", [])
            if gene:
                score += 1.0  # in a gene
                annotations.append(f"gene:{gene}")
            if pathways:
                score += len(pathways) * 0.5  # in pathway(s)
                for pw in pathways:
                    if pw:
                        annotations.append(f"pathway:{pw}")

        z_func[i] = np.log1p(score)  # log(1 + score) for diminishing returns
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
        print(f"\n  Top 5:")
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
                          f"primary" if is_primary else "secondary"],
            n_ld_neighbors=int(np.sum(R[i] > r2_threshold)) - 1,
        ))

    candidates.sort(key=lambda c: -c.pip)

    if verbose:
        cs_size = sum(1 for c in candidates if c.in_credible_set)
        print(f"  Credible set ({credible_set_coverage*100:.0f}%): {cs_size} variants")
        print(f"\n  Clusters:")
        for c, info in sorted(cluster_info.items()):
            best = variants[info["best_z"]]
            print(f"    Cluster {c}: {info['n_members']} variants, "
                  f"top={best['variantId'][:40]} z={z_stats[info['best_z']]:.2f}")
        print(f"\n  Top 5 candidates:")
        print(f"  {'Variant':<50s} {'z_stat':>7s} {'PIP':>6s} {'Cluster':>7s} {'CS':>3s}")
        print(f"  {'-'*75}")
        for c in candidates[:5]:
            cluster = [a for a in c.annotations if a.startswith("cluster:")][0]
            cs_mark = "***" if c.in_credible_set else ""
            print(f"  {c.variant_id:<50s} {c.z_stat:>7.2f} {c.pip:>6.4f} {cluster:>7s} {cs_mark:>3s}")

    return candidates


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
