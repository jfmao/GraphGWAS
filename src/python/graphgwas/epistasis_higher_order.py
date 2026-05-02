"""Higher-order epistasis detection on the variant-gene bipartite graph.

Implements **M5 Random Walk with Restart (RWR)** from
docs/EPISTASIS_V2_DESIGN.md.  The bipartite graph has variant nodes on
one side and gene nodes on the other, with edges drawn from the
graph_cache annotation (each variant linked to the genes whose body or
±2-kb flanking region it overlaps).

For a seed variant v_seed, the RWR stationary distribution over the
graph is the long-run probability of being at each node when, at every
step, the walker either follows a random outgoing edge (with prob 1-α)
or teleports back to v_seed (with prob α).  Variants with high RWR
score under v_seed are those biologically linked to v_seed through
shared genes — the natural pool from which to enumerate epistatic
partners.

Pair candidate ranking: ``mutual_rwr(v_i, v_j) = p_{v_i}[v_j] * p_{v_j}[v_i]``,
the probability that walks rooted at each end visit the other end.
This is the "joint RWR" score from EPISTASIS_V2_DESIGN.md §M5.

Higher-order (k ≥ 3) extension: triplet score is the product of three
pairwise mutual_rwr values.  Quadruplets and beyond follow the same
construction at exponentially increasing pair-counts.

This module is the structurally-unique-claim kernel for paper #2
(`§Y.5`): pairwise BOOST/MAPIT/MDR baselines cannot formulate k ≥ 3
interactions in their canonical form, while RWR on the typed bipartite
graph extends naturally.

NOT YET IMPLEMENTED (deferred, see TODOs):
  - M6 GNN-attention extraction (needs PyTorch Geometric infrastructure)
  - M7 persistent-homology / topological epistasis (needs Gudhi or Giotto-TDA)
"""
from __future__ import annotations

import json
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np


# ===========================================================================
# Bipartite graph from graph_cache
# ===========================================================================

def build_bipartite_adjacency(
    graph_cache_path: Path,
    variant_ids_to_keep: list[str] | None = None,
) -> tuple[np.ndarray, list[str], list[str]]:
    """Load variant→gene bipartite adjacency from a graph_cache JSON.

    Args:
        graph_cache_path: path to a graph_cache JSON whose top-level keys
            are variant IDs (``chrN:pos:REF:ALT`` format) and values are
            dicts with a ``"genes"`` list.
        variant_ids_to_keep: optional whitelist; if given, only variants
            in this list are included.  Useful for restricting to MAF-
            filtered variants.

    Returns:
        (A, variant_ids, gene_ids) where:
            A : adjacency matrix, shape (n_variants, n_genes), values in {0, 1}
            variant_ids : list of variant ID strings, length n_variants
            gene_ids    : list of gene symbols, length n_genes
    """
    cache = json.loads(graph_cache_path.read_text())
    if variant_ids_to_keep is not None:
        keep = set(variant_ids_to_keep)
        cache = {k: v for k, v in cache.items() if k in keep}

    gene_set: set[str] = set()
    for entry in cache.values():
        gene_set.update(entry.get("genes", []))
    gene_ids = sorted(gene_set)
    gene_to_col = {g: i for i, g in enumerate(gene_ids)}

    variant_ids = list(cache.keys())
    n_var = len(variant_ids)
    n_gene = len(gene_ids)

    A = np.zeros((n_var, n_gene), dtype=np.float64)
    for vi, vid in enumerate(variant_ids):
        for gene in cache[vid].get("genes", []):
            A[vi, gene_to_col[gene]] = 1.0
    return A, variant_ids, gene_ids


# ===========================================================================
# Random Walk with Restart on the bipartite graph
# ===========================================================================

def rwr_stationary(
    A: np.ndarray,
    seed_idx: int,
    alpha: float = 0.15,
    n_iter: int = 50,
    tol: float = 1e-6,
) -> np.ndarray:
    """Compute the RWR stationary distribution from one seed variant.

    The walk alternates variant→gene→variant on the bipartite graph; at
    each step it teleports back to the seed with probability α.  We
    iterate the matrix recurrence

        p_{t+1} = (1 - α) · M · p_t + α · e_seed

    until ‖p_{t+1} - p_t‖_∞ < tol or n_iter is reached.

    Args:
        A: variant→gene adjacency, shape (n_var, n_gene), 0/1 entries.
        seed_idx: index of the seed variant.
        alpha: restart probability in (0, 1].
        n_iter: maximum iterations.
        tol: convergence tolerance.

    Returns:
        p: stationary distribution over variants, shape (n_var,) summing to 1.
    """
    n_var, n_gene = A.shape
    if not (0.0 < alpha <= 1.0):
        raise ValueError(f"alpha must be in (0, 1]; got {alpha}")
    if n_var == 0:
        raise ValueError("Empty adjacency matrix")

    # Row-normalize variant→gene (each variant's outgoing weight sums to 1)
    var_deg = A.sum(axis=1, keepdims=True).clip(min=1e-12)
    P_vg = A / var_deg
    # Column-normalize gene→variant (each gene distributes uniformly over its variants)
    gene_deg = A.sum(axis=0, keepdims=True).clip(min=1e-12)
    P_gv = (A / gene_deg).T

    # Two-step transition matrix on variants (variant→gene→variant)
    M_vv = P_vg @ P_gv  # shape (n_var, n_var)

    # Seed vector
    e = np.zeros(n_var, dtype=np.float64)
    e[seed_idx] = 1.0

    p = e.copy()
    for _ in range(n_iter):
        p_new = (1.0 - alpha) * (M_vv @ p) + alpha * e
        if np.max(np.abs(p_new - p)) < tol:
            p = p_new
            break
        p = p_new

    # Re-normalize (matrix multiply may drift slightly)
    p = p / max(p.sum(), 1e-12)
    return p


# ===========================================================================
# Pair and triplet candidate generation
# ===========================================================================

def mutual_rwr_pair_scores(
    A: np.ndarray,
    seed_indices: list[int] | np.ndarray,
    alpha: float = 0.15,
    top_k: int = 1000,
    n_iter: int = 50,
) -> list[tuple[int, int, float]]:
    """Compute joint-RWR scores for variant pairs via the seed pool.

    For each (v_i, v_j) where both are in seed_indices, the score is::

        score(i, j) = p_{v_i}[v_j] * p_{v_j}[v_i]

    The top-k pairs are returned.

    Args:
        A: bipartite adjacency, (n_var, n_gene).
        seed_indices: variant indices to seed RWR from.  Restrict to a
            tractable set (e.g., variants with annotated genes) to avoid
            O(n_var^2) pair enumeration.
        alpha: RWR restart probability.
        top_k: number of top-scoring pairs to return.
        n_iter: max RWR iterations per seed.

    Returns:
        List of (i, j, score) tuples, sorted by score descending.
    """
    seeds = list(seed_indices)
    n_seed = len(seeds)
    print(f"  computing RWR for {n_seed} seeds...", flush=True)

    # Cache p_{v_i} for each seed
    P = np.zeros((n_seed, A.shape[0]), dtype=np.float64)
    for k, si in enumerate(seeds):
        P[k] = rwr_stationary(A, si, alpha=alpha, n_iter=n_iter)
        if (k + 1) % max(1, n_seed // 10) == 0:
            print(f"    seed {k+1}/{n_seed} done", flush=True)

    # Score every pair
    print(f"  scoring {n_seed*(n_seed-1)//2:,} pairs ...", flush=True)
    pair_scores: list[tuple[int, int, float]] = []
    for ki, kj in combinations(range(n_seed), 2):
        i, j = seeds[ki], seeds[kj]
        score = P[ki, j] * P[kj, i]
        if score > 0:
            pair_scores.append((i, j, float(score)))

    pair_scores.sort(key=lambda x: x[2], reverse=True)
    return pair_scores[:top_k]


def gene_pair_rwr_scores(
    A: np.ndarray,
    gene_seed_indices: list[int] | np.ndarray,
    alpha: float = 0.15,
    top_gene_pairs: int = 200,
    n_iter: int = 30,
) -> list[tuple[int, int, float]]:
    """Run RWR on the *gene* side of the bipartite graph and rank gene pairs.

    Symmetric counterpart to :func:`mutual_rwr_pair_scores`, but with seeds
    on genes rather than variants.  Walks gene→variant→gene; mutual score
    between two seeded genes (g_a, g_b) is

        score(a, b) = p_{g_a}[g_b] * p_{g_b}[g_a]

    Use when the biological hypothesis is a *gene-pair* relationship (e.g.,
    BCY1 × TPK1) and the variant-priority seed pool is too sparse to put
    both genes' variants in the top-K (the §Y.4 yeast failure mode).
    Combine with :func:`expand_gene_pairs_to_variant_pairs` to recover all
    variant-variant cross products for testing.

    Args:
        A: variant→gene adjacency, (n_var, n_gene).
        gene_seed_indices: gene indices to seed RWR from.  Pass all gene
            indices for full enumeration, or a subset for targeted scans.
        alpha: restart probability.
        top_gene_pairs: number of top gene pairs to return.
        n_iter: max iterations.

    Returns:
        List of (gene_i, gene_j, score) tuples, sorted by score descending.
    """
    n_var, n_gene = A.shape
    seeds = list(gene_seed_indices)
    n_seed = len(seeds)
    if n_seed == 0:
        return []
    if not (0.0 < alpha <= 1.0):
        raise ValueError(f"alpha must be in (0, 1]; got {alpha}")

    var_deg = A.sum(axis=1, keepdims=True).clip(min=1e-12)
    gene_deg = A.sum(axis=0, keepdims=True).clip(min=1e-12)
    P_vg = A / var_deg                            # variant→gene
    P_gv = (A / gene_deg).T                       # gene→variant
    # Two-step gene→variant→gene transition (n_gene × n_gene)
    M_gg = P_gv @ P_vg

    E = np.zeros((n_gene, n_seed), dtype=np.float64)
    for k, gi in enumerate(seeds):
        E[gi, k] = 1.0
    P = E.copy()
    for _ in range(n_iter):
        P = (1.0 - alpha) * (M_gg @ P) + alpha * E
    P = P / np.clip(P.sum(axis=0, keepdims=True), 1e-12, None)

    seed_arr = np.asarray(seeds)
    P_seed = P[seed_arr, :]                       # (n_seed, n_seed)
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


def expand_gene_pairs_to_variant_pairs(
    gene_pairs: list[tuple[int, int, float]],
    A: np.ndarray,
    variant_ids: list[str],
    *,
    max_variants_per_gene: int = 50,
) -> list[tuple[int, int, float]]:
    """Cross-product expansion: each gene pair → all (v ∈ g_a) × (v ∈ g_b) variant pairs.

    The variant indices come from ``A`` (the bipartite adjacency).  Pairs
    are deduplicated when the same variant appears in multiple gene pairs;
    the kept score is the maximum gene-pair score so far.

    Args:
        gene_pairs: output of :func:`gene_pair_rwr_scores`.
        A: bipartite adjacency, (n_var, n_gene).  Same matrix used to seed
            ``gene_pair_rwr_scores``.
        variant_ids: list aligned with rows of A (used only for logging).
        max_variants_per_gene: cap on per-gene variant fanout (avoids
            quadratic blow-up on hub genes).  When exceeded, randomly
            sample without replacement.

    Returns:
        List of (var_i, var_j, score) tuples — variant column indices into
        the dosage matrix that built ``A``.  Sorted by score descending.
    """
    rng = np.random.default_rng(0)  # deterministic per-pair sampling
    pair_score: dict[tuple[int, int], float] = {}
    for gene_a, gene_b, s in gene_pairs:
        vars_a = np.where(A[:, gene_a] > 0)[0]
        vars_b = np.where(A[:, gene_b] > 0)[0]
        if len(vars_a) > max_variants_per_gene:
            vars_a = rng.choice(vars_a, size=max_variants_per_gene, replace=False)
        if len(vars_b) > max_variants_per_gene:
            vars_b = rng.choice(vars_b, size=max_variants_per_gene, replace=False)
        for vi in vars_a:
            for vj in vars_b:
                if vi == vj:
                    continue
                key = (int(min(vi, vj)), int(max(vi, vj)))
                # Keep best score for duplicate variant pairs
                if key not in pair_score or pair_score[key] < s:
                    pair_score[key] = s
    out = [(i, j, s) for (i, j), s in pair_score.items()]
    out.sort(key=lambda t: t[2], reverse=True)
    return out


def mutual_rwr_triplet_scores(
    A: np.ndarray,
    seed_indices: list[int] | np.ndarray,
    alpha: float = 0.15,
    top_k: int = 200,
    pair_top_k: int = 5000,
    n_iter: int = 50,
) -> list[tuple[int, int, int, float]]:
    """Compute joint-RWR scores for variant triplets (k=3).

    Strategy: first generate the top-pair candidates via
    :func:`mutual_rwr_pair_scores`; then for each triplet (i, j, k) where
    at least two of {i,j}, {i,k}, {j,k} are in the top-pair list, compute
    the score as the product of the three pairwise scores.

    Args:
        A: bipartite adjacency.
        seed_indices: seeds for RWR.
        alpha: restart probability.
        top_k: number of top-scoring triplets to return.
        pair_top_k: how many top pairs to use as the seed for triplet
            enumeration.
        n_iter: max RWR iterations per seed.

    Returns:
        List of (i, j, k, score) tuples, sorted descending.
    """
    pairs = mutual_rwr_pair_scores(A, seed_indices, alpha=alpha,
                                    top_k=pair_top_k, n_iter=n_iter)
    pair_score = {(min(i, j), max(i, j)): s for i, j, s in pairs}

    # Build adjacency from top pairs to enumerate triangles
    adj: dict[int, set[int]] = defaultdict(set)
    for i, j, _ in pairs:
        adj[i].add(j)
        adj[j].add(i)

    print(f"  enumerating triangles in top-{pair_top_k} pair graph ...",
          flush=True)
    triplets: list[tuple[int, int, int, float]] = []
    seen: set[tuple[int, int, int]] = set()
    for i, neigh_i in adj.items():
        for j in neigh_i:
            if j <= i:
                continue
            common = neigh_i & adj[j]
            for k in common:
                if k <= j:
                    continue
                key = (i, j, k)
                if key in seen:
                    continue
                seen.add(key)
                s = (pair_score[(i, j)] *
                     pair_score[(min(i, k), max(i, k))] *
                     pair_score[(min(j, k), max(j, k))])
                triplets.append((i, j, k, s))

    triplets.sort(key=lambda x: x[3], reverse=True)
    return triplets[:top_k]


# ===========================================================================
# Bridge to graphgwas.bias: build Z from M5 candidate pairs
# ===========================================================================

def m5_interaction_matrix(
    pair_scores: list[tuple[int, int, float]],
    dosages: np.ndarray,
    *,
    top_k: int | None = None,
) -> np.ndarray:
    """Build interaction-feature matrix Z from M5 candidate pairs.

    Equivalent to :func:`graphgwas.epistasis_v2.motif_interaction_matrix`
    but seeded by M5 RWR scores instead of M2 motif memberships.  The
    resulting Z is a graph-typed interaction subspace of dimension
    ``min(top_k, len(pair_scores))``, suitable for ``rho_max(g, Z)``.

    Args:
        pair_scores: output of :func:`mutual_rwr_pair_scores`.
        dosages: shape (n_samples, n_variants).  Indices in ``pair_scores``
            must reference columns of this matrix.
        top_k: keep only the top-k pairs (None ⇒ keep all).

    Returns:
        Z, shape (n_samples, n_pairs_kept), centred per-column.
    """
    pairs = pair_scores if top_k is None else pair_scores[:top_k]
    if not pairs:
        return np.empty((dosages.shape[0], 0), dtype=np.float64)
    cols = []
    for i, j, _ in pairs:
        prod = dosages[:, i] * dosages[:, j]
        cols.append(prod - prod.mean())
    return np.column_stack(cols)


# ===========================================================================
# Higher-order placeholders (M6, M7)
# ===========================================================================

def gnn_attention_pairs(*args, **kwargs):
    """M6 GNN attention-based interaction extraction.

    [TODO paper2-Y5-M6]: train a heterogeneous GNN over the
    variant–gene–pathway–PPI graph, then derive pairwise interaction
    scores from cross-attention weights or integrated gradients.
    Requires PyTorch Geometric (already an optional dep via
    ``graphgwas[gnn]``).
    """
    raise NotImplementedError(
        "M6 GNN-attention extraction not yet implemented. "
        "Install: pip install graphgwas[gnn]; "
        "training infrastructure design in docs/EPISTASIS_V2_DESIGN.md §M6."
    )


def topological_epistasis(*args, **kwargs):
    """M7 Persistent-homology / topological epistasis.

    [TODO paper2-Y5-M7]: build the case-only co-occurrence simplicial
    complex; compute persistent 1-cycles via Gudhi or Giotto-TDA;
    rank multi-way interactions by birth-death persistence.
    Requires gudhi or giotto-tda (not currently a dep).
    """
    raise NotImplementedError(
        "M7 persistent homology not yet implemented. "
        "Install: pip install gudhi  OR  pip install giotto-tda; "
        "design in docs/EPISTASIS_V2_DESIGN.md §M7."
    )
