"""Weighted gene-gene substrate for the M5★ RWR kernel.

Replaces the binary `1[(g_a, g_b) ∈ canonical_PPI]` gene-gene transition
matrix used in §Y.7 with a data-driven weighted matrix that consumes
multiple evidence sources (Billmann HAP1 qGI scores, DepMap
co-essentiality, STRING / canonical PPI as fallback). Each edge can be
backed by multiple sources; we max-pool across sources to avoid
double-counting.

Edge file format (TSV, one row per (gene_a, gene_b, source) triple,
gene order is canonical alphabetical so each undirected edge appears once):

    gene_a    gene_b    weight    source

Where:
  - gene_a, gene_b: gene IDs (HGNC symbols for human; locus IDs for
    plants; systematic names for yeast).
  - weight: float in [0, 1]. Already-normalised by the loader; raw qGI
    or correlation values are mapped via clip(|x|/scale, 0, 1).
  - source: short string label (e.g. "billmann_qgi", "depmap_ed",
    "string_ppi", "canonical_ppi").

The module is backwards-compatible: callers that don't pass weighted
edges get the binary cache.ppi behavior unchanged. The §Y.7 M5★ kernel
becomes a special case where weight=1 for every canonical-PPI edge.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_weighted_edges(path: Path | str) -> dict[tuple[str, str], dict[str, float]]:
    """Load a weighted edge TSV.

    Returns: {(g_a, g_b): {source: weight, ...}}, with (g_a, g_b)
    canonicalised so g_a < g_b (alphabetical). Each undirected edge
    can have multiple sources; the consumer can max-pool, sum, or
    weight-by-source as needed.
    """
    out: dict[tuple[str, str], dict[str, float]] = defaultdict(dict)
    path = Path(path)
    if not path.exists():
        return {}
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            g_a, g_b = row["gene_a"], row["gene_b"]
            if g_a == g_b:
                continue
            if g_a > g_b:
                g_a, g_b = g_b, g_a
            try:
                w = float(row["weight"])
            except (KeyError, ValueError):
                continue
            if not np.isfinite(w):
                continue
            src = row.get("source", "unknown")
            # Keep max if multiple rows for same (pair, source)
            if src in out[(g_a, g_b)]:
                out[(g_a, g_b)][src] = max(out[(g_a, g_b)][src], w)
            else:
                out[(g_a, g_b)][src] = w
    return dict(out)


def edges_from_billmann_qgi_tsv(path: Path | str,
                                  use_stringent_only: bool = False,
                                  qgi_scale: float = 1.0) -> list[dict]:
    """Convert billmann_hc_gi_edges.tsv → weighted-edge rows.

    qGI scores typically in [-3, 3]; we map weight = clip(|qGI|/qgi_scale, 0, 1).
    Default qgi_scale=1.0 means a |qGI| of 1.0 → weight 1.0 (saturated);
    smaller scale = sharper saturation.
    """
    rows: list[dict] = []
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if use_stringent_only and not row.get("GI_stringent"):
                continue
            if not row.get("GI_standard"):
                continue
            try:
                qgi = float(row["qGI_score"])
            except (KeyError, ValueError):
                continue
            w = float(np.clip(abs(qgi) / qgi_scale, 0.0, 1.0))
            g_a, g_b = row["library_gene"], row["query_gene"]
            if g_a == g_b:
                continue
            if g_a > g_b:
                g_a, g_b = g_b, g_a
            rows.append({"gene_a": g_a, "gene_b": g_b,
                          "weight": f"{w:.4f}",
                          "source": "billmann_qgi"})
    return rows


def edges_from_depmap_ed_tsv(path: Path | str,
                              ed_scale: float = 0.5) -> list[dict]:
    """Convert billmann_depmap_ed_edges.tsv → weighted-edge rows.

    File_S21 reports two ED scores (lib-vs-query, query-vs-lib). We take
    the max(|ed_lib_v_query|, |ed_query_v_lib|) and map
    weight = clip(max_ed / ed_scale, 0, 1). Default ed_scale=0.5 because
    DepMap correlations rarely exceed 0.5 in absolute value.
    """
    rows: list[dict] = []
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            try:
                ed1 = float(row.get("ED_score_lib_v_query") or 0.0)
                ed2 = float(row.get("ED_score_query_v_lib") or 0.0)
            except ValueError:
                continue
            ed_max = max(abs(ed1), abs(ed2))
            w = float(np.clip(ed_max / ed_scale, 0.0, 1.0))
            g_a, g_b = row["library_gene"], row["query_gene"]
            if g_a == g_b:
                continue
            if g_a > g_b:
                g_a, g_b = g_b, g_a
            rows.append({"gene_a": g_a, "gene_b": g_b,
                          "weight": f"{w:.4f}",
                          "source": "depmap_ed"})
    return rows


def edges_from_canonical_ppi(canonical_edges: Iterable[tuple[str, str]],
                              weight: float = 1.0) -> list[dict]:
    """Convert a canonical-PPI tuple list → weighted-edge rows."""
    out: list[dict] = []
    for g_a, g_b in canonical_edges:
        if g_a == g_b:
            continue
        if g_a > g_b:
            g_a, g_b = g_b, g_a
        out.append({"gene_a": g_a, "gene_b": g_b,
                     "weight": f"{weight:.4f}",
                     "source": "canonical_ppi"})
    return out


def write_edges_tsv(rows: Iterable[dict], path: Path | str) -> None:
    """Persist a weighted-edge row list to TSV."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        w = csv.DictWriter(f, fieldnames=["gene_a", "gene_b", "weight", "source"],
                            delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


# ---------------------------------------------------------------------------
# Substrate construction
# ---------------------------------------------------------------------------

def build_weighted_gene_gene_matrix(
    gene_ids: list[str],
    weighted_edges: dict[tuple[str, str], dict[str, float]],
    *,
    source_weights: dict[str, float] | None = None,
    pool: str = "max",
    include_self_loops: bool = True,
) -> np.ndarray:
    """Construct a (n_gene, n_gene) gene-gene adjacency from weighted edges.

    For each undirected (g_a, g_b) edge in `weighted_edges`, the entry
    G[a, b] = G[b, a] is set as follows:

      pool="max":   G[a,b] = max_s(source_weights[s] * w_s)
      pool="sum":   G[a,b] = sum_s(source_weights[s] * w_s)
      pool="prod":  G[a,b] = 1 - prod_s(1 - source_weights[s] * w_s)  (noisy-OR)

    Self-loops can be added (default), giving each gene a chance to stay
    where it is during the random walk.

    Args:
      gene_ids: ordered list of gene IDs that defines the matrix index
        space. Edges referencing genes outside this list are dropped.
      weighted_edges: output of `load_weighted_edges`.
      source_weights: optional per-source multiplier (e.g.
        {"billmann_qgi": 1.0, "depmap_ed": 0.6, "canonical_ppi": 1.0}).
        Defaults to all 1.0.
      pool: aggregation across sources for the same edge.
      include_self_loops: add 1.0 on the diagonal.

    Returns:
      G: shape (n_gene, n_gene), float64. Not yet row-normalised — the
      caller is responsible for that (matches the existing §Y.7 kernel).
    """
    n = len(gene_ids)
    gene_to_idx = {g: i for i, g in enumerate(gene_ids)}
    G = np.zeros((n, n), dtype=np.float64)
    if include_self_loops:
        np.fill_diagonal(G, 1.0)
    if source_weights is None:
        source_weights = {}

    n_edges_applied = 0
    n_edges_dropped = 0
    for (g_a, g_b), source_to_w in weighted_edges.items():
        if g_a not in gene_to_idx or g_b not in gene_to_idx:
            n_edges_dropped += 1
            continue
        ia, ib = gene_to_idx[g_a], gene_to_idx[g_b]
        # Aggregate across sources
        contrib = []
        for src, w in source_to_w.items():
            mult = source_weights.get(src, 1.0)
            contrib.append(mult * w)
        if not contrib:
            continue
        if pool == "max":
            edge_w = max(contrib)
        elif pool == "sum":
            edge_w = sum(contrib)
        elif pool == "prod":
            # Noisy-OR
            p = 1.0
            for c in contrib:
                p *= (1.0 - min(c, 1.0))
            edge_w = 1.0 - p
        else:
            raise ValueError(f"unknown pool: {pool}")
        G[ia, ib] = edge_w
        G[ib, ia] = edge_w
        n_edges_applied += 1

    return G


def row_normalise(G: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """Convert a non-negative adjacency to a row-stochastic transition matrix."""
    row_sum = G.sum(axis=1, keepdims=True).clip(min=eps)
    return G / row_sum


def build_weighted_p_gg(
    gene_ids: list[str],
    weighted_edges: dict[tuple[str, str], dict[str, float]],
    *,
    source_weights: dict[str, float] | None = None,
    pool: str = "max",
) -> np.ndarray:
    """Convenience wrapper: weighted edges → row-normalised P_gg.

    Drop-in replacement for the binary `P_gg` constructed inline in
    `tests/m5_variants_benchmark.compute_rwr_matrix(..., use_ppi_bridge=True)`.
    """
    G = build_weighted_gene_gene_matrix(
        gene_ids, weighted_edges,
        source_weights=source_weights, pool=pool,
        include_self_loops=True,
    )
    return row_normalise(G)


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

def edge_diagnostics(
    weighted_edges: dict[tuple[str, str], dict[str, float]],
    gene_ids: list[str] | None = None,
) -> dict:
    """Quick summary of an edge dict (counts, weight distribution, source mix)."""
    n_edges = len(weighted_edges)
    if n_edges == 0:
        return {"n_edges": 0}
    src_counts: dict[str, int] = defaultdict(int)
    weights: list[float] = []
    for source_to_w in weighted_edges.values():
        for src, w in source_to_w.items():
            src_counts[src] += 1
            weights.append(w)
    arr = np.asarray(weights)
    out = {
        "n_edges": n_edges,
        "n_source_observations": len(weights),
        "source_counts": dict(src_counts),
        "weight_min": float(arr.min()),
        "weight_max": float(arr.max()),
        "weight_median": float(np.median(arr)),
        "weight_mean": float(arr.mean()),
    }
    if gene_ids is not None:
        gset = set(gene_ids)
        in_gset = sum(1 for (a, b) in weighted_edges
                       if a in gset and b in gset)
        out["n_edges_with_both_genes_in_panel"] = in_gset
    return out
