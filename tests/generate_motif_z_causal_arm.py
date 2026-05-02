"""Causal-SNP arm of paper #2 §Y.1 — does motif-typed Z suppress ρ_max for
biology-aligned targets relative to random Z?

The null arm (``generate_motif_z_comparison.py``) showed that for null target
SNPs (sampled to be disjoint from Z's column-set), motif-Z and random-Z give
near-identical ρ_max distributions on 1KG chr22 (mean 0.45 vs 0.46).  This
is the "no penalty on null" claim.

This script tests the second-half claim: when target SNPs are themselves
biology-aligned (i.e., they sit in at least one same-gene motif pair, so
the *true* interaction direction lies in the motif span by construction),
the random-Z ρ_max is more inflated than the motif-Z ρ_max.  The economic
reason: random Z spans an arbitrary 100-dimensional subspace that happens
to overlap with the target SNP through generic Gram-matrix noise, while
motif Z spans the small biology-aligned subspace that contains the actual
interaction direction without the off-axis noise.

Output:
  paper/epistasis_v1/figures/fig3_motif_z_causal_arm.{png,pdf,json}

Usage:
  python tests/generate_motif_z_causal_arm.py  [--quick]
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from scipy import stats as sp_stats  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src" / "python"))

from graphgwas.bgen_reader import BgenReader  # noqa: E402
from graphgwas.bias import (  # noqa: E402
    conservativeness_ratio,
    mu_shift,
    rho_max,
)


BGEN_DIR = REPO_ROOT / "tests" / "data" / "human" / "1kGP_bgen"
CACHE_PATH = REPO_ROOT / "data" / "annotations" / "human_graph_cache_v2_chr22.json"
OUT_DIR = REPO_ROOT / "paper" / "epistasis_v1" / "figures"

CHR22_WINDOW = (30_000_000, 35_000_000)
MAF_MIN = 0.05


# ===========================================================================
# Data loading (same as fig2 script)
# ===========================================================================

def load_chr22(reader: BgenReader) -> tuple[np.ndarray, list[str]]:
    print(f"  loading chr22:{CHR22_WINDOW[0]}-{CHR22_WINDOW[1]} ...", flush=True)
    df, dos = reader.load_locus("22", *CHR22_WINDOW, format="dosage")
    af = np.nanmean(dos, axis=0) / 2.0
    maf = np.minimum(af, 1.0 - af)
    keep = maf >= MAF_MIN
    dos = dos[:, keep].astype(np.float64)
    df = df[keep].reset_index(drop=True)
    vids = [f"{r.chr}:{int(r.pos)}:{r.a2}:{r.a1}"
            for r in df.itertuples()]
    print(f"    n_samples={dos.shape[0]}, n_variants={dos.shape[1]}", flush=True)
    return dos, vids


def load_motif_pairs(vids: list[str]) -> tuple[
    list[tuple[int, int]],
    set[int],   # variants with at least one same-gene partner (CAUSAL pool)
    set[int],   # variants with no same-gene partner (NULL pool)
    dict,
]:
    """Build motif pairs and the causal/null target pools simultaneously.

    Returns:
        pairs:        list of (i,j) variant index tuples (motif P1 same-gene)
        causal_pool:  variant indices that appear in at least one motif pair
        null_pool:    variant indices that appear in NO motif pair
        stats:        bookkeeping dict
    """
    print(f"  loading graph cache: {CACHE_PATH.name}")
    cache = json.loads(CACHE_PATH.read_text())
    print(f"    {len(cache):,} variants in cache")

    vid_to_idx = {v: i for i, v in enumerate(vids)}
    gene_to_vars: dict[str, list[int]] = defaultdict(list)
    for vid, idx in vid_to_idx.items():
        entry = cache.get(vid)
        if entry is None:
            continue
        for gene in entry.get("genes", []):
            gene_to_vars[gene].append(idx)

    MAX_PAIRS_PER_GENE = 500
    rng = np.random.default_rng(0)
    pairs: list[tuple[int, int]] = []
    n_gene_capped = 0
    for gene, var_idx in gene_to_vars.items():
        if len(var_idx) < 2:
            continue
        n_full = len(var_idx) * (len(var_idx) - 1) // 2
        if n_full <= MAX_PAIRS_PER_GENE:
            for i, j in itertools.combinations(var_idx, 2):
                pairs.append((i, j))
        else:
            n_gene_capped += 1
            for _ in range(MAX_PAIRS_PER_GENE):
                a, b = rng.choice(var_idx, size=2, replace=False)
                pairs.append((int(a), int(b)))

    pairs_set = set()
    deduped: list[tuple[int, int]] = []
    for i, j in pairs:
        key = (min(i, j), max(i, j))
        if key in pairs_set:
            continue
        pairs_set.add(key)
        deduped.append(key)

    # Causal pool = union of all variants that participate in any motif pair
    causal_pool: set[int] = set()
    for i, j in deduped:
        causal_pool.add(i)
        causal_pool.add(j)
    # Null pool = all other variants (not in any motif pair)
    all_idx = set(range(len(vids)))
    null_pool = all_idx - causal_pool

    stats = {
        "n_variants_total": len(vids),
        "n_variants_in_cache": sum(1 for v in vids if v in cache),
        "n_genes_with_2plus_variants": sum(
            1 for vs in gene_to_vars.values() if len(vs) >= 2
        ),
        "n_genes_capped": n_gene_capped,
        "max_pairs_per_gene": MAX_PAIRS_PER_GENE,
        "n_motif_pairs_unique": len(deduped),
        "n_causal_pool": len(causal_pool),
        "n_null_pool": len(null_pool),
    }
    print(f"    {stats['n_motif_pairs_unique']:,} unique same-gene pairs")
    print(f"    causal pool (≥1 same-gene partner): "
          f"{stats['n_causal_pool']:,} variants")
    print(f"    null pool   (no same-gene partner): "
          f"{stats['n_null_pool']:,} variants")
    return deduped, causal_pool, null_pool, stats


# ===========================================================================
# 2x2 grid: target_type × Z_type
# ===========================================================================

def sweep_one_cell(target_dos: np.ndarray,
                   target_indices_pool: np.ndarray,
                   z_builder,                 # callable: rng → Z
                   n_z: int,
                   n_target: int,
                   rng: np.random.Generator,
                   name: str) -> np.ndarray:
    """Run one cell: build n_z Z matrices, draw n_target targets each, compute ρ_max.

    target_indices_pool: variant column indices to sample targets from.  We
        ALWAYS exclude variants that appear in the just-built Z columns
        (strict no-path null) so the "biology-aligned" claim is geometric,
        not artifactual from g being literally in col(Z).
    """
    rho = np.empty(n_z * n_target, dtype=np.float64)
    k = 0
    for zi in range(n_z):
        Z, used_var_idx = z_builder(rng)
        Z = Z - Z.mean(axis=0, keepdims=True)
        # Exclude any target SNP that appears as a factor in Z
        candidate = np.setdiff1d(target_indices_pool,
                                 np.fromiter(used_var_idx, dtype=int))
        if len(candidate) < n_target:
            raise RuntimeError(
                f"[{name}] target pool too small after Z exclusion: "
                f"{len(candidate)} < {n_target}"
            )
        sel = rng.choice(candidate, size=n_target, replace=False)
        for ti in sel:
            rho[k] = rho_max(target_dos[:, ti], Z, preprocess_inputs=True)
            k += 1
        if (zi + 1) % max(1, n_z // 5) == 0:
            print(f"    [{name}] {zi+1}/{n_z} done", flush=True)
    return rho


def make_random_z_builder(dos: np.ndarray, m_int: int):
    """Returns callable(rng) → (Z, used_var_idx_set)."""
    n_var = dos.shape[1]
    def builder(rng):
        a = rng.integers(0, n_var, size=m_int)
        b = rng.integers(0, n_var, size=m_int)
        Z = dos[:, a] * dos[:, b]
        used = set(a.tolist()) | set(b.tolist())
        return Z, used
    return builder


def make_motif_z_builder(dos: np.ndarray,
                          motif_pairs: list[tuple[int, int]],
                          m_int: int):
    """Returns callable(rng) → (Z, used_var_idx_set)."""
    pair_arr = np.asarray(motif_pairs)  # (P, 2)
    n_pairs = len(motif_pairs)
    def builder(rng):
        sel = rng.choice(n_pairs, size=min(m_int, n_pairs), replace=False)
        a = pair_arr[sel, 0]
        b = pair_arr[sel, 1]
        Z = dos[:, a] * dos[:, b]
        used = set(a.tolist()) | set(b.tolist())
        return Z, used
    return builder


# ===========================================================================
# Plot
# ===========================================================================

def _summarise(arr: np.ndarray) -> dict:
    return {
        "n": int(arr.size),
        "min": float(arr.min()),
        "mean": float(arr.mean()),
        "median": float(np.median(arr)),
        "p95": float(np.percentile(arr, 95)),
        "max": float(arr.max()),
    }


def plot(cells: dict[str, np.ndarray],
         out_dir: Path,
         params: dict,
         motif_stats: dict) -> None:
    """2×2 grid: rows = target type (null, causal); cols = Z type (random, motif)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.4),
                             sharex=True, sharey=True)
    bins = np.linspace(0, 1.0, 80)

    layout = [
        ("Null target × Random Z",   cells["null_random"],   "#7f7f7f", axes[0, 0]),
        ("Null target × Motif Z",    cells["null_motif"],    "#1f77b4", axes[0, 1]),
        ("Causal target × Random Z", cells["causal_random"], "#7f7f7f", axes[1, 0]),
        ("Causal target × Motif Z",  cells["causal_motif"],  "#1f77b4", axes[1, 1]),
    ]
    for title, rho, color, ax in layout:
        ax.hist(rho, bins=bins, color=color, edgecolor="white", linewidth=0.3)
        ax.axvline(rho.mean(), color="black", linestyle="--", linewidth=0.8)
        ax.set_yscale("log")
        ax.set_title(
            f"{title}\nmean={rho.mean():.3f}, p95={np.percentile(rho, 95):.3f}, "
            f"max={rho.max():.3f}",
            fontsize=10,
        )
    for ax in axes[1]:
        ax.set_xlabel(r"$\rho_{\max}$")
    for ax in axes[:, 0]:
        ax.set_ylabel("Counts (log scale)")

    fig.suptitle(
        f"§Y.1 Causal arm: motif-typed Z vs random Z, biology-aligned vs "
        f"biology-isolated targets\n"
        f"1KG chr22, n={params['n_samples']}, "
        f"motif pool: {motif_stats['n_motif_pairs_unique']:,} same-gene pairs",
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_dir / "fig3_motif_z_causal_arm.png", dpi=200,
                bbox_inches="tight")
    fig.savefig(out_dir / "fig3_motif_z_causal_arm.pdf",
                bbox_inches="tight")
    plt.close(fig)

    # KS tests + bias-derived headline numbers
    n_samp = params["n_samples"]
    lambda_var = 0.1   # mid-range from Yelmen et al. 2026 (their Fig 1 reports λ ~ 0.1)
    x_thresh = 5.45    # genome-wide significance threshold

    def headline_from_rho(arr: np.ndarray) -> dict:
        """At the mean rho, what is mu shift and R(5.45) under λ=0.1?"""
        rho_mean = float(arr.mean())
        return {
            "rho_mean": rho_mean,
            "mu_shift_at_rho_mean": mu_shift(rho_mean, lambda_var, n_samp),
            "R_at_x5.45_lambda0.1_rho_mean": conservativeness_ratio(
                x_thresh, n_samp, lambda_var, rho_mean
            ),
        }

    null_ks = sp_stats.ks_2samp(cells["null_random"], cells["null_motif"])
    causal_ks = sp_stats.ks_2samp(cells["causal_random"], cells["causal_motif"])

    sidecar = {
        "params": params,
        "motif_stats": motif_stats,
        "lambda_var_assumed": lambda_var,
        "x_threshold": x_thresh,
        "cells": {
            "null_random":   _summarise(cells["null_random"]),
            "null_motif":    _summarise(cells["null_motif"]),
            "causal_random": _summarise(cells["causal_random"]),
            "causal_motif":  _summarise(cells["causal_motif"]),
        },
        "headline_at_mean_rho": {
            "null_random":   headline_from_rho(cells["null_random"]),
            "null_motif":    headline_from_rho(cells["null_motif"]),
            "causal_random": headline_from_rho(cells["causal_random"]),
            "causal_motif":  headline_from_rho(cells["causal_motif"]),
        },
        "ks_tests": {
            "null_random_vs_null_motif":     {
                "statistic": float(null_ks.statistic),
                "p_value":   float(null_ks.pvalue),
            },
            "causal_random_vs_causal_motif": {
                "statistic": float(causal_ks.statistic),
                "p_value":   float(causal_ks.pvalue),
            },
        },
    }
    (out_dir / "fig3_motif_z_causal_arm.json").write_text(
        json.dumps(sidecar, indent=2)
    )
    print(f"\nWrote {out_dir / 'fig3_motif_z_causal_arm'}.{{png,pdf,json}}")


# ===========================================================================
# Main
# ===========================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()

    if args.quick:
        n_z, m_int, n_target = 10, 50, 30
    else:
        n_z, m_int, n_target = 30, 100, 100

    print("=== §Y.1 Causal-SNP arm: motif-Z vs random-Z, causal vs null targets ===")
    print(f"BGEN dir   : {BGEN_DIR}")
    print(f"Cache file : {CACHE_PATH}")
    print(f"Output dir : {OUT_DIR}")
    print(f"Sweep      : n_z={n_z}, m_int={m_int}, n_target_per_z={n_target}")
    print()

    if not BGEN_DIR.exists():
        sys.exit(f"BGEN dir not found: {BGEN_DIR}")
    if not CACHE_PATH.exists():
        sys.exit(f"Graph cache not found: {CACHE_PATH}")

    rng = np.random.default_rng(args.seed)
    with BgenReader(BGEN_DIR) as reader:
        dos, vids = load_chr22(reader)
    motif_pairs, causal_pool, null_pool, motif_stats = load_motif_pairs(vids)
    if len(motif_pairs) < n_z * m_int:
        sys.exit(f"Too few motif pairs ({len(motif_pairs)}) for "
                 f"n_z*m_int = {n_z * m_int}")
    if len(null_pool) < n_target * 2:
        sys.exit(f"Null pool too small ({len(null_pool)}) for "
                 f"n_target = {n_target}")

    causal_arr = np.fromiter(causal_pool, dtype=int)
    null_arr = np.fromiter(null_pool, dtype=int)

    rand_z = make_random_z_builder(dos, m_int)
    motif_z = make_motif_z_builder(dos, motif_pairs, m_int)

    print("\n[1/4] null target × random Z")
    rho_null_rand = sweep_one_cell(dos, null_arr, rand_z, n_z, n_target,
                                    rng, "null·random")
    print("\n[2/4] null target × motif Z")
    rho_null_motif = sweep_one_cell(dos, null_arr, motif_z, n_z, n_target,
                                     rng, "null·motif")
    print("\n[3/4] causal target × random Z")
    rho_causal_rand = sweep_one_cell(dos, causal_arr, rand_z, n_z, n_target,
                                      rng, "causal·random")
    print("\n[4/4] causal target × motif Z")
    rho_causal_motif = sweep_one_cell(dos, causal_arr, motif_z, n_z, n_target,
                                       rng, "causal·motif")
    print()

    cells = {
        "null_random":   rho_null_rand,
        "null_motif":    rho_null_motif,
        "causal_random": rho_causal_rand,
        "causal_motif":  rho_causal_motif,
    }
    plot(cells, OUT_DIR, params={
        "n_samples": int(dos.shape[0]),
        "n_z": int(n_z),
        "m_int": int(m_int),
        "n_target": int(n_target),
        "maf_min": float(MAF_MIN),
        "chr22_window": list(CHR22_WINDOW),
        "seed": int(args.seed),
    }, motif_stats=motif_stats)

    print("\nSummary (ρ_max means):")
    print(f"  null   target × random Z : {rho_null_rand.mean():.4f}")
    print(f"  null   target × motif  Z : {rho_null_motif.mean():.4f}")
    print(f"  causal target × random Z : {rho_causal_rand.mean():.4f}")
    print(f"  causal target × motif  Z : {rho_causal_motif.mean():.4f}")
    delta_null = rho_null_rand.mean() - rho_null_motif.mean()
    delta_causal = rho_causal_rand.mean() - rho_causal_motif.mean()
    print(f"  Δ(random − motif) on null   : {delta_null:+.4f}")
    print(f"  Δ(random − motif) on causal : {delta_causal:+.4f}")
    print(f"  → motif-Z bias suppression on causal targets: "
          f"{delta_causal - delta_null:+.4f}")


if __name__ == "__main__":
    main()
