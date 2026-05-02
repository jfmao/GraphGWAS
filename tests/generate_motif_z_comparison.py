"""Compare random-Z vs motif-typed-Z rho_max distributions on 1KG chr22.

The novel paper-#2 contribution beyond Yelmen et al. 2026.  Their analysis
uses a *random* interaction-feature matrix Z (each column a product of two
randomly chosen SNPs).  We instead build Z from **biology-typed motif pairs**:
pairs of SNPs that share at least one annotated gene (motif P1, same-gene
cis interaction).

Hypothesis: the motif-typed Z spans a smaller, biology-aligned subspace.
For null target SNPs (truly uncorrelated with the trait), the motif-Z and
random-Z rho_max distributions should be similar (the motif typing does
not penalise null SNPs).  For *causal* target SNPs (correlated with the
true interaction signal), motif-Z should suppress rho_max relative to
random-Z (because the motif subspace is aligned with biology, not with
the arbitrary interaction direction the random Z spans).

This script tests the first half (null SNPs ⇒ similar distributions); the
causal-SNP arm requires phenotype simulation under a non-null DGP and is
deferred to the REGENIE rescue demo.

Output:
  paper/epistasis_v1/figures/fig2_motif_z_comparison.{png,pdf,json}

Usage:
  python tests/generate_motif_z_comparison.py  [--quick]
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

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src" / "python"))

from graphgwas.bgen_reader import BgenReader  # noqa: E402
from graphgwas.bias import rho_max  # noqa: E402


BGEN_DIR = REPO_ROOT / "tests" / "data" / "human" / "1kGP_bgen"
CACHE_PATH = REPO_ROOT / "data" / "annotations" / "human_graph_cache_v2_chr22.json"
OUT_DIR = REPO_ROOT / "paper" / "epistasis_v1" / "figures"

CHR22_WINDOW = (30_000_000, 35_000_000)
MAF_MIN = 0.05


# ===========================================================================
# Data loading
# ===========================================================================

def load_chr22(reader: BgenReader) -> tuple[np.ndarray, list[str]]:
    """Load chr22 dosages + variant IDs in 'chrN:pos:a1:a2' format."""
    print(f"  loading chr22:{CHR22_WINDOW[0]}-{CHR22_WINDOW[1]} ...", flush=True)
    df, dos = reader.load_locus("22", *CHR22_WINDOW, format="dosage")
    af = np.nanmean(dos, axis=0) / 2.0
    maf = np.minimum(af, 1.0 - af)
    keep = maf >= MAF_MIN
    dos = dos[:, keep].astype(np.float64)
    df = df[keep].reset_index(drop=True)
    # Build variant IDs in the 'chr22:pos:REF:ALT' format used by graph_cache.
    # BGEN reader's `a1, a2` are (ALT, REF) per plink2's VCF→BGEN convention,
    # so we swap to match the cache's VCF-standard REF:ALT key format.
    vids = [f"{r.chr}:{int(r.pos)}:{r.a2}:{r.a1}"
            for r in df.itertuples()]
    print(f"    n_samples={dos.shape[0]}, n_variants={dos.shape[1]} "
          f"(after MAF≥{MAF_MIN})", flush=True)
    return dos, vids


def load_motif_pairs(vids: list[str]) -> tuple[list[tuple[int, int]], dict]:
    """Build P1 (same-gene) motif pairs from human_graph_cache_v2_chr22.json.

    Returns:
        (pairs, stats) where pairs is a list of (i, j) index tuples into vids
        and stats is a dict of bookkeeping numbers.
    """
    print(f"  loading graph cache: {CACHE_PATH.name}")
    cache = json.loads(CACHE_PATH.read_text())
    print(f"    {len(cache)} variants in cache")

    # Build vid -> column index for fast lookup
    vid_to_idx = {v: i for i, v in enumerate(vids)}

    # Build gene -> [variant indices] for variants present in our dosage matrix
    gene_to_vars: dict[str, list[int]] = defaultdict(list)
    n_in_cache = 0
    for vid, idx in vid_to_idx.items():
        entry = cache.get(vid)
        if entry is None:
            continue
        n_in_cache += 1
        for gene in entry.get("genes", []):
            gene_to_vars[gene].append(idx)

    # Generate same-gene pairs (motif P1)
    # Cap per-gene pairs to avoid O(n^2) explosion on large genes; we sample
    # up to MAX_PAIRS_PER_GENE distinct pairs from each gene's variant list.
    MAX_PAIRS_PER_GENE = 500
    rng = np.random.default_rng(0)  # deterministic motif sampling
    pairs: list[tuple[int, int]] = []
    n_pairwise_blowup = 0
    for gene, var_idx in gene_to_vars.items():
        if len(var_idx) < 2:
            continue
        n_full = len(var_idx) * (len(var_idx) - 1) // 2
        if n_full <= MAX_PAIRS_PER_GENE:
            for i, j in itertools.combinations(var_idx, 2):
                pairs.append((i, j))
        else:
            # Random pair sampling
            n_pairwise_blowup += 1
            for _ in range(MAX_PAIRS_PER_GENE):
                a, b = rng.choice(var_idx, size=2, replace=False)
                pairs.append((int(a), int(b)))

    # Deduplicate
    pairs_set = set()
    deduped = []
    for i, j in pairs:
        key = (min(i, j), max(i, j))
        if key in pairs_set:
            continue
        pairs_set.add(key)
        deduped.append(key)

    stats = {
        "n_variants_total": len(vids),
        "n_variants_in_cache": n_in_cache,
        "n_genes_with_2plus_variants": sum(
            1 for vs in gene_to_vars.values() if len(vs) >= 2
        ),
        "n_genes_capped": n_pairwise_blowup,
        "max_pairs_per_gene": MAX_PAIRS_PER_GENE,
        "n_motif_pairs_unique": len(deduped),
    }
    print(f"    {n_in_cache}/{len(vids)} variants have cache entry")
    print(f"    {stats['n_genes_with_2plus_variants']} genes have ≥2 variants "
          f"({stats['n_genes_capped']} genes hit per-gene pair cap)")
    print(f"    {len(deduped):,} unique same-gene pairs (motif P1)")
    return deduped, stats


# ===========================================================================
# rho_max sweeps
# ===========================================================================

def sweep_random_z(target_dos: np.ndarray,
                   interaction_dos: np.ndarray,
                   n_z: int,
                   m_int: int,
                   n_target: int,
                   rng: np.random.Generator,
                   name: str) -> np.ndarray:
    """Random-Z sweep: Z columns are products of two randomly-chosen SNPs."""
    n_t, n_i = target_dos.shape[1], interaction_dos.shape[1]
    rho = np.empty(n_z * n_target, dtype=np.float64)
    k = 0
    for zi in range(n_z):
        a_idx = rng.integers(0, n_i, size=m_int)
        b_idx = rng.integers(0, n_i, size=m_int)
        Z = interaction_dos[:, a_idx] * interaction_dos[:, b_idx]
        Z = Z - Z.mean(axis=0, keepdims=True)
        for ti in rng.choice(n_t, size=n_target, replace=False):
            rho[k] = rho_max(target_dos[:, ti], Z, preprocess_inputs=True)
            k += 1
        if (zi + 1) % max(1, n_z // 5) == 0:
            print(f"    [{name}] {zi+1}/{n_z} done", flush=True)
    return rho


def sweep_motif_z(dosages: np.ndarray,
                  motif_pairs: list[tuple[int, int]],
                  n_z: int,
                  m_int: int,
                  n_target: int,
                  rng: np.random.Generator,
                  name: str) -> np.ndarray:
    """Motif-Z sweep: Z columns are products of *same-gene* pairs.

    Crucially, when sampling target SNPs we exclude any variant that appears
    in any of the m_int sampled motif pairs (strict no-path null).
    """
    n_pairs = len(motif_pairs)
    n_var = dosages.shape[1]
    pair_arr = np.asarray(motif_pairs)  # (n_pairs, 2)
    rho = np.empty(n_z * n_target, dtype=np.float64)
    k = 0
    for zi in range(n_z):
        # Sample m_int motif pairs without replacement
        pair_idx = rng.choice(n_pairs, size=min(m_int, n_pairs), replace=False)
        used_vars = set(pair_arr[pair_idx].ravel().tolist())
        # Build Z from sampled pairs
        a_cols = pair_arr[pair_idx, 0]
        b_cols = pair_arr[pair_idx, 1]
        Z = dosages[:, a_cols] * dosages[:, b_cols]
        Z = Z - Z.mean(axis=0, keepdims=True)
        # Sample target SNPs from variants NOT in used_vars (strict no-path null)
        candidate_targets = np.setdiff1d(np.arange(n_var),
                                         np.fromiter(used_vars, dtype=int))
        if len(candidate_targets) < n_target:
            raise RuntimeError("Too many variants used in Z; reduce m_int")
        for ti in rng.choice(candidate_targets, size=n_target, replace=False):
            rho[k] = rho_max(dosages[:, ti], Z, preprocess_inputs=True)
            k += 1
        if (zi + 1) % max(1, n_z // 5) == 0:
            print(f"    [{name}] {zi+1}/{n_z} done", flush=True)
    return rho


# ===========================================================================
# Plot
# ===========================================================================

def plot(rho_random: np.ndarray,
         rho_motif: np.ndarray,
         out_dir: Path,
         params: dict,
         motif_stats: dict) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, (ax_r, ax_m) = plt.subplots(1, 2, figsize=(11, 4.4),
                                     sharex=True, sharey=True)

    bins = np.linspace(0, 1.0, 80)
    for ax, rho, color, design in [(ax_r, rho_random, "#7f7f7f", "Random Z (Yelmen baseline)"),
                                   (ax_m, rho_motif, "#1f77b4", "Motif-typed Z (P1 same-gene, GraphGWAS)")]:
        ax.hist(rho, bins=bins, color=color, edgecolor="white", linewidth=0.3)
        ax.axvline(rho.mean(), color="black", linestyle="--", linewidth=0.8)
        ax.set_yscale("log")
        ax.set_xlabel(r"$\rho_{\max}$")
        ax.set_title(
            f"{design}\n"
            f"min={rho.min():.3f}, mean={rho.mean():.3f}, "
            f"max={rho.max():.3f}"
        )
    ax_r.set_ylabel("Counts (log scale)")

    # Overlay second on first as a thin step for visual contrast
    fig.suptitle(
        f"Motif-typed Z vs random Z on 1KG chr22 (n={params['n_samples']}, "
        f"null target SNPs only)\n"
        f"motif pool: {motif_stats['n_motif_pairs_unique']:,} same-gene pairs "
        f"from {motif_stats['n_genes_with_2plus_variants']:,} multi-variant genes",
        fontsize=10,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.93])

    fig.savefig(out_dir / "fig2_motif_z_comparison.png", dpi=200,
                bbox_inches="tight")
    fig.savefig(out_dir / "fig2_motif_z_comparison.pdf",
                bbox_inches="tight")
    plt.close(fig)

    (out_dir / "fig2_motif_z_comparison.json").write_text(json.dumps({
        "params": params,
        "motif_stats": motif_stats,
        "random_z": _summarise(rho_random),
        "motif_z": _summarise(rho_motif),
        "ks_test_note": (
            "If random_z and motif_z distributions are statistically similar "
            "for null target SNPs (which they should be — this is the "
            "no-penalty-on-null claim), the next step is the causal-SNP arm: "
            "simulate phenotypes under non-null DGP and verify motif-Z "
            "rho_max is suppressed for the causal target."
        ),
    }, indent=2))
    print(f"\nWrote {out_dir / 'fig2_motif_z_comparison'}.{{png,pdf,json}}")


def _summarise(arr: np.ndarray) -> dict:
    return {
        "n": int(arr.size),
        "min": float(arr.min()),
        "mean": float(arr.mean()),
        "median": float(np.median(arr)),
        "p95": float(np.percentile(arr, 95)),
        "max": float(arr.max()),
    }


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

    print("=== Motif-typed Z comparison ===")
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
    motif_pairs, motif_stats = load_motif_pairs(vids)
    if len(motif_pairs) < n_z * m_int:
        sys.exit(f"Too few motif pairs ({len(motif_pairs)}) for "
                 f"n_z={n_z}*m_int={m_int}={n_z*m_int}")

    print("\nSweep 1/2: random Z")
    rho_rand = sweep_random_z(dos, dos, n_z, m_int, n_target, rng, "random-Z")

    print("\nSweep 2/2: motif-typed Z (P1 same-gene)")
    rho_motif = sweep_motif_z(dos, motif_pairs, n_z, m_int, n_target, rng,
                              "motif-Z")
    print()

    plot(rho_rand, rho_motif, OUT_DIR, params={
        "n_samples": int(dos.shape[0]),
        "n_z": int(n_z),
        "m_int": int(m_int),
        "n_target": int(n_target),
        "maf_min": float(MAF_MIN),
        "chr22_window": list(CHR22_WINDOW),
        "seed": int(args.seed),
    }, motif_stats=motif_stats)

    print("\nSummary:")
    print(f"  random Z   rho_max  mean={rho_rand.mean():.4f}  "
          f"max={rho_rand.max():.4f}")
    print(f"  motif Z    rho_max  mean={rho_motif.mean():.4f}  "
          f"max={rho_motif.max():.4f}")


if __name__ == "__main__":
    main()
