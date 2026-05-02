"""Reproduce Yelmen et~al.\ 2026 Fig~1 on 1KG chr21+chr22.

Validates the ``graphgwas.bias.rho_max`` implementation against real
human genotypes.  Two designs:

  - **same-chromosome**:  target SNPs and interaction SNPs both on chr22
  - **different-chromosome**:  target SNPs on chr22, interaction SNPs on chr21

For each design we generate ``N_Z`` random interaction-feature matrices
(each: ``M_INT`` columns of two-way SNP×SNP products), and compute
``rho_max`` for ``N_TARGET`` random target SNPs against each Z.

The resulting histograms reproduce Yelmen Fig~1's qualitative claim:
the different-chromosome distribution is tight near zero while the
same-chromosome distribution has a heavy upper tail.

Output:
  paper/epistasis_v1/figures/fig1_rho_max_distribution.{png,pdf}

Usage:
  python tests/generate_yelmen_fig1.py  [--quick]

The ``--quick`` flag uses smaller ``N_Z``, ``M_INT``, ``N_TARGET`` for
fast iteration (~30 s instead of ~5 min).
"""
from __future__ import annotations

import argparse
import json
import sys
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
OUT_DIR = REPO_ROOT / "paper" / "epistasis_v1" / "figures"


# ===========================================================================
# Locus selection
# ===========================================================================

# 5-Mb windows in the body of each chromosome; large enough to give plenty of
# common variants without hammering the BGEN-load cost.
CHR21_WINDOW = (30_000_000, 35_000_000)
CHR22_WINDOW = (30_000_000, 35_000_000)

MAF_MIN = 0.05  # Yelmen et al. use the GSA bi-allelic SNP framework; we
                # filter 1KG dosages to roughly the same regime.


def load_chrom(reader: BgenReader, chrom: str, window: tuple[int, int]) -> np.ndarray:
    """Load + MAF-filter dosage matrix for one chromosome window.

    Returns:
        Dosage matrix, shape (n_samples, n_variants_after_filter).
    """
    print(f"  loading {chrom}:{window[0]}-{window[1]} ...", flush=True)
    _, dos = reader.load_locus(chrom, window[0], window[1], format="dosage")
    if dos.size == 0:
        raise RuntimeError(f"No variants in {chrom}:{window}")
    # MAF filter
    af = np.nanmean(dos, axis=0) / 2.0
    maf = np.minimum(af, 1.0 - af)
    keep = maf >= MAF_MIN
    n_before = dos.shape[1]
    dos = dos[:, keep].astype(np.float64)
    print(f"    n_samples={dos.shape[0]}, n_variants_before_maf={n_before}, "
          f"after_maf={dos.shape[1]}", flush=True)
    return dos


# ===========================================================================
# rho_max sweep
# ===========================================================================

def sweep_design(target_dos: np.ndarray,
                 interaction_dos: np.ndarray,
                 n_z: int,
                 m_int: int,
                 n_target: int,
                 rng: np.random.Generator,
                 design_name: str) -> np.ndarray:
    """Run the (n_z × n_target) sweep for one design.

    Args:
        target_dos: dosages from which target SNPs are drawn (n_samples, n_t).
        interaction_dos: dosages from which the columns of Z are built
            (n_samples, n_i).  Z is constructed by sampling random pairs and
            taking element-wise products.
        n_z: number of random Z matrices.
        m_int: number of interaction-feature columns per Z.
        n_target: number of target SNPs to test per Z matrix.
        rng: numpy Generator.
        design_name: for progress logging.
    """
    n_t = target_dos.shape[1]
    n_i = interaction_dos.shape[1]
    print(f"  [{design_name}] n_t={n_t} target SNPs, n_i={n_i} pool, "
          f"n_z={n_z}, m_int={m_int}, n_target_per_z={n_target}", flush=True)

    rho_values = np.empty(n_z * n_target, dtype=np.float64)
    k = 0
    for zi in range(n_z):
        # Build Z: m_int columns, each a centred product of two random
        # interaction SNPs (sampled with replacement; collisions are fine
        # — Yelmen et al. don't enforce uniqueness either).
        a_idx = rng.integers(0, n_i, size=m_int)
        b_idx = rng.integers(0, n_i, size=m_int)
        Z = interaction_dos[:, a_idx] * interaction_dos[:, b_idx]
        Z = Z - Z.mean(axis=0, keepdims=True)

        # Random target SNPs (preprocessed via centring inside rho_max())
        t_idx = rng.choice(n_t, size=n_target, replace=False)
        for ti in t_idx:
            g = target_dos[:, ti]
            rho_values[k] = rho_max(g, Z, preprocess_inputs=True)
            k += 1
        if (zi + 1) % max(1, n_z // 10) == 0:
            print(f"    [{design_name}] {zi+1}/{n_z} Z matrices done", flush=True)
    return rho_values


# ===========================================================================
# Plotting
# ===========================================================================

def plot(rho_diff: np.ndarray,
         rho_same: np.ndarray,
         out_dir: Path,
         params: dict) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, (ax_d, ax_s) = plt.subplots(1, 2, figsize=(11, 4.4))

    # Yelmen et al. plot |rho_max| on log-counts; we mirror their styling
    bins_d = np.linspace(0, max(0.05, float(np.percentile(rho_diff, 99.5)) * 1.05), 60)
    ax_d.hist(rho_diff, bins=bins_d, color="#1f77b4", edgecolor="white",
              linewidth=0.3)
    ax_d.axvline(rho_diff.mean(), color="black", linestyle="--", linewidth=0.8)
    ax_d.set_yscale("log")
    ax_d.set_xlabel(r"$\rho_{\max}$")
    ax_d.set_ylabel("Counts (log scale)")
    ax_d.set_title(
        f"Different chromosomes (target=chr22, Z=chr21)\n"
        f"min={rho_diff.min():.3f}, mean={rho_diff.mean():.3f}, "
        f"max={rho_diff.max():.3f}"
    )

    bins_s = np.linspace(0, 1.0, 80)
    ax_s.hist(rho_same, bins=bins_s, color="#1f77b4", edgecolor="white",
              linewidth=0.3)
    ax_s.axvline(rho_same.mean(), color="black", linestyle="--", linewidth=0.8)
    ax_s.set_yscale("log")
    ax_s.set_xlabel(r"$\rho_{\max}$")
    ax_s.set_ylabel("Counts (log scale)")
    ax_s.set_title(
        f"Same chromosome (target=chr22, Z=chr22)\n"
        f"min={rho_same.min():.3f}, mean={rho_same.mean():.3f}, "
        f"max={rho_same.max():.3f}"
    )

    fig.suptitle(
        f"Reproduction of Yelmen et al. 2026 Fig 1 on 1KG chr21+22 "
        f"(n={params['n_samples']}, "
        f"n_z={params['n_z']}, m_int={params['m_int']}, "
        f"n_target_per_z={params['n_target']})",
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    fig.savefig(out_dir / "fig1_rho_max_distribution.png", dpi=200,
                bbox_inches="tight")
    fig.savefig(out_dir / "fig1_rho_max_distribution.pdf",
                bbox_inches="tight")
    plt.close(fig)

    # Sidecar JSON for reproducibility audits
    (out_dir / "fig1_rho_max_distribution.json").write_text(json.dumps({
        "params": params,
        "different_chr": {
            "n_values": int(rho_diff.size),
            "min": float(rho_diff.min()),
            "mean": float(rho_diff.mean()),
            "median": float(np.median(rho_diff)),
            "p95": float(np.percentile(rho_diff, 95)),
            "max": float(rho_diff.max()),
        },
        "same_chr": {
            "n_values": int(rho_same.size),
            "min": float(rho_same.min()),
            "mean": float(rho_same.mean()),
            "median": float(np.median(rho_same)),
            "p95": float(np.percentile(rho_same, 95)),
            "max": float(rho_same.max()),
        },
        "yelmen_2026_reference": {
            "different_chr_mean": 0.032,
            "different_chr_max": 0.042,
            "same_chr_mean": 0.035,
            "same_chr_max": 0.849,
            "note": "Reference values from Yelmen et al. 2026 Fig 1 caption "
                    "(Estonian Biobank n=210,145, 8170 GSA SNPs on chr21+22).",
        },
    }, indent=2))
    print(f"\nWrote {out_dir / 'fig1_rho_max_distribution'}.{{png,pdf,json}}")


# ===========================================================================
# Main
# ===========================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--quick", action="store_true",
                        help="Smaller sweep for fast iteration")
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()

    if args.quick:
        n_z, m_int, n_target = 10, 50, 30
    else:
        n_z, m_int, n_target = 30, 100, 100

    print("=== Yelmen Fig 1 reproduction ===")
    print(f"BGEN dir : {BGEN_DIR}")
    print(f"Output   : {OUT_DIR}")
    print(f"Sweep    : n_z={n_z}, m_int={m_int}, n_target_per_z={n_target}")
    print()

    if not BGEN_DIR.exists():
        sys.exit(f"BGEN directory not found: {BGEN_DIR}")

    rng = np.random.default_rng(args.seed)
    with BgenReader(BGEN_DIR) as reader:
        chr22 = load_chrom(reader, "22", CHR22_WINDOW)
        chr21 = load_chrom(reader, "21", CHR21_WINDOW)

    n_samples = chr22.shape[0]
    assert chr21.shape[0] == n_samples, "Sample counts must match across chromosomes"
    print()

    # Different-chromosome design: target on chr22, Z built from chr21
    print("Sweep 1/2: different-chromosome design")
    rho_diff = sweep_design(chr22, chr21, n_z, m_int, n_target, rng,
                            "diff-chr")
    print()

    # Same-chromosome design: target on chr22, Z built from chr22
    # Important: when sampling Z columns, we should not let the target SNP
    # itself appear in Z (Yelmen's strict no-path null).  We approximate this
    # by sampling Z columns from a disjoint half of chr22 variants.
    print("Sweep 2/2: same-chromosome design (target/interaction SNPs disjoint)")
    n_chr22 = chr22.shape[1]
    half = n_chr22 // 2
    target_pool = chr22[:, :half]
    int_pool = chr22[:, half:]
    rho_same = sweep_design(target_pool, int_pool, n_z, m_int, n_target, rng,
                            "same-chr")
    print()

    plot(rho_diff, rho_same, OUT_DIR, params={
        "n_samples": int(n_samples),
        "n_z": int(n_z),
        "m_int": int(m_int),
        "n_target": int(n_target),
        "maf_min": float(MAF_MIN),
        "chr21_window": list(CHR21_WINDOW),
        "chr22_window": list(CHR22_WINDOW),
        "seed": int(args.seed),
    })

    print("\nSummary:")
    print(f"  diff-chr  rho_max  mean={rho_diff.mean():.4f}  "
          f"max={rho_diff.max():.4f}")
    print(f"  same-chr  rho_max  mean={rho_same.mean():.4f}  "
          f"max={rho_same.max():.4f}")
    print(f"  Yelmen 2026 ref:  diff-chr mean=0.032, max=0.042; "
          f"same-chr mean=0.035, max=0.849")


if __name__ == "__main__":
    main()
