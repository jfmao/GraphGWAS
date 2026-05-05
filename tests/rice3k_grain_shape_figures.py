"""Generate figures for the 3kRG grain weight + shape analysis.

Reads the GWAS sumstats + fine-mapping summary + recovery scorecard, and emits
paper-ready figures.

Outputs (paper/finemapping_v1/figures/):
  rice_grain_manhattan.{png,pdf}     — 4-trait Manhattan grid + lambda_GC inset
  rice_grain_qq.{png,pdf}            — 4-trait QQ panel
  rice_grain_recovery_scorecard.{png,pdf}  — bar chart, 5 methods × 3 tiers
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

DATA = Path("/mnt/data/GraphGWAS/data/rice_3k")
RES = DATA / "results"
GWAS_OUT = RES / "grain_gwas"
FM_OUT = RES / "grain_finemap"
FIG = Path("/mnt/data/GraphGWAS/paper/finemapping_v1/figures")
FIG.mkdir(parents=True, exist_ok=True)

TRAITS = ["TGW", "GL", "GW", "RLW"]
TRAIT_COLORS = {"TGW": "#1f77b4", "GL": "#2ca02c", "GW": "#d62728", "RLW": "#9467bd"}
GW_THRESHOLD = 5e-8
SUG_THRESHOLD = 1e-5


def _read_glm(trait):
    fp = GWAS_OUT / f"gwas_grain.{trait}.glm.linear"
    if not fp.exists():
        return None
    df = pd.read_csv(fp, sep="\t", low_memory=False, usecols=["#CHROM", "POS", "P", "TEST", "ERRCODE"])
    df.columns = [c.lstrip("#") for c in df.columns]
    df = df[(df["TEST"] == "ADD") & (df["ERRCODE"] == ".")]
    df["P"] = pd.to_numeric(df["P"], errors="coerce")
    df = df.dropna(subset=["P"])
    df["CHROM"] = df["CHROM"].astype(str).apply(
        lambda c: c if c.startswith("Chr") else f"Chr{c}"
    )
    return df[["CHROM", "POS", "P"]]


def _lambda_gc(P):
    chi2 = stats.norm.isf(P / 2) ** 2
    return float(np.median(chi2) / stats.chi2.ppf(0.5, df=1))


def _chrom_offsets(df):
    """Cumulative position offsets so the Manhattan plots concatenate chromosomes."""
    chroms = sorted(df["CHROM"].unique(), key=lambda c: int(c.replace("Chr", "")))
    offsets, cum = {}, 0
    for c in chroms:
        offsets[c] = cum
        cum += int(df.loc[df["CHROM"] == c, "POS"].max()) + 1_000_000
    return offsets, chroms, cum


def manhattan_grid():
    fig, axes = plt.subplots(4, 1, figsize=(12, 9), sharex=True)
    for i, trait in enumerate(TRAITS):
        df = _read_glm(trait)
        if df is None or df.empty:
            continue
        # downsample non-significant points to keep the figure light
        sig_mask = df["P"] < 1e-3
        df_sig = df[sig_mask]
        df_bg = df[~sig_mask].sample(n=min(50_000, (~sig_mask).sum()), random_state=0)
        plot_df = pd.concat([df_bg, df_sig], ignore_index=True)
        offsets, chroms, total = _chrom_offsets(df)
        plot_df["x"] = plot_df.apply(lambda r: r["POS"] + offsets[r["CHROM"]], axis=1)
        plot_df["color"] = plot_df["CHROM"].apply(
            lambda c: "#9aa6b2" if int(c.replace("Chr", "")) % 2 else "#cdd5dd"
        )
        # Color significant points by trait
        ax = axes[i]
        ax.scatter(plot_df["x"], -np.log10(plot_df["P"]),
                   c=plot_df["color"], s=2, alpha=0.6, rasterized=True)
        # Highlight significant
        s = plot_df[plot_df["P"] < GW_THRESHOLD]
        ax.scatter(s["x"], -np.log10(s["P"]), c=TRAIT_COLORS[trait], s=6, alpha=0.9, rasterized=True)
        ax.axhline(-np.log10(GW_THRESHOLD), color="r", lw=0.7, ls="--")
        ax.axhline(-np.log10(SUG_THRESHOLD), color="orange", lw=0.5, ls=":")
        lam = _lambda_gc(df["P"].values)
        ax.set_ylabel(f"$-\\log_{{10}}(P)$\n{trait}\n($\\lambda_{{GC}}$={lam:.2f})")
        ax.set_ylim(0, max(15, np.ceil(-np.log10(df["P"].min()))))
        if i == len(TRAITS) - 1:
            mids = [offsets[c] + (max(df.loc[df["CHROM"] == c, "POS"]) / 2) for c in chroms]
            ax.set_xticks(mids)
            ax.set_xticklabels([c.replace("Chr", "") for c in chroms], fontsize=9)
            ax.set_xlabel("Chromosome")
    fig.suptitle("3kRG grain weight + shape — multi-trait GWAS", y=0.995, fontsize=12)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = FIG / f"rice_grain_manhattan.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"  wrote {out.name}")
    plt.close(fig)


def qq_panel():
    """Standard GWAS Q-Q plot with log-rank subsampling.

    Naive uniform-rank subsampling misrepresents the curve: keeping the top
    K = 5000 ranks plus a linspace over [0, n-1] gives dense coverage in the
    upper tail (expected ~ 3.7-7.4) and very sparse coverage in the bulk
    (expected ~ 0-3), so the rendered curve looks like a steep ramp from
    origin instead of the canonical "diagonal in the bulk, deflect upward
    at the tail" shape. We use log-spaced ranks instead so points are
    distributed evenly along expected -log10(p), which is the natural axis
    for visual interpretation.
    """
    fig, axes = plt.subplots(2, 2, figsize=(8.5, 8.0))
    axes_flat = axes.flatten()
    for i, trait in enumerate(TRAITS):
        df = _read_glm(trait)
        if df is None or df.empty:
            continue
        P = np.sort(df["P"].values)
        n = len(P)
        expected = -np.log10((np.arange(1, n + 1)) / (n + 1))
        observed = -np.log10(P)

        # Log-rank subsampling: dense at extremes, log-spaced for bulk.
        # All top-K most-significant points + log-spaced indices for the rest.
        top_k = min(2000, n)
        log_n = max(top_k, 5000)
        log_idx = np.unique(np.round(
            np.logspace(np.log10(top_k), np.log10(n - 1), log_n)
        ).astype(int))
        idx = np.unique(np.r_[np.arange(top_k), log_idx])
        idx = idx[idx < n]

        ax = axes_flat[i]
        ax.scatter(expected[idx], observed[idx], s=6, c=TRAIT_COLORS[trait],
                   alpha=0.7, rasterized=True, edgecolor="none")
        # Diagonal y = x reference (the null expectation).
        # Cap axis limits so the diagonal stays visible regardless of inflation.
        x_max = float(expected.max()) * 1.05
        y_max = max(x_max, float(observed.max()) * 1.05)
        ax.plot([0, x_max], [0, x_max], "k--", lw=0.8, alpha=0.6,
                label="y = x (null)")
        # λ_GC reference line (slope = √λ_GC under chi-square inflation).
        lam = _lambda_gc(df["P"].values)
        ax.plot([0, x_max], [0, x_max * np.sqrt(lam)], color="#888888",
                lw=0.8, alpha=0.7, ls=":",
                label=f"slope = √λ_GC = {np.sqrt(lam):.2f}")
        ax.set_xlim(0, x_max)
        ax.set_ylim(0, y_max)
        ax.set_xlabel(r"Expected $-\log_{10}(P)$")
        ax.set_ylabel(r"Observed $-\log_{10}(P)$")
        ax.set_title(f"{trait}  ($\\lambda_{{GC}}$={lam:.2f})")
        ax.legend(fontsize=7, loc="upper left")
    fig.suptitle("3kRG grain — Q–Q plots", y=0.995)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = FIG / f"rice_grain_qq.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"  wrote {out.name}")
    plt.close(fig)


def recovery_scorecard():
    sc_path = FM_OUT / "grain_recovery_scorecard.tsv"
    if not sc_path.exists():
        print(f"  skip (no {sc_path.name})")
        return
    sc = pd.read_csv(sc_path, sep="\t")
    methods = sc["method"].tolist()
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    x = np.arange(len(methods))
    width = 0.27
    # Tier 1 — Niu 2021 21 QTNs
    t1_pct = sc["tier1_recovery_pct"].fillna(0)
    t2_pct = sc["tier2_recovery_pct"].fillna(0)
    # Tier 3 — normalize as fraction of max across methods (it's an absolute count)
    t3_total = sc["tier3_ren_genes_recovered"].fillna(0)
    t3_pct = 100 * t3_total / max(t3_total.max(), 1)

    ax.bar(x - width, t1_pct, width, label="Tier 1: Niu 2021 21 QTNs (CS)", color="#1f77b4")
    ax.bar(x, t2_pct, width, label="Tier 2: Niu 2021 7 NEW genes (locus)", color="#2ca02c")
    ax.bar(x + width, t3_pct, width,
           label=f"Tier 3: Ren 2023 catalogue (normalised; max = {int(t3_total.max())} genes)",
           color="#d62728")
    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=15, ha="right")
    ax.set_ylabel("Recovery (%)")
    ax.set_title("3kRG grain — fine-mapping recovery vs. Niu 2021 + Ren 2023 ground truth")
    ax.legend(fontsize=8, loc="upper right")
    ax.set_ylim(0, 105)
    ax.grid(axis="y", alpha=0.3, ls=":")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = FIG / f"rice_grain_recovery_scorecard.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"  wrote {out.name}")
    plt.close(fig)


def main():
    print("Manhattan grid:")
    manhattan_grid()
    print("Q-Q panel:")
    qq_panel()
    print("Recovery scorecard:")
    recovery_scorecard()


if __name__ == "__main__":
    main()
