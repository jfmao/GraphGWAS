"""Generate supplementary Pan-UKB cross-ancestry fine-mapping figure.

Reads results/panukb/summary.json and produces a 2x2 panel figure
showing per-ancestry credible-set size, lead-variant PIP, and a
heatmap view.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    fig_dir = Path("/mnt/data/GraphGWAS/results/benchmark_v2/paper_figures")
    fig_dir.mkdir(parents=True, exist_ok=True)
    results_path = Path("/mnt/data/GraphGWAS/results/panukb/summary.json")
    with open(results_path) as f:
        data = json.load(f)

    loci = list(data.keys())
    ancestries = ["EUR", "CSA", "AFR", "EAS"]

    pip_matrix = np.zeros((len(loci), len(ancestries)))
    cs_matrix = np.zeros((len(loci), len(ancestries)), dtype=int)
    n_matrix = np.zeros((len(loci), len(ancestries)))
    for i, locus in enumerate(loci):
        for j, anc in enumerate(ancestries):
            rec = data[locus]["ancestries"].get(anc, {})
            pip_matrix[i, j] = rec.get("top_hbp_pip", np.nan) or np.nan
            cs_matrix[i, j] = rec.get("n_cs_hbp", 0) or 0
            n_matrix[i, j] = rec.get("N_sumstats", 0) or 0

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Panel A: top-variant PIP heatmap
    ax = axes[0, 0]
    im = ax.imshow(pip_matrix, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(len(ancestries)), ancestries)
    ax.set_yticks(range(len(loci)), [f"{l}\n({data[l]['trait']})" for l in loci])
    for i in range(len(loci)):
        for j in range(len(ancestries)):
            v = pip_matrix[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    color="white" if v < 0.5 else "black", fontsize=9)
    ax.set_title("(a) Top-variant HBP PIP")
    plt.colorbar(im, ax=ax, label="PIP")

    # Panel B: credible-set size
    ax = axes[0, 1]
    im = ax.imshow(np.log10(cs_matrix + 1), aspect="auto", cmap="magma_r")
    ax.set_xticks(range(len(ancestries)), ancestries)
    ax.set_yticks(range(len(loci)), [f"{l}\n({data[l]['trait']})" for l in loci])
    for i in range(len(loci)):
        for j in range(len(ancestries)):
            ax.text(j, i, f"{cs_matrix[i,j]}", ha="center", va="center",
                    color="white" if cs_matrix[i, j] > 50 else "black",
                    fontsize=9)
    ax.set_title("(b) 95% credible-set size")
    plt.colorbar(im, ax=ax, label="log₁₀(CS+1)")

    # Panel C: CS size vs sample size (log-log)
    ax = axes[1, 0]
    colors = {"EUR": "#1f77b4", "CSA": "#ff7f0e", "AFR": "#2ca02c", "EAS": "#d62728"}
    for j, anc in enumerate(ancestries):
        ns = n_matrix[:, j]
        css = cs_matrix[:, j]
        ax.scatter(ns, css, s=80, c=colors[anc], label=anc, alpha=0.75,
                   edgecolors="k", linewidth=0.5)
        for i, locus in enumerate(loci):
            ax.annotate(locus, (ns[i], css[i]), fontsize=7,
                        xytext=(3, 3), textcoords="offset points")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Pan-UKB N (per ancestry)")
    ax.set_ylabel("HBP credible-set size")
    ax.set_title("(c) Credible-set size scales with sample size")
    ax.legend(title="Ancestry")
    ax.grid(True, alpha=0.3)

    # Panel D: PIP across ancestries per locus (grouped bars)
    ax = axes[1, 1]
    width = 0.2
    x = np.arange(len(loci))
    for j, anc in enumerate(ancestries):
        vals = pip_matrix[:, j]
        ax.bar(x + (j - 1.5) * width, vals, width, label=anc,
               color=colors[anc], edgecolor="k", linewidth=0.5)
    ax.set_xticks(x, [f"{l}\n{data[l]['trait']}" for l in loci])
    ax.set_ylabel("Top-variant HBP PIP")
    ax.set_title("(d) Top-variant PIP by locus × ancestry")
    ax.axhline(0.95, ls="--", c="grey", lw=1, alpha=0.5, label="PIP=0.95")
    ax.legend(title="Ancestry", loc="upper right")
    ax.set_ylim(0, 1.05)

    fig.suptitle(
        "Supplementary Fig. S2 — Pan-UKB cross-ancestry fine-mapping preview\n"
        "HBP + 1KG-matched LD on 4 canonical loci × 4 ancestries (16 cells).\n"
        "EUR N ≈ 420 K; CSA 8.9 K; AFR 6.6 K; EAS 2.7 K. Sumstats via tabix.",
        fontsize=11,
    )
    fig.tight_layout()

    for ext in ("png", "pdf"):
        fig.savefig(fig_dir / f"figS2_panukb_cross_ancestry.{ext}", dpi=200,
                    bbox_inches="tight")
    print(f"Wrote {fig_dir}/figS2_panukb_cross_ancestry.{{png,pdf}}")


if __name__ == "__main__":
    main()
