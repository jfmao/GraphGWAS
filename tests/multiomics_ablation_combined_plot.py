"""Combined multi-species ablation plot + summary.

Loads rice and human ablation TSVs, produces a 2-row / 3-col figure
that separates "resolvable" loci (baseline causal PIP ≥ threshold) from
"unresolvable" ones. This is the figure to include in the paper.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RICE_TSV = Path("/mnt/data/GraphGWAS/data/rice_3k/results/ablation_rice_combined.tsv")
HUMAN_TSV = Path("/mnt/data/GraphGWAS/results/ablation_human_chr22/ablation_human_combined.tsv")
OUT_DIR = Path("/mnt/data/GraphGWAS/data/rice_3k/results")

PIP_RESOLVABLE_THRESHOLD = 5e-3  # baseline causal_pip at 100% coverage


def load_species(tsv: Path, species: str) -> pd.DataFrame:
    if not tsv.exists():
        print(f"[warn] missing {tsv}")
        return pd.DataFrame()
    d = pd.read_csv(tsv, sep="\t")
    d["species"] = species
    return d


def classify_locus(df: pd.DataFrame) -> dict[str, str]:
    """Return dict {locus: 'resolvable' or 'unresolvable'} based on PIP at 100%."""
    at_100 = df[(df["method"]=="HBP") & (df["fraction"]==1.0)]
    cls = {}
    for locus, sub in at_100.groupby("locus"):
        pip = sub["causal_pip"].mean()
        cls[locus] = "resolvable" if pip >= PIP_RESOLVABLE_THRESHOLD else "unresolvable"
    return cls


def main():
    rice = load_species(RICE_TSV, "rice")
    human = load_species(HUMAN_TSV, "human")
    all_df = pd.concat([rice, human], ignore_index=True)
    if all_df.empty:
        print("No data.")
        return

    cls = classify_locus(all_df)
    all_df["class"] = all_df["locus"].map(cls)

    print("Locus classification:")
    for k, v in sorted(cls.items()):
        print(f"  {k:<30s} {v}")

    # Aggregate per (species, locus, fraction)
    agg = (all_df[all_df["method"]=="HBP"]
           .groupby(["species","locus","fraction","class"])
           .agg(cs_mean=("cs_size","mean"), cs_std=("cs_size","std"),
                pip_mean=("causal_pip","mean"), pip_std=("causal_pip","std"),
                rank_mean=("causal_rank","mean"))
           .reset_index())

    # Figure: 2 rows (resolvable, unresolvable) × 3 cols (CS, PIP, rank)
    fig, axes = plt.subplots(2, 3, figsize=(14, 8.5), sharex=True)
    row_titles = ["Resolvable loci (baseline PIP ≥ 0.005)",
                  "Unresolvable loci (baseline PIP < 0.005)"]

    colormap = {}
    palette = plt.cm.tab10(np.linspace(0, 1, 10))
    for i, locus in enumerate(sorted(agg["locus"].unique())):
        colormap[locus] = palette[i % 10]

    for r, cls_name in enumerate(["resolvable", "unresolvable"]):
        sub = agg[agg["class"]==cls_name]
        a1, a2, a3 = axes[r]
        for locus in sub["locus"].unique():
            s = sub[sub["locus"]==locus]
            sp = s["species"].iloc[0]
            label = f"{sp}:{locus}"
            c = colormap[locus]
            a1.errorbar(s["fraction"]*100, s["cs_mean"], yerr=s["cs_std"],
                        marker="o", label=label, capsize=3, color=c)
            a2.errorbar(s["fraction"]*100, s["pip_mean"], yerr=s["pip_std"],
                        marker="o", label=label, capsize=3, color=c)
            a3.plot(s["fraction"]*100, s["rank_mean"], marker="o",
                    label=label, color=c)
        a1.set_yscale("log"); a1.set_ylabel("95% CS size")
        a2.set_ylabel("Causal PIP")
        a3.set_yscale("log"); a3.set_ylabel("Causal rank")
        for a in (a1,a2,a3):
            a.grid(True, alpha=0.3)
            if r == 1: a.set_xlabel("Multi-omics edges retained (%)")
        a1.set_title(f"{row_titles[r]} — CS")
        a2.set_title(f"{row_titles[r]} — causal PIP")
        a3.set_title(f"{row_titles[r]} — causal rank")
        if len(sub) > 0:
            axes[r, 0].legend(loc="best", fontsize=7)

    plt.tight_layout()
    plt.savefig(OUT_DIR / "ablation_multi_species.png", dpi=150, bbox_inches="tight")
    plt.savefig(OUT_DIR / "ablation_multi_species.pdf", bbox_inches="tight")
    print(f"\nWrote {OUT_DIR}/ablation_multi_species.{{png,pdf}}")

    # Summary table: rank improvement per locus
    pivot = (agg.pivot_table(index=["species","locus","class"],
                             columns="fraction",
                             values="rank_mean")
             .reset_index())
    pivot["delta_rank"] = pivot[1.0] - pivot[0.1]
    pivot["fold_improvement"] = np.where(
        pivot[1.0] > 0,
        pivot[0.1] / pivot[1.0],
        np.nan,
    )
    print("\nRank-improvement summary:")
    print(pivot[["species","locus","class",0.1,1.0,"delta_rank","fold_improvement"]]
          .to_string(index=False))
    pivot.to_csv(OUT_DIR / "ablation_multi_species_summary.tsv",
                 sep="\t", index=False)


if __name__ == "__main__":
    main()
