"""Paper #2 Fig 7 (was placeholder) — 4-species §Y.8 catalogue grid.

Bar chart:
  - x-axis: species (yeast, Arabidopsis, rice, human)
  - y-axis: count of M5★ rank-1 / tested
  - 2 bars per species: M5★ vs M2 rank-1 counts
  - For human: a 3rd cluster showing weighted-substrate (M5★ binary +1)

Sources:
  - results/paper2_epistasis/multi_pair_catalogue_benchmark_4species.json
  - results/paper2_epistasis/multi_pair_catalogue_benchmark_human_binary.json
  - results/paper2_epistasis/multi_pair_catalogue_benchmark_human_weighted.json

Output:
  paper/epistasis_v1/figures/fig7_y8_4species_grid.{pdf,png,json}
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results" / "paper2_epistasis"
OUT_DIR = REPO / "paper" / "epistasis_v1" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Load
all4 = json.loads((RESULTS / "multi_pair_catalogue_benchmark_4species.json").read_text())
hum_w = json.loads((RESULTS / "multi_pair_catalogue_benchmark_human_weighted.json").read_text())


def count_pair_status(pairs):
    """Returns (n_total, n_tested, n_m5_rank1, n_m5_failed, n_m2_rank1)."""
    n_total = len(pairs)
    n_tested = 0
    n_m5_rank1 = 0
    n_m5_failed = 0
    n_m2_rank1 = 0
    for p in pairs:
        if p.get("status") != "TESTED":
            continue
        n_tested += 1
        m5r = p.get("M5_star", {}).get("ground_truth_rank")
        if m5r == 1:
            n_m5_rank1 += 1
        if m5r is None:
            n_m5_failed += 1
        m2r = p.get("M2", {}).get("ground_truth_rank")
        if m2r == 1:
            n_m2_rank1 += 1
    return n_total, n_tested, n_m5_rank1, n_m5_failed, n_m2_rank1


species_data = {}
for sp in ["yeast", "arabidopsis", "rice"]:
    pairs = all4[sp]
    species_data[sp] = count_pair_status(pairs)
# Human binary (from 4species file)
species_data["human (binary)"] = count_pair_status(all4["human"])
# Human weighted
species_data["human (weighted)"] = count_pair_status(hum_w["human"])

# --------- Plotting ---------
fig, ax = plt.subplots(1, 1, figsize=(9.5, 4.8))

species_order = ["yeast", "arabidopsis", "rice",
                 "human (binary)", "human (weighted)"]
labels = ["S. cerevisiae",
          "Arabidopsis thaliana",
          "Oryza sativa",
          "Homo sapiens\n(binary STRING-PPI)",
          "Homo sapiens\n(weighted +Billmann+DepMap)"]

x = np.arange(len(species_order))
bar_w = 0.30

# Per-species bars: M5★ rank-1, M2 rank-1, M5★ failed (NF/0)
m5_rank1 = [species_data[sp][2] for sp in species_order]
m2_rank1 = [species_data[sp][4] for sp in species_order]
m5_failed = [species_data[sp][3] for sp in species_order]
n_tested = [species_data[sp][1] for sp in species_order]
n_total = [species_data[sp][0] for sp in species_order]

ax.bar(x - bar_w, m5_rank1, bar_w,
        label=r"M5$^\star$ rank-1", color="#2E7D32", edgecolor="black", linewidth=0.5)
ax.bar(x, m5_failed, bar_w,
        label=r"M5$^\star$ failed (NF/0)", color="#C62828", edgecolor="black", linewidth=0.5)
ax.bar(x + bar_w, m2_rank1, bar_w,
        label="M2 rank-1", color="#1565C0", edgecolor="black", linewidth=0.5)

# Annotate: N tested label above each species cluster
for i, sp in enumerate(species_order):
    ax.text(i, max(m5_rank1[i], m5_failed[i], m2_rank1[i]) + 1.2,
            f"n={n_tested[i]}/{n_total[i]}",
            ha="center", va="bottom", fontsize=9, color="#444")
    # Highlight rescue: human weighted has +1 vs binary
    if sp == "human (weighted)":
        binary_rank1 = species_data["human (binary)"][2]
        diff = m5_rank1[i] - binary_rank1
        if diff > 0:
            ax.annotate(f"+{diff} rescue\n(NCAPD3×PARP1)",
                         xy=(i - bar_w, m5_rank1[i]),
                         xytext=(i - bar_w, m5_rank1[i] + 8),
                         ha="center", fontsize=8, color="#1B5E20",
                         arrowprops=dict(arrowstyle="->", color="#1B5E20", lw=1))

ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel("Count of testable catalogue pairs", fontsize=11)
ax.set_title(r"§Y.8 catalogue grid: M5$^\star$ vs M2 across 4 species "
              r"(94 testable pairs of 107 catalogue entries)",
              fontsize=11, pad=12)
ax.legend(loc="upper left", fontsize=9, framealpha=0.95)
ax.grid(axis="y", linestyle=":", alpha=0.4)
ax.set_ylim(0, max(max(m5_rank1), max(m5_failed), max(m2_rank1)) * 1.32)

plt.tight_layout()
fig.savefig(OUT_DIR / "fig7_y8_4species_grid.pdf", bbox_inches="tight")
fig.savefig(OUT_DIR / "fig7_y8_4species_grid.png", dpi=150, bbox_inches="tight")
print(f"Wrote {OUT_DIR / 'fig7_y8_4species_grid.pdf'}")
print(f"Wrote {OUT_DIR / 'fig7_y8_4species_grid.png'}")

# Save JSON metadata
metadata = {
    "figure": "Fig 7 — §Y.8 catalogue grid",
    "species_data": {sp: {"n_total": species_data[sp][0],
                           "n_tested": species_data[sp][1],
                           "M5_star_rank1": species_data[sp][2],
                           "M5_star_failed": species_data[sp][3],
                           "M2_rank1": species_data[sp][4]}
                      for sp in species_order},
    "totals": {
        "tested": sum(species_data[sp][1] for sp in ["yeast","arabidopsis","rice","human (binary)"]),
        "M5_star_rank1_binary": sum(species_data[sp][2] for sp in ["yeast","arabidopsis","rice","human (binary)"]),
        "M5_star_rank1_weighted": species_data["yeast"][2] + species_data["arabidopsis"][2] +
                                    species_data["rice"][2] + species_data["human (weighted)"][2],
    },
    "rescue_pairs": [
        {"pair_id": "human.billmann.Condensin_I-PARP-1-XRCC1_complex",
         "binary_M5_star": "NF/0", "weighted_M5_star": "1/21730",
         "g1": "NCAPD3", "g2": "PARP1"}
    ],
}
(OUT_DIR / "fig7_y8_4species_grid.json").write_text(json.dumps(metadata, indent=2))
print(f"Wrote {OUT_DIR / 'fig7_y8_4species_grid.json'}")

print()
print("Figure 7 summary:")
for sp in species_order:
    n_t, n_te, n_r1, n_f, n_m2 = species_data[sp]
    print(f"  {sp:<28}  total={n_t:>3}  tested={n_te:>3}  M5★ rank-1={n_r1:>3}  M5★ NF/0={n_f:>3}  M2 rank-1={n_m2:>3}")
