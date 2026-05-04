"""Paper #2 Fig 8 — binary STRING-PPI vs weighted Billmann+DepMap+STRING
substrate, head-to-head on human (§Y.8 + §Y.9).

Two-panel figure:
  (A) Per-pair scatter: binary M5★ best-pair p-value (x) vs weighted (y),
      with the +1 rescue pair highlighted.
  (B) Bar chart of M5★ rank-1 / M5★ NF/0 counts before/after the
      weighted substrate (binary 50/20, weighted 51/19).

Sources:
  results/paper2_epistasis/multi_pair_catalogue_benchmark_human_binary.json
  results/paper2_epistasis/multi_pair_catalogue_benchmark_human_weighted.json

Output:
  paper/epistasis_v1/figures/fig8_y8_binary_vs_weighted.{pdf,png,json}
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

binary = json.loads((RESULTS / "multi_pair_catalogue_benchmark_human_binary.json").read_text())
weighted = json.loads((RESULTS / "multi_pair_catalogue_benchmark_human_weighted.json").read_text())

bp = {r["pair_id"]: r for r in binary["human"]}
wp = {r["pair_id"]: r for r in weighted["human"]}

xs, ys, labels, colors = [], [], [], []
rescue_pair = None
for pid, b_r in bp.items():
    w_r = wp.get(pid)
    if not w_r or b_r.get("status") != "TESTED":
        continue
    bp_p = b_r.get("M5_star", {}).get("ground_truth_p")
    wp_p = w_r.get("M5_star", {}).get("ground_truth_p")
    if bp_p is None and wp_p is None:
        continue  # both NF; not informative
    if bp_p is None:
        # Rescue: NF→hit. Plot at edge with annotation
        xs.append(1e-1)  # placeholder x for "missing"
        ys.append(wp_p)
        labels.append(pid)
        colors.append("#1B5E20")
        rescue_pair = (pid, wp_p,
                        b_r.get("g1_symbol"), b_r.get("g2_symbol"))
        continue
    if wp_p is None:
        # Regression: hit→NF (we don't expect any)
        xs.append(bp_p); ys.append(1e-1); labels.append(pid)
        colors.append("#C62828")
        continue
    xs.append(bp_p); ys.append(wp_p); labels.append(pid)
    colors.append("#37474F")

xs = np.asarray(xs); ys = np.asarray(ys)
print(f"  pairs plotted: {len(xs)}")
print(f"  rescue pair: {rescue_pair}")

# --------- Plot ---------
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7),
                           gridspec_kw={"width_ratios": [1.4, 1]})

# Panel A: scatter
ax = axes[0]
ax.scatter(xs, ys, c=colors, s=22, alpha=0.6, edgecolors="black", linewidths=0.4)
# Identity line
lim = max(np.nanmax(xs), np.nanmax(ys))
ax.plot([1e-300, 1.0], [1e-300, 1.0], "k--", alpha=0.4, lw=0.8, label="y = x")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(1e-12, 2)
ax.set_ylim(1e-12, 2)
ax.set_xlabel("M5$^\\star$ best-pair $p$ (binary STRING-PPI)", fontsize=10)
ax.set_ylabel("M5$^\\star$ best-pair $p$ (weighted +Billmann+DepMap)", fontsize=10)
ax.set_title("(A) Per-pair $p_{\\min}$: binary vs weighted substrate",
              fontsize=11)
ax.grid(True, which="both", linestyle=":", alpha=0.3)
# Annotate the rescue pair
if rescue_pair:
    pid, wp_p, g1, g2 = rescue_pair
    ax.annotate(f"+1 rescue: {g1}×{g2}\n(NF/0 → rank 1, p={wp_p:.1e})",
                 xy=(0.1, wp_p),
                 xytext=(0.001, wp_p * 100),
                 fontsize=8, color="#1B5E20",
                 arrowprops=dict(arrowstyle="->", color="#1B5E20", lw=1))
# Mark GWAS thresholds
for thr, lbl in [(5e-8, "GWAS 5e-8"), (0.05, "α=0.05")]:
    ax.axvline(thr, color="gray", lw=0.4, ls=":", alpha=0.5)
    ax.axhline(thr, color="gray", lw=0.4, ls=":", alpha=0.5)
ax.legend(loc="lower right", fontsize=8)

# Panel B: rank-1 / NF count comparison
ax = axes[1]
b_rank1 = sum(1 for r in binary["human"]
                if r.get("status")=="TESTED"
                and r.get("M5_star", {}).get("ground_truth_rank")==1)
w_rank1 = sum(1 for r in weighted["human"]
                if r.get("status")=="TESTED"
                and r.get("M5_star", {}).get("ground_truth_rank")==1)
b_failed = sum(1 for r in binary["human"]
                if r.get("status")=="TESTED"
                and r.get("M5_star", {}).get("ground_truth_rank") is None)
w_failed = sum(1 for r in weighted["human"]
                if r.get("status")=="TESTED"
                and r.get("M5_star", {}).get("ground_truth_rank") is None)
labels_b = ["Binary STRING", "Weighted +Billmann\n+DepMap+STRING"]
x = np.arange(2)
bw = 0.35
ax.bar(x - bw/2, [b_rank1, w_rank1], bw, label=r"M5$^\star$ rank-1",
        color="#2E7D32", edgecolor="black", linewidth=0.5)
ax.bar(x + bw/2, [b_failed, w_failed], bw, label=r"M5$^\star$ NF/0",
        color="#C62828", edgecolor="black", linewidth=0.5)
for i, (r1, fl) in enumerate([(b_rank1, b_failed), (w_rank1, w_failed)]):
    ax.text(i - bw/2, r1 + 1, str(r1), ha="center", fontsize=10, fontweight="bold")
    ax.text(i + bw/2, fl + 1, str(fl), ha="center", fontsize=10, fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels(labels_b, fontsize=9)
ax.set_ylabel("Count of testable human pairs (n=70)", fontsize=10)
ax.set_title(f"(B) M5$^\\star$ recovery: +{w_rank1-b_rank1} rescue, "
              f"−{b_failed-w_failed} NF/0",
              fontsize=11)
ax.legend(loc="upper right", fontsize=9)
ax.grid(axis="y", linestyle=":", alpha=0.3)
ax.set_ylim(0, max(b_rank1, w_rank1) * 1.18)

plt.tight_layout()
fig.savefig(OUT_DIR / "fig8_y8_binary_vs_weighted.pdf", bbox_inches="tight")
fig.savefig(OUT_DIR / "fig8_y8_binary_vs_weighted.png", dpi=150, bbox_inches="tight")
print(f"Wrote {OUT_DIR / 'fig8_y8_binary_vs_weighted.pdf'}")
print(f"Wrote {OUT_DIR / 'fig8_y8_binary_vs_weighted.png'}")

# Save JSON metadata
metadata = {
    "figure": "Fig 8 — binary STRING-PPI vs weighted Billmann+DepMap substrate",
    "binary": {"M5_star_rank1": b_rank1, "M5_star_NF": b_failed},
    "weighted": {"M5_star_rank1": w_rank1, "M5_star_NF": w_failed},
    "rescue_pair": {"id": rescue_pair[0], "weighted_p": rescue_pair[1],
                     "g1": rescue_pair[2], "g2": rescue_pair[3]} if rescue_pair else None,
    "n_pairs_plotted_panel_A": int(len(xs)),
    "n_pairs_with_both_p": int(sum(1 for c in colors if c == "#37474F")),
}
(OUT_DIR / "fig8_y8_binary_vs_weighted.json").write_text(json.dumps(metadata, indent=2))
print(f"Wrote {OUT_DIR / 'fig8_y8_binary_vs_weighted.json'}")
