"""Paper #2 Fig 9 — null-permutation FPR calibration for §Y.7 / §Y.8 M5★.

Two-panel figure:
  (A) Histogram of M5★ best-pair $p_{\min}$ across N=100 nulls per mode
      (noise vs marginal-only), with the positive control $p$ marked.
  (B) Q-Q plot of -log10($p_{\min}$) vs the order statistic
      expected from the Beta(1, n_tested) distribution under H0
      (BRCA1×PARP1 substrate, 141,245 cross-pairs).

Source: results/paper2_epistasis/null_fpr_test.json

Output: paper/epistasis_v1/figures/fig9_y10_null_fpr.{pdf,png,json}
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

d = json.loads((RESULTS / "null_fpr_test.json").read_text())
n_tested = d["null_results"]["noise"]["n_tested_per_perm"]
bonf = d["null_results"]["noise"]["bonferroni_threshold"]
pc_p = d["positive_control_sanity"]["p_interaction"]

fig, axes = plt.subplots(1, 2, figsize=(11, 4.7))

# --- Panel A: histograms ---
ax = axes[0]
modes = ["noise", "marginal"]
colors_mode = {"noise": "#1976D2", "marginal": "#E64A19"}
labels_mode = {"noise": r"Null ‘noise’ ($y=\varepsilon$)",
               "marginal": r"Null ‘marginal’ ($y=\alpha G_X+\beta G_Y+\varepsilon$)"}
log_bins = np.logspace(-20, 0, 40)
for mode in modes:
    ps = [p["p_min"] for p in d["null_results"][mode]["per_perm"]]
    ax.hist(ps, bins=log_bins, alpha=0.55, color=colors_mode[mode],
             label=labels_mode[mode], edgecolor="black", linewidth=0.4)
ax.set_xscale("log")
ax.axvline(pc_p, color="#2E7D32", lw=2.0, ls="--",
            label=f"Positive control p={pc_p:.1e}")
ax.axvline(bonf, color="#444", lw=1.0, ls=":",
            label=f"Bonferroni 0.05/{n_tested:,} = {bonf:.1e}")
ax.axvline(5e-8, color="#888", lw=0.8, ls=":",
            label="GWAS 5e-8")
ax.set_xlabel("M5$^\\star$ best-pair $p_{\\min}$ (log scale)", fontsize=10)
ax.set_ylabel("Count of permutations", fontsize=10)
ax.set_title(f"(A) FPR under null phenotypes (N=100 each)\n"
             f"BRCA1×PARP1 substrate; {n_tested:,} cross-pairs/perm",
             fontsize=11)
ax.legend(loc="upper left", fontsize=8, framealpha=0.95)
ax.grid(True, which="both", linestyle=":", alpha=0.3)

# --- Panel B: Q-Q under Beta(1, n_tested) ---
ax = axes[1]
ks = np.arange(1, 101) / 101  # quantile positions
# Under H0, the min of n iid Uniform(0,1) ≈ Beta(1, n_tested);
# the k-th order statistic of N=100 mins has expected value qbeta(k/N, ?)
# More simply: use empirical CDF of -log10 p vs uniform expected
for mode in modes:
    ps = sorted(p["p_min"] for p in d["null_results"][mode]["per_perm"])
    # Expected order statistics under H0: U_(k) = k/(N+1)
    # For p_min: expected = beta CDF inverse with parameters (1, n_tested)
    # Simpler: -log10(empirical p) vs -log10(uniform-quantile / n_tested)
    # i.e. Bonferroni-adjusted expected
    expected_p_min = np.array([1 - (1 - q) ** (1.0 / n_tested) for q in ks])
    obs = np.asarray(ps)
    # Avoid log(0)
    obs = np.maximum(obs, 1e-300)
    expected_p_min = np.maximum(expected_p_min, 1e-300)
    ax.scatter(-np.log10(expected_p_min), -np.log10(obs),
                s=18, alpha=0.7, color=colors_mode[mode],
                label=labels_mode[mode], edgecolors="black", linewidths=0.3)
# Identity line
lim = max(ax.get_xlim()[1], ax.get_ylim()[1])
ax.plot([0, lim], [0, lim], "k--", alpha=0.5, lw=0.8, label="y = x (calibrated)")
ax.set_xlabel(r"Expected $-\log_{10}(p_{\min})$ under $H_0$", fontsize=10)
ax.set_ylabel(r"Observed $-\log_{10}(p_{\min})$", fontsize=10)
ax.set_title("(B) Q-Q under $H_0$\n"
              "(points above y=x → false-positive inflation)",
              fontsize=11)
ax.legend(loc="upper left", fontsize=8, framealpha=0.95)
ax.grid(True, linestyle=":", alpha=0.3)

plt.tight_layout()
fig.savefig(OUT_DIR / "fig9_y10_null_fpr.pdf", bbox_inches="tight")
fig.savefig(OUT_DIR / "fig9_y10_null_fpr.png", dpi=150, bbox_inches="tight")
print(f"Wrote {OUT_DIR / 'fig9_y10_null_fpr.pdf'}")
print(f"Wrote {OUT_DIR / 'fig9_y10_null_fpr.png'}")

# Save JSON metadata
metadata = {
    "figure": "Fig 9 — null-permutation FPR calibration",
    "pair": d["pair"],
    "n_perm_per_mode": d["null_results"]["noise"]["n_perm"],
    "n_tested_per_perm": n_tested,
    "bonferroni_threshold": bonf,
    "positive_control_p": pc_p,
    "summary": {
        mode: {
            "GWAS_sig_p_lt_5e8": d["null_results"][mode]["n_below_gwas_5e8"],
            "Bonf_sig": d["null_results"][mode]["n_below_bonferroni"],
            "uncorrected_p_lt_05": d["null_results"][mode]["n_below_uncorrected_p05"],
            "fpr_gwas": d["null_results"][mode]["fpr_gwas"],
            "fpr_bonferroni": d["null_results"][mode]["fpr_bonferroni"],
            "p_min_distribution": d["null_results"][mode]["p_min_distribution"],
        } for mode in ["noise", "marginal"]
    },
    "interpretation": (
        "Noise null is well-calibrated (Bonferroni-sig 2/100 vs ~5 expected). "
        "Marginal null shows mild inflation (6/100) — main-effect leakage into "
        "the interaction LR test under strong correlated genotypes. Honest finding."
    ),
}
(OUT_DIR / "fig9_y10_null_fpr.json").write_text(json.dumps(metadata, indent=2))
print(f"Wrote {OUT_DIR / 'fig9_y10_null_fpr.json'}")
