"""v0.1.5 hyperparameter sensitivity: vary mixture π and γ on chr22 sims.

Tests robustness of GAFM-MX to perturbations of the SBayesRC default
mixture weights and effect-size scales. Reports rank-1 rate and mean PIP
at causal across 30 reps × 5 hyperparameter perturbations × 2 h² levels.

Default: π=(0.005, 0.003, 0.001, 0.001), γ=(0.001, 0.01, 0.1, 1.0)

Perturbations:
  - default
  - π_uniform: π=(1, 1, 1, 1)/4 (equal weights)
  - π_concentrated: π=(0.001, 0.001, 0.001, 0.005) (more weight on largest variance)
  - γ_shifted_left: γ=(0.0001, 0.001, 0.01, 0.1) (everything 10× smaller)
  - γ_shifted_right: γ=(0.01, 0.1, 1.0, 10.0) (everything 10× larger)

Output: results/benchmark_v2/v15_hyperparam_sensitivity/sensitivity.{json,tsv}
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src" / "python"))
sys.path.insert(0, "/mnt/data/GraphGWAS/src/python")

from graphgwas.bgen_reader import BgenReader  # noqa: E402
from graphgwas.finemapping_v2 import (  # noqa: E402
    DEFAULT_MIXTURE_GAMMA,
    DEFAULT_MIXTURE_PI,
    gafm_mx_from_sumstats,
)

OUT = Path("/mnt/data/GraphGWAS/results/benchmark_v2/v15_hyperparam_sensitivity")
OUT.mkdir(parents=True, exist_ok=True)

N_REPS = int(os.environ.get("N_REPS", "30"))
H2_LEVELS = [float(x) for x in os.environ.get("H2_LEVELS", "0.05,0.10").split(",")]
WINDOW_BP = 50000
COVERAGE = 0.95


PERTURBATIONS = {
    "default": {"pi": DEFAULT_MIXTURE_PI, "gamma": DEFAULT_MIXTURE_GAMMA},
    "pi_uniform": {"pi": (0.25, 0.25, 0.25, 0.25), "gamma": DEFAULT_MIXTURE_GAMMA},
    "pi_concentrated": {"pi": (0.001, 0.001, 0.001, 0.005), "gamma": DEFAULT_MIXTURE_GAMMA},
    "gamma_left": {"pi": DEFAULT_MIXTURE_PI, "gamma": (0.0001, 0.001, 0.01, 0.1)},
    "gamma_right": {"pi": DEFAULT_MIXTURE_PI, "gamma": (0.01, 0.1, 1.0, 10.0)},
}


def _simulate_locus(reader, rng, h2):
    positions = reader._pos_array("22")
    pos_min, pos_max = positions[0] + WINDOW_BP, positions[-1] - WINDOW_BP
    for _ in range(10):
        center = int(rng.integers(pos_min, pos_max))
        vdf, dosage = reader.load_locus(
            "22", center - WINDOW_BP // 2, center + WINDOW_BP // 2,
            format="dosage",
        )
        if len(vdf) == 0:
            continue
        col_mean = np.nanmean(dosage, axis=0)
        dosage = np.where(np.isnan(dosage), col_mean, dosage)
        af = col_mean / 2.0
        maf = np.minimum(af, 1 - af)
        keep = (maf > 0.05) & np.isfinite(col_mean)
        if keep.sum() < 30:
            continue
        dosage = dosage[:, keep]
        vdf = vdf.loc[keep].reset_index(drop=True)
        af_keep = af[keep]
        cands = np.where((af_keep > 0.1) & (af_keep < 0.4))[0]
        if len(cands) == 0:
            continue
        causal_idx = int(rng.choice(cands))
        g = dosage[:, causal_idx]
        g_c = g - g.mean()
        var_g = g_c.var()
        beta = np.sqrt(h2 / var_g) if var_g > 1e-8 else 0
        eps = rng.standard_normal(dosage.shape[0]) * np.sqrt(1 - h2)
        y = beta * g_c + eps
        return vdf, dosage, y, causal_idx
    return None


def _summary_inputs(X, y):
    n, p = X.shape
    Xc = X - X.mean(axis=0, keepdims=True)
    yc = y - y.mean()
    sx = Xc.std(axis=0)
    sx = np.where(sx < 1e-10, 1.0, sx)
    Xs = Xc / sx
    sy = yc.std()
    if sy < 1e-10:
        sy = 1.0
    ys = yc / sy
    z = (Xs.T @ ys) / np.sqrt(n)
    R = (Xs.T @ Xs) / n
    R = np.clip(R, -1.0, 1.0)
    np.fill_diagonal(R, 1.0)
    R_sq = np.clip(R ** 2, 0, 1)
    np.fill_diagonal(R_sq, 1.0)
    return z, R, R_sq, n


def _build_variants(vdf, X):
    return [
        {
            "variantId": f"chr22:{int(vdf.iloc[j]['pos'])}:{vdf.iloc[j]['a1']}:{vdf.iloc[j]['a2']}",
            "chr": "chr22",
            "pos": int(vdf.iloc[j]["pos"]),
            "ref": vdf.iloc[j]["a1"],
            "alt": vdf.iloc[j]["a2"],
            "af_total": float(X[:, j].mean() / 2),
        }
        for j in range(X.shape[1])
    ]


def main():
    reader = BgenReader("/mnt/data/GraphGWAS/tests/data/human/1kGP_bgen")
    print(f"v0.1.5 hyperparameter sensitivity — N_REPS={N_REPS} × {len(H2_LEVELS)} h²"
          f" × {len(PERTURBATIONS)} perturbations", flush=True)

    rows = []
    for h2 in H2_LEVELS:
        for i in range(N_REPS):
            rng = np.random.default_rng(7000 + int(h2 * 1000) + i)
            sim = _simulate_locus(reader, rng, h2)
            if sim is None:
                continue
            vdf, X, y, causal_idx = sim
            z, R, R_sq, n = _summary_inputs(X, y)
            variants = _build_variants(vdf, X)

            for pert_name, params in PERTURBATIONS.items():
                t0 = time.time()
                try:
                    out = gafm_mx_from_sumstats(
                        variants, z, R_sq, n_samples=n, alpha=0.7,
                        credible_set_coverage=COVERAGE, chr_name="chr22",
                        mixture_pi=params["pi"], mixture_gamma=params["gamma"],
                    )
                    by_id = {c.variant_id: c.pip for c in out}
                    pips = np.array([by_id[v["variantId"]] for v in variants])
                    pip_c = float(pips[causal_idx])
                    rank = int((pips >= pip_c).sum())
                except Exception:
                    pip_c = 0.0
                    rank = -1
                rows.append({
                    "h2": h2, "rep": i, "perturbation": pert_name,
                    "rank": rank, "pip_at_causal": pip_c,
                    "time_s": time.time() - t0,
                })
        print(f"  h2={h2} done", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "sensitivity.tsv", sep="\t", index=False)

    print(f"\n{'h2':>5} {'perturbation':<18} {'rank-1':>8} {'mean_pip':>9} {'med_t (s)':>10}")
    summary = {}
    for h2 in H2_LEVELS:
        for p in PERTURBATIONS.keys():
            sub = df[(df["h2"] == h2) & (df["perturbation"] == p)]
            r1 = int((sub["rank"] == 1).sum())
            mp = float(sub["pip_at_causal"].mean())
            mt = float(sub["time_s"].median())
            print(f"{h2:>5.2f} {p:<18} {r1:>4}/{len(sub):<3} {mp:>9.3f} {mt:>10.4f}")
            summary[f"{h2}_{p}"] = {
                "h2": h2, "perturbation": p, "n": len(sub),
                "rank_1": r1, "mean_pip": mp, "median_time_s": mt,
            }
    json.dump({"config": {"n_reps": N_REPS, "h2_levels": H2_LEVELS,
                          "perturbations": list(PERTURBATIONS.keys())},
               "summary": summary}, open(OUT / "sensitivity.json", "w"), indent=2)
    print(f"\nWrote {OUT}/sensitivity.{{json,tsv}}")


if __name__ == "__main__":
    main()
