"""v0.1.5 PIP calibration + null FPR benchmark on chr22 sims.

Two tests:
  1. Calibration (under H1 with single causal): bin all variant PIPs across
     50 reps × 3 h² levels and compute empirical TDR = fraction in bin that
     are the causal variant. Well-calibrated method → TDR ≈ midpoint(bin).
  2. Null FPR (under H0): simulate y = ε (no causal), count fraction of
     loci with any variant PIP ≥ {0.5, 0.9}. Well-behaved method → FPR low.

Methods: GAFM, HBP, GAFM-MX, HBP-MX, ENS, SuSiE, SuSiE-inf, FINEMAP-inf.

Output: results/benchmark_v2/v15_calibration_null/{calibration,null}.{json,tsv}
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src" / "python"))
sys.path.insert(0, "/mnt/data/GraphGWAS/src/python")

import susieinf  # noqa: E402
import finemapinf  # noqa: E402

from graphgwas.bgen_reader import BgenReader  # noqa: E402
from graphgwas.finemapping_v2 import (  # noqa: E402
    ensemble_from_sumstats,
    gafm_mx_from_sumstats,
    hbp_finemap_from_sumstats,
    hbp_mx_from_sumstats,
    l1_finemap_from_sumstats,
)

OUT = Path("/mnt/data/GraphGWAS/results/benchmark_v2/v15_calibration_null")
OUT.mkdir(parents=True, exist_ok=True)

N_REPS_H1 = int(os.environ.get("N_REPS_H1", "50"))
N_REPS_H0 = int(os.environ.get("N_REPS_H0", "100"))
H2_LEVELS = [float(x) for x in os.environ.get("H2_LEVELS", "0.02,0.05,0.10").split(",")]
WINDOW_BP = int(os.environ.get("WINDOW_BP", "50000"))
L_EFFECTS = 5
COVERAGE = 0.95


def _simulate_locus_h1(reader, rng, h2, window_bp=WINDOW_BP):
    """Single-causal simulation."""
    positions = reader._pos_array("22")
    pos_min, pos_max = positions[0] + window_bp, positions[-1] - window_bp
    for _ in range(10):
        center = int(rng.integers(pos_min, pos_max))
        vdf, dosage = reader.load_locus(
            "22", center - window_bp // 2, center + window_bp // 2,
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


def _simulate_locus_h0(reader, rng, window_bp=WINDOW_BP):
    """Null simulation: y = ε, no causal."""
    positions = reader._pos_array("22")
    pos_min, pos_max = positions[0] + window_bp, positions[-1] - window_bp
    for _ in range(10):
        center = int(rng.integers(pos_min, pos_max))
        vdf, dosage = reader.load_locus(
            "22", center - window_bp // 2, center + window_bp // 2,
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
        y = rng.standard_normal(dosage.shape[0])
        return vdf, dosage, y
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


def _pip_to_array(out, variants):
    by_id = {c.variant_id: c.pip for c in out}
    return np.array([by_id[v["variantId"]] for v in variants])


def _run_susie_rss(z, R, n):
    tmp = tempfile.mkdtemp(prefix="susie_rss_")
    np.savetxt(f"{tmp}/z.csv", z)
    np.savetxt(f"{tmp}/R.csv", R, delimiter=",")
    r_script = f'''
    library(susieR)
    z <- scan("{tmp}/z.csv")
    R <- as.matrix(read.csv("{tmp}/R.csv", header=FALSE))
    fit <- tryCatch(
        susie_rss(z=z, R=R, n={n}, L={L_EFFECTS}, coverage={COVERAGE}, verbose=FALSE),
        error=function(e) NULL
    )
    if(!is.null(fit)) write.table(data.frame(idx=seq_along(fit$pip), pip=fit$pip),
        file="{tmp}/pips.tsv", sep="\\t", row.names=FALSE, quote=FALSE)
    '''
    proc = subprocess.run(["R", "--no-save", "--no-restore", "-e", r_script],
                          capture_output=True, text=True, timeout=300)
    if not Path(f"{tmp}/pips.tsv").exists():
        return None
    return pd.read_csv(f"{tmp}/pips.tsv", sep="\t")["pip"].values


def _run_susie_inf(z, R, n):
    try:
        fit = susieinf.susie(z=z, meansq=1.0, n=n, L=L_EFFECTS, LD=R,
                              method="MLE", verbose=False, maxiter=50)
        raw_pip = fit["PIP"] if isinstance(fit, dict) else fit.PIP
        pips = (1 - np.prod(1 - np.nan_to_num(raw_pip, nan=0.0), axis=1)
                if raw_pip.ndim == 2 else np.nan_to_num(raw_pip, nan=0.0))
        return np.clip(pips, 0, 1)
    except Exception:
        return None


def _run_finemap_inf(z, R, n):
    try:
        fit = finemapinf.finemap(z=z, meansq=1.0, n=n, L=L_EFFECTS, LD=R,
                                  verbose=0, sched_sss=[50, 50, 500])
        raw_pip = fit["PIP"] if isinstance(fit, dict) else fit.PIP
        pips = (1 - np.prod(1 - np.nan_to_num(raw_pip, nan=0.0), axis=1)
                if raw_pip.ndim == 2 else np.nan_to_num(raw_pip, nan=0.0))
        return np.clip(pips, 0, 1)
    except Exception:
        return None


def _all_method_pips(variants, z, R, R_sq, n_samples):
    """Return dict {method_label: pips_array_in_variants_order}."""
    out = {}
    out["GAFM"] = _pip_to_array(l1_finemap_from_sumstats(
        variants, z, R_sq, alpha=0.7, credible_set_coverage=COVERAGE,
        chr_name="chr22"), variants)
    out["HBP"] = _pip_to_array(hbp_finemap_from_sumstats(
        variants, z, R_sq, graph_cache={}, credible_set_coverage=COVERAGE,
        chr_name="chr22"), variants)
    out["GAFM-MX"] = _pip_to_array(gafm_mx_from_sumstats(
        variants, z, R_sq, n_samples=n_samples, alpha=0.7,
        credible_set_coverage=COVERAGE, chr_name="chr22"), variants)
    out["HBP-MX"] = _pip_to_array(hbp_mx_from_sumstats(
        variants, z, R_sq, n_samples=n_samples,
        credible_set_coverage=COVERAGE, chr_name="chr22"), variants)
    out["ENS"] = _pip_to_array(ensemble_from_sumstats(
        variants, z, R_sq, n_samples=n_samples, alpha=0.7,
        credible_set_coverage=COVERAGE, chr_name="chr22"), variants)
    p = _run_susie_rss(z, R, n_samples)
    if p is not None:
        out["SuSiE"] = p
    p = _run_susie_inf(z, R, n_samples)
    if p is not None:
        out["SuSiE-inf"] = p
    p = _run_finemap_inf(z, R, n_samples)
    if p is not None:
        out["FINEMAP-inf"] = p
    return out


def calibration_run(reader):
    """Per-method calibration: bin all variant PIPs and compute empirical TDR."""
    bin_edges = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.01])
    bin_mid = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    counts = {}  # method -> (n_bins, count_total, count_causal)

    seed_offset = 5000
    for h2 in H2_LEVELS:
        for i in range(N_REPS_H1):
            rng = np.random.default_rng(seed_offset + int(h2 * 1000) * 100 + i)
            sim = _simulate_locus_h1(reader, rng, h2)
            if sim is None:
                continue
            vdf, X, y, causal_idx = sim
            z, R, R_sq, n = _summary_inputs(X, y)
            variants = _build_variants(vdf, X)
            try:
                method_pips = _all_method_pips(variants, z, R, R_sq, n)
            except Exception as e:
                print(f"  h2={h2} rep {i}: skipping ({e})", flush=True)
                continue
            for method, pips in method_pips.items():
                if method not in counts:
                    counts[method] = {
                        "total": np.zeros(len(bin_mid), dtype=int),
                        "causal": np.zeros(len(bin_mid), dtype=int),
                    }
                bin_idx = np.digitize(pips, bin_edges) - 1
                bin_idx = np.clip(bin_idx, 0, len(bin_mid) - 1)
                for j in range(len(pips)):
                    b = bin_idx[j]
                    counts[method]["total"][b] += 1
                    if j == causal_idx:
                        counts[method]["causal"][b] += 1
            if (i + 1) % 10 == 0:
                print(f"  h2={h2} {i+1}/{N_REPS_H1} reps done", flush=True)

    cal = {"bin_edges": bin_edges.tolist(), "bin_mid": bin_mid.tolist(),
           "h2_levels": H2_LEVELS, "n_reps_per_h2": N_REPS_H1, "methods": {}}
    for method, c in counts.items():
        tdr = np.where(c["total"] > 0, c["causal"] / np.maximum(c["total"], 1), np.nan)
        cal["methods"][method] = {
            "total": c["total"].tolist(),
            "causal": c["causal"].tolist(),
            "tdr": tdr.tolist(),
        }
    return cal


def null_fpr_run(reader):
    """Per-method null FPR: count loci with any variant PIP ≥ {0.5, 0.9}."""
    thresholds = [0.5, 0.9]
    summary = {}  # method -> {threshold: count_loci_with_any_above}

    for i in range(N_REPS_H0):
        rng = np.random.default_rng(9000 + i)
        sim = _simulate_locus_h0(reader, rng)
        if sim is None:
            continue
        vdf, X, y = sim
        z, R, R_sq, n = _summary_inputs(X, y)
        variants = _build_variants(vdf, X)
        try:
            method_pips = _all_method_pips(variants, z, R, R_sq, n)
        except Exception as e:
            print(f"  null rep {i}: skipping ({e})", flush=True)
            continue
        for method, pips in method_pips.items():
            if method not in summary:
                summary[method] = {"n_loci": 0, "thresholds": {t: 0 for t in thresholds}}
            summary[method]["n_loci"] += 1
            for t in thresholds:
                if (pips >= t).any():
                    summary[method]["thresholds"][t] += 1
        if (i + 1) % 10 == 0:
            print(f"  null {i+1}/{N_REPS_H0} reps done", flush=True)

    out = {"thresholds": thresholds, "n_reps": N_REPS_H0, "methods": {}}
    for method, c in summary.items():
        out["methods"][method] = {
            "n_loci": c["n_loci"],
            "fpr": {str(t): c["thresholds"][t] / max(c["n_loci"], 1) for t in thresholds},
        }
    return out


def main():
    reader = BgenReader("/mnt/data/GraphGWAS/tests/data/human/1kGP_bgen")
    print(f"v0.1.5 calibration + null FPR — N_H1={N_REPS_H1} × {len(H2_LEVELS)} h²; "
          f"N_H0={N_REPS_H0}", flush=True)

    print("\n[Calibration: H1 single causal across h² levels]", flush=True)
    cal = calibration_run(reader)
    json.dump(cal, open(OUT / "calibration.json", "w"), indent=2)
    rows = []
    for method, c in cal["methods"].items():
        for j, mid in enumerate(cal["bin_mid"]):
            rows.append({
                "method": method, "bin_mid": mid, "total": c["total"][j],
                "causal": c["causal"][j], "tdr": c["tdr"][j],
            })
    pd.DataFrame(rows).to_csv(OUT / "calibration.tsv", sep="\t", index=False)
    print("\nCalibration (TDR per PIP bin):", flush=True)
    print(f"  {'method':<14} " + " ".join(f"{m:>5.2f}" for m in cal["bin_mid"]))
    for method, c in cal["methods"].items():
        tdr_str = " ".join(f"{t:>5.3f}" if not np.isnan(t) else "  nan" for t in c["tdr"])
        print(f"  {method:<14} {tdr_str}")

    print("\n[Null FPR: H0 pure noise]", flush=True)
    null = null_fpr_run(reader)
    json.dump(null, open(OUT / "null.json", "w"), indent=2)
    rows = []
    for method, c in null["methods"].items():
        for t in null["thresholds"]:
            rows.append({"method": method, "threshold": t, "fpr": c["fpr"][str(t)],
                         "n_loci": c["n_loci"]})
    pd.DataFrame(rows).to_csv(OUT / "null.tsv", sep="\t", index=False)
    print("\nNull FPR (fraction of null loci with any PIP ≥ threshold):", flush=True)
    print(f"  {'method':<14} {'PIP≥0.5':>8} {'PIP≥0.9':>8}")
    for method, c in null["methods"].items():
        print(f"  {method:<14} {c['fpr']['0.5']:>8.3f} {c['fpr']['0.9']:>8.3f}")


if __name__ == "__main__":
    main()
