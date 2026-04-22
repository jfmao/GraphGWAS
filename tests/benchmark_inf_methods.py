"""Benchmark vs Cui et al. 2024 infinitesimal-effects baselines.

Adds SuSiE-inf and FINEMAP-inf (FinucaneLab/fine-mapping-inf) to the
existing SuSiE / FINEMAP / L1 / HBP comparison on 1KG chr22 F1
single-causal-in-LD simulations.

These baselines model a polygenic background term that handles non-
sparse architectures — Wu et al. 2026 (Fig. 4b) shows they are the
state-of-art region-specific fine-mappers.

Output: results/benchmark_v2/inf_methods/inf_methods.{json,tsv}
"""
from __future__ import annotations

import bisect
import json
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

import susieinf
import finemapinf

from graphgwas.bgen_reader import BgenReader
from graphgwas.finemapping_v2 import fast_hbp_finemap

OUT = Path("/mnt/data/GraphGWAS/results/benchmark_v2/inf_methods")
OUT.mkdir(parents=True, exist_ok=True)


def _simulate_locus(reader: BgenReader, rng, h2: float, window_bp: int = 50_000):
    """Pick a 50 kb chr22 window and simulate a single-causal phenotype."""
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
        # Impute NaN with column mean, MAF filter
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
        # Pick causal at AF 0.1-0.4
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


def _run_susie_r(X: np.ndarray, y: np.ndarray, causal_idx: int) -> dict:
    """Vanilla SuSiE via R subprocess."""
    tmp = tempfile.mkdtemp()
    np.savetxt(f"{tmp}/X.csv", X, delimiter=",")
    np.savetxt(f"{tmp}/y.csv", y)
    r_script = f'''
    library(susieR)
    X <- as.matrix(read.csv("{tmp}/X.csv", header=FALSE))
    y <- scan("{tmp}/y.csv")
    fit <- tryCatch(susie(X, y, L=5, verbose=FALSE), error=function(e) NULL)
    if(!is.null(fit)) {{
        pips <- fit$pip
        ci <- {causal_idx + 1}
        cat(sprintf("PIP=%.6f\\nRANK=%d\\n", pips[ci], sum(pips >= pips[ci])))
    }} else cat("FAILED\\n")
    '''
    t0 = time.time()
    proc = subprocess.run(["R", "--no-save", "-e", r_script],
                          capture_output=True, text=True, timeout=60)
    dt = time.time() - t0
    rank = pip = None
    for line in proc.stdout.split("\n"):
        if line.startswith("PIP="):
            pip = float(line.split("=")[1])
        elif line.startswith("RANK="):
            rank = int(line.split("=")[1])
    return {"rank": rank, "pip": pip, "time": dt}


def _summary_inputs(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float, np.ndarray]:
    """Compute (z, meansq, LD) from raw genotype matrix and phenotype."""
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
    meansq = (ys ** 2).mean()
    LD = (Xs.T @ Xs) / n
    return z, float(meansq), LD


def _run_susie_inf(X: np.ndarray, y: np.ndarray, causal_idx: int, L: int = 5) -> dict:
    z, meansq, LD = _summary_inputs(X, y)
    n = X.shape[0]
    t0 = time.time()
    try:
        fit = susieinf.susie(z=z, meansq=meansq, n=n, L=L, LD=LD,
                              method="moments", verbose=False, maxiter=50)
    except Exception as e:
        return {"rank": None, "pip": 0, "time": time.time() - t0, "error": str(e)}
    dt = time.time() - t0
    raw_pip = fit["PIP"] if isinstance(fit, dict) else fit.PIP
    # SuSiE returns alpha per (variant, effect); per-variant PIP = 1 - prod_l (1-alpha_il)
    if raw_pip.ndim == 2:
        pips = 1 - np.prod(1 - raw_pip, axis=1)
    else:
        pips = raw_pip
    pip_c = float(pips[causal_idx])
    rank = int((pips >= pip_c).sum())
    return {"rank": rank, "pip": pip_c, "time": dt}


def _run_finemap_inf(X: np.ndarray, y: np.ndarray, causal_idx: int, L: int = 5) -> dict:
    z, meansq, LD = _summary_inputs(X, y)
    n = X.shape[0]
    t0 = time.time()
    try:
        fit = finemapinf.finemap(z=z, meansq=meansq, n=n, L=L, LD=LD,
                                  verbose=0, sched_sss=[50, 50, 1000])
    except Exception as e:
        return {"rank": None, "pip": 0, "time": time.time() - t0, "error": str(e)}
    dt = time.time() - t0
    raw_pip = fit["PIP"] if isinstance(fit, dict) else fit.PIP
    if raw_pip.ndim == 2:
        pips = 1 - np.prod(1 - raw_pip, axis=1)
    else:
        pips = raw_pip
    pip_c = float(pips[causal_idx])
    rank = int((pips >= pip_c).sum())
    return {"rank": rank, "pip": pip_c, "time": dt}


def _run_l1_or_hbp(variants_df, X, y, causal_idx, fast_hbp: bool = True):
    """Run HBP via fast_hbp_finemap (no annotations → pure stat L1 / HBP equiv)."""
    variants = []
    for j in range(X.shape[1]):
        row = variants_df.iloc[j]
        variants.append({
            "variantId": f"chr22:{int(row['pos'])}:{row['a1']}:{row['a2']}",
            "chr": "chr22",
            "pos": int(row["pos"]),
            "ref": row["a1"],
            "alt": row["a2"],
            "af_total": float(X[:, j].mean() / 2),
            "dosage": X[:, j],
        })
    t0 = time.time()
    cands = fast_hbp_finemap(variants, y, graph_cache={}, n_rounds=3, chr_name="chr22")
    dt = time.time() - t0
    causal_pos = variants[causal_idx]["pos"]
    cands_sorted = sorted(cands, key=lambda c: -c.pip)
    rank = pip = None
    for r, c in enumerate(cands_sorted, 1):
        if c.pos == causal_pos:
            rank = r
            pip = float(c.pip)
            break
    return {"rank": rank, "pip": pip, "time": dt}


def main():
    reader = BgenReader("/mnt/data/GraphGWAS/tests/data/human/1kGP_bgen")
    n_reps = 30
    h2 = 0.05  # Moderate weak signal
    results = []
    t_start = time.time()
    for i in range(n_reps):
        rng = np.random.default_rng(4000 + i)
        sim = _simulate_locus(reader, rng, h2)
        if sim is None:
            continue
        vdf, X, y, causal_idx = sim
        n_var = X.shape[1]
        causal_pos = int(vdf.iloc[causal_idx]["pos"])

        susie = _run_susie_r(X, y, causal_idx)
        susie_inf = _run_susie_inf(X, y, causal_idx)
        finemap_inf = _run_finemap_inf(X, y, causal_idx)
        hbp = _run_l1_or_hbp(vdf, X, y, causal_idx)

        row = {
            "rep": i,
            "n_var": n_var,
            "causal_pos": causal_pos,
            "causal_af": float(X[:, causal_idx].mean() / 2),
            "susie_rank": susie["rank"], "susie_pip": susie["pip"], "susie_time": susie["time"],
            "susie_inf_rank": susie_inf["rank"], "susie_inf_pip": susie_inf["pip"], "susie_inf_time": susie_inf["time"],
            "finemap_inf_rank": finemap_inf["rank"], "finemap_inf_pip": finemap_inf["pip"], "finemap_inf_time": finemap_inf["time"],
            "hbp_rank": hbp["rank"], "hbp_pip": hbp["pip"], "hbp_time": hbp["time"],
        }
        results.append(row)
        print(f"rep {i:2d}: nv={n_var:4d} | SuSiE=#{susie['rank']} "
              f"SuSiE-inf=#{susie_inf['rank']} FINEMAP-inf=#{finemap_inf['rank']} "
              f"HBP=#{hbp['rank']}")

    elapsed = time.time() - t_start

    # Summary
    summary = {}
    for label, key in [("SuSiE", "susie"), ("SuSiE-inf", "susie_inf"),
                        ("FINEMAP-inf", "finemap_inf"), ("HBP", "hbp")]:
        ranks = [r[f"{key}_rank"] for r in results if r[f"{key}_rank"] is not None]
        pips = [r[f"{key}_pip"] for r in results if r[f"{key}_pip"] is not None]
        times = [r[f"{key}_time"] for r in results if r[f"{key}_time"] is not None]
        summary[label] = {
            "n": len(ranks),
            "rank_1": int(sum(1 for r in ranks if r == 1)),
            "mean_rank": float(np.mean(ranks)) if ranks else None,
            "median_rank": float(np.median(ranks)) if ranks else None,
            "mean_pip": float(np.mean(pips)) if pips else None,
            "mean_time_s": float(np.mean(times)) if times else None,
        }

    print("\n=== Cui et al. 2024 infinitesimal-effects benchmark ===")
    print(f"({n_reps} reps, 1KG chr22, h²={h2}, weak signal)")
    print(json.dumps(summary, indent=2))
    print(f"Total: {elapsed:.1f}s")

    (OUT / "inf_methods.json").write_text(json.dumps(
        {"summary": summary, "results": results,
         "params": {"n_reps": n_reps, "h2": h2}}, indent=2))
    with (OUT / "inf_methods.tsv").open("w") as fh:
        if results:
            cols = list(results[0].keys())
            fh.write("\t".join(cols) + "\n")
            for r in results:
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"\nWrote {OUT}/inf_methods.{{json,tsv}}")


if __name__ == "__main__":
    main()
