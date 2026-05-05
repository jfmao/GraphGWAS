"""v0.1.5 9-method head-to-head benchmark on chr22 single-causal F1 simulations.

Extends the v0.1.4 6-method benchmark (tests/benchmark_inf_methods.py) by
adding GAFM-MX, HBP-MX, and ENS. Reports per-method:
  - Rank-1 rate (causal variant is the single highest-PIP variant)
  - Top-5 rate
  - Mean PIP at the causal variant
  - Mean wall-clock time per locus

All methods receive the same per-variant z-scores and signed-correlation LD
matrix, computed once from the raw 1000G chr22 BGEN. Methods that take R²
(GAFM, HBP, GAFM-MX, HBP-MX, ENS) get the squared matrix.

Output: results/benchmark_v2/v15_chr22/v15_chr22.{json,tsv}

Default config: 30 reps, h²=0.05, single-causal F1 simulation (Wu et al. 2026).
Per-rep wall-clock: ~30-60 sec depending on locus density. Total: ~30 min.

Usage:
  python tests/benchmark_v15_chr22.py                     # 30 reps, h2=0.05
  N_REPS=10 H2=0.10 python tests/benchmark_v15_chr22.py   # custom
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

import susieinf  # noqa: E402
import finemapinf  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402

from graphgwas.bgen_reader import BgenReader  # noqa: E402
from graphgwas.finemapping_v2 import (  # noqa: E402
    apply_mixture_posterior,
    credible_set_from_pips,
    deflate_z_for_lambda_gc,
    ensemble_from_sumstats,
    gafm_mx_from_sumstats,
    hbp_finemap_from_sumstats,
    hbp_mx_from_sumstats,
    l1_finemap_from_sumstats,
)

OUT = Path("/mnt/data/GraphGWAS/results/benchmark_v2/v15_chr22")
OUT.mkdir(parents=True, exist_ok=True)

N_REPS = int(os.environ.get("N_REPS", "30"))
H2 = float(os.environ.get("H2", "0.05"))
WINDOW_BP = int(os.environ.get("WINDOW_BP", "50000"))
L_EFFECTS = 5
COVERAGE = 0.95


def _simulate_locus(reader: BgenReader, rng, h2: float, window_bp: int = WINDOW_BP):
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


def _summary_inputs(X: np.ndarray, y: np.ndarray):
    """Returns (z, R, R_sq, n_samples) — z = per-variant z-score, R = signed corr,
    R_sq = squared corr, n = sample size."""
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
    # Marginal z = sqrt(n) * cor(X_j, y)
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


def _rank_pip_at_causal(pips, causal_idx):
    pips = np.asarray(pips, dtype=float)
    pip_c = float(pips[causal_idx])
    rank = int((pips >= pip_c).sum())
    return rank, pip_c


def _run_susie_rss(z, R, n, causal_idx):
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
    t0 = time.time()
    proc = subprocess.run(["R", "--no-save", "--no-restore", "-e", r_script],
                           capture_output=True, text=True, timeout=300)
    dt = time.time() - t0
    if not Path(f"{tmp}/pips.tsv").exists():
        return {"rank": None, "pip": 0, "time": dt, "error": (proc.stderr or "")[-200:]}
    pips = pd.read_csv(f"{tmp}/pips.tsv", sep="\t")["pip"].values
    rank, pip = _rank_pip_at_causal(pips, causal_idx)
    return {"rank": rank, "pip": pip, "time": dt}


def _run_susie_inf(z, R, n, causal_idx):
    t0 = time.time()
    try:
        fit = susieinf.susie(z=z, meansq=1.0, n=n, L=L_EFFECTS, LD=R,
                              method="MLE", verbose=False, maxiter=50)
        raw_pip = fit["PIP"] if isinstance(fit, dict) else fit.PIP
        pips = (1 - np.prod(1 - np.nan_to_num(raw_pip, nan=0.0), axis=1)
                if raw_pip.ndim == 2 else np.nan_to_num(raw_pip, nan=0.0))
        pips = np.clip(pips, 0, 1)
    except Exception as e:
        return {"rank": None, "pip": 0, "time": time.time() - t0, "error": str(e)}
    rank, pip = _rank_pip_at_causal(pips, causal_idx)
    return {"rank": rank, "pip": pip, "time": time.time() - t0}


def _run_finemap_inf(z, R, n, causal_idx):
    t0 = time.time()
    try:
        fit = finemapinf.finemap(z=z, meansq=1.0, n=n, L=L_EFFECTS, LD=R,
                                  verbose=0, sched_sss=[50, 50, 500])
        raw_pip = fit["PIP"] if isinstance(fit, dict) else fit.PIP
        pips = (1 - np.prod(1 - np.nan_to_num(raw_pip, nan=0.0), axis=1)
                if raw_pip.ndim == 2 else np.nan_to_num(raw_pip, nan=0.0))
        pips = np.clip(pips, 0, 1)
    except Exception as e:
        return {"rank": None, "pip": 0, "time": time.time() - t0, "error": str(e)}
    rank, pip = _rank_pip_at_causal(pips, causal_idx)
    return {"rank": rank, "pip": pip, "time": time.time() - t0}


def _run_gafm(variants, z, R_sq, causal_idx):
    t0 = time.time()
    out = l1_finemap_from_sumstats(
        variants, z, R_sq, alpha=0.7, credible_set_coverage=COVERAGE,
        chr_name="chr22",
    )
    pips = np.array([c.pip for c in out])
    # the output is sorted by pip; reorder to original variant order
    by_id = {c.variant_id: c.pip for c in out}
    pips = np.array([by_id[v["variantId"]] for v in variants])
    rank, pip = _rank_pip_at_causal(pips, causal_idx)
    return {"rank": rank, "pip": pip, "time": time.time() - t0}


def _run_hbp(variants, z, R_sq, causal_idx):
    t0 = time.time()
    out = hbp_finemap_from_sumstats(
        variants, z, R_sq, graph_cache={},
        credible_set_coverage=COVERAGE, chr_name="chr22",
    )
    by_id = {c.variant_id: c.pip for c in out}
    pips = np.array([by_id[v["variantId"]] for v in variants])
    rank, pip = _rank_pip_at_causal(pips, causal_idx)
    return {"rank": rank, "pip": pip, "time": time.time() - t0}


def _run_gafm_mx(variants, z, R_sq, n_samples, causal_idx):
    t0 = time.time()
    out = gafm_mx_from_sumstats(
        variants, z, R_sq, n_samples=n_samples, alpha=0.7,
        credible_set_coverage=COVERAGE, chr_name="chr22",
    )
    by_id = {c.variant_id: c.pip for c in out}
    pips = np.array([by_id[v["variantId"]] for v in variants])
    rank, pip = _rank_pip_at_causal(pips, causal_idx)
    return {"rank": rank, "pip": pip, "time": time.time() - t0}


def _run_hbp_mx(variants, z, R_sq, n_samples, causal_idx):
    t0 = time.time()
    out = hbp_mx_from_sumstats(
        variants, z, R_sq, n_samples=n_samples,
        credible_set_coverage=COVERAGE, chr_name="chr22",
    )
    by_id = {c.variant_id: c.pip for c in out}
    pips = np.array([by_id[v["variantId"]] for v in variants])
    rank, pip = _rank_pip_at_causal(pips, causal_idx)
    return {"rank": rank, "pip": pip, "time": time.time() - t0}


def _run_ens(variants, z, R_sq, n_samples, causal_idx):
    t0 = time.time()
    out = ensemble_from_sumstats(
        variants, z, R_sq, n_samples=n_samples, alpha=0.7,
        credible_set_coverage=COVERAGE, chr_name="chr22",
    )
    by_id = {c.variant_id: c.pip for c in out}
    pips = np.array([by_id[v["variantId"]] for v in variants])
    rank, pip = _rank_pip_at_causal(pips, causal_idx)
    return {"rank": rank, "pip": pip, "time": time.time() - t0}


METHODS = [
    ("GAFM",        "gafm",        _run_gafm),
    ("HBP",         "hbp",         _run_hbp),
    ("GAFM-MX",     "gafm_mx",     _run_gafm_mx),
    ("HBP-MX",      "hbp_mx",      _run_hbp_mx),
    ("ENS",         "ens",         _run_ens),
    ("SuSiE",       "susie",       _run_susie_rss),
    ("SuSiE-inf",   "susie_inf",   _run_susie_inf),
    ("FINEMAP-inf", "finemap_inf", _run_finemap_inf),
]


def main():
    reader = BgenReader("/mnt/data/GraphGWAS/tests/data/human/1kGP_bgen")
    print(f"v0.1.5 9-method chr22 benchmark — N_REPS={N_REPS} H2={H2} WINDOW_BP={WINDOW_BP}",
          flush=True)
    results = []
    t_start = time.time()
    for i in range(N_REPS):
        rng = np.random.default_rng(4000 + i)
        sim = _simulate_locus(reader, rng, H2)
        if sim is None:
            continue
        vdf, X, y, causal_idx = sim
        n_var, n_samples = X.shape[1], X.shape[0]
        causal_pos = int(vdf.iloc[causal_idx]["pos"])
        z, R, R_sq, n = _summary_inputs(X, y)
        variants = _build_variants(vdf, X)

        row = {"rep": i, "n_var": n_var, "n_samples": n_samples, "causal_pos": causal_pos}
        for label, key, runner in METHODS:
            if key == "susie":
                r = runner(z, R, n, causal_idx)
            elif key in ("susie_inf", "finemap_inf"):
                r = runner(z, R, n, causal_idx)
            elif key in ("gafm_mx", "hbp_mx", "ens"):
                r = runner(variants, z, R_sq, n_samples, causal_idx)
            else:
                r = runner(variants, z, R_sq, causal_idx)
            row[f"{key}_rank"] = r["rank"]
            row[f"{key}_pip"] = r["pip"]
            row[f"{key}_time"] = r["time"]
        results.append(row)
        ranks_str = " ".join(f"{label}#{row[f'{key}_rank']}" for label, key, _ in METHODS)
        print(f"  rep {i:2d} (nv={n_var:4d}): {ranks_str}", flush=True)

    elapsed = time.time() - t_start

    summary = {}
    for label, key, _ in METHODS:
        ranks = [r[f"{key}_rank"] for r in results if r[f"{key}_rank"] is not None]
        pips = [r[f"{key}_pip"] for r in results if r[f"{key}_pip"] is not None]
        times = [r[f"{key}_time"] for r in results if r[f"{key}_time"] is not None]
        summary[label] = {
            "n": len(ranks),
            "rank_1": int(sum(1 for r in ranks if r == 1)),
            "rank_le5": int(sum(1 for r in ranks if r <= 5)),
            "mean_rank": float(np.mean(ranks)) if ranks else None,
            "median_rank": float(np.median(ranks)) if ranks else None,
            "mean_pip": float(np.mean(pips)) if pips else None,
            "median_time_s": float(np.median(times)) if times else None,
            "mean_time_s": float(np.mean(times)) if times else None,
        }

    out = {
        "config": {"n_reps": N_REPS, "h2": H2, "window_bp": WINDOW_BP, "L": L_EFFECTS},
        "elapsed_s": elapsed,
        "summary": summary,
        "per_rep": results,
    }
    h2_tag = f"h{int(round(H2 * 100)):02d}"
    json.dump(out, open(OUT / f"v15_chr22_{h2_tag}.json", "w"), indent=2)
    pd.DataFrame(results).to_csv(OUT / f"v15_chr22_{h2_tag}.tsv", sep="\t", index=False)
    print(f"\nelapsed: {elapsed/60:.1f} min")
    print(f"\n{'method':<14} {'rank-1':>8} {'rank≤5':>8} {'mean rank':>10} "
          f"{'mean PIP':>10} {'median t (s)':>14}")
    for label, _, _ in METHODS:
        s = summary[label]
        print(f"{label:<14} {s['rank_1']:>4d}/{s['n']:<3d} "
              f"{s['rank_le5']:>4d}/{s['n']:<3d} {s['mean_rank']:>10.2f} "
              f"{s['mean_pip']:>10.3f} {s['median_time_s']:>14.4f}")


if __name__ == "__main__":
    main()
