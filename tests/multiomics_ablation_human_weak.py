"""Human chr22 ablation — WEAK-signal regime (λ=3.5).

Same scaffold as multiomics_ablation_human.py, but with λ=3.5 instead
of λ=6.0. In this regime LD deconvolution cannot cleanly resolve the
causal on its own, so the multi-omics prior has room to influence
ranking. This matches the regime where BG1/TGW6 live in rice.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/mnt/data/GraphGWAS/tests")
from multiomics_ablation_test import run_ablation

VCF = "/mnt/data/GraphPop/data/raw/1000g/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz"
CACHE = Path("/mnt/data/GraphGWAS/data/annotations/human_graph_cache_chr22.json")
OUT = Path("/mnt/data/GraphGWAS/results/ablation_human_chr22_weak")
OUT.mkdir(parents=True, exist_ok=True)

LOCI = [
    ("APOL1",  "chr22", 36661906, 500_000),
    ("TIMP3",  "chr22", 32858635, 500_000),
    ("BCR",    "chr22", 23523082, 500_000),
    ("PANX2",  "chr22", 50608284, 500_000),
    ("SNRPD3", "chr22", 24283168, 500_000),
]

LAMBDA = 3.5   # WEAK
SIGMA = 1.0


def _parse_gt(g):
    if g in ("0|0","0/0"): return 0.0
    if g in ("1|0","0|1","1/0","0/1"): return 1.0
    if g in ("1|1","1/1"): return 2.0
    return np.nan


def load_locus_dosage(chrom, start, end):
    q = subprocess.run(
        ["bcftools","query","-r", f"{chrom}:{start}-{end}",
         "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", VCF],
        capture_output=True, text=True, timeout=600,
    )
    rows, mats = [], []
    for line in q.stdout.strip().split("\n"):
        if not line: continue
        f = line.split("\t")
        if len(f) < 5: continue
        d = np.array([_parse_gt(g) for g in f[4:]], dtype=np.float32)
        mu = np.nanmean(d) if np.isfinite(d).any() else 0.0
        d = np.where(np.isnan(d), mu, d)
        af = mu / 2.0 if 0 < mu < 2 else 0.0
        if af < 0.01 or af > 0.99: continue
        if "," in f[3]: continue
        c = f[0] if f[0].startswith("chr") else f"chr{f[0]}"
        rows.append((c, int(f[1]), f[2], f[3], af))
        mats.append(d)
    if not rows: return pd.DataFrame(), np.zeros((0,0))
    return pd.DataFrame(rows, columns=["chr","pos","ref","alt","af"]), np.vstack(mats)


def main():
    print(f"Loading graph cache ...", flush=True)
    full_cache = json.load(open(CACHE))
    all_results = []

    for gene_name, chrom, lead_pos, half_win in LOCI:
        lo, hi = lead_pos - half_win, lead_pos + half_win
        print(f"\n=== {gene_name} @ {chrom}:{lead_pos:,} (λ={LAMBDA}) ===", flush=True)
        var_df, D = load_locus_dosage(chrom, lo, hi)
        if len(var_df) < 20:
            print(f"  skip — {len(var_df)} variants"); continue
        mu = D.mean(axis=1, keepdims=True)
        std = np.where(D.std(axis=1, keepdims=True) < 1e-10, 1.0,
                       D.std(axis=1, keepdims=True))
        Z = (D - mu) / std
        R = (Z @ Z.T) / Z.shape[1]
        R_sq = np.clip(R**2, 0, 1).astype(np.float64)
        np.fill_diagonal(R_sq, 1.0)

        vid = (var_df["chr"] + ":" + var_df["pos"].astype(str) + ":"
               + var_df["ref"] + ":" + var_df["alt"])
        var_df["variant_id"] = vid.values

        causal_idx, best_q = None, -1
        for i, v in enumerate(var_df["variant_id"].values):
            info = full_cache.get(v)
            if not info: continue
            q = len(info.get("genes",[])) + len(info.get("pathways",[])) \
                + min(len(info.get("ppi",[])), 10)
            if q > best_q:
                best_q = q; causal_idx = i
        if causal_idx is None:
            print("  skip — no cached causal"); continue
        causal_id = var_df["variant_id"].iloc[causal_idx]
        print(f"  causal: {causal_id} quality={best_q}")

        np.random.seed(42)
        z = R[causal_idx, :] * LAMBDA + np.random.normal(0, SIGMA, len(var_df))
        z[causal_idx] = LAMBDA

        variants = [
            {"variantId": v["variant_id"], "chr": v["chr"], "pos": int(v["pos"]),
             "ref": v["ref"], "alt": v["alt"], "af_total": float(v["af"])}
            for _, v in var_df.iterrows()
        ]
        window_cache = {k: full_cache[k] for k in var_df["variant_id"].values
                         if k in full_cache}
        print(f"  locus: {len(variants)} variants, {len(window_cache)} annotated",
              flush=True)

        df = run_ablation(
            species="human_weak", locus_name=gene_name,
            variants=variants, z=z, R_sq=R_sq,
            full_cache=window_cache, causal_variant_id=causal_id,
            fractions=(0.1, 0.25, 0.5, 0.75, 1.0),
            n_seeds=10, output_dir=OUT,
        )
        all_results.append(df)

    if not all_results: return
    combined = pd.concat(all_results, ignore_index=True)
    out_tsv = OUT / "ablation_human_weak_combined.tsv"
    combined.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nCombined: {out_tsv}")

    agg = (combined[combined["method"]=="HBP"]
           .groupby(["locus","fraction"])
           .agg(cs_mean=("cs_size","mean"), cs_std=("cs_size","std"),
                pip_mean=("causal_pip","mean"), pip_std=("causal_pip","std"),
                rank_mean=("causal_rank","mean")).reset_index())
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for locus in agg["locus"].unique():
        s = agg[agg["locus"]==locus]
        axes[0].errorbar(s["fraction"]*100, s["cs_mean"], yerr=s["cs_std"],
                         marker="o", label=locus, capsize=3)
        axes[1].errorbar(s["fraction"]*100, s["pip_mean"], yerr=s["pip_std"],
                         marker="o", label=locus, capsize=3)
        axes[2].plot(s["fraction"]*100, s["rank_mean"], marker="o", label=locus)
    axes[0].set_ylabel("95% CS size"); axes[0].set_yscale("log")
    axes[0].set_title(f"CS — human weak-signal (λ={LAMBDA})")
    axes[1].set_ylabel("Causal PIP"); axes[1].set_title("Causal PIP")
    axes[2].set_ylabel("Causal rank"); axes[2].set_yscale("log")
    axes[2].set_title("Causal rank")
    for a in axes:
        a.set_xlabel("Multi-omics edges retained (%)")
        a.grid(True, alpha=0.3)
    axes[0].legend(loc="best", fontsize=7)
    plt.tight_layout()
    plt.savefig(OUT / "ablation_human_weak.png", dpi=150, bbox_inches="tight")
    plt.savefig(OUT / "ablation_human_weak.pdf", bbox_inches="tight")
    print(f"Wrote {OUT}/ablation_human_weak.png")


if __name__ == "__main__":
    main()
