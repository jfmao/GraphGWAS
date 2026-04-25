"""Yeast multi-omics ablation on the 5 best-resolved loci.

Picks the 5 loci with smallest L1 CS from yeast_finemap_summary.tsv,
then re-runs HBP at 5 cache-coverage fractions × 10 seeds.
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

CACHE = Path("/mnt/data/GraphGWAS/data/yeast/yeast_graph_cache.json")
SUMMARY = Path("/mnt/data/GraphGWAS/results/yeast_finemap/yeast_finemap_summary.tsv")
LD_VCF = "/mnt/data/GraphGWAS/tests/data/yeast/1011_import_ready.vcf.gz"
OUT_DIR = Path("/mnt/data/GraphGWAS/results/yeast_finemap/ablation")
OUT_DIR.mkdir(parents=True, exist_ok=True)
GWAS_DIR = Path("/mnt/data/GraphGWAS/results/grammar_corrected")
WINDOW_BP = 15_000


def _parse_gt(g):
    if g in ("0|0","0/0"): return 0.0
    if g in ("1|0","0|1","1/0","0/1"): return 1.0
    if g in ("1|1","1/1"): return 2.0
    return np.nan


def load_locus(chrom, lo, hi):
    q = subprocess.run(
        ["bcftools","query","-r", f"{chrom}:{lo}-{hi}",
         "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", LD_VCF],
        capture_output=True, text=True, timeout=120,
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
        rows.append((f[0], int(f[1]), f[2], f[3], af))
        mats.append(d)
    if not rows:
        return pd.DataFrame(), np.zeros((0, 0))
    return pd.DataFrame(rows, columns=["chr","pos","ref","alt","af"]), np.vstack(mats)


def main():
    print("Loading yeast graph cache ...", flush=True)
    full_cache = json.load(CACHE.open())
    print(f"  {len(full_cache):,} variants in cache")

    df = pd.read_csv(SUMMARY, sep="\t")
    # Pick top 5 best-resolved loci (smallest L1 CS, with PIP > 0.3 for meaningful prior effect)
    df_sorted = df.sort_values(["l1_cs_size", "lead_p"]).head(5)
    print(f"\nSelected 5 loci for ablation:")
    print(df_sorted[["trait","lead_chr","lead_pos","lead_p","l1_cs_size","l1_top_pip"]].to_string(index=False))

    all_results = []
    for _, row in df_sorted.iterrows():
        trait = row["trait"]; chrom = row["lead_chr"]; pos = int(row["lead_pos"])
        causal_id = row["l1_top_variant"]
        lo, hi = pos - WINDOW_BP, pos + WINDOW_BP

        # Load sumstats
        ss_file = GWAS_DIR / f"gwas_{trait}_grammar.tsv"
        ss = pd.read_csv(ss_file, sep="\t", low_memory=False)
        ss = ss.dropna(subset=["p_value","beta","se"])
        ss = ss[(ss["chr"]==chrom) & (ss["pos"].between(lo, hi))].copy()
        ss["z"] = ss["beta"] / ss["se"]
        ss["key"] = ss["variantId"]

        # Load LD
        var_df, D = load_locus(chrom, lo, hi)
        if len(var_df) < 10:
            print(f"  [{trait} {chrom}:{pos}] insufficient LD")
            continue
        mu = D.mean(axis=1, keepdims=True)
        std = np.where(D.std(axis=1, keepdims=True) < 1e-10, 1.0,
                       D.std(axis=1, keepdims=True))
        Z = (D - mu) / std
        R_sq = np.clip(((Z @ Z.T) / Z.shape[1])**2, 0, 1).astype(np.float64)
        np.fill_diagonal(R_sq, 1.0)
        vid = (var_df["chr"] + ":" + var_df["pos"].astype(str) + ":"
               + var_df["ref"] + ":" + var_df["alt"])
        var_df["variant_id"] = vid.values
        z_map = dict(zip(ss["key"].values, ss["z"].values))
        z_flip = {f"{r['chr']}:{r['pos']}:{r['alt']}:{r['ref']}": -r["z"]
                  for _, r in ss.iterrows()}
        z = np.full(len(var_df), np.nan)
        for i, k in enumerate(var_df["variant_id"].values):
            if k in z_map: z[i] = z_map[k]
            elif k in z_flip: z[i] = z_flip[k]
        keep = np.isfinite(z)
        if causal_id not in var_df["variant_id"].values:
            print(f"  [{trait} {chrom}:{pos}] causal {causal_id} not in LD ref"); continue
        if keep.sum() < 10: continue
        idx = np.where(keep)[0]
        var_df = var_df.loc[keep].reset_index(drop=True)
        z = z[keep]; R_sq = R_sq[np.ix_(idx, idx)]
        variants = [
            {"variantId": v["variant_id"], "chr": v["chr"], "pos": int(v["pos"]),
             "ref": v["ref"], "alt": v["alt"], "af_total": float(v["af"])}
            for _, v in var_df.iterrows()
        ]
        window_cache = {k: full_cache[k] for k in var_df["variant_id"].values
                        if k in full_cache}
        print(f"\n=== {trait} {chrom}:{pos:,} (causal {causal_id}) ===", flush=True)
        print(f"  variants: {len(variants)}, cached: {len(window_cache)} "
              f"({100*len(window_cache)/len(variants):.0f}%)")
        out = run_ablation(
            species="yeast", locus_name=f"{trait}_{chrom}_{pos}",
            variants=variants, z=z, R_sq=R_sq,
            full_cache=window_cache, causal_variant_id=causal_id,
            fractions=(0.1, 0.25, 0.5, 0.75, 1.0),
            n_seeds=10, output_dir=OUT_DIR,
        )
        all_results.append(out)

    if not all_results: return
    combined = pd.concat(all_results, ignore_index=True)
    combined.to_csv(OUT_DIR.parent / "ablation_yeast_combined.tsv", sep="\t",
                    index=False)
    print(f"\nWrote ablation_yeast_combined.tsv ({len(combined)} rows)")

    # Plot
    agg = (combined[combined["method"]=="HBP"]
           .groupby(["locus","fraction"])
           .agg(rank_mean=("causal_rank","mean"),
                rank_std=("causal_rank","std"),
                pip_mean=("causal_pip","mean"))
           .reset_index())
    fig, ax = plt.subplots(figsize=(7,4.5))
    for locus in agg["locus"].unique():
        s = agg[agg["locus"]==locus]
        baseline_pip = s[s["fraction"]==1.0]["pip_mean"].iloc[0]
        ax.plot(s["fraction"]*100, s["rank_mean"], marker="o",
                label=f"{locus} (PIP={baseline_pip:.2f})", linewidth=2)
    ax.set_xlabel("Multi-omics edges retained (%)")
    ax.set_ylabel("Causal rank by HBP PIP (lower is better)")
    ax.set_yscale("log")
    ax.set_title("Yeast — multi-omics coverage ablation\n(5 best-resolved loci)")
    ax.grid(True, alpha=0.3); ax.legend(loc="best", fontsize=8)
    plt.tight_layout()
    plt.savefig(OUT_DIR.parent / "ablation_yeast.png", dpi=150, bbox_inches="tight")
    plt.savefig(OUT_DIR.parent / "ablation_yeast.pdf", bbox_inches="tight")
    print(f"Wrote ablation_yeast.png/pdf")


if __name__ == "__main__":
    main()
