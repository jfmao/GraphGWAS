"""Rice ablation test driver.

Runs the multi-omics ablation test on rice for each of the 5 simulated
causal loci (where we have ground truth) using the GWAS sumstats +
rice multi-omics graph cache we already built.

Produces one long-format TSV per locus + a combined one, plus a
matplotlib figure showing CS size and causal PIP vs coverage.
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

DATA_DIR = Path("/mnt/data/GraphGWAS/data/rice_3k")
RES = DATA_DIR / "results"
OUT_DIR = RES / "ablation"
OUT_DIR.mkdir(parents=True, exist_ok=True)

GWAS = RES / "rice_gwas.GRAIN_SIZE.glm.linear"
VCF = "/mnt/data/GraphPop/data/raw/3kRG_data/NB_final_snp.vcf.gz"
SIM = pd.read_csv(DATA_DIR / "pheno" / "causal_variants.tsv", sep="\t")

FM_WINDOW = 250_000


def _parse_gt(g):
    if g in ("0|0","0/0"): return 0.0
    if g in ("1|0","0|1","1/0","0/1"): return 1.0
    if g in ("1|1","1/1"): return 2.0
    return np.nan


def load_locus(chrom, lo, hi):
    q = subprocess.run(
        ["bcftools","query","-r",f"{chrom}:{lo}-{hi}",
         "-f","%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", VCF],
        capture_output=True, text=True, timeout=300,
    )
    rows, mats = [], []
    for line in q.stdout.strip().split("\n"):
        f = line.split("\t")
        if len(f) < 5: continue
        gts = f[4:]
        d = np.array([_parse_gt(g) for g in gts], dtype=np.float32)
        mu = np.nanmean(d) if np.isfinite(d).any() else 0.0
        d = np.where(np.isnan(d), mu, d)
        af = mu / 2.0 if 0 < mu < 2 else 0.0
        if af < 0.01 or af > 0.99: continue
        rows.append((f[0], int(f[1]), f[2], f[3], af))
        mats.append(d)
    if not rows:
        return pd.DataFrame(), np.zeros((0,0))
    return pd.DataFrame(rows, columns=["chr","pos","ref","alt","af"]), np.vstack(mats)


def main():
    # Load GWAS sumstats (same normalisation as finemap driver)
    gw = pd.read_csv(GWAS, sep="\t", low_memory=False)
    gw.columns = [c.lstrip("#") for c in gw.columns]
    gw = gw[(gw["TEST"]=="ADD") & (gw["ERRCODE"]==".")].copy()
    gw["CHROM"] = gw["CHROM"].astype(str).apply(
        lambda c: c if c.startswith("Chr") else f"Chr{c}")
    for c in ("BETA","SE","P"):
        gw[c] = pd.to_numeric(gw[c], errors="coerce")
    gw = gw.dropna(subset=["P","BETA","SE"])
    gw["z"] = gw["BETA"] / gw["SE"]

    all_results = []
    # Track which chromosomes we've already loaded (cache is big)
    loaded: dict[str, dict] = {}
    def _get_chr_cache(chrom):
        if chrom not in loaded:
            p = DATA_DIR / "annotations" / f"rice_graph_cache_{chrom}.json"
            print(f"  loading {p.name}...", flush=True)
            loaded[chrom] = json.load(open(p))
            print(f"    {len(loaded[chrom]):,} variants cached", flush=True)
        return loaded[chrom]

    for _, sc in SIM.iterrows():
        chrom = sc["chr"]; cpos = int(sc["pos"])
        causal_id = f"{chrom}:{cpos}:{sc['ref']}:{sc['alt']}"
        lo = cpos - FM_WINDOW; hi = cpos + FM_WINDOW
        print(f"\n=== {sc['gene']} @ {chrom}:{cpos:,} ===", flush=True)
        var_df, D = load_locus(chrom, lo, hi)
        if len(var_df) < 10:
            print(f"  skip — only {len(var_df)} variants")
            continue

        # Build LD + z per variant
        mu = D.mean(axis=1, keepdims=True)
        std = np.where(D.std(axis=1, keepdims=True) < 1e-10, 1.0,
                       D.std(axis=1, keepdims=True))
        Z = (D - mu) / std
        R2 = np.clip(((Z @ Z.T) / Z.shape[1])**2, 0, 1).astype(np.float64)
        np.fill_diagonal(R2, 1.0)
        vid = (var_df["chr"] + ":" + var_df["pos"].astype(str) + ":"
               + var_df["ref"] + ":" + var_df["alt"])
        var_df["variant_id"] = vid.values

        # Map sumstats onto window variants (flipping allele-swapped z)
        gw_sub = gw[(gw["CHROM"]==chrom) & (gw["POS"].between(lo,hi))].copy()
        gw_sub["key"] = (gw_sub["CHROM"].astype(str)+":"+gw_sub["POS"].astype(str)
                         +":"+gw_sub["REF"]+":"+gw_sub["ALT"])
        z_map = dict(zip(gw_sub["key"].values, gw_sub["z"].values))
        z_flip = {f"{r['CHROM']}:{r['POS']}:{r['ALT']}:{r['REF']}": -r["z"]
                  for _, r in gw_sub.iterrows()}
        z = np.full(len(var_df), np.nan)
        for i, k in enumerate(var_df["variant_id"].values):
            if k in z_map: z[i] = z_map[k]
            elif k in z_flip: z[i] = z_flip[k]
        keep = np.isfinite(z)
        if keep.sum() < 10 or causal_id not in var_df["variant_id"].values:
            print(f"  skip — causal {causal_id} not in locus window "
                  f"or <10 sumstats matches")
            continue
        idx = np.where(keep)[0]
        var_df = var_df.loc[keep].reset_index(drop=True)
        z = z[keep]
        R2 = R2[np.ix_(idx, idx)]

        variants = [
            {"variantId": v["variant_id"], "chr": v["chr"], "pos": int(v["pos"]),
             "ref": v["ref"], "alt": v["alt"], "af_total": float(v["af"])}
            for _, v in var_df.iterrows()
        ]

        # Subset the chromosome cache to this locus window only
        full_chr = _get_chr_cache(chrom)
        window_cache = {k: full_chr[k] for k in var_df["variant_id"].values if k in full_chr}
        print(f"  locus window: {len(variants)} variants, "
              f"{len(window_cache)} with annotations "
              f"({100*len(window_cache)/len(variants):.0f}%)", flush=True)

        df = run_ablation(
            species="rice",
            locus_name=sc["gene"].replace("/", "_"),
            variants=variants, z=z, R_sq=R2,
            full_cache=window_cache,
            causal_variant_id=causal_id,
            fractions=(0.1, 0.25, 0.5, 0.75, 1.0),
            n_seeds=10,
            output_dir=OUT_DIR,
        )
        all_results.append(df)

    if not all_results:
        print("\nNo loci produced results.")
        return

    combined = pd.concat(all_results, ignore_index=True)
    out = RES / "ablation_rice_combined.tsv"
    combined.to_csv(out, sep="\t", index=False)
    print(f"\n\nCombined results: {out}")

    # --- Plot: CS size and causal PIP vs coverage ---
    agg = (combined[combined["method"]=="HBP"]
           .groupby(["locus", "fraction"])
           .agg(cs_mean=("cs_size", "mean"),
                cs_std=("cs_size", "std"),
                pip_mean=("causal_pip", "mean"),
                pip_std=("causal_pip", "std"),
                rank_mean=("causal_rank", "mean"))
           .reset_index())
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4.5))
    for locus in agg["locus"].unique():
        sub = agg[agg["locus"]==locus]
        ax1.errorbar(sub["fraction"]*100, sub["cs_mean"], yerr=sub["cs_std"],
                     marker="o", label=locus, capsize=3)
        ax2.errorbar(sub["fraction"]*100, sub["pip_mean"], yerr=sub["pip_std"],
                     marker="o", label=locus, capsize=3)
        ax3.plot(sub["fraction"]*100, sub["rank_mean"], marker="o", label=locus)
    ax1.set_xlabel("Multi-omics edges retained (%)")
    ax1.set_ylabel("95% credible-set size")
    ax1.set_title("CS size — rice (HBP, mean ± SD, 10 seeds)")
    ax1.set_yscale("log")
    ax1.grid(True, alpha=0.3)
    ax2.set_xlabel("Multi-omics edges retained (%)")
    ax2.set_ylabel("Causal variant PIP")
    ax2.set_title("Causal PIP — rice (HBP)")
    ax2.grid(True, alpha=0.3)
    ax3.set_xlabel("Multi-omics edges retained (%)")
    ax3.set_ylabel("Rank of causal variant by PIP")
    ax3.set_title("Causal rank — rice (HBP)")
    ax3.grid(True, alpha=0.3)
    ax3.set_yscale("log")
    ax1.legend(loc="best", fontsize=7)
    plt.tight_layout()
    plt.savefig(RES / "ablation_rice.png", dpi=150, bbox_inches="tight")
    plt.savefig(RES / "ablation_rice.pdf", bbox_inches="tight")
    print(f"Wrote {RES}/ablation_rice.{{png,pdf}}")


if __name__ == "__main__":
    main()
