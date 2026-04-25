"""Human chr22 multi-omics ablation on the 5 best-resolved Pan-UKB loci."""
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
from pyliftover import LiftOver

sys.path.insert(0, "/mnt/data/GraphGWAS/tests")
from multiomics_ablation_test import run_ablation
sys.path.insert(0, "/mnt/data/GraphGWAS/src/python")
from graphgwas.panukb import fetch_sumstats_locus

CACHE = Path("/mnt/data/GraphGWAS/data/annotations/human_graph_cache_chr22.json")
LD_VCF = "/mnt/data/GraphPop/data/raw/1000g/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz"
SUMMARY = Path("/mnt/data/GraphGWAS/results/human_finemap/human_finemap_summary.tsv")
OUT = Path("/mnt/data/GraphGWAS/results/human_finemap/ablation")
OUT.mkdir(parents=True, exist_ok=True)
WINDOW_BP = 500_000

PHENO_MAP = {
    "BMI":   {"phenocode":"21001","trait_type":"continuous","modifier":"irnt"},
    "Height":{"phenocode":"50",   "trait_type":"continuous","modifier":"irnt"},
    "LDL":   {"phenocode":"30780","trait_type":"biomarkers","modifier":"irnt"},
    "TG":    {"phenocode":"30870","trait_type":"biomarkers","modifier":"irnt"},
}
_LIFT = LiftOver("hg19", "hg38")


def _b37_to_b38(pos):
    res = _LIFT.convert_coordinate("chr22", pos)
    return int(res[0][1]) if res else None


def _parse_gt(g):
    if g in ("0|0","0/0"): return 0.0
    if g in ("1|0","0|1","1/0","0/1"): return 1.0
    if g in ("1|1","1/1"): return 2.0
    return np.nan


def load_locus(start, end):
    q = subprocess.run(
        ["bcftools","query","-r", f"chr22:{start}-{end}",
         "-f","%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", LD_VCF],
        capture_output=True, text=True, timeout=300,
    )
    rows, mats = [], []
    for line in q.stdout.strip().split("\n"):
        if not line: continue
        f = line.split("\t")
        if len(f) < 5: continue
        d = np.array([_parse_gt(g) for g in f[4:]], dtype=np.float32)
        mu = np.nanmean(d) if np.isfinite(d).any() else 0.0
        d = np.where(np.isnan(d), mu, d)
        af = mu/2.0 if 0 < mu < 2 else 0.0
        if af < 0.01 or af > 0.99: continue
        if "," in f[3]: continue
        c = f[0] if f[0].startswith("chr") else f"chr{f[0]}"
        rows.append((c, int(f[1]), f[2], f[3], af))
        mats.append(d)
    if not rows:
        return pd.DataFrame(), np.zeros((0,0))
    return pd.DataFrame(rows, columns=["chr","pos","ref","alt","af"]), np.vstack(mats)


def main():
    print("Loading human chr22 cache + finemap summary ...", flush=True)
    full_cache = json.load(CACHE.open())
    df = pd.read_csv(SUMMARY, sep="\t")
    # Pick 5 best-resolved (smallest L1 CS, prefer with PIP > 0.05)
    df_sorted = df.sort_values(["l1_cs_size","lead_p"]).head(5)
    print("\n5 ablation loci:")
    print(df_sorted[["trait","lead_chr","lead_pos","lead_p","l1_cs_size","l1_top_pip"]].to_string(index=False))

    all_results = []
    # Cache fetched sumstats per (trait, window) to avoid double tabix
    fetched = {}

    for _, row in df_sorted.iterrows():
        trait = row["trait"]; pos = int(row["lead_pos"])
        causal_id = row["l1_top_variant"]
        # Sumstats: tabix the locus window in b37 coords.
        # Need to convert pos b38 → b37 for fetch
        # Reverse liftover via chain hg38→hg19
        rev = LiftOver("hg38","hg19")
        b37 = rev.convert_coordinate("chr22", pos)
        if not b37:
            print(f"  [{trait} chr22:{pos}] reverse liftover failed"); continue
        b37_pos = int(b37[0][1])
        meta = PHENO_MAP[trait]
        key = (trait, b37_pos)
        if key not in fetched:
            print(f"\n  Fetching {trait} chr22 b37:{b37_pos:,} ±{WINDOW_BP//2}", flush=True)
            res = fetch_sumstats_locus(
                phenocode=meta["phenocode"], chr="22",
                start=max(1, b37_pos - WINDOW_BP), end=b37_pos + WINDOW_BP,
                ancestries=["EUR"],
                trait_type=meta["trait_type"], modifier=meta["modifier"],
            )
            ss = res["EUR"].variants.copy()
            ss["pos_b37"] = ss["pos"].astype(int)
            ss["pos_b38"] = ss["pos_b37"].apply(_b37_to_b38)
            ss = ss.dropna(subset=["pos_b38"])
            ss["pos_b38"] = ss["pos_b38"].astype(int)
            ss["pos"] = ss["pos_b38"]
            ss["chr"] = "chr22"
            ss["z"] = ss["beta"] / ss["se"]
            ss["key"] = "chr22:" + ss["pos"].astype(str) + ":" + ss["ref"] + ":" + ss["alt"]
            fetched[key] = ss
        ss = fetched[key]

        # Locus window in b38
        lo = pos - WINDOW_BP//2; hi = pos + WINDOW_BP//2
        sub = ss[ss["pos"].between(lo, hi)].copy()
        if len(sub) < 20:
            print(f"  [{trait} chr22:{pos}] insufficient sumstats"); continue

        var_df, D = load_locus(lo, hi)
        if len(var_df) < 20:
            print(f"  [{trait} chr22:{pos}] insufficient LD"); continue
        mu = D.mean(axis=1, keepdims=True)
        std = np.where(D.std(axis=1, keepdims=True) < 1e-10, 1.0,
                       D.std(axis=1, keepdims=True))
        Z = (D - mu) / std
        R_sq = np.clip(((Z @ Z.T) / Z.shape[1])**2, 0, 1).astype(np.float64)
        np.fill_diagonal(R_sq, 1.0)
        vid = (var_df["chr"] + ":" + var_df["pos"].astype(str) + ":"
               + var_df["ref"] + ":" + var_df["alt"])
        var_df["variant_id"] = vid.values
        z_map = dict(zip(sub["key"].values, sub["z"].values))
        z_flip = {f"chr22:{r['pos']}:{r['alt']}:{r['ref']}": -r["z"]
                  for _, r in sub.iterrows()}
        z = np.full(len(var_df), np.nan)
        for i, k in enumerate(var_df["variant_id"].values):
            if k in z_map: z[i] = z_map[k]
            elif k in z_flip: z[i] = z_flip[k]
        keep = np.isfinite(z)
        if causal_id not in var_df["variant_id"].values:
            print(f"  [{trait} chr22:{pos}] causal {causal_id} not in LD ref"); continue
        if keep.sum() < 20: continue
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
        print(f"\n=== {trait} chr22:{pos:,} (causal {causal_id}) ===", flush=True)
        print(f"  variants: {len(variants)}, cached: {len(window_cache)} "
              f"({100*len(window_cache)/len(variants):.0f}%)")
        out = run_ablation(
            species="human", locus_name=f"{trait}_chr22_{pos}",
            variants=variants, z=z, R_sq=R_sq,
            full_cache=window_cache, causal_variant_id=causal_id,
            fractions=(0.1, 0.25, 0.5, 0.75, 1.0),
            n_seeds=10, output_dir=OUT,
        )
        all_results.append(out)

    if not all_results: return
    combined = pd.concat(all_results, ignore_index=True)
    combined.to_csv(OUT.parent / "ablation_human_combined.tsv", sep="\t", index=False)
    print(f"\nWrote ablation_human_combined.tsv ({len(combined)} rows)")

    agg = (combined[combined["method"]=="HBP"]
           .groupby(["locus","fraction"])
           .agg(rank_mean=("causal_rank","mean"),
                rank_std=("causal_rank","std"),
                pip_mean=("causal_pip","mean"))
           .reset_index())
    fig, ax = plt.subplots(figsize=(7,4.5))
    for locus in agg["locus"].unique():
        s = agg[agg["locus"]==locus]
        bp = s[s["fraction"]==1.0]["pip_mean"].iloc[0]
        ax.plot(s["fraction"]*100, s["rank_mean"], marker="o",
                label=f"{locus} (PIP={bp:.2f})", linewidth=2)
    ax.set_xlabel("Multi-omics edges retained (%)")
    ax.set_ylabel("Causal rank by HBP PIP (lower is better)")
    ax.set_yscale("log")
    ax.set_title("Human chr22 — multi-omics coverage ablation\n(5 best-resolved Pan-UKB loci)")
    ax.grid(True, alpha=0.3); ax.legend(loc="best", fontsize=8)
    plt.tight_layout()
    plt.savefig(OUT.parent / "ablation_human.png", dpi=150, bbox_inches="tight")
    plt.savefig(OUT.parent / "ablation_human.pdf", bbox_inches="tight")
    print(f"Wrote ablation_human.png/pdf")


if __name__ == "__main__":
    main()
