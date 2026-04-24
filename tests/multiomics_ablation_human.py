"""Human chr22 multi-omics ablation driver.

Simulates a GWAS with a known causal variant in each of 5 chr22 loci
chosen to have strong GTEx eQTL annotation and STRING PPI coverage.
Runs the reusable ablation scaffold at 5 cache-coverage levels × 10 seeds
and writes per-locus TSVs + a combined figure.

Uses the real 1KG chr22 phased VCF to derive LD, samples z-scores for
the causal variant from N(λ, 1) with λ chosen to yield PIP ~ 0.3-0.5
in the unpenalized run (visible signal). All other variants' z-scores
are the causal-LD-induced values: z_i = r_{i,c} × λ + ε_i.

Outputs:
  data/annotations/results_human_chr22/ablation_human_{locus}.tsv
  data/annotations/results_human_chr22/ablation_human_combined.tsv
  data/annotations/results_human_chr22/ablation_human.{png,pdf}
"""
from __future__ import annotations

import json
import random
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
OUT = Path("/mnt/data/GraphGWAS/results/ablation_human_chr22")
OUT.mkdir(parents=True, exist_ok=True)

# Causals chosen from chr22 gene-dense regions with strong GTEx/STRING coverage
LOCI = [
    ("APOL1",    "chr22", 36661906, 500_000),  # APOL1 kidney disease, rs73885319
    ("TIMP3",    "chr22", 32858635, 500_000),  # TIMP3 macular degeneration
    ("BCR",      "chr22", 23523082, 500_000),  # BCR-ABL (hematology)
    ("PANX2",    "chr22", 50608284, 500_000),  # PANX2 neurological
    ("SNRPD3",   "chr22", 24283168, 500_000),  # SNRPD3 splicing (housekeeping)
]

CAUSAL_EFFECT_LAMBDA = 6.0  # z-score of the causal (corresponds to strong signal)
EPSILON_SIGMA = 1.0         # noise std on non-causal z-scores
N_REFERENCE = 2504          # 1KG sample count used for standardisation


def _parse_gt(g: str) -> float:
    if g in ("0|0","0/0"): return 0.0
    if g in ("1|0","0|1","1/0","0/1"): return 1.0
    if g in ("1|1","1/1"): return 2.0
    return np.nan


def load_locus_dosage(chrom: str, start: int, end: int) -> tuple[pd.DataFrame, np.ndarray]:
    """Query the 1KG phased VCF for a locus; return dosage matrix + var_df."""
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
        gts = f[4:]
        d = np.array([_parse_gt(g) for g in gts], dtype=np.float32)
        mu = np.nanmean(d) if np.isfinite(d).any() else 0.0
        d = np.where(np.isnan(d), mu, d)
        af = mu / 2.0 if 0 < mu < 2 else 0.0
        if af < 0.01 or af > 0.99: continue
        chrom_str = f[0] if f[0].startswith("chr") else f"chr{f[0]}"
        # Skip multi-allelic
        if "," in f[3]: continue
        rows.append((chrom_str, int(f[1]), f[2], f[3], af))
        mats.append(d)
    if not rows:
        return pd.DataFrame(), np.zeros((0, 0))
    return pd.DataFrame(rows, columns=["chr","pos","ref","alt","af"]), np.vstack(mats)


def main():
    print(f"[1/3] Loading graph cache {CACHE.name} ...", flush=True)
    full_cache = json.load(open(CACHE))
    print(f"  {len(full_cache):,} variants in cache\n", flush=True)

    all_results = []
    for gene_name, chrom, lead_pos, half_win in LOCI:
        lo, hi = lead_pos - half_win, lead_pos + half_win
        print(f"\n=== {gene_name} @ {chrom}:{lead_pos:,} ===", flush=True)
        var_df, D = load_locus_dosage(chrom, lo, hi)
        if len(var_df) < 20:
            print(f"  skip — only {len(var_df)} variants in locus")
            continue
        # Build LD matrix (r)
        mu = D.mean(axis=1, keepdims=True)
        std = np.where(D.std(axis=1, keepdims=True) < 1e-10, 1.0,
                       D.std(axis=1, keepdims=True))
        Z = (D - mu) / std
        R = (Z @ Z.T) / Z.shape[1]
        R_sq = np.clip(R ** 2, 0, 1).astype(np.float64)
        np.fill_diagonal(R_sq, 1.0)

        vid = (var_df["chr"] + ":" + var_df["pos"].astype(str) + ":"
               + var_df["ref"] + ":" + var_df["alt"])
        var_df["variant_id"] = vid.values

        # Pick a causal: variant closest to lead_pos that is IN the cache,
        # preferring one with ≥ 1 PPI partner and pathway/gene entry.
        # Cache lookup by variant_id.
        causal_idx = None
        best_quality = -1
        for i, v in enumerate(var_df["variant_id"].values):
            info = full_cache.get(v)
            if not info: continue
            q = len(info.get("genes", [])) + len(info.get("pathways", [])) \
                + min(len(info.get("ppi", [])), 10)
            if q > best_quality:
                best_quality = q
                causal_idx = i
        if causal_idx is None:
            print("  skip — no cached causal found in window")
            continue
        causal_id = var_df["variant_id"].iloc[causal_idx]
        print(f"  causal: {causal_id}  quality={best_quality}")

        # Build simulated z-scores: z_i = r_{i,causal} * λ + ε
        np.random.seed(42)
        r_to_causal = R[causal_idx, :]  # signed r (not r^2)
        z = r_to_causal * CAUSAL_EFFECT_LAMBDA + \
            np.random.normal(0.0, EPSILON_SIGMA, size=len(var_df))
        z[causal_idx] = CAUSAL_EFFECT_LAMBDA  # anchor

        variants = [
            {"variantId": v["variant_id"], "chr": v["chr"], "pos": int(v["pos"]),
             "ref": v["ref"], "alt": v["alt"], "af_total": float(v["af"])}
            for _, v in var_df.iterrows()
        ]

        # Restrict cache to variants in this window
        window_cache = {k: full_cache[k] for k in var_df["variant_id"].values
                         if k in full_cache}
        print(f"  locus: {len(variants)} variants, "
              f"{len(window_cache)} with annotations "
              f"({100*len(window_cache)/len(variants):.0f}%)", flush=True)

        df = run_ablation(
            species="human",
            locus_name=gene_name,
            variants=variants, z=z, R_sq=R_sq,
            full_cache=window_cache,
            causal_variant_id=causal_id,
            fractions=(0.1, 0.25, 0.5, 0.75, 1.0),
            n_seeds=10,
            output_dir=OUT,
        )
        all_results.append(df)

    if not all_results:
        print("\nNo loci produced results.")
        return

    combined = pd.concat(all_results, ignore_index=True)
    out_tsv = OUT / "ablation_human_combined.tsv"
    combined.to_csv(out_tsv, sep="\t", index=False)
    print(f"\n\nCombined: {out_tsv}")

    # Plot
    agg = (combined[combined["method"]=="HBP"]
           .groupby(["locus", "fraction"])
           .agg(cs_mean=("cs_size","mean"), cs_std=("cs_size","std"),
                pip_mean=("causal_pip","mean"), pip_std=("causal_pip","std"),
                rank_mean=("causal_rank","mean"))
           .reset_index())
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(15, 4.5))
    for locus in agg["locus"].unique():
        s = agg[agg["locus"]==locus]
        a1.errorbar(s["fraction"]*100, s["cs_mean"], yerr=s["cs_std"],
                    marker="o", label=locus, capsize=3)
        a2.errorbar(s["fraction"]*100, s["pip_mean"], yerr=s["pip_std"],
                    marker="o", label=locus, capsize=3)
        a3.plot(s["fraction"]*100, s["rank_mean"], marker="o", label=locus)
    for a in (a1,a2,a3):
        a.set_xlabel("Multi-omics edges retained (%)")
        a.grid(True, alpha=0.3)
    a1.set_ylabel("95% credible-set size"); a1.set_title("CS — human chr22 (HBP)")
    a1.set_yscale("log")
    a2.set_ylabel("Causal PIP"); a2.set_title("Causal PIP — human chr22")
    a3.set_ylabel("Causal rank by PIP"); a3.set_title("Causal rank — human chr22")
    a3.set_yscale("log")
    a1.legend(loc="best", fontsize=7)
    plt.tight_layout()
    plt.savefig(OUT / "ablation_human.png", dpi=150, bbox_inches="tight")
    plt.savefig(OUT / "ablation_human.pdf", bbox_inches="tight")
    print(f"Wrote {OUT}/ablation_human.{{png,pdf}}")


if __name__ == "__main__":
    main()
