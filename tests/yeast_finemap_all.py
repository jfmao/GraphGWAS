"""Fine-map all 35 yeast grammar-corrected GWAS sumstats files.

For each trait:
  - Cluster GW-sig + suggestive leads at 30 kb window (yeast LD scale)
  - Cap at top-3 GW + top-1 suggestive per trait
  - For each lead:
      - Extract sumstats in ±30 kb window
      - Build LD matrix from 1011_snps_maf05.vcf.gz (LD ref)
      - Run L1 + HBP fine-mapping using yeast multi-omics cache
  - Write per-trait .tsv + combined summary

Parallelism: traits are chunked across N parallel workers.
"""
from __future__ import annotations

import argparse
import glob
import json
import multiprocessing as mp
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd

sys.path.insert(0, "/mnt/data/GraphGWAS/src/python")
from graphgwas.finemapping_v2 import (
    hbp_finemap_from_sumstats,
    l1_finemap_from_sumstats,
)

YEAST_DATA = Path("/mnt/data/GraphGWAS/tests/data/yeast")
RES_DIR = Path("/mnt/data/GraphGWAS/results/grammar_corrected")
CACHE_VERSION = os.environ.get("CACHE_VERSION", "v1")
_cache_name = "yeast_graph_cache.json" if CACHE_VERSION == "v1" else f"yeast_graph_cache_{CACHE_VERSION}.json"
CACHE_PATH = Path(f"/mnt/data/GraphGWAS/data/yeast/{_cache_name}")
OUT_DIR = Path(f"/mnt/data/GraphGWAS/results/yeast_finemap_{CACHE_VERSION}") if CACHE_VERSION != "v1" \
          else Path("/mnt/data/GraphGWAS/results/yeast_finemap")
OUT_DIR.mkdir(parents=True, exist_ok=True)
LD_VCF = YEAST_DATA / "1011_import_ready.vcf.gz"  # has .csi index

WINDOW_BP = 15_000        # ±15 kb (yeast LD decays fast)
GW_THRESHOLD = 5e-8
SUG_THRESHOLD = 1e-5
MAX_GW_PER_TRAIT = 5
MAX_SUG_PER_TRAIT = 2
N_WORKERS = 8


def _parse_gt(g: str) -> float:
    if g in ("0|0","0/0"): return 0.0
    if g in ("1|0","0|1","1/0","0/1"): return 1.0
    if g in ("1|1","1/1"): return 2.0
    return np.nan


def load_locus_dosage(chrom: str, lo: int, hi: int) -> Tuple[pd.DataFrame, np.ndarray]:
    # yeast VCF uses chromosome1..chromosome16 names
    q = subprocess.run(
        ["bcftools","query","-r", f"{chrom}:{lo}-{hi}",
         "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", str(LD_VCF)],
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
        return pd.DataFrame(), np.zeros((0,0))
    return pd.DataFrame(rows, columns=["chr","pos","ref","alt","af"]), np.vstack(mats)


def cluster_leads(df: pd.DataFrame, p_thresh: float, window: int) -> list[dict]:
    d = df[df["p_value"] <= p_thresh].copy()
    if d.empty: return []
    d = d.sort_values(["chr", "p_value"])
    leads = []
    for chrom, sub in d.groupby("chr", sort=False):
        taken = []
        for _, r in sub.iterrows():
            if any(abs(r["pos"] - t) < window for t in taken):
                continue
            taken.append(r["pos"])
            leads.append({
                "chrom": chrom, "pos": int(r["pos"]),
                "p": float(r["p_value"]), "beta": float(r["beta"]),
                "ref": r["ref"], "alt": r["alt"],
                "variant_id": r["variantId"],
            })
    return leads


def finemap_locus(trait: str, lead: dict, sumstats: pd.DataFrame,
                  full_cache: dict) -> dict | None:
    chrom = lead["chrom"]; pos = lead["pos"]
    lo, hi = pos - WINDOW_BP, pos + WINDOW_BP
    sub = sumstats[(sumstats["chr"]==chrom) & (sumstats["pos"].between(lo, hi))].copy()
    if len(sub) < 10:
        return None
    sub["z"] = sub["beta"] / sub["se"]
    sub["key"] = sub["variantId"]

    var_df, D = load_locus_dosage(chrom, lo, hi)
    if len(var_df) < 10:
        return None
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
    z_flip = {f"{r['chr']}:{r['pos']}:{r['alt']}:{r['ref']}": -r["z"]
              for _, r in sub.iterrows()}
    z = np.full(len(var_df), np.nan)
    for i, k in enumerate(var_df["variant_id"].values):
        if k in z_map: z[i] = z_map[k]
        elif k in z_flip: z[i] = z_flip[k]
    keep = np.isfinite(z)
    if keep.sum() < 10:
        return None
    idx = np.where(keep)[0]
    var_df = var_df.loc[keep].reset_index(drop=True)
    z = z[keep]
    R_sq = R_sq[np.ix_(idx, idx)]
    variants = [
        {"variantId": v["variant_id"], "chr": v["chr"], "pos": int(v["pos"]),
         "ref": v["ref"], "alt": v["alt"], "af_total": float(v["af"])}
        for _, v in var_df.iterrows()
    ]
    window_cache = {k: full_cache[k] for k in var_df["variant_id"].values
                    if k in full_cache} if full_cache else {}

    l1 = l1_finemap_from_sumstats(
        variants, z, R_sq, z_func=None, alpha=0.7,
        r2_smooth=0.3, credible_set_coverage=0.95, chr_name=chrom,
        graph_cache=window_cache,
    )
    hbp = hbp_finemap_from_sumstats(
        variants, z, R_sq, graph_cache=window_cache,
        r2_smooth=0.3, credible_set_coverage=0.95, chr_name=chrom,
    )
    l1_cs = [c for c in l1 if c.in_credible_set]
    hbp_cs = [c for c in hbp if c.in_credible_set]
    l1_top = sorted(l1, key=lambda c: -c.pip)[:3]
    hbp_top = sorted(hbp, key=lambda c: -c.pip)[:3]
    return {
        "trait": trait, "lead_chr": chrom, "lead_pos": pos,
        "lead_p": lead["p"], "sig_level": lead["sig_level"],
        "n_variants": len(variants), "cache_coverage": len(window_cache),
        "l1_cs_size": len(l1_cs), "hbp_cs_size": len(hbp_cs),
        "l1_top_variant": l1_top[0].variant_id if l1_top else None,
        "l1_top_pip": l1_top[0].pip if l1_top else None,
        "hbp_top_variant": hbp_top[0].variant_id if hbp_top else None,
        "hbp_top_pip": hbp_top[0].pip if hbp_top else None,
    }


def process_trait(args):
    trait_file, full_cache = args
    trait = Path(trait_file).stem.replace("gwas_","").replace("_grammar","")
    print(f"  [{trait}] loading sumstats...", flush=True)
    df = pd.read_csv(trait_file, sep="\t", low_memory=False)
    df = df.dropna(subset=["p_value","beta","se"])
    df = df[df["se"] > 0]
    n_var = len(df)
    n_gw = int((df["p_value"] <= GW_THRESHOLD).sum())
    n_sug = int((df["p_value"] <= SUG_THRESHOLD).sum())

    leads_gw = cluster_leads(df, GW_THRESHOLD, WINDOW_BP)
    leads_sug = cluster_leads(df, SUG_THRESHOLD, WINDOW_BP)
    for l in leads_gw: l["sig_level"] = "GW"
    for l in leads_sug:
        if any(l2["chrom"]==l["chrom"] and abs(l2["pos"]-l["pos"]) < WINDOW_BP
               for l2 in leads_gw):
            continue
        l["sig_level"] = "suggestive"
    leads_gw = sorted(leads_gw, key=lambda x: x["p"])[:MAX_GW_PER_TRAIT]
    leads_sug = sorted([l for l in leads_sug if l.get("sig_level")=="suggestive"],
                       key=lambda x: x["p"])[:MAX_SUG_PER_TRAIT]
    leads = leads_gw + leads_sug

    print(f"  [{trait}] {n_var:,} vars, {n_gw} GW, {n_sug} sug; "
          f"fine-mapping {len(leads)} leads", flush=True)

    rows = []
    for lead in leads:
        try:
            r = finemap_locus(trait, lead, df, full_cache)
            if r:
                rows.append(r)
        except Exception as e:
            print(f"    [{trait} {lead['chrom']}:{lead['pos']}] error {e}", flush=True)
    out = OUT_DIR / f"yeast_finemap_{trait}.tsv"
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    return {"trait": trait, "n_var": n_var, "n_gw_indep": len(leads_gw),
            "n_sug_indep": len(leads_sug), "n_finemapped": len(rows),
            "min_p": float(df["p_value"].min())}, rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=N_WORKERS)
    parser.add_argument("--limit", type=int, default=0,
                        help="limit number of traits (0=all)")
    args = parser.parse_args()

    print("Loading yeast multi-omics graph cache ...", flush=True)
    t0 = time.time()
    full_cache = json.load(CACHE_PATH.open())
    print(f"  {len(full_cache):,} variants loaded in {time.time()-t0:.0f}s", flush=True)

    trait_files = sorted(glob.glob(str(RES_DIR / "gwas_*_grammar.tsv")))
    if args.limit:
        trait_files = trait_files[:args.limit]
    print(f"\n{len(trait_files)} trait files to process with {args.workers} workers")

    summary_rows, all_loci = [], []
    if args.workers > 1:
        with mp.Pool(args.workers) as pool:
            for s, rows in pool.imap_unordered(
                    process_trait, [(f, full_cache) for f in trait_files]):
                summary_rows.append(s)
                all_loci.extend(rows)
                print(f"  done {s['trait']} ({len(all_loci)} loci so far)",
                      flush=True)
    else:
        for f in trait_files:
            s, rows = process_trait((f, full_cache))
            summary_rows.append(s)
            all_loci.extend(rows)

    summary = pd.DataFrame(summary_rows).sort_values("min_p")
    summary.to_csv(OUT_DIR / "yeast_gwas_summary.tsv", sep="\t", index=False)
    pd.DataFrame(all_loci).to_csv(OUT_DIR / "yeast_finemap_summary.tsv",
                                   sep="\t", index=False)
    print(f"\nWrote {OUT_DIR}/yeast_finemap_summary.tsv  ({len(all_loci)} loci)")
    print(f"Wrote {OUT_DIR}/yeast_gwas_summary.tsv  ({len(summary)} traits)")


if __name__ == "__main__":
    main()
