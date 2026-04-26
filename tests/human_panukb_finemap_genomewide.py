"""Genome-wide human Pan-UKB fine-mapping for 4 phenotypes.

Per phenotype, fetches sumstats per chromosome via tabix-over-HTTPS,
clusters GW-sig leads at 500 kb, fine-maps top-N leads using the
genome-wide GTEx + STRING graph cache.

Parallelism: chromosomes processed in parallel (multiprocessing.Pool)
to fully use the 32-core machine.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from pyliftover import LiftOver

sys.path.insert(0, "/mnt/data/GraphGWAS/src/python")
from graphgwas.finemapping_v2 import (
    hbp_finemap_from_sumstats,
    l1_finemap_from_sumstats,
)
from graphgwas.panukb import fetch_sumstats_locus

ANN = Path("/mnt/data/GraphGWAS/data/annotations")
LD_VCF_DIR = Path("/mnt/data/GraphPop/data/raw/1000g")
# CACHE_VERSION env var: "v1" = GTEx+STRING (default), "v2" = + cCRE+prior_score
CACHE_VERSION = os.environ.get("CACHE_VERSION", "v1")
CACHE_PREFIX = "human_graph_cache_v2_chr" if CACHE_VERSION == "v2" \
               else "human_graph_cache_chr"
OUT_DIR_NAME = f"human_finemap_gw_{CACHE_VERSION}"
OUT_DIR = Path(f"/mnt/data/GraphGWAS/results/{OUT_DIR_NAME}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

_ALL_PHENOS = [
    {"phenocode":"21001","trait_type":"continuous","modifier":"irnt","label":"BMI"},
    {"phenocode":"50",   "trait_type":"continuous","modifier":"irnt","label":"Height"},
    {"phenocode":"30780","trait_type":"biomarkers","modifier":"irnt","label":"LDL"},
    {"phenocode":"30870","trait_type":"biomarkers","modifier":"irnt","label":"TG"},
]
# Skip Height: Pan-UKB EUR has all-low_confidence for chr22 + most chrs after liftover
PHENOS = [p for p in _ALL_PHENOS if p["label"] != "Height"]
GW_THRESHOLD = 5e-8
WINDOW_BP = 500_000
MAX_LEADS_PER_CHROM = 5

CHROM_LENGTHS_B37 = {
    "1": 249_250_621, "2": 243_199_373, "3": 198_022_430, "4": 191_154_276,
    "5": 180_915_260, "6": 171_115_067, "7": 159_138_663, "8": 146_364_022,
    "9": 141_213_431, "10":135_534_747, "11":135_006_516, "12":133_851_895,
    "13":115_169_878, "14":107_349_540, "15":102_531_392, "16":90_354_753,
    "17":81_195_210,  "18":78_077_248,  "19":59_128_983,  "20":63_025_520,
    "21":48_129_895,  "22":51_304_566,
}

_LIFT_38 = LiftOver("hg19", "hg38")
_LIFT_37 = LiftOver("hg38", "hg19")


def _b37_to_b38(chrom: str, pos: int):
    res = _LIFT_38.convert_coordinate(f"chr{chrom}", pos)
    return int(res[0][1]) if res else None


def _parse_gt(g):
    if g in ("0|0","0/0"): return 0.0
    if g in ("1|0","0|1","1/0","0/1"): return 1.0
    if g in ("1|1","1/1"): return 2.0
    return np.nan


def get_ld_vcf(chrom: str) -> Path:
    """1KG NYGC chr22 lives in GraphPop; chr1-21 + chr22 in GraphMana 1kGP_high_coverage."""
    # Primary: per-chrom 1KG high-coverage from GraphMana
    p = Path(f"/mnt/data/GraphMana/data/1000g/vcf/1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
    if p.exists():
        return p
    # Fallback: original CCDG chr22-only
    candidates = list(LD_VCF_DIR.glob(f"CCDG_*chr{chrom}.filtered*.vcf.gz"))
    return candidates[0] if candidates else None


def load_locus_dosage(chrom_b38: str, start: int, end: int):
    vcf = get_ld_vcf(chrom_b38.replace("chr",""))
    if not vcf or not vcf.exists():
        return pd.DataFrame(), np.zeros((0, 0))
    q = subprocess.run(
        ["bcftools","query","-r", f"{chrom_b38}:{start}-{end}",
         "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", str(vcf)],
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
        af = mu / 2.0 if 0 < mu < 2 else 0.0
        if af < 0.01 or af > 0.99: continue
        if "," in f[3]: continue
        c = f[0] if f[0].startswith("chr") else f"chr{f[0]}"
        rows.append((c, int(f[1]), f[2], f[3], af))
        mats.append(d)
    if not rows: return pd.DataFrame(), np.zeros((0,0))
    return pd.DataFrame(rows, columns=["chr","pos","ref","alt","af"]), np.vstack(mats)


def fetch_chr_sumstats(p: dict, chrom: str) -> pd.DataFrame | None:
    """Fetch a whole chromosome of Pan-UKB sumstats."""
    end_b37 = CHROM_LENGTHS_B37.get(chrom, 250_000_000)
    try:
        res = fetch_sumstats_locus(
            phenocode=p["phenocode"], chr=chrom, start=1, end=end_b37,
            ancestries=["EUR"],
            trait_type=p["trait_type"], modifier=p["modifier"],
        )
    except Exception as e:
        return None
    if "EUR" not in res:
        return None
    df = res["EUR"].variants.copy()
    df["pos_b37"] = df["pos"].astype(int)
    df["pos_b38"] = df["pos_b37"].apply(lambda x: _b37_to_b38(chrom, x))
    df = df.dropna(subset=["pos_b38"])
    df["pos_b38"] = df["pos_b38"].astype(int)
    df["pos"] = df["pos_b38"]
    df["chr"] = f"chr{chrom}"
    df["p"] = 10 ** (-df["log10p"])
    return df


def cluster_leads(df, p_thr, w):
    if "low_confidence" in df.columns:
        d = df[(df["p"] <= p_thr) & (~df["low_confidence"].astype(bool))].copy()
    else:
        d = df[df["p"] <= p_thr].copy()
    if d.empty: return []
    d = d.sort_values("p")
    leads, taken = [], []
    for _, r in d.iterrows():
        if any(abs(r["pos"] - t) < w for t in taken): continue
        taken.append(r["pos"])
        leads.append({"chrom": r["chr"], "pos": int(r["pos"]),
                      "p": float(r["p"]), "ref": r["ref"], "alt": r["alt"]})
    return leads


def fm_lead(label, lead, ss, full_cache):
    pos = lead["pos"]
    chrom = lead["chrom"]
    lo, hi = pos - WINDOW_BP//2, pos + WINDOW_BP//2
    sub = ss[ss["pos"].between(lo, hi)].copy()
    if len(sub) < 20: return None
    sub["z"] = sub["beta"] / sub["se"]
    sub["key"] = chrom + ":" + sub["pos"].astype(str) + ":" + sub["ref"] + ":" + sub["alt"]

    var_df, D = load_locus_dosage(chrom, lo, hi)
    if len(var_df) < 20: return None
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
    z_flip = {f"{chrom}:{r['pos']}:{r['alt']}:{r['ref']}": -r["z"]
              for _, r in sub.iterrows()}
    z = np.full(len(var_df), np.nan)
    for i, k in enumerate(var_df["variant_id"].values):
        if k in z_map: z[i] = z_map[k]
        elif k in z_flip: z[i] = z_flip[k]
    keep = np.isfinite(z)
    if keep.sum() < 20: return None
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
    l1_top = sorted(l1, key=lambda c: -c.pip)[0] if l1 else None
    hbp_top = sorted(hbp, key=lambda c: -c.pip)[0] if hbp else None
    return {
        "trait": label, "lead_chr": chrom, "lead_pos": pos,
        "lead_p": lead["p"], "n_variants": len(variants),
        "cache_coverage": len(window_cache),
        "l1_cs_size": len(l1_cs), "hbp_cs_size": len(hbp_cs),
        "l1_top_variant": l1_top.variant_id if l1_top else None,
        "l1_top_pip": float(l1_top.pip) if l1_top else None,
        "hbp_top_variant": hbp_top.variant_id if hbp_top else None,
        "hbp_top_pip": float(hbp_top.pip) if hbp_top else None,
    }


def process_chrom_phenotype(args):
    p, chrom, cache_path = args
    label = p["label"]
    print(f"  [{label} chr{chrom}] fetching sumstats ...", flush=True)
    ss = fetch_chr_sumstats(p, chrom)
    if ss is None or ss.empty:
        return label, chrom, [], None
    leads = cluster_leads(ss, GW_THRESHOLD, WINDOW_BP)
    leads = leads[:MAX_LEADS_PER_CHROM]
    summary = {"trait": label, "chrom": chrom, "n_var": len(ss),
               "n_gw": int((ss["p"] <= GW_THRESHOLD).sum()),
               "n_indep_leads": len(leads),
               "min_p": float(ss["p"].min()) if len(ss)>0 else 1.0}
    if not leads:
        return label, chrom, [], summary

    # Load chrom-specific cache only (saves memory)
    if Path(cache_path).exists():
        cache = json.load(open(cache_path))
    else:
        cache = {}
    print(f"  [{label} chr{chrom}] {len(leads)} leads, cache={len(cache):,}",
          flush=True)
    rows = []
    for lead in leads:
        try:
            r = fm_lead(label, lead, ss, cache)
            if r:
                rows.append(r)
        except Exception as e:
            print(f"    [{label} chr{chrom}:{lead['pos']}] error: {e}")
    return label, chrom, rows, summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--chroms", type=str, default="1-22")
    args = ap.parse_args()

    if "-" in args.chroms:
        a, b = args.chroms.split("-")
        chroms = [str(c) for c in range(int(a), int(b)+1)]
    else:
        chroms = args.chroms.split(",")

    tasks = []
    for p in PHENOS:
        for chrom in chroms:
            cache_path = ANN / f"{CACHE_PREFIX}{chrom}.json"
            tasks.append((p, chrom, str(cache_path)))
    print(f"Launching {len(tasks)} (phenotype × chromosome) tasks "
          f"with {args.workers} workers")

    all_loci, all_summary = [], []
    with mp.Pool(args.workers) as pool:
        for label, chrom, rows, summary in pool.imap_unordered(
                process_chrom_phenotype, tasks):
            all_loci.extend(rows)
            if summary:
                all_summary.append(summary)
            print(f"  done {label} chr{chrom}: {len(rows)} loci "
                  f"({len(all_loci)} total so far)", flush=True)

    pd.DataFrame(all_summary).to_csv(OUT_DIR / "human_gwas_summary_gw.tsv",
                                      sep="\t", index=False)
    pd.DataFrame(all_loci).to_csv(OUT_DIR / "human_finemap_summary_gw.tsv",
                                   sep="\t", index=False)
    print(f"\nWrote {OUT_DIR}/human_finemap_summary_gw.tsv ({len(all_loci)} loci)")


if __name__ == "__main__":
    main()
