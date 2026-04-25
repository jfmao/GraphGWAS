"""Human Pan-UKB chr22 fine-mapping for 4 phenotypes (BMI, Height, LDL, TG).

For each phenotype:
  - Use tabix-over-HTTPS to fetch chr22 sumstats from Pan-UKB
  - Cluster GW-sig leads (p ≤ 5e-8) at 500 kb window
  - Fine-map top-N leads with L1 + HBP using human chr22 graph cache
  - Write per-phenotype TSV + combined summary

Outputs:
  results/human_finemap/human_finemap_<TRAIT>.tsv
  results/human_finemap/human_finemap_summary.tsv
  results/human_finemap/human_chr22_leads_<TRAIT>.tsv  (full lead list)
"""
from __future__ import annotations

import json
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

# Pan-UKB is GRCh37; our 1KG LD ref + graph cache are GRCh38.
_LIFT = LiftOver("hg19", "hg38")


def _b37_to_b38(pos: int) -> int | None:
    res = _LIFT.convert_coordinate("chr22", pos)
    return int(res[0][1]) if res else None

CACHE = Path("/mnt/data/GraphGWAS/data/annotations/human_graph_cache_chr22.json")
LD_VCF = "/mnt/data/GraphPop/data/raw/1000g/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz"
OUT_DIR = Path("/mnt/data/GraphGWAS/results/human_finemap")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# 4 phenotypes
PHENOS = [
    {"phenocode":"21001", "trait_type":"continuous", "modifier":"irnt",
     "label":"BMI", "n_cases":439590},
    {"phenocode":"50",    "trait_type":"continuous", "modifier":"irnt",
     "label":"Height", "n_cases":438478},
    {"phenocode":"30780", "trait_type":"biomarkers", "modifier":"irnt",
     "label":"LDL", "n_cases":419831},
    {"phenocode":"30870", "trait_type":"biomarkers", "modifier":"irnt",
     "label":"TG", "n_cases":420271},
]
ANCESTRIES = ("EUR",)
GW_THRESHOLD = 5e-8
WINDOW_BP = 500_000   # 500 kb (longer LD in humans)
MAX_LEADS_PER_TRAIT = 5


def _parse_gt(g: str) -> float:
    if g in ("0|0","0/0"): return 0.0
    if g in ("1|0","0|1","1/0","0/1"): return 1.0
    if g in ("1|1","1/1"): return 2.0
    return np.nan


def load_locus_dosage(start: int, end: int):
    q = subprocess.run(
        ["bcftools","query","-r", f"chr22:{start}-{end}",
         "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", LD_VCF],
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
    if not rows:
        return pd.DataFrame(), np.zeros((0, 0))
    return pd.DataFrame(rows, columns=["chr","pos","ref","alt","af"]), np.vstack(mats)


def fetch_chr22_sumstats(p: dict) -> pd.DataFrame:
    """Fetch entire chr22 in one tabix call. Liftover b37→b38 inline."""
    print(f"  [{p['label']}] fetching chr22 sumstats from Pan-UKB ...", flush=True)
    res = fetch_sumstats_locus(
        phenocode=p["phenocode"], chr="22", start=1, end=51_304_566,
        ancestries=ANCESTRIES,
        trait_type=p["trait_type"], modifier=p["modifier"],
    )
    if "EUR" not in res:
        return pd.DataFrame()
    df = res["EUR"].variants.copy()  # GRCh37 chr22:pos:ref:alt
    df["pos_b37"] = df["pos"].astype(int)
    print(f"    {len(df):,} variants on chr22 b37; lifting over to b38 ...", flush=True)
    df["pos_b38"] = df["pos_b37"].apply(_b37_to_b38)
    df = df.dropna(subset=["pos_b38"])
    df["pos_b38"] = df["pos_b38"].astype(int)
    df = df.rename(columns={"pos": "pos_orig"})
    df["pos"] = df["pos_b38"]
    df["p"] = 10 ** (-df["log10p"])
    df["chr"] = "chr22"
    print(f"    {len(df):,} after liftover", flush=True)
    return df


def cluster_leads(df: pd.DataFrame, p_thresh: float, window: int) -> list[dict]:
    d = df[df["p"] <= p_thresh].copy()
    if "low_confidence" in d.columns:
        d = d[~d["low_confidence"].astype(bool)].copy()
    if d.empty: return []
    d = d.sort_values("p")
    leads, taken = [], []
    for _, r in d.iterrows():
        if any(abs(r["pos"] - t) < window for t in taken):
            continue
        taken.append(r["pos"])
        leads.append({"chrom":"chr22", "pos":int(r["pos"]),
                      "p":float(r["p"]), "beta":float(r["beta"]),
                      "ref":r["ref"], "alt":r["alt"]})
    return leads


def finemap_lead(label, lead, ss_df, full_cache):
    pos = lead["pos"]
    lo, hi = pos - WINDOW_BP//2, pos + WINDOW_BP//2
    sub = ss_df[ss_df["pos"].between(lo, hi)].copy()
    if len(sub) < 20: return None
    sub["z"] = sub["beta"] / sub["se"]
    sub["key"] = "chr22:" + sub["pos"].astype(str) + ":" + sub["ref"] + ":" + sub["alt"]

    var_df, D = load_locus_dosage(lo, hi)
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
    z_flip = {f"chr22:{r['pos']}:{r['alt']}:{r['ref']}": -r["z"]
              for _, r in sub.iterrows()}
    z = np.full(len(var_df), np.nan)
    for i, k in enumerate(var_df["variant_id"].values):
        if k in z_map: z[i] = z_map[k]
        elif k in z_flip: z[i] = z_flip[k]
    keep = np.isfinite(z)
    if keep.sum() < 20: return None
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
                    if k in full_cache}

    l1 = l1_finemap_from_sumstats(
        variants, z, R_sq, z_func=None, alpha=1.0,
        r2_smooth=0.3, credible_set_coverage=0.95, chr_name="chr22",
    )
    hbp = hbp_finemap_from_sumstats(
        variants, z, R_sq, graph_cache=window_cache,
        r2_smooth=0.3, credible_set_coverage=0.95, chr_name="chr22",
    )
    l1_cs = [c for c in l1 if c.in_credible_set]
    hbp_cs = [c for c in hbp if c.in_credible_set]
    l1_top = sorted(l1, key=lambda c: -c.pip)[:1]
    hbp_top = sorted(hbp, key=lambda c: -c.pip)[:1]
    return {
        "trait": label, "lead_chr": "chr22", "lead_pos": pos,
        "lead_p": lead["p"], "n_variants": len(variants),
        "cache_coverage": len(window_cache),
        "l1_cs_size": len(l1_cs), "hbp_cs_size": len(hbp_cs),
        "l1_top_variant": l1_top[0].variant_id if l1_top else None,
        "l1_top_pip": l1_top[0].pip if l1_top else None,
        "hbp_top_variant": hbp_top[0].variant_id if hbp_top else None,
        "hbp_top_pip": hbp_top[0].pip if hbp_top else None,
    }


def main():
    print("Loading human chr22 graph cache ...", flush=True)
    full_cache = json.load(CACHE.open())
    print(f"  {len(full_cache):,} variants in cache")

    summary = []
    all_leads = []
    finemap_rows = []
    for p in PHENOS:
        print(f"\n=== {p['label']} (Pan-UKB pheno {p['phenocode']}) ===", flush=True)
        try:
            ss = fetch_chr22_sumstats(p)
        except Exception as e:
            print(f"  fetch failed: {e}")
            continue
        if ss.empty:
            print("  empty sumstats")
            continue
        n_gw = int((ss["p"] <= GW_THRESHOLD).sum())
        leads = cluster_leads(ss, GW_THRESHOLD, WINDOW_BP)
        leads = leads[:MAX_LEADS_PER_TRAIT]
        print(f"  {n_gw} GW-sig variants, {len(leads)} indep leads (top {MAX_LEADS_PER_TRAIT})")
        all_leads.extend([{**l, "trait": p["label"]} for l in leads])
        summary.append({
            "trait": p["label"], "phenocode": p["phenocode"],
            "n_variants_chr22": len(ss),
            "n_gw_chr22": n_gw, "n_indep_leads_chr22": len(leads),
            "min_p": float(ss["p"].min()) if len(ss)>0 else 1.0,
        })

        for lead in leads:
            try:
                r = finemap_lead(p["label"], lead, ss, full_cache)
                if r:
                    finemap_rows.append(r)
                    print(f"    [{r['lead_chr']}:{r['lead_pos']}] "
                          f"L1 CS={r['l1_cs_size']} pip={r['l1_top_pip']:.3f}  "
                          f"HBP CS={r['hbp_cs_size']} pip={r['hbp_top_pip']:.3f}",
                          flush=True)
            except Exception as e:
                print(f"    [{lead['chrom']}:{lead['pos']}] error: {e}")

    pd.DataFrame(summary).to_csv(OUT_DIR / "human_gwas_summary.tsv", sep="\t",
                                  index=False)
    pd.DataFrame(all_leads).to_csv(OUT_DIR / "human_chr22_leads.tsv", sep="\t",
                                    index=False)
    pd.DataFrame(finemap_rows).to_csv(OUT_DIR / "human_finemap_summary.tsv",
                                       sep="\t", index=False)
    print(f"\nWrote {OUT_DIR}/human_gwas_summary.tsv ({len(summary)} traits)")
    print(f"Wrote {OUT_DIR}/human_finemap_summary.tsv ({len(finemap_rows)} loci)")


if __name__ == "__main__":
    main()
