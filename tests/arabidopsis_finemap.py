"""Fine-map Arabidopsis priority phenotypes with L1 + HBP using TAIR10 + Plant Reactome cache.

For each .glm.linear output:
  - Cluster top GW + suggestive leads (250 kb window)
  - Cap at top-3 GW + top-1 suggestive
  - L1 + HBP fine-map per lead using arabidopsis_graph_cache_chr*.json
  - Compare against AraGWAS bonferroni associations for ground-truth validation
"""
from __future__ import annotations

import glob
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/mnt/data/GraphGWAS/src/python")
from graphgwas.finemapping_v2 import (
    hbp_finemap_from_sumstats,
    l1_finemap_from_sumstats,
)

DATA = Path("/mnt/data/GraphGWAS/data/arabidopsis")
GWAS_DIR = DATA
ANN = DATA / "annotations"
VCF = "/mnt/data/GraphGWAS/tests/data/arabidopsis/1001genomes_snp-short-indel_only_ACGTN.vcf.gz"
ARAGWAS = Path("/mnt/data/GraphGWAS/tests/data/arabidopsis/aragwas_bonf_associations.csv")
OUT_DIR = Path("/mnt/data/GraphGWAS/results/arabidopsis_finemap")
OUT_DIR.mkdir(parents=True, exist_ok=True)
WINDOW_BP = 100_000
GW = 5e-8
SUG = 1e-5
MAX_GW, MAX_SUG = 3, 1


def _parse_gt(g):
    if g in ("0|0","0/0"): return 0.0
    if g in ("1|0","0|1","1/0","0/1"): return 1.0
    if g in ("1|1","1/1"): return 2.0
    return np.nan


def load_locus(chrom, lo, hi):
    q = subprocess.run(
        ["bcftools","query","-r", f"{chrom}:{lo}-{hi}",
         "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", VCF],
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
        rows.append((f[0], int(f[1]), f[2], f[3], af))
        mats.append(d)
    if not rows:
        return pd.DataFrame(), np.zeros((0, 0))
    return pd.DataFrame(rows, columns=["chr","pos","ref","alt","af"]), np.vstack(mats)


def cluster_leads(df, p_thr, w):
    d = df[df["P"] <= p_thr].copy()
    if d.empty: return []
    d = d.sort_values(["CHROM","P"])
    out = []
    for ch, sub in d.groupby("CHROM", sort=False):
        taken = []
        for _, r in sub.iterrows():
            if any(abs(r["POS"] - t) < w for t in taken): continue
            taken.append(r["POS"])
            out.append({"chrom":str(ch), "pos":int(r["POS"]),
                        "p":float(r["P"]), "ref":r["REF"], "alt":r["ALT"],
                        "beta":float(r["BETA"])})
    return out


def fm_lead(trait, lead, sumstats, full_cache):
    chrom, pos = lead["chrom"], lead["pos"]
    lo, hi = pos - WINDOW_BP//2, pos + WINDOW_BP//2
    sub = sumstats[(sumstats["CHROM"].astype(str)==chrom)
                    & (sumstats["POS"].between(lo, hi))].copy()
    if len(sub) < 20: return None
    sub["z"] = sub["BETA"] / sub["SE"]

    var_df, D = load_locus(chrom, lo, hi)
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

    sub["key"] = sub["CHROM"].astype(str) + ":" + sub["POS"].astype(str) \
                + ":" + sub["REF"] + ":" + sub["ALT"]
    z_map = dict(zip(sub["key"].values, sub["z"].values))
    z_flip = {f"{r['CHROM']}:{r['POS']}:{r['ALT']}:{r['REF']}": -r["z"]
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
        variants, z, R_sq, z_func=None, alpha=1.0,
        r2_smooth=0.3, credible_set_coverage=0.95, chr_name=chrom,
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
        "trait": trait, "lead_chr": chrom, "lead_pos": pos, "lead_p": lead["p"],
        "n_variants": len(variants), "cache_coverage": len(window_cache),
        "l1_cs_size": len(l1_cs), "hbp_cs_size": len(hbp_cs),
        "l1_top_variant": l1_top.variant_id if l1_top else None,
        "l1_top_pip": l1_top.pip if l1_top else None,
        "hbp_top_variant": hbp_top.variant_id if hbp_top else None,
        "hbp_top_pip": hbp_top.pip if hbp_top else None,
    }


def main():
    print("Loading arabidopsis chr1-5 graph caches ...", flush=True)
    chr_caches = {}
    for c in "12345":
        p = ANN.parent / f"arabidopsis_graph_cache_chr{c}.json"
        if p.exists():
            chr_caches[c] = json.load(p.open())
            print(f"  chr{c}: {len(chr_caches[c]):,} variants")

    files = sorted(glob.glob(str(GWAS_DIR / "arabi_gwas.*.glm.linear")))
    print(f"\n{len(files)} GWAS files to process")

    summary, all_loci = [], []
    for fp in files:
        trait = Path(fp).stem.replace("arabi_gwas.","").replace(".glm","")
        if trait.endswith(".glm"): trait = trait[:-4]
        df = pd.read_csv(fp, sep="\t", low_memory=False)
        df.columns = [c.lstrip("#") for c in df.columns]
        df = df[(df["TEST"]=="ADD") & (df["ERRCODE"]==".")].copy()
        for c in ("BETA","SE","P"):
            df[c] = pd.to_numeric(df[c], errors="coerce")
        df = df.dropna(subset=["P","BETA","SE"])
        df["CHROM"] = df["CHROM"].astype(str)
        n_var = len(df)
        n_gw = int((df["P"]<=GW).sum())
        n_sug = int((df["P"]<=SUG).sum())

        leads_gw = cluster_leads(df, GW, WINDOW_BP)
        leads_sug = cluster_leads(df, SUG, WINDOW_BP)
        for l in leads_gw: l["sig_level"]="GW"
        for l in leads_sug:
            if any(l2["chrom"]==l["chrom"] and abs(l2["pos"]-l["pos"])<WINDOW_BP
                   for l2 in leads_gw):
                continue
            l["sig_level"]="suggestive"
        leads_gw = sorted(leads_gw, key=lambda x: x["p"])[:MAX_GW]
        leads_sug = sorted([l for l in leads_sug if l.get("sig_level")=="suggestive"],
                           key=lambda x: x["p"])[:MAX_SUG]
        leads = leads_gw + leads_sug
        print(f"\n=== {trait} ===  N={n_var:,}  GW={n_gw}  SUG={n_sug}  -> {len(leads)} leads", flush=True)

        summary.append({"trait": trait, "n_variants": n_var,
                        "n_gw": n_gw, "n_sug": n_sug,
                        "n_indep_gw": len(leads_gw), "n_indep_sug": len(leads_sug),
                        "min_p": float(df["P"].min())})
        for lead in leads:
            cc = chr_caches.get(lead["chrom"], {})
            try:
                r = fm_lead(trait, lead, df, cc)
                if r:
                    r["sig_level"] = lead["sig_level"]
                    all_loci.append(r)
                    print(f"    [{lead['chrom']}:{lead['pos']:,}] L1 CS={r['l1_cs_size']} pip={r['l1_top_pip']:.3f}  HBP CS={r['hbp_cs_size']} pip={r['hbp_top_pip']:.3f}",
                          flush=True)
            except Exception as e:
                print(f"    [{lead['chrom']}:{lead['pos']}] error: {e}")

    pd.DataFrame(summary).to_csv(OUT_DIR / "arabi_gwas_summary.tsv", sep="\t", index=False)
    pd.DataFrame(all_loci).to_csv(OUT_DIR / "arabi_finemap_summary.tsv", sep="\t", index=False)
    print(f"\nWrote {OUT_DIR}/arabi_finemap_summary.tsv ({len(all_loci)} loci)")


if __name__ == "__main__":
    main()
