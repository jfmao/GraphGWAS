"""Fine-map IRRI GWAS lead loci with L1 + HBP using the rice multi-omics
graph cache. For each (trait, lead locus) pair in irri_lead_loci.tsv,
loads sumstats, builds LD from the 3kRG VCF, runs both methods, and
writes the credible sets + PIP-ranked variants.

Emits a summary listing the top PIP variant per locus along with any
overlap with Ren 2023 ground-truth grain-quality genes (for traits that
are grain-related).
"""
from __future__ import annotations

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

DATA = Path("/mnt/data/GraphGWAS/data/rice_3k")
RES = DATA / "results"
LEADS = RES / "irri_gwas" / "irri_lead_loci.tsv"
GWAS_OUT = RES / "irri_gwas"
FM_OUT = RES / "irri_finemap"
FM_OUT.mkdir(parents=True, exist_ok=True)

VCF = "/mnt/data/GraphPop/data/raw/3kRG_data/NB_final_snp.vcf.gz"
GT_GENES = pd.read_csv(DATA / "ground_truth" / "grain_quality_causal_genes.tsv",
                        sep="\t")
FM_WINDOW = 250_000
MAX_LEADS_PER_TRAIT = 5   # keep it tractable
MAX_SUGGESTIVE_PER_TRAIT = 3


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
        rows.append((f[0], int(f[1]), f[2], f[3], af))
        mats.append(d)
    if not rows:
        return pd.DataFrame(), np.zeros((0,0))
    return pd.DataFrame(rows, columns=["chr","pos","ref","alt","af"]), np.vstack(mats)


def run_finemap_for_lead(trait, lead_chrom, lead_pos, gw_path, chr_cache):
    """Fine-map one locus; return dict of results + top PIP variants."""
    lo, hi = lead_pos - FM_WINDOW, lead_pos + FM_WINDOW
    # Load sumstats in window
    gw = pd.read_csv(gw_path, sep="\t", low_memory=False)
    gw.columns = [c.lstrip("#") for c in gw.columns]
    gw = gw[(gw["TEST"]=="ADD") & (gw["ERRCODE"]==".")].copy()
    gw["CHROM"] = gw["CHROM"].astype(str).apply(
        lambda c: c if c.startswith("Chr") else f"Chr{c}")
    for c in ("BETA","SE","P"):
        gw[c] = pd.to_numeric(gw[c], errors="coerce")
    gw = gw.dropna(subset=["P","BETA","SE"])
    gw["z"] = gw["BETA"] / gw["SE"]
    gw_sub = gw[(gw["CHROM"]==lead_chrom) & (gw["POS"].between(lo,hi))].copy()
    if len(gw_sub) < 20:
        return None

    # Load VCF locus and build LD
    var_df, D = load_locus(lead_chrom, lo, hi)
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
    window_cache = {k: chr_cache[k] for k in var_df["variant_id"].values
                    if k in chr_cache} if chr_cache else {}

    # Run L1 + HBP
    l1 = l1_finemap_from_sumstats(
        variants, z, R_sq, z_func=None, alpha=1.0,
        r2_smooth=0.3, credible_set_coverage=0.95, chr_name=lead_chrom,
    )
    hbp = hbp_finemap_from_sumstats(
        variants, z, R_sq, graph_cache=window_cache,
        r2_smooth=0.3, credible_set_coverage=0.95, chr_name=lead_chrom,
    )
    l1_cs = [c for c in l1 if c.in_credible_set]
    hbp_cs = [c for c in hbp if c.in_credible_set]
    l1_top = sorted(l1, key=lambda c: -c.pip)[:3]
    hbp_top = sorted(hbp, key=lambda c: -c.pip)[:3]
    return {
        "trait": trait, "lead_chr": lead_chrom, "lead_pos": lead_pos,
        "n_variants": len(variants),
        "cache_coverage": len(window_cache),
        "l1_cs_size": len(l1_cs), "hbp_cs_size": len(hbp_cs),
        "l1_top_variants": [(c.variant_id, c.pip) for c in l1_top],
        "hbp_top_variants": [(c.variant_id, c.pip) for c in hbp_top],
    }


def main():
    if not LEADS.exists():
        print(f"No lead loci at {LEADS}")
        return
    leads = pd.read_csv(LEADS, sep="\t")
    print(f"Loaded {len(leads)} total lead loci across {leads['trait'].nunique()} traits")

    # Cap per trait to keep tractable
    def pick(g):
        gw = g[g["sig_level"]=="GW"].head(MAX_LEADS_PER_TRAIT)
        sug = g[g["sig_level"]=="suggestive"].head(MAX_SUGGESTIVE_PER_TRAIT)
        return pd.concat([gw, sug])
    leads = leads.sort_values("p").groupby("trait", group_keys=False).apply(pick).reset_index(drop=True)
    print(f"After cap: {len(leads)} leads selected")

    # Cache loading
    chr_caches = {}
    def get_chr_cache(chrom):
        if chrom not in chr_caches:
            p = DATA / "annotations" / f"rice_graph_cache_{chrom}.json"
            if p.exists():
                print(f"  loading {p.name}...", flush=True)
                chr_caches[chrom] = json.load(open(p))
            else:
                chr_caches[chrom] = {}
        return chr_caches[chrom]

    all_results = []
    for _, lead in leads.iterrows():
        trait = lead["trait"]; chrom = lead["chrom"]; pos = int(lead["pos"])
        sig = lead["sig_level"]
        print(f"\n--- {trait} {chrom}:{pos:,} [{sig}, p={lead['p']:.2e}] ---",
              flush=True)
        # Find the GWAS sumstats file for this trait
        gwf = list(GWAS_OUT.glob(f"gwas_{trait}.*.glm.linear"))
        if not gwf:
            print(f"  no GWAS output for {trait}"); continue
        gw_path = gwf[0]
        cc = get_chr_cache(chrom)
        try:
            r = run_finemap_for_lead(trait, chrom, pos, gw_path, cc)
        except Exception as e:
            print(f"  ERROR: {e}")
            continue
        if r is None:
            print("  skipped (insufficient data)")
            continue
        r["sig_level"] = sig
        r["lead_p"] = float(lead["p"])
        all_results.append(r)
        print(f"  n_var={r['n_variants']}  cache={r['cache_coverage']}  "
              f"L1 CS={r['l1_cs_size']}  HBP CS={r['hbp_cs_size']}")
        print(f"    L1 top: {r['l1_top_variants']}")
        print(f"    HBP top: {r['hbp_top_variants']}")

    # Compile summary
    rows = []
    for r in all_results:
        l1_top = r["l1_top_variants"][0] if r["l1_top_variants"] else (None, None)
        hbp_top = r["hbp_top_variants"][0] if r["hbp_top_variants"] else (None, None)
        rows.append({
            "trait": r["trait"], "lead_chr": r["lead_chr"], "lead_pos": r["lead_pos"],
            "sig_level": r["sig_level"], "lead_p": r["lead_p"],
            "n_variants": r["n_variants"], "cache_coverage": r["cache_coverage"],
            "l1_cs_size": r["l1_cs_size"], "hbp_cs_size": r["hbp_cs_size"],
            "l1_top_variant": l1_top[0], "l1_top_pip": l1_top[1],
            "hbp_top_variant": hbp_top[0], "hbp_top_pip": hbp_top[1],
        })
    df = pd.DataFrame(rows)
    out = FM_OUT / "irri_finemap_summary.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"\nWrote {out}  ({len(df)} loci fine-mapped)")
    if len(df) > 0:
        print(df.to_string(index=False))


if __name__ == "__main__":
    main()
