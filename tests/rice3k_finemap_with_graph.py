"""Re-run rice 3K fine-mapping with the rice multi-omics graph_cache
loaded. Compare against the baseline (empty graph_cache) run from
tests/rice3k_finemap_all_loci.py.

Writes:
  results/finemap_with_graph/{locus}_{l1|hbp}.tsv
  results/rice3k_graph_vs_baseline.tsv  — per-locus comparison
  results/rice3k_graph_vs_baseline.md   — summary report
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

DATA_DIR = Path("/mnt/data/GraphGWAS/data/rice_3k")
RES = DATA_DIR / "results"
FM_BASE = RES / "finemap"          # baseline (empty graph_cache)
FM_GRAPH = RES / "finemap_with_graph"
FM_GRAPH.mkdir(parents=True, exist_ok=True)

GWAS = RES / "rice_gwas.GRAIN_SIZE.glm.linear"
LEADS = RES / "lead_loci.tsv"
CACHE_DIR = DATA_DIR / "annotations"


def _load_chr_cache(chrom: str) -> dict:
    """Lazy loader for per-chromosome cache JSON."""
    p = CACHE_DIR / f"rice_graph_cache_{chrom}.json"
    if not p.exists():
        return {}
    with open(p) as f:
        return json.load(f)
VCF = "/mnt/data/GraphPop/data/raw/3kRG_data/NB_final_snp.vcf.gz"
SIMULATED_CAUSALS = DATA_DIR / "pheno" / "causal_variants.tsv"

FM_WINDOW = 250_000
R2_SMOOTH = 0.3


def load_sumstats() -> pd.DataFrame:
    df = pd.read_csv(GWAS, sep="\t", low_memory=False)
    df.columns = [c.lstrip("#") for c in df.columns]
    df = df[(df["TEST"] == "ADD") & (df["ERRCODE"] == ".")].copy()
    df["CHROM"] = df["CHROM"].astype(str).apply(
        lambda c: c if c.startswith("Chr") else f"Chr{c}")
    for c in ("BETA", "SE", "P"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["P", "BETA", "SE"])
    df["z"] = df["BETA"] / df["SE"]
    return df


def _parse_gt(g: str) -> float:
    if g in ("0|0", "0/0"): return 0.0
    if g in ("1|0", "0|1", "1/0", "0/1"): return 1.0
    if g in ("1|1", "1/1"): return 2.0
    return np.nan


def load_locus(chrom: str, start: int, end: int, min_af: float = 0.01):
    region = f"{chrom}:{start}-{end}"
    q = subprocess.run(
        ["bcftools", "query", "-r", region,
         "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", VCF],
        capture_output=True, text=True, timeout=300,
    )
    if q.returncode != 0 or not q.stdout:
        return pd.DataFrame(), np.zeros((0, 0))
    rows, mats = [], []
    for line in q.stdout.strip().split("\n"):
        f = line.split("\t")
        gts = f[4:]
        d = np.array([_parse_gt(g) for g in gts], dtype=np.float32)
        mu = np.nanmean(d) if np.isfinite(d).any() else 0.0
        d = np.where(np.isnan(d), mu, d)
        af = mu / 2.0 if 0 < mu < 2 else 0.0
        if af < min_af or af > 1 - min_af:
            continue
        rows.append((f[0], int(f[1]), f[2], f[3], af))
        mats.append(d)
    if not rows:
        return pd.DataFrame(), np.zeros((0, 0))
    return pd.DataFrame(rows, columns=["chr","pos","ref","alt","af"]), np.vstack(mats)


def finemap_locus(lead_row, gwas, cache):
    chrom = lead_row["CHROM"]
    lead_pos = int(lead_row["POS"])
    lo, hi = lead_pos - FM_WINDOW, lead_pos + FM_WINDOW
    var_df, D = load_locus(chrom, lo, hi)
    if len(var_df) < 10:
        return [], [], 0, 0
    # LD matrix
    mu = D.mean(axis=1, keepdims=True)
    std = np.where(D.std(axis=1, keepdims=True) < 1e-10, 1.0,
                   D.std(axis=1, keepdims=True))
    Z = (D - mu) / std
    R2 = np.clip(((Z @ Z.T) / Z.shape[1]) ** 2, 0, 1).astype(np.float64)
    np.fill_diagonal(R2, 1.0)
    vid = (var_df["chr"] + ":" + var_df["pos"].astype(str) + ":"
           + var_df["ref"] + ":" + var_df["alt"])
    var_df["variant_id"] = vid.values

    gw = gwas[(gwas["CHROM"]==chrom) & (gwas["POS"].between(lo, hi))].copy()
    gw["key"] = gw["CHROM"].astype(str)+":"+gw["POS"].astype(str)+":"+gw["REF"]+":"+gw["ALT"]
    z_map = dict(zip(gw["key"].values, gw["z"].values))
    z_flip = {f"{r['CHROM']}:{r['POS']}:{r['ALT']}:{r['REF']}": -r["z"]
              for _, r in gw.iterrows()}
    z = np.full(len(var_df), np.nan)
    for i, k in enumerate(var_df["variant_id"].values):
        if k in z_map:
            z[i] = z_map[k]
        elif k in z_flip:
            z[i] = z_flip[k]
    keep = np.isfinite(z)
    if keep.sum() < 10:
        return [], [], 0, 0
    idx = np.where(keep)[0]
    var_df = var_df.loc[keep].reset_index(drop=True)
    z = z[keep]
    R2 = R2[np.ix_(idx, idx)]

    variants = [
        {"variantId": v["variant_id"], "chr": v["chr"], "pos": int(v["pos"]),
         "ref": v["ref"], "alt": v["alt"], "af_total": float(v["af"])}
        for _, v in var_df.iterrows()
    ]

    # Build per-variant functional score for L1 from the graph cache
    gene_hits = 0
    pathway_hits = 0
    z_func = np.zeros(len(variants))
    for i, v in enumerate(variants):
        info = cache.get(v["variantId"])
        if info is None:
            continue
        # Weighted functional score: genes + 0.5×pathways + 0.3×ppi partners
        score = 1.0 * len(info.get("genes", [])) \
            + 0.5 * len(info.get("pathways", [])) \
            + 0.3 * min(len(info.get("ppi", [])), 10)
        z_func[i] = np.log1p(score)
        if info.get("genes"):
            gene_hits += 1
        if info.get("pathways") or info.get("ppi"):
            pathway_hits += 1
    if z_func.max() > 0:
        z_func = z_func / z_func.max() * np.max(np.abs(z))

    l1 = l1_finemap_from_sumstats(
        variants, z, R2, z_func=z_func, alpha=0.5,
        r2_smooth=R2_SMOOTH, chr_name=chrom,
    )
    hbp = hbp_finemap_from_sumstats(
        variants, z, R2, graph_cache=cache,
        r2_smooth=R2_SMOOTH, chr_name=chrom,
    )
    return l1, hbp, gene_hits, pathway_hits


def compare_one(locus_id, l1, hbp, sim_causals):
    """Summarise CS size, top-PIP, and simulated-causal recovery."""
    def _stats(cand):
        cs = [c for c in cand if c.in_credible_set]
        if not cs:
            return None, 0, 0.0
        top = max(cs, key=lambda c: c.pip)
        # Does any CS variant fall within 100 kb of any simulated causal?
        within = None
        for c in cs:
            chrom = c.chr
            for _, s in sim_causals.iterrows():
                if chrom == s["chr"] and abs(c.pos - s["pos"]) <= 100_000:
                    within = (s["gene"], abs(c.pos - s["pos"]))
                    break
            if within is not None:
                break
        return within, len(cs), float(top.pip)

    return {"L1": _stats(l1), "HBP": _stats(hbp)}


def main():
    print("[1/3] Loading GWAS sumstats + leads + simulated causals...")
    gwas = load_sumstats()
    leads = pd.read_csv(LEADS, sep="\t")
    sim = pd.read_csv(SIMULATED_CAUSALS, sep="\t")
    print(f"  GWAS: {len(gwas):,} rows; leads: {len(leads)}")
    print(f"  Simulated causals: {len(sim)}")

    # Lazy-load chromosome caches on demand
    chr_caches: dict[str, dict] = {}
    def _get_cache(chrom: str) -> dict:
        if chrom not in chr_caches:
            print(f"    loading graph cache for {chrom}...", flush=True)
            chr_caches[chrom] = _load_chr_cache(chrom)
            print(f"      {len(chr_caches[chrom]):,} variants with annotations", flush=True)
        return chr_caches[chrom]

    print("\n[2/3] Re-fine-mapping 21 loci WITH graph cache (lazy per-chrom)...")
    rows = []
    for i, lead in leads.iterrows():
        lid = f"{lead['CHROM']}_{int(lead['POS'])}"
        cache = _get_cache(lead["CHROM"])
        l1, hbp, gene_hits, pw_hits = finemap_locus(lead, gwas, cache)
        if not l1:
            continue
        # Save outputs
        def _save(cands, method):
            df = pd.DataFrame([{
                "variant_id": c.variant_id, "chr": c.chr, "pos": c.pos,
                "pip": c.pip, "in_cs": c.in_credible_set,
                "combined_score": c.combined_score,
            } for c in cands])
            df.to_csv(FM_GRAPH / f"{lid}_{method}.tsv", sep="\t", index=False)
        _save(l1, "l1")
        _save(hbp, "hbp")

        stats = compare_one(lid, l1, hbp, sim)
        rows.append({
            "locus_id": lid,
            "lead_chr": lead["CHROM"],
            "lead_pos": int(lead["POS"]),
            "lead_p": lead["P"],
            "var_with_gene_info": gene_hits,
            "var_with_pathway_or_ppi": pw_hits,
            "L1_cs_size_with_graph": stats["L1"][1],
            "L1_top_pip_with_graph": stats["L1"][2],
            "L1_sim_hit": stats["L1"][0][0] if stats["L1"][0] else None,
            "L1_sim_dist": stats["L1"][0][1] if stats["L1"][0] else None,
            "HBP_cs_size_with_graph": stats["HBP"][1],
            "HBP_top_pip_with_graph": stats["HBP"][2],
            "HBP_sim_hit": stats["HBP"][0][0] if stats["HBP"][0] else None,
            "HBP_sim_dist": stats["HBP"][0][1] if stats["HBP"][0] else None,
        })
        print(f"  [{i+1}/{len(leads)}] {lid}  "
              f"gene-rich:{gene_hits:>4}  prior-rich:{pw_hits:>3}  "
              f"L1 CS={stats['L1'][1]:>5} pip={stats['L1'][2]:.4f}  "
              f"HBP CS={stats['HBP'][1]:>5} pip={stats['HBP'][2]:.4f}",
              flush=True)

    out = pd.DataFrame(rows)
    with_graph = RES / "rice3k_graph_vs_baseline.tsv"
    out.to_csv(with_graph, sep="\t", index=False)
    print(f"\n  Wrote {with_graph}")

    print("\n[3/3] Computing baseline-vs-graph deltas...")
    # Baseline stats from the first run (files in FM_BASE)
    base_stats = []
    for f in FM_BASE.glob("*_l1.tsv"):
        lid = f.stem.replace("_l1", "")
        l1 = pd.read_csv(f, sep="\t")
        hbp = pd.read_csv(FM_BASE / f"{lid}_hbp.tsv", sep="\t")
        base_stats.append({
            "locus_id": lid,
            "L1_cs_baseline": int((l1["in_cs"]).sum()),
            "L1_pip_baseline": float(l1[l1["in_cs"]]["pip"].max()) if (l1["in_cs"]).any() else 0.0,
            "HBP_cs_baseline": int((hbp["in_cs"]).sum()),
            "HBP_pip_baseline": float(hbp[hbp["in_cs"]]["pip"].max()) if (hbp["in_cs"]).any() else 0.0,
        })
    base = pd.DataFrame(base_stats)
    merged = out.merge(base, on="locus_id", how="left")
    merged["L1_cs_reduction"] = merged["L1_cs_baseline"] - merged["L1_cs_size_with_graph"]
    merged["L1_cs_reduction_pct"] = (100 * merged["L1_cs_reduction"] / merged["L1_cs_baseline"]).round(1)
    merged["HBP_cs_reduction"] = merged["HBP_cs_baseline"] - merged["HBP_cs_size_with_graph"]
    merged["HBP_cs_reduction_pct"] = (100 * merged["HBP_cs_reduction"] / merged["HBP_cs_baseline"]).round(1)
    merged["L1_pip_boost"] = merged["L1_top_pip_with_graph"] - merged["L1_pip_baseline"]
    merged["HBP_pip_boost"] = merged["HBP_top_pip_with_graph"] - merged["HBP_pip_baseline"]
    merged.to_csv(with_graph, sep="\t", index=False)

    # Markdown summary
    lines = [
        "# Rice 3K — effect of loading the multi-omics graph cache",
        "",
        "Compares fine-mapping of the **same 21 GW-significant leads** "
        "with and without the rice multi-omics graph prior loaded.",
        "",
        f"- Graph cache variants: {len(cache):,}",
        f"- Variants cached with pathway or PPI info (the 'prior-rich' subset): "
        f"{sum(1 for v in cache.values() if v['pathways'] or v['ppi']):,}",
        "",
        "## Per-locus change (CS size, top PIP)",
        "",
        "| Locus | Prior-rich vars | L1 CS before→after | L1 PIP before→after | HBP CS before→after | HBP PIP before→after |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for _, r in merged.iterrows():
        lines.append(
            f"| {r['locus_id']} | {r['var_with_pathway_or_ppi']} | "
            f"{r['L1_cs_baseline']} → **{r['L1_cs_size_with_graph']}** "
            f"({r['L1_cs_reduction_pct']:+.0f}%) | "
            f"{r['L1_pip_baseline']:.4f} → **{r['L1_top_pip_with_graph']:.4f}** "
            f"({r['L1_pip_boost']:+.4f}) | "
            f"{r['HBP_cs_baseline']} → **{r['HBP_cs_size_with_graph']}** "
            f"({r['HBP_cs_reduction_pct']:+.0f}%) | "
            f"{r['HBP_pip_baseline']:.4f} → **{r['HBP_top_pip_with_graph']:.4f}** "
            f"({r['HBP_pip_boost']:+.4f}) |"
        )
    lines.extend([
        "",
        "## Aggregate",
        "",
        f"- Mean L1 CS reduction: **{merged['L1_cs_reduction_pct'].mean():.1f}%**",
        f"- Mean HBP CS reduction: **{merged['HBP_cs_reduction_pct'].mean():.1f}%**",
        f"- Mean L1 top-PIP boost: **{merged['L1_pip_boost'].mean():+.4f}**",
        f"- Mean HBP top-PIP boost: **{merged['HBP_pip_boost'].mean():+.4f}**",
    ])
    md = RES / "rice3k_graph_vs_baseline.md"
    md.write_text("\n".join(lines))
    print(f"  Wrote {md}")


if __name__ == "__main__":
    main()
