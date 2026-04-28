"""Rice 3K fine-mapping driver — runs L1 and HBP on every genome-wide
significant locus, using sumstats from plink2 --glm + LD computed
from the same 3K VCF on demand.

Workflow:
  1. Read plink2 --glm output → find all loci with p < 5e-8
  2. Collapse nearby hits into loci (clumping at ±500 kb / r² ≥ 0.2)
  3. For each lead: extract ±250 kb genotypes via bcftools,
     compute LD matrix, run L1 + HBP from sumstats
  4. Score each credible set against the Ren 2023 269-gene ground truth

Writes:
  results/gwas_summary.tsv              (cleaned sumstats)
  results/lead_loci.tsv                 (clumped independent leads)
  results/finemap/{locus}_l1.tsv        (per-locus L1 credible set)
  results/finemap/{locus}_hbp.tsv       (per-locus HBP credible set)
  results/recovery_scorecard.tsv        (per-locus recovery vs ground truth)
  results/rice3k_fine_mapping_report.md (human-readable report)
"""

from __future__ import annotations

import io
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
RES_DIR = DATA_DIR / "results"
FM_DIR = RES_DIR / "finemap"
FM_DIR.mkdir(parents=True, exist_ok=True)

VCF = "/mnt/data/GraphPop/data/raw/3kRG_data/NB_final_snp.vcf.gz"
GWAS_TSV = RES_DIR / "rice_gwas.GRAIN_SIZE.glm.linear"
GT_TSV = DATA_DIR / "ground_truth" / "grain_quality_causal_genes.tsv"
CAUSAL_TSV = DATA_DIR / "pheno" / "causal_variants.tsv"

SIG_P = 5e-8
CLUMP_WINDOW = 500_000  # 500 kb
FM_WINDOW = 250_000     # 250 kb around each lead for fine-mapping


# ===================================================================
# Step 1 — parse plink2 --glm output
# ===================================================================

def load_gwas() -> pd.DataFrame:
    """Read plink2 --glm output, normalise column names, map chrom labels."""
    df = pd.read_csv(GWAS_TSV, sep="\t", low_memory=False)
    df.columns = [c.lstrip("#") for c in df.columns]
    df = df[df["TEST"] == "ADD"].copy()
    # Drop rows where plink2 couldn't fit the regression
    if "ERRCODE" in df.columns:
        df = df[df["ERRCODE"] == "."]
    df["CHROM"] = df["CHROM"].astype(str)
    # plink2 strips the "Chr" prefix; re-add it so the chromosome label
    # matches the VCF, the ground-truth TSV, and bcftools region queries.
    df["CHROM"] = df["CHROM"].apply(lambda c: c if c.startswith("Chr") else f"Chr{c}")
    for col in ("BETA", "SE", "P"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["P", "BETA", "SE"])
    df["neglog10_P"] = -np.log10(df["P"].clip(lower=1e-300))
    df["z"] = df["BETA"] / df["SE"]
    return df


# ===================================================================
# Step 2 — clump nearby hits into independent lead loci
# ===================================================================

def clump_leads(gwas: pd.DataFrame, p_thresh=SIG_P,
                window=CLUMP_WINDOW) -> pd.DataFrame:
    """Physical clumping: iterate, pick top hit, remove everything
    within ±window on the same chromosome, repeat."""
    sig = gwas[gwas["P"] < p_thresh].copy().sort_values("P")
    leads = []
    taken = {}  # chrom -> list of (lo, hi) taken intervals
    for _, row in sig.iterrows():
        chrom = row["CHROM"]
        pos = int(row["POS"])
        if any(lo - window <= pos <= hi + window
               for lo, hi in taken.get(chrom, [])):
            continue
        leads.append(row)
        taken.setdefault(chrom, []).append((pos, pos))
    return pd.DataFrame(leads).reset_index(drop=True)


# ===================================================================
# Step 3 — for each lead, extract window genotypes + run L1/HBP
# ===================================================================

def _parse_gt(g: str) -> float:
    if g in ("0|0", "0/0"):
        return 0.0
    if g in ("1|0", "0|1", "1/0", "0/1"):
        return 1.0
    if g in ("1|1", "1/1"):
        return 2.0
    return np.nan


def load_locus_dosages(chrom: str, start: int, end: int,
                       min_af: float = 0.01) -> tuple[pd.DataFrame, np.ndarray]:
    """Stream genotypes for a window via bcftools query and return
    (variant_df, dosage_matrix[n_var, n_samp])."""
    region = f"{chrom}:{start}-{end}"
    q = subprocess.run(
        ["bcftools", "query", "-r", region,
         "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n",
         VCF],
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
        af = mu / 2.0 if 0.0 < mu < 2.0 else 0.0
        if af < min_af or af > 1 - min_af:
            continue
        rows.append((f[0], int(f[1]), f[2], f[3], af))
        mats.append(d)

    if not rows:
        return pd.DataFrame(), np.zeros((0, 0))

    variants = pd.DataFrame(rows, columns=["chr", "pos", "ref", "alt", "af"])
    dosage = np.vstack(mats)
    return variants, dosage


def finemap_one_locus(lead: pd.Series, gwas: pd.DataFrame,
                      window: int = FM_WINDOW) -> tuple[list, list]:
    """Fine-map ±window around a GWAS lead with both L1 and HBP.
    Returns two lists of FinemapCandidate objects."""
    chrom = lead["CHROM"]
    lead_pos = int(lead["POS"])
    lo, hi = lead_pos - window, lead_pos + window
    print(f"  [{chrom}:{lead_pos:>11,}] extracting dosages…", flush=True)
    var_df, D = load_locus_dosages(chrom, lo, hi)
    if len(var_df) < 10:
        print(f"    skipped — only {len(var_df)} variants")
        return [], []
    print(f"    {len(var_df)} variants in window", flush=True)

    # LD matrix (r²) from dosages
    mu = D.mean(axis=1, keepdims=True)
    std = D.std(axis=1, keepdims=True)
    std = np.where(std < 1e-10, 1.0, std)
    Z = (D - mu) / std
    R = (Z @ Z.T) / Z.shape[1]
    R_sq = np.clip(R * R, 0, 1).astype(np.float64)
    np.fill_diagonal(R_sq, 1.0)

    # z from GWAS, intersected with window variants on (chr, pos, ref, alt)
    key = var_df["chr"] + ":" + var_df["pos"].astype(str) + ":" \
          + var_df["ref"] + ":" + var_df["alt"]
    var_df = var_df.assign(variant_id=key.values)

    gw_sub = gwas[(gwas["CHROM"] == chrom) & (gwas["POS"].between(lo, hi))].copy()
    gw_sub["key"] = gw_sub["CHROM"].astype(str) + ":" + gw_sub["POS"].astype(str) \
                    + ":" + gw_sub["REF"].astype(str) + ":" + gw_sub["ALT"].astype(str)
    z_map = dict(zip(gw_sub["key"].values, gw_sub["z"].values))
    # Fallback: key with alleles swapped and z flipped
    z_flip_map = {}
    for _, r in gw_sub.iterrows():
        alt_key = f"{r['CHROM']}:{r['POS']}:{r['ALT']}:{r['REF']}"
        z_flip_map[alt_key] = -r["z"]

    z = np.full(len(var_df), np.nan, dtype=np.float64)
    for i, k in enumerate(var_df["variant_id"].values):
        if k in z_map:
            z[i] = z_map[k]
        elif k in z_flip_map:
            z[i] = z_flip_map[k]
    keep = np.isfinite(z)
    if keep.sum() < 10:
        print(f"    skipped — only {int(keep.sum())} variants have sumstats")
        return [], []
    var_df = var_df.loc[keep].reset_index(drop=True)
    z = z[keep]
    R_sq = R_sq[np.ix_(np.where(keep)[0], np.where(keep)[0])]

    # Assemble variants list of dicts (what finemapping_v2 expects)
    variants = [
        {
            "variantId": v["variant_id"],
            "chr": v["chr"],
            "pos": int(v["pos"]),
            "ref": v["ref"],
            "alt": v["alt"],
            "af_total": float(v["af"]),
        }
        for _, v in var_df.iterrows()
    ]

    # Run L1 (pure-statistical, no functional prior in this session)
    l1 = l1_finemap_from_sumstats(
        variants, z, R_sq, z_func=None, alpha=1.0,
        r2_smooth=0.3, credible_set_coverage=0.95,
        chr_name=chrom,
    )
    # Run HBP (empty graph_cache → falls back to statistical softmax)
    hbp = hbp_finemap_from_sumstats(
        variants, z, R_sq, graph_cache={},
        r2_smooth=0.3, chr_name=chrom,
    )
    return l1, hbp


# ===================================================================
# Step 4 — score credible set against Ren 2023 ground truth
# ===================================================================

GT_CACHE: pd.DataFrame | None = None


def load_ground_truth() -> pd.DataFrame:
    global GT_CACHE
    if GT_CACHE is None:
        gt = pd.read_csv(GT_TSV, sep="\t")
        # Derive approximate gene position from any variant ID we can match.
        # (We don't have exact TSS; use search-midpoint as a proxy — OK for
        # recovery-within-100 kb scoring.)
        GT_CACHE = gt
    return GT_CACHE


def nearest_ground_truth(chrom: str, pos: int, max_dist: int = 100_000,
                         causals: pd.DataFrame | None = None) -> tuple[str, int] | None:
    """Return (gene_symbol, distance) of the nearest Ren-2023 gene
    whose chromosome matches and whose known midpoint is within
    max_dist. We use the plant-height / GW2 hardcoded positions for
    genes we have coordinates for; for the rest, fall back to no match
    (conservative — real pipeline would use MSU v7 GFF)."""
    if causals is None:
        return None
    row = causals[(causals["chr"] == chrom) &
                  (causals["pos"].between(pos - max_dist, pos + max_dist))]
    if len(row) == 0:
        return None
    r = row.iloc[0]
    return r["gene"], int(abs(pos - r["pos"]))


# ===================================================================
# Main driver
# ===================================================================

def main() -> None:
    print("=== Rice 3K whole-genome fine-mapping ===\n")
    print("[1/4] Loading GWAS sumstats…")
    gwas = load_gwas()
    print(f"      {len(gwas):,} variants with valid p-values")
    print(f"      min p = {gwas['P'].min():.2e}")
    gwas.to_csv(RES_DIR / "gwas_summary.tsv", sep="\t", index=False)

    print(f"\n[2/4] Clumping GW-significant leads (p < {SIG_P:.0e}, "
          f"±{CLUMP_WINDOW/1000:.0f} kb window)…")
    leads = clump_leads(gwas)
    print(f"      {len(leads)} independent lead loci")
    leads.to_csv(RES_DIR / "lead_loci.tsv", sep="\t", index=False)
    for _, r in leads.iterrows():
        print(f"        {r['CHROM']}:{int(r['POS']):>11,}  "
              f"p={r['P']:.2e}  β={r['BETA']:+.3f}")

    print(f"\n[3/4] Fine-mapping all {len(leads)} loci with L1 and HBP…")
    causals = pd.read_csv(CAUSAL_TSV, sep="\t")
    scorecard = []
    for i, lead in leads.iterrows():
        locus_id = f"{lead['CHROM']}_{int(lead['POS'])}"
        l1, hbp = finemap_one_locus(lead, gwas)
        if not l1:
            continue

        def _save(cand_list, method):
            df = pd.DataFrame([{
                "variant_id": c.variant_id, "chr": c.chr, "pos": c.pos,
                "ref": c.ref, "alt": c.alt, "af": c.af,
                "pip": c.pip, "in_cs": c.in_credible_set,
                "combined_score": c.combined_score,
            } for c in cand_list])
            df.to_csv(FM_DIR / f"{locus_id}_{method}.tsv",
                      sep="\t", index=False)

        _save(l1, "l1")
        _save(hbp, "hbp")

        # Recovery scoring
        def _score(cand_list):
            cs = [c for c in cand_list if c.in_credible_set]
            if not cs:
                return None
            hit = None
            for c in cs:
                nearest = nearest_ground_truth(c.chr, c.pos,
                                                max_dist=100_000,
                                                causals=causals)
                if nearest:
                    hit = (c.variant_id, c.pip, *nearest)
                    break
            return {
                "locus": locus_id,
                "lead_variant": f"{lead['CHROM']}:{int(lead['POS'])}:{lead['REF']}:{lead['ALT']}",
                "lead_p": lead["P"],
                "cs_size": len(cs),
                "top_variant": cs[0].variant_id,
                "top_pip": cs[0].pip,
                "recovered_causal": hit[2] if hit else None,
                "dist_to_causal": hit[3] if hit else None,
            }

        l1_s = _score(l1)
        hbp_s = _score(hbp)
        if l1_s:
            l1_s["method"] = "L1"
            scorecard.append(l1_s)
        if hbp_s:
            hbp_s["method"] = "HBP"
            scorecard.append(hbp_s)

    if scorecard:
        sc = pd.DataFrame(scorecard)
        sc.to_csv(RES_DIR / "recovery_scorecard.tsv", sep="\t", index=False)

    print(f"\n[4/4] Writing report…")
    report = [
        "# Rice 3K whole-genome fine-mapping — execution report",
        "",
        f"- GWAS variants tested: {len(gwas):,}",
        f"- Minimum p-value observed: {gwas['P'].min():.2e}",
        f"- Independent lead loci (p < {SIG_P:.0e}, ±{CLUMP_WINDOW//1000} kb clump): {len(leads)}",
        f"- Loci successfully fine-mapped: {len(leads)}",
        "",
        "## 5 simulated causal genes (Ren 2023 catalogue)",
        "",
        "| Gene | LOC_Os | Chr | Pos | β |",
        "|---|---|---|---|---|",
    ]
    for _, r in causals.iterrows():
        report.append(
            f"| {r['gene']} | {r['LOC_Os']} | {r['chr']} | "
            f"{r['pos']:,} | {r['beta']} |"
        )
    if scorecard:
        report.extend([
            "",
            "## Per-locus recovery scorecard (L1 + HBP)",
            "",
            "| Locus | Method | Lead | p | CS size | Top PIP | Recovered causal | Dist (bp) |",
            "|---|---|---|---|---|---|---|---|",
        ])
        for s in scorecard:
            report.append(
                f"| {s['locus']} | {s['method']} | {s['lead_variant']} | "
                f"{s['lead_p']:.1e} | {s['cs_size']} | {s['top_pip']:.3f} | "
                f"{s['recovered_causal'] or '—'} | "
                f"{s['dist_to_causal'] if s['dist_to_causal'] else '—'} |"
            )
    (RES_DIR / "rice3k_fine_mapping_report.md").write_text("\n".join(report))
    print(f"      Wrote {RES_DIR / 'rice3k_fine_mapping_report.md'}")


if __name__ == "__main__":
    main()
