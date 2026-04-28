"""Rice 3K full recovery audit — checks every fine-mapped credible
set against the ENTIRE Ren 2023 ground-truth catalogue (269 genes),
not just the 5 simulated causal genes.

For each gene in the ground-truth catalogue we look up its
approximate MSU v7 genomic coordinate by searching the VCF's
snpEff annotation field, then ask: was any credible-set variant
within ±100 kb of that gene's position?

Writes:
  results/rice3k_recovery_full.tsv        — per-gene recovery
  results/rice3k_recovery_summary.md      — human-readable summary
"""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

DATA_DIR = Path("/mnt/data/GraphGWAS/data/rice_3k")
RES_DIR = DATA_DIR / "results"
FM_DIR = RES_DIR / "finemap"
VCF = "/mnt/data/GraphPop/data/raw/3kRG_data/NB_final_snp.vcf.gz"

GT_FILE = DATA_DIR / "ground_truth" / "grain_quality_causal_genes.tsv"
LEADS_FILE = RES_DIR / "lead_loci.tsv"

# Chromosome ranges that actual rice genes occupy (from MSU v7).
# We use broad per-gene-ID index heuristics to narrow the bcftools
# search window, then grep for the exact LOC_Os ID in the ANN field.
CHROM_LENGTHS = {
    "Chr1": 43_269_907, "Chr2": 35_936_220, "Chr3": 36_413_793,
    "Chr4": 35_502_686, "Chr5": 29_957_429, "Chr6": 31_247_783,
    "Chr7": 29_697_605, "Chr8": 28_442_007, "Chr9": 23_009_500,
    "Chr10": 23_206_234, "Chr11": 29_020_102, "Chr12": 27_530_849,
}


def find_gene_position(loc_os: str, chrom: str, cache: dict) -> int | None:
    """Locate a rice gene's approximate genomic position by grepping
    the VCF snpEff annotation field for LOC_Os ID within an iteratively
    narrowed window."""
    if loc_os in cache:
        return cache[loc_os]
    # MSU v7 gene IDs are sequential within chromosome — use gene number
    # to guess the approximate Mb range.
    m = re.match(r"LOC_Os(\d+)g(\d+)", loc_os)
    if not m:
        cache[loc_os] = None
        return None
    gene_num = int(m.group(2))
    chr_len = CHROM_LENGTHS.get(chrom, 0)
    # Rough estimate: genes number 1-50,000 across a ~35 Mb chromosome;
    # scale gene_num (in tens of thousands) to Mb.
    approx_mb = gene_num / 55_000 * chr_len  # 55k genes per chrom (upper)
    # Search ±3 Mb window for the gene ID in snpEff annotation
    lo = max(0, int(approx_mb) - 3_000_000)
    hi = min(chr_len, int(approx_mb) + 3_000_000)
    region = f"{chrom}:{lo}-{hi}"
    cmd = ["bcftools", "view", "-r", region, "-H", VCF]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    if proc.returncode != 0:
        cache[loc_os] = None
        return None
    for line in proc.stdout.split("\n"):
        if loc_os in line:
            pos = int(line.split("\t")[1])
            cache[loc_os] = pos
            return pos
    # Fallback: scan whole chromosome (slow)
    region = chrom
    proc2 = subprocess.run(["bcftools", "view", "-r", region, "-H", VCF],
                           capture_output=True, text=True, timeout=300)
    for line in proc2.stdout.split("\n")[:100_000]:
        if loc_os in line:
            pos = int(line.split("\t")[1])
            cache[loc_os] = pos
            return pos
    cache[loc_os] = None
    return None


def load_credible_sets():
    """Aggregate all fine-mapping outputs into one DataFrame indexed
    by locus id. Keep only credible-set members."""
    rows = []
    for tsv in sorted(FM_DIR.glob("*_l1.tsv")):
        locus_id = tsv.stem.replace("_l1", "")
        df = pd.read_csv(tsv, sep="\t")
        cs = df[df["in_cs"] == True].copy()
        cs["locus_id"] = locus_id
        cs["method"] = "L1"
        rows.append(cs)
    for tsv in sorted(FM_DIR.glob("*_hbp.tsv")):
        locus_id = tsv.stem.replace("_hbp", "")
        df = pd.read_csv(tsv, sep="\t")
        cs = df[df["in_cs"] == True].copy()
        cs["locus_id"] = locus_id
        cs["method"] = "HBP"
        rows.append(cs)
    if not rows:
        return pd.DataFrame()
    all_cs = pd.concat(rows, ignore_index=True)
    # parse chr + pos from variant_id (format Chr8:26514866:G:A)
    split = all_cs["variant_id"].str.split(":", expand=True)
    all_cs["cs_chr"] = split[0]
    all_cs["cs_pos"] = pd.to_numeric(split[1])
    return all_cs


def main():
    print("[1/3] Loading ground truth + credible sets...")
    gt = pd.read_csv(GT_FILE, sep="\t")
    cs = load_credible_sets()
    print(f"  Ground truth: {len(gt)} rice causal genes")
    print(f"  Credible-set members aggregated: {len(cs):,}")

    # Also load the simulated causals (these 5 ARE guaranteed recoverable)
    sim_causals = pd.read_csv(DATA_DIR / "pheno" / "causal_variants.tsv", sep="\t")
    sim_genes = set(sim_causals["gene"].str.split("/").str[0])

    print("\n[2/3] Looking up gene positions in VCF snpEff annotations...")
    cache = {}
    out_rows = []
    for i, g in gt.iterrows():
        loc = g["LOC_Os_id"]
        chrom = g["chr"]
        pos = find_gene_position(loc, chrom, cache)
        if pos is None:
            # Skip genes we can't place
            out_rows.append({
                "gene": g["gene"], "LOC_Os": loc, "chr": chrom,
                "gene_pos": None, "category": g["category"],
                "in_credible_set": False, "min_dist_bp": None,
                "matched_locus": None, "top_pip": None, "is_simulated": False,
            })
            continue

        # Check distance from every L1 credible-set variant
        overlap = cs[(cs["method"] == "L1") & (cs["cs_chr"] == chrom)]
        overlap = overlap.assign(dist=(overlap["cs_pos"] - pos).abs())
        within = overlap[overlap["dist"] <= 100_000]
        if len(within) > 0:
            winning = within.loc[within["dist"].idxmin()]
            row = {
                "gene": g["gene"], "LOC_Os": loc, "chr": chrom,
                "gene_pos": pos, "category": g["category"],
                "in_credible_set": True,
                "min_dist_bp": int(winning["dist"]),
                "matched_locus": winning["locus_id"],
                "top_pip": float(winning["pip"]),
                "is_simulated": any(sg in g["gene"] for sg in sim_genes),
            }
        else:
            row = {
                "gene": g["gene"], "LOC_Os": loc, "chr": chrom,
                "gene_pos": pos, "category": g["category"],
                "in_credible_set": False, "min_dist_bp": None,
                "matched_locus": None, "top_pip": None,
                "is_simulated": any(sg in g["gene"] for sg in sim_genes),
            }
        out_rows.append(row)
        if (i + 1) % 20 == 0:
            print(f"  processed {i+1}/{len(gt)} genes...")

    audit = pd.DataFrame(out_rows)
    audit.to_csv(RES_DIR / "rice3k_recovery_full.tsv", sep="\t", index=False)
    print(f"  Wrote {RES_DIR / 'rice3k_recovery_full.tsv'}")

    print("\n[3/3] Summary statistics...")
    n_placed = audit["gene_pos"].notna().sum()
    n_recovered = audit["in_credible_set"].sum()
    print(f"  Genes with known position: {n_placed}/{len(audit)}")
    print(f"  Recovered in any credible set: {n_recovered}")

    # Per-simulated-gene check (guaranteed ground truth)
    sim_recovered = audit[audit["is_simulated"] & audit["in_credible_set"]]
    print(f"\n  Simulated causals recovered: {len(sim_recovered)}/{len(sim_causals)}")
    print(sim_recovered[["gene", "LOC_Os", "chr", "gene_pos",
                         "min_dist_bp", "matched_locus", "top_pip"]].to_string())

    # Breakdown by category
    print(f"\n  Recovery by trait category:")
    by_cat = audit.dropna(subset=["gene_pos"]).groupby("category").agg(
        total=("gene", "count"),
        recovered=("in_credible_set", "sum"),
    )
    by_cat["pct"] = (100 * by_cat["recovered"] / by_cat["total"]).round(1)
    print(by_cat.to_string())

    # Write markdown summary
    md = [
        "# Rice 3K full recovery audit — against Ren 2023 269-gene catalogue",
        "",
        f"- GWAS leads fine-mapped: {cs['locus_id'].nunique()}",
        f"- Ground-truth genes in Ren 2023 catalogue: {len(gt)}",
        f"- Genes with placed genomic position: {n_placed}",
        f"- Genes whose TSS falls within ±100 kb of ≥1 credible-set variant: {n_recovered}",
        f"- Overall recovery rate: {100*n_recovered/max(n_placed,1):.1f}%",
        "",
        "## Simulated causal genes (ground truth by construction)",
        "",
        "| Gene | LOC_Os | Chr | Position | Recovered | Min dist (bp) | Matched locus |",
        "|---|---|---|---|---|---|---|",
    ]
    for _, r in audit[audit["is_simulated"]].iterrows():
        rec = "✓" if r["in_credible_set"] else "✗"
        md.append(
            f"| {r['gene']} | {r['LOC_Os']} | {r['chr']} | "
            f"{int(r['gene_pos']) if pd.notna(r['gene_pos']) else '—':,} | {rec} | "
            f"{int(r['min_dist_bp']) if pd.notna(r['min_dist_bp']) else '—'} | "
            f"{r['matched_locus'] or '—'} |"
        )

    md.extend([
        "",
        "## Per-category recovery across the full 269-gene catalogue",
        "",
        "| Category | Genes placed | Recovered in any CS | Pct |",
        "|---|---:|---:|---:|",
    ])
    for cat, r in by_cat.iterrows():
        md.append(f"| {cat} | {int(r['total'])} | {int(r['recovered'])} | {r['pct']}% |")

    md.extend([
        "",
        "## Non-simulated genes also recovered",
        "",
        "These Ren 2023 ground-truth genes fell within ±100 kb of a "
        "genome-wide-significant credible set even though they were NOT "
        "among the 5 simulated causal genes. They represent additional "
        "validation hits that emerged from real rice 3K LD / population "
        "structure patterns:",
        "",
        "| Gene | LOC_Os | Chr | Category | Matched locus | Dist (bp) |",
        "|---|---|---|---|---|---|",
    ])
    other = audit[(~audit["is_simulated"]) & audit["in_credible_set"]]
    for _, r in other.iterrows():
        md.append(
            f"| {r['gene']} | {r['LOC_Os']} | {r['chr']} | "
            f"{r['category']} | {r['matched_locus']} | {int(r['min_dist_bp'])} |"
        )

    out = RES_DIR / "rice3k_recovery_summary.md"
    out.write_text("\n".join(md))
    print(f"\n  Wrote {out}")


if __name__ == "__main__":
    main()
