"""Rice 3K recovery scoring — against the Ren 2023 269-gene catalogue.

Uses the precomputed gene-position index from the VCF's snpEff ANN
field, then asks for each Ren 2023 ground-truth gene: does any
credible-set variant (from our 21 fine-mapped loci) fall within
±100 kb of the gene's position?

Writes:
  results/rice3k_recovery_full.tsv     — one row per gene
  results/rice3k_recovery_summary.md   — human-readable report
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

DATA_DIR = Path("/mnt/data/GraphGWAS/data/rice_3k")
RES_DIR = DATA_DIR / "results"
FM_DIR = RES_DIR / "finemap"

GT = pd.read_csv(DATA_DIR / "ground_truth" / "grain_quality_causal_genes.tsv", sep="\t")
IDX = pd.read_csv(DATA_DIR / "ground_truth" / "gene_position_index.tsv", sep="\t")
SIM = pd.read_csv(DATA_DIR / "pheno" / "causal_variants.tsv", sep="\t")

print(f"Ground-truth genes (Ren 2023):      {len(GT)}")
print(f"Gene-position index from VCF:       {len(IDX):,}")
print(f"Simulated causal genes:              {len(SIM)}")

# Join ground truth with index to get positions
gt_placed = GT.merge(IDX, on="LOC_Os_id", how="left", suffixes=("", "_idx"))
n_placed = gt_placed["pos"].notna().sum()
print(f"Ren 2023 genes with placed position: {n_placed}/{len(GT)}")


# --- Load all credible-set members from the 21 fine-mapped loci ---
cs_rows = []
for tsv in sorted(FM_DIR.glob("*_l1.tsv")):
    locus_id = tsv.stem.replace("_l1", "")
    df = pd.read_csv(tsv, sep="\t")
    cs = df[df["in_cs"]].copy()
    cs["locus_id"] = locus_id
    cs["method"] = "L1"
    cs_rows.append(cs)
for tsv in sorted(FM_DIR.glob("*_hbp.tsv")):
    locus_id = tsv.stem.replace("_hbp", "")
    df = pd.read_csv(tsv, sep="\t")
    cs = df[df["in_cs"]].copy()
    cs["locus_id"] = locus_id
    cs["method"] = "HBP"
    cs_rows.append(cs)
all_cs = pd.concat(cs_rows, ignore_index=True)
split = all_cs["variant_id"].str.split(":", expand=True)
all_cs["cs_chr"] = split[0]
all_cs["cs_pos"] = pd.to_numeric(split[1])
print(f"Total credible-set variant-entries: {len(all_cs):,} "
      f"(L1: {(all_cs['method']=='L1').sum():,}  "
      f"HBP: {(all_cs['method']=='HBP').sum():,})")


# --- Score recovery for each Ren 2023 gene ---
# Match simulated causals by LOC_Os ID (robust; avoids
# separator and substring pitfalls with gene symbols).
sim_loc_ids = set(SIM["LOC_Os"])

rows = []
for _, g in gt_placed.iterrows():
    loc = g["LOC_Os_id"]
    chrom = g["chr"]
    pos = g["pos"]
    is_sim = loc in sim_loc_ids
    if pd.isna(pos):
        rows.append(dict(gene=g["gene"], LOC_Os=loc, chr=chrom,
                         gene_pos=None, category=g["category"],
                         in_CS_L1=False, min_dist_L1=None,
                         in_CS_HBP=False, min_dist_HBP=None,
                         locus_matched=None, is_simulated=is_sim))
        continue
    chrom_match = all_cs[all_cs["cs_chr"] == chrom].copy()
    chrom_match["d"] = (chrom_match["cs_pos"] - int(pos)).abs()
    near = chrom_match[chrom_match["d"] <= 100_000]
    l1_near = near[near["method"] == "L1"]
    hbp_near = near[near["method"] == "HBP"]
    min_l1 = int(l1_near["d"].min()) if len(l1_near) else None
    min_hbp = int(hbp_near["d"].min()) if len(hbp_near) else None
    matched_locus = None
    if len(near):
        matched_locus = near.sort_values("d").iloc[0]["locus_id"]
    rows.append(dict(
        gene=g["gene"], LOC_Os=loc, chr=chrom, gene_pos=int(pos),
        category=g["category"],
        in_CS_L1=bool(min_l1 is not None), min_dist_L1=min_l1,
        in_CS_HBP=bool(min_hbp is not None), min_dist_HBP=min_hbp,
        locus_matched=matched_locus, is_simulated=is_sim,
    ))

out = pd.DataFrame(rows)
out.to_csv(RES_DIR / "rice3k_recovery_full.tsv", sep="\t", index=False)
print(f"\nWrote {RES_DIR / 'rice3k_recovery_full.tsv'}")


# --- Summary ---
placed = out[out["gene_pos"].notna()].copy()
n_l1 = placed["in_CS_L1"].sum()
n_hbp = placed["in_CS_HBP"].sum()
print(f"\nOverall recovery (any of 21 GW-significant fine-mapped loci):")
print(f"  L1   in-CS hits within ±100 kb: {n_l1}/{len(placed)}"
      f" ({100*n_l1/max(len(placed),1):.1f}%)")
print(f"  HBP  in-CS hits within ±100 kb: {n_hbp}/{len(placed)}"
      f" ({100*n_hbp/max(len(placed),1):.1f}%)")

print(f"\nSimulated causals recovered:")
sim = placed[placed["is_simulated"]]
for _, r in sim.iterrows():
    print(f"  {r['gene']:20} chr={r['chr']:5} pos={r['gene_pos']:>10,}  "
          f"L1_dist={r['min_dist_L1']}   HBP_dist={r['min_dist_HBP']}   "
          f"locus={r['locus_matched']}")

print(f"\nBy trait category:")
by_cat = placed.groupby("category").agg(
    total=("gene", "count"),
    in_L1=("in_CS_L1", "sum"),
    in_HBP=("in_CS_HBP", "sum"),
)
by_cat["L1_pct"] = (100 * by_cat["in_L1"] / by_cat["total"]).round(1)
by_cat["HBP_pct"] = (100 * by_cat["in_HBP"] / by_cat["total"]).round(1)
print(by_cat.to_string())

# --- Markdown summary ---
md = [
    "# Rice 3K full recovery audit vs Ren 2023 269-gene catalogue",
    "",
    f"- Fine-mapped GW-significant loci: {all_cs['locus_id'].nunique()}",
    f"- Ground-truth genes (Ren 2023): {len(GT)}",
    f"- Genes placed by VCF gene-position index: {n_placed}",
    f"- L1 recovery (any CS within ±100 kb): **{n_l1}/{n_placed} "
    f"({100*n_l1/max(n_placed,1):.1f}%)**",
    f"- HBP recovery (any CS within ±100 kb): **{n_hbp}/{n_placed} "
    f"({100*n_hbp/max(n_placed,1):.1f}%)**",
    "",
    "## 5 simulated causals (guaranteed ground truth)",
    "",
    "| Gene | LOC_Os | Chr | TSS pos | L1 min-dist | HBP min-dist | Locus |",
    "|---|---|---|---:|---:|---:|---|",
]
for _, r in sim.iterrows():
    md.append(
        f"| {r['gene']} | {r['LOC_Os']} | {r['chr']} | "
        f"{r['gene_pos']:,} | "
        f"{(str(r['min_dist_L1'])+' bp') if r['min_dist_L1'] is not None else '—'} | "
        f"{(str(r['min_dist_HBP'])+' bp') if r['min_dist_HBP'] is not None else '—'} | "
        f"{r['locus_matched'] or '—'} |"
    )
md.extend([
    "",
    "## Per-category recovery across the full 269-gene catalogue",
    "",
    "| Category | Genes placed | In L1 CS | In HBP CS | L1 % | HBP % |",
    "|---|---:|---:|---:|---:|---:|",
])
for cat, r in by_cat.iterrows():
    md.append(f"| {cat} | {int(r['total'])} | "
              f"{int(r['in_L1'])} | {int(r['in_HBP'])} | "
              f"{r['L1_pct']}% | {r['HBP_pct']}% |")

md.extend([
    "",
    "## Non-simulated Ren 2023 genes that fell into a credible set",
    "",
    "These are genes the scan did NOT have as simulated causals but "
    "whose genomic position is within ±100 kb of a GW-significant "
    "fine-mapped lead. They represent additional co-located rice "
    "grain genes picked up incidentally by the real 3K LD and "
    "subpopulation-structure signal.",
    "",
    "| Gene | LOC_Os | Chr | Position | Category | L1 dist | HBP dist | Locus |",
    "|---|---|---|---:|---|---:|---:|---|",
])
other = placed[(~placed["is_simulated"]) & placed["in_CS_L1"]].sort_values("min_dist_L1")
for _, r in other.iterrows():
    md.append(
        f"| {r['gene']} | {r['LOC_Os']} | {r['chr']} | {r['gene_pos']:,} | "
        f"{r['category']} | {r['min_dist_L1']} | {r['min_dist_HBP']} | "
        f"{r['locus_matched']} |"
    )

(RES_DIR / "rice3k_recovery_summary.md").write_text("\n".join(md))
print(f"\nWrote {RES_DIR / 'rice3k_recovery_summary.md'}")
