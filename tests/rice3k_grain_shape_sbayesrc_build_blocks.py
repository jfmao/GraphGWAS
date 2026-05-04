"""Build a SBayesRC block-reference file targeting the 41 fine-mapped grain loci.

SBayesRC's LDstep1 needs a block-ref TSV whose rows are non-overlapping
[chr, start_bp, end_bp] intervals. We don't need the full 752-block
genome-wide reference — only the blocks containing our 41 fine-mapped
leads. Collapsing 41 leads (across 4 traits) into ±250 kb intervals and
merging overlaps gives a much smaller block set tractable for the rice rerun.

Outputs:
  data/rice_3k/sbayesrc_rice/rice_targeted_blocks.txt
      (block-ref format: chr_int start_bp end_bp, one row per merged block)
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

DATA = Path("/mnt/data/GraphGWAS/data/rice_3k")
SBAYES_DIR = DATA / "sbayesrc_rice"
SBAYES_DIR.mkdir(parents=True, exist_ok=True)
WINDOW_BP = 250_000  # ±250 kb around each lead

LEADS = []
fm_dir = DATA / "results" / "grain_finemap"
for jf in sorted(fm_dir.glob("*_Chr*_*.json")):
    parts = jf.stem.split("_")
    if len(parts) < 3:
        continue
    chrom = parts[1]
    pos = int(parts[2])
    LEADS.append({"chrom": chrom, "pos": pos})

leads = pd.DataFrame(LEADS).drop_duplicates(subset=["chrom", "pos"])
print(f"Unique leads (across traits): {len(leads)}")

# Merge overlapping intervals per chromosome
rows = []
for chrom, sub in leads.groupby("chrom"):
    chrom_n = int(chrom.replace("Chr", ""))
    intervals = sorted(
        (max(1, p - WINDOW_BP), p + WINDOW_BP) for p in sub["pos"]
    )
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    for s, e in merged:
        rows.append({"chr": chrom_n, "start": s, "end": e})

blocks = pd.DataFrame(rows).sort_values(["chr", "start"]).reset_index(drop=True)
print(f"Merged into {len(blocks)} non-overlapping blocks")
print(blocks.to_string(index=False))

# Write block-ref. SBayesRC's blockRef format: tab-separated chr start end
out = SBAYES_DIR / "rice_targeted_blocks.txt"
blocks.to_csv(out, sep="\t", index=False, header=False)
print(f"\nWrote {out}")
