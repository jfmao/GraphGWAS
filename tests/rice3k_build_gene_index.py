"""Build a LOC_Os → (chrom, pos) index from the VCF's snpEff ANN field.

Single streaming pass through the gzipped VCF body. For each LOC_Os
ID we keep the first (smallest-POS) variant annotated in/near that
gene — an approximation of the gene's genomic start.

Writes:
  data/rice_3k/ground_truth/gene_position_index.tsv
"""
from __future__ import annotations

import gzip
import re
from pathlib import Path

VCF = Path("/mnt/data/GraphPop/data/raw/3kRG_data/NB_bialSNP_pseudo_canonical_ALL.vcf.gz")
# Note: this file has snpEff ANN annotations; chromosome labels are
# numeric 1..12, not Chr1..Chr12 — we normalise to the `Chr` prefix so
# the index matches our ground-truth TSV (which uses Chr1..Chr12).
OUT = Path("/mnt/data/GraphGWAS/data/rice_3k/ground_truth/gene_position_index.tsv")

gene_rx = re.compile(rb"LOC_Os\d+g\d+")
first_pos: dict[str, tuple[str, int]] = {}

n = 0
n_ann = 0
print(f"Streaming {VCF} (gzip)...", flush=True)

with gzip.open(VCF, "rb") as fh:
    for line in fh:
        if line.startswith(b"#"):
            continue
        n += 1
        if n % 2_000_000 == 0:
            print(f"  {n/1e6:.1f}M variants · {len(first_pos):,} genes indexed",
                  flush=True)

        # Very fast split — only need chrom, pos, info
        parts = line.split(b"\t", 8)
        if len(parts) < 8:
            continue
        chrom_raw = parts[0].decode()
        chrom = chrom_raw if chrom_raw.startswith("Chr") else f"Chr{chrom_raw}"
        pos = int(parts[1])
        info = parts[7]

        matches = gene_rx.findall(info)
        if not matches:
            continue
        n_ann += 1
        seen = set()
        for mb in matches:
            loc = mb.decode()
            if loc in seen:
                continue
            seen.add(loc)
            prev = first_pos.get(loc)
            if prev is None or pos < prev[1]:
                first_pos[loc] = (chrom, pos)

print(f"\nTotal variants: {n:,}")
print(f"Variants with gene annotation: {n_ann:,}")
print(f"Distinct LOC_Os IDs: {len(first_pos):,}")

OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w") as f:
    f.write("LOC_Os_id\tchr\tpos\n")
    for loc, (c, p) in sorted(first_pos.items()):
        f.write(f"{loc}\t{c}\t{p}\n")
print(f"Wrote {OUT}")
