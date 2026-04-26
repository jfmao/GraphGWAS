"""Phase 3 — augment Arabidopsis graph_cache with prior_score + DAP-seq peak overlap.

Inputs:
  - arabidopsis_graph_cache_chr<N>.json (v2 with proper PPI from UniProt fix)
  - DAP-seq peaks (O'Malley 2016) if downloaded; else skip
  - TAIR10 GFF (for coding/promoter prior_score)

Adds:
  - prior_score = 1.0 if in DAP-seq peak (any TF) + 0.5 base if in coding region
                  + 0.3 if in pathway-annotated gene's promoter window (±2 kb)
  - tf_binding (list of TF names whose DAP-seq peak overlaps the variant) —
    optional layer when DAP-seq data is available

Outputs:
  arabidopsis_graph_cache_v3_chr<N>.json
"""
from __future__ import annotations

import bisect
import glob
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path

DATA = Path("/mnt/data/GraphGWAS/data/arabidopsis")
ANN = DATA / "annotations"
TAIR_GFF = ANN / "tair10.gff3.gz"
DAP_DIR = ANN / "dap_seq_peaks"   # if present, BED files per TF
NCBI_TO_CHROM = {
    "NC_003070.9": "1", "NC_003071.7": "2", "NC_003074.8": "3",
    "NC_003075.7": "4", "NC_003076.8": "5",
}


def parse_gff_genes_with_promoters() -> dict[str, list[tuple[int,int,str,str]]]:
    """{chrom: [(start_with_promoter, end, AT_id, ftype)] sorted}.
    Promoter = ±2 kb upstream of gene start.
    """
    AT_RE = re.compile(r"AT[1-5MC]G\d+")
    out: dict[str, list] = defaultdict(list)
    with gzip.open(TAIR_GFF, "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9: continue
            chrom_n = f[0]
            if chrom_n not in NCBI_TO_CHROM: continue
            ftype = f[2]
            if ftype not in ("gene", "exon", "CDS"): continue
            try:
                start, end = int(f[3]), int(f[4])
            except ValueError:
                continue
            m = AT_RE.search(f[8])
            if not m: continue
            chrom = NCBI_TO_CHROM[chrom_n]
            # Add 2 kb promoter window for "gene" rows
            ext_start = max(1, start - 2000) if ftype == "gene" else start
            out[chrom].append((ext_start, end, m.group(0), ftype))
    for c in out:
        out[c].sort()
    return dict(out)


def load_dap_peaks() -> dict[str, list[tuple[int,int,str]]]:
    """{chrom: [(start, end, TF_id)] sorted}. Empty if no DAP-seq downloaded."""
    out: dict[str, list] = defaultdict(list)
    if not DAP_DIR.exists():
        return {}
    n = 0
    for fp in sorted(DAP_DIR.glob("*.bed")):
        tf = fp.stem.split("_")[0]  # convention: <TF>_<replicate>.bed
        with fp.open() as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 3: continue
                chrom = f[0].lstrip("Cc").lstrip("h")
                if chrom not in {"1","2","3","4","5"}: continue
                try:
                    start, end = int(f[1]), int(f[2])
                except ValueError:
                    continue
                out[chrom].append((start, end, tf))
                n += 1
    for c in out:
        out[c].sort()
    print(f"  DAP-seq: {n:,} peaks across {len(out)} chroms")
    return dict(out)


def main():
    print("Parsing TAIR10 GFF for gene + promoter ranges ...", flush=True)
    chrom_genes = parse_gff_genes_with_promoters()
    for c in sorted(chrom_genes):
        print(f"  chr{c}: {len(chrom_genes[c]):,} feature ranges")

    print("\nLoading DAP-seq peaks (optional) ...", flush=True)
    dap = load_dap_peaks()

    print("\nAugmenting per-chrom caches ...", flush=True)
    for chrom in "12345":
        in_path = DATA / f"arabidopsis_graph_cache_chr{chrom}.json"
        if not in_path.exists():
            print(f"  [chr{chrom}] missing v2 cache, skip"); continue
        cache = json.load(in_path.open())
        gs = chrom_genes.get(chrom, [])
        starts = [g[0] for g in gs]
        dap_chr = dap.get(chrom, [])
        dap_starts = [p[0] for p in dap_chr]
        n_with_prior = 0
        n_with_tf = 0
        for vid, info in cache.items():
            try:
                pos = int(vid.split(":", 2)[1])
            except (IndexError, ValueError):
                continue
            # prior_score: max over feature overlap
            score = 0.0
            i0 = bisect.bisect_right(starts, pos) - 1
            for j in range(max(0, i0-3), min(len(gs), i0+2)):
                s, e, gid, ft = gs[j]
                if s <= pos <= e:
                    if ft == "CDS": score = max(score, 1.0)
                    elif ft == "exon": score = max(score, 0.8)
                    elif ft == "gene": score = max(score, 0.5)
            # DAP-seq overlap → boost score and add tf_binding
            if dap_chr:
                k0 = bisect.bisect_right(dap_starts, pos) - 1
                tfs = []
                for j in range(max(0, k0-3), min(len(dap_chr), k0+2)):
                    s, e, tf = dap_chr[j]
                    if s <= pos <= e:
                        tfs.append(tf)
                if tfs:
                    score = max(score, 1.5)
                    info["tf_binding"] = list(set(tfs))[:10]
                    n_with_tf += 1
            if score > 0:
                info["prior_score"] = float(score)
                n_with_prior += 1
        out = DATA / f"arabidopsis_graph_cache_v3_chr{chrom}.json"
        with out.open("w") as fh:
            json.dump(cache, fh)
        size_mb = out.stat().st_size / 1e6
        n_total = len(cache)
        print(f"  chr{chrom}: {n_with_prior:,}/{n_total:,} prior_score "
              f"({100*n_with_prior/n_total:.0f}%); {n_with_tf:,} tf_binding; "
              f"wrote {out.name} ({size_mb:.0f} MB)")


if __name__ == "__main__":
    main()
