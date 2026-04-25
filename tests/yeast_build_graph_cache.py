"""Build a yeast multi-omics graph_cache (variant→gene→GO_slim→PPI).

Inputs (all already on disk under tests/data/yeast/):
  - SGD_features.tab            (gene + chromosome + start/end + strand)
  - go_slim_mapping.tab         (gene → GO_slim term, used as pathway proxy)
  - BIOGRID downloaded if needed (gene↔gene PPI; optional)
  - One yeast GWAS sumstats file (gives variant id format)

Output:
  data/yeast/yeast_graph_cache.json
  - {variantId: {"genes": [...], "pathways": [GO_slim_terms], "ppi": [partner_orfs]}}
"""
from __future__ import annotations

import gzip
import json
import re
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path
from urllib.request import urlretrieve

YEAST_DATA = Path("/mnt/data/GraphGWAS/tests/data/yeast")
OUT_DIR = Path("/mnt/data/GraphGWAS/data/yeast")
OUT_DIR.mkdir(parents=True, exist_ok=True)
GWAS_REF = Path("/mnt/data/GraphGWAS/results/grammar_corrected/gwas_YPETHANOL_grammar.tsv")
OUT_CACHE = OUT_DIR / "yeast_graph_cache.json"

SGD_FEATURES = YEAST_DATA / "SGD_features.tab"
GO_SLIM = YEAST_DATA / "go_slim_mapping.tab"
BIOGRID_OUT = OUT_DIR / "biogrid_yeast_559292.tab3.txt"
BIOGRID_URL = "https://thebiogrid.org/downloads/archives/Latest%20Release/BIOGRID-ORGANISM-LATEST.tab3.zip"

# yeast chromosome name normaliser: SGD says "chromosome 4" or "4"; GWAS uses "chromosome4"
def _norm_chr(s: str) -> str:
    s = str(s).strip()
    if s.startswith("chromosome "):
        s = "chromosome" + s.split(" ", 1)[1]
    elif s.isdigit():
        s = f"chromosome{s}"
    return s


def parse_sgd_features() -> tuple[dict, dict]:
    """Returns:
      gene_pos: {gene_orf: (chrom, start, end, strand)}
      pos_index: {chrom: [(start, end, gene_orf)] sorted}
    """
    gene_pos = {}
    pos_idx = defaultdict(list)
    seen_orfs = set()
    with SGD_FEATURES.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12: continue
            sgd_id, ftype, status, orf = f[0], f[1], f[2], f[3]
            chrom_raw = f[8]
            try:
                start = int(f[9]); end = int(f[10])
            except ValueError:
                continue
            if not orf or not chrom_raw: continue
            # ORF + tRNA + ncRNA are reasonable; skip CDS subfeatures
            if ftype not in ("ORF","tRNA_gene","ncRNA_gene","rRNA_gene",
                             "snRNA_gene","snoRNA_gene","pseudogene","transposable_element_gene"):
                continue
            if orf in seen_orfs:
                continue
            seen_orfs.add(orf)
            chrom = f"chromosome{chrom_raw.strip()}"
            lo, hi = min(start, end), max(start, end)
            gene_pos[orf] = (chrom, lo, hi, "+" if start <= end else "-")
            pos_idx[chrom].append((lo, hi, orf))
    for c in pos_idx:
        pos_idx[c].sort()
    return gene_pos, dict(pos_idx)


def parse_go_slim() -> dict[str, list[str]]:
    """Returns {gene_orf: [GO_slim_terms]}"""
    out = defaultdict(list)
    with GO_SLIM.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7: continue
            orf = f[0]
            term = f[4]
            ftype = f[6]
            if not orf.startswith("Y") or not term: continue
            if ftype not in ("ORF|Verified", "ORF|Uncharacterized"): continue
            if term not in out[orf]:
                out[orf].append(term)
    return dict(out)


def download_biogrid_ppi() -> dict[str, set[str]]:
    """Download BIOGRID, extract S. cerevisiae 559292 interactions."""
    sce_tab = OUT_DIR / "BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-LATEST.tab3.txt"
    if not sce_tab.exists():
        zip_path = OUT_DIR / "biogrid.zip"
        if not zip_path.exists():
            print(f"  Downloading BIOGRID ({BIOGRID_URL}) — may take a few min...",
                  flush=True)
            try:
                urlretrieve(BIOGRID_URL, str(zip_path))
            except Exception as e:
                print(f"  Download failed: {e}; skipping PPI layer.", flush=True)
                return {}
        # Extract just the cerevisiae S288c file
        try:
            r = subprocess.run(
                ["unzip", "-o", str(zip_path),
                 "BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-*.tab3.txt",
                 "-d", str(OUT_DIR)],
                capture_output=True, text=True, timeout=600,
            )
            # Move/symlink the extracted file
            import glob
            cands = glob.glob(str(OUT_DIR / "BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-*.tab3.txt"))
            if cands:
                Path(cands[0]).rename(sce_tab)
        except Exception as e:
            print(f"  Unzip failed: {e}; skipping PPI layer.", flush=True)
            return {}
    if not sce_tab.exists():
        return {}
    # Parse tab3 — columns 6,7 are systematic name (yeast ORF) of A and B
    ppi = defaultdict(set)
    with sce_tab.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            iA = header.index("Systematic Name Interactor A")
            iB = header.index("Systematic Name Interactor B")
            iEv = header.index("Experimental System Type")
        except ValueError:
            print("  BIOGRID header missing expected columns; skipping.")
            return {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(iA, iB, iEv): continue
            a, b, ev = f[iA], f[iB], f[iEv]
            if not (a.startswith("Y") and b.startswith("Y")): continue
            if ev != "physical": continue
            ppi[a].add(b)
            ppi[b].add(a)
    print(f"  PPI: {len(ppi):,} genes with ≥1 partner", flush=True)
    return dict(ppi)


def main():
    t0 = time.time()
    print("[1/4] Parsing SGD features ...", flush=True)
    gene_pos, pos_idx = parse_sgd_features()
    print(f"  {len(gene_pos):,} genes parsed across {len(pos_idx)} chromosomes")

    print("\n[2/4] Parsing GO slim mapping (pathway proxy) ...", flush=True)
    gene_pathways = parse_go_slim()
    n_with = sum(1 for v in gene_pathways.values() if v)
    print(f"  {n_with:,} genes with GO_slim annotations")

    print("\n[3/4] Loading BIOGRID PPI for S. cerevisiae 559292 ...", flush=True)
    gene_ppi = download_biogrid_ppi()

    # Variant -> overlapping gene(s)
    print("\n[4/4] Mapping variants in yeast GWAS sumstats to genes ...", flush=True)
    cache = {}
    n = 0
    n_annot = 0
    with GWAS_REF.open() as fh:
        next(fh)  # header
        for line in fh:
            f = line.split("\t")
            if len(f) < 4: continue
            n += 1
            if n % 200_000 == 0:
                print(f"  {n/1e6:.1f}M scanned, {len(cache):,} cached", flush=True)
            vid = f[0]
            chrom = f[1]
            try:
                pos = int(f[2])
            except ValueError:
                continue
            # Find overlapping genes (linear scan; yeast genes are sparse enough)
            overlaps = []
            for lo, hi, orf in pos_idx.get(chrom, []):
                if lo <= pos <= hi:
                    overlaps.append(orf)
                if lo > pos:
                    break  # sorted
            if not overlaps:
                continue
            n_annot += 1
            paths = []
            ppi = []
            for g in overlaps[:4]:
                for p in gene_pathways.get(g, []):
                    if p not in paths:
                        paths.append(p)
                if g in gene_ppi:
                    for partner in list(gene_ppi[g])[:20]:
                        if partner not in ppi:
                            ppi.append(partner)
            if not paths and not ppi:
                # Skip variants with only gene info — cache bloat for HBP
                continue
            cache[vid] = {"genes": overlaps[:4], "pathways": paths, "ppi": ppi}

    elapsed = time.time() - t0
    print(f"\nScanned {n:,} variants; {n_annot:,} hit a gene; "
          f"{len(cache):,} kept (with pathway or PPI)")
    print(f"Time: {elapsed:.0f}s")
    print(f"\nWriting {OUT_CACHE} ...", flush=True)
    with OUT_CACHE.open("w") as fh:
        json.dump(cache, fh)
    print(f"  Size: {OUT_CACHE.stat().st_size/1e6:.1f} MB")
    print("Done.")


if __name__ == "__main__":
    main()
