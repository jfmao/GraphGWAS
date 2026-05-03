"""Build an Arabidopsis multi-omics graph_cache.

Inputs:
  - tair10.gff3.gz (NCBI TAIR10.1) — variant→gene by overlap
  - Plant Reactome (Arabidopsis subset) — gene→pathway
  - STRING v12 species 3702 — gene↔gene PPI (combined_score≥700)

Output:
  data/arabidopsis/arabidopsis_graph_cache_chr{1..5}.json (per-chrom for lazy loading)
"""
from __future__ import annotations

import gzip
import json
import re
import time
from collections import defaultdict
from pathlib import Path

ANN_DIR = Path("/mnt/data/GraphGWAS/data/arabidopsis/annotations")
TAIR_GFF = ANN_DIR / "tair10.gff3.gz"
STRING_LINKS = ANN_DIR / "3702.protein.links.v12.0.txt.gz"
STRING_INFO = ANN_DIR / "3702.protein.info.v12.0.txt.gz"
PLANT_REACTOME = Path("/mnt/data/GraphPop/data/raw/plant_reactome/gene_ids_by_pathway_and_species.tab")
VCF = Path("/mnt/data/GraphGWAS/tests/data/arabidopsis/1001genomes_snp-short-indel_only_ACGTN.vcf.gz")
OUT_DIR = Path("/mnt/data/GraphGWAS/data/arabidopsis")
OUT_DIR.mkdir(parents=True, exist_ok=True)
STRING_THRESHOLD = 700

# NCBI chromosome IDs → 1KG-style numeric chrom names used in the VCF
NCBI_TO_CHROM = {
    "NC_003070.9": "1",  # Chr1
    "NC_003071.7": "2",  # Chr2
    "NC_003074.8": "3",  # Chr3
    "NC_003075.7": "4",  # Chr4
    "NC_003076.8": "5",  # Chr5
}

# AT-gene ID regex
AT_RE = re.compile(r"AT[1-5MC]G\d+")


def parse_gff_genes() -> dict[str, list[tuple[int, int, str]]]:
    """Return {chrom_str: [(start, end, AT_id)] sorted by start}."""
    out = defaultdict(list)
    with gzip.open(TAIR_GFF, "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9: continue
            chrom_ncbi = f[0]
            if chrom_ncbi not in NCBI_TO_CHROM: continue
            # Allow `pseudogene` too — some catalogue genes (e.g. AT4G03060
            # / AOP2) are marked pseudogene in COL-0 but functional in
            # other ecotypes (Cvi). See scripts/patch_caches_for_missing_genes.py.
            if f[2] not in ("gene", "pseudogene"): continue
            try:
                start, end = int(f[3]), int(f[4])
            except ValueError:
                continue
            attrs = f[8]
            m = AT_RE.search(attrs)
            if not m: continue
            chrom = NCBI_TO_CHROM[chrom_ncbi]
            out[chrom].append((start, end, m.group(0)))
    for c in out:
        out[c].sort()
    return dict(out)


def parse_plant_reactome() -> dict[str, list[str]]:
    """Return {AT_id: [pathway_names]} for Arabidopsis."""
    out = defaultdict(list)
    with PLANT_REACTOME.open() as fh:
        # Format: pathway_id\tpathway_name\tspecies\tgene_id
        next(fh, None)  # skip header if present (try)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4: continue
            if "Arabidopsis thaliana" not in f[2]: continue
            pname, gid = f[1], f[3].upper()
            if not gid.startswith("AT"): continue
            if pname not in out[gid]:
                out[gid].append(pname)
    return dict(out)


def parse_string_ath_ppi() -> dict[str, set[str]]:
    """Return {AT_id: {AT_id partners}} via UniProt→AT mapping (idmapping.dat).
    Arabidopsis STRING IDs are UniProt accessions (3702.A0A0A7EPL0 etc.),
    not AT-IDs. We use the UniProt knowledgebase idmapping file
    (Gene_OrderedLocusName field) to map UniProt → ATxGyyyyy.
    """
    UNIPROT_DAT = ANN_DIR / "uniprot_3702_idmapping.tab.gz"
    print(f"  Parsing UniProt → AT mapping ({UNIPROT_DAT.name}) ...", flush=True)
    uniprot_to_at: dict[str, str] = {}
    with gzip.open(UNIPROT_DAT, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3: continue
            if f[1] != "Gene_OrderedLocusName": continue
            uid, at = f[0], f[2].upper()
            m = AT_RE.search(at)
            if m:
                uniprot_to_at[uid] = m.group(0)
    print(f"    {len(uniprot_to_at):,} UniProt→AT mappings")

    # Pull symbol→UniProt from STRING info as a fallback path
    print(f"  Parsing STRING info (3702 species) ...", flush=True)
    ensp_to_at: dict[str, str] = {}
    with gzip.open(STRING_INFO, "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2: continue
            ensp = f[0].split(".", 1)[1] if "." in f[0] else f[0]
            # 1) UniProt → AT lookup (primary)
            at = uniprot_to_at.get(ensp)
            # 2) AT-regex on preferred_name (fallback)
            if not at:
                m = AT_RE.search(f[1])
                if m: at = m.group(0)
            if at:
                ensp_to_at[ensp] = at
    print(f"    {len(ensp_to_at):,} STRING-protein → AT mappings")

    print(f"  Parsing STRING links (combined_score≥{STRING_THRESHOLD}) ...", flush=True)
    ppi: dict[str, set[str]] = defaultdict(set)
    n_total = 0; n_kept = 0
    with gzip.open(STRING_LINKS, "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split(" ")
            if len(f) < 3: continue
            n_total += 1
            try:
                s = int(f[2])
            except ValueError:
                continue
            if s < STRING_THRESHOLD: continue
            n_kept += 1
            a = f[0].split(".", 1)[1] if "." in f[0] else f[0]
            b = f[1].split(".", 1)[1] if "." in f[1] else f[1]
            ag = ensp_to_at.get(a) or (AT_RE.search(a).group(0) if AT_RE.search(a) else None)
            bg = ensp_to_at.get(b) or (AT_RE.search(b).group(0) if AT_RE.search(b) else None)
            if not ag or not bg: continue
            ppi[ag].add(bg); ppi[bg].add(ag)
    print(f"    {n_total:,} edges total, {n_kept:,} kept; "
          f"{len(ppi):,} genes with ≥1 partner")
    return dict(ppi)


def main():
    t0 = time.time()
    print("[1/4] Parsing TAIR10 GFF (gene → chrom + range) ...", flush=True)
    chrom_genes = parse_gff_genes()
    print(f"  Genes per chrom: " +
          ", ".join(f"chr{c}={len(chrom_genes[c])}" for c in sorted(chrom_genes)))

    print("\n[2/4] Parsing Plant Reactome (Arabidopsis gene→pathway) ...", flush=True)
    gene_paths = parse_plant_reactome()
    print(f"  {len(gene_paths)} genes with pathway annotations")

    print("\n[3/4] Parsing STRING v12 (species 3702) PPI ...", flush=True)
    gene_ppi = parse_string_ath_ppi()

    print("\n[4/4] Streaming VCF and building per-chrom cache ...", flush=True)
    import subprocess
    out_caches = {c: {} for c in NCBI_TO_CHROM.values()}
    cmd = ["bcftools","query","-f","%CHROM\t%POS\t%REF\t%ALT\n", str(VCF)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1<<20)
    n = 0; n_annot = 0
    for line in proc.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 4: continue
        n += 1
        if n % 1_000_000 == 0:
            print(f"  {n/1e6:.1f}M variants scanned, "
                  f"{sum(len(c) for c in out_caches.values()):,} cached",
                  flush=True)
        chrom, pos, ref, alt = f[0], int(f[1]), f[2], f[3]
        if chrom not in out_caches: continue
        # binary search in chrom_genes
        gs = chrom_genes.get(chrom, [])
        # linear works since avg ~5500 genes/chrom; but use binary search via bisect
        import bisect
        # gs is sorted by start; find first gene with start <= pos and end >= pos
        starts = [g[0] for g in gs]  # could precompute outside loop, but tolerable
        # quicker: iterate over a small window
        i0 = bisect.bisect_right(starts, pos) - 1
        overlaps = []
        for j in range(max(0, i0-3), min(len(gs), i0+2)):
            s, e, gid = gs[j]
            if s <= pos <= e:
                overlaps.append(gid)
        if not overlaps: continue
        n_annot += 1
        paths = []
        ppi = []
        for g in overlaps[:4]:
            for p in gene_paths.get(g, []):
                if p not in paths: paths.append(p)
            if g in gene_ppi:
                for partner in list(gene_ppi[g])[:20]:
                    if partner not in ppi: ppi.append(partner)
        if not paths and not ppi: continue
        vid = f"{chrom}:{pos}:{ref}:{alt}"
        out_caches[chrom][vid] = {
            "genes": overlaps[:4], "pathways": paths, "ppi": ppi,
        }
    proc.wait()
    print(f"\n  {n:,} total variants scanned; {n_annot:,} hit a gene; "
          f"{sum(len(c) for c in out_caches.values()):,} kept (with pathway/PPI)")

    for chrom, cache in out_caches.items():
        out = OUT_DIR / f"arabidopsis_graph_cache_chr{chrom}.json"
        with out.open("w") as fh:
            json.dump(cache, fh)
        print(f"  Wrote {out.name}: {len(cache):,} variants, "
              f"{out.stat().st_size/1e6:.1f} MB")
    print(f"\nTotal time: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
