"""Build a chr22 multi-omics graph_cache for human ablation test.

Schema matches the rice cache (variant → {genes, pathways, ppi}):
  - Variant→Gene: from GTEx chr22 eQTL cache (ENSG IDs)
  - Gene→Pathway: empty for now (Reactome not pre-parsed on this host)
  - Gene↔Gene PPI: STRING v12 9606 links, combined_score ≥ 700 (high-conf)

Writes:
  data/annotations/human_graph_cache_chr22.json
"""
from __future__ import annotations

import gzip
import json
from collections import defaultdict
from pathlib import Path

ANNO = Path("/mnt/data/GraphGWAS/data/annotations")
EQTL_PATH = ANNO / "gtex_chr22_eqtl_cache.json"
STRING_LINKS = ANNO / "9606.protein.links.v12.0.txt.gz"
STRING_INFO = ANNO / "9606.protein.info.v12.0.txt.gz"
OUT = ANNO / "human_graph_cache_chr22.json"
STRING_THRESHOLD = 700  # STRING combined_score (high confidence)

# ---------------------------------------------------------
# 1. ENSP -> ENSG (approximate via preferred_name only if needed)
# ---------------------------------------------------------
# STRING uses ENSP (protein) IDs; GTEx uses ENSG (gene) IDs.
# We strip the "9606." prefix and use ENSP directly as the PPI node.
# Then we need an ENSG -> ENSP map, which we derive from a simpler
# gene-name pivot: a gene with any eQTL variant gets its ENSP via
# a permissive fallback — if unknown, include the ENSG itself as a
# node that only connects to other ENSG nodes with known ENSP.
# Simpler: use ENSG ↔ ENSG PPI by resolving STRING ENSP → gene name.
print("[1/3] Parsing STRING info: ENSP → gene symbol ...", flush=True)
ensp_to_sym: dict[str, str] = {}
with gzip.open(STRING_INFO, "rt") as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 2: continue
        ensp = f[0].split(".", 1)[1] if "." in f[0] else f[0]
        sym = f[1]
        ensp_to_sym[ensp] = sym
print(f"  {len(ensp_to_sym):,} ENSP→symbol mappings")

# ---------------------------------------------------------
# 2. High-confidence PPI (ENSP space, then map to symbols)
# ---------------------------------------------------------
print(f"\n[2/3] Parsing STRING links (combined_score ≥ {STRING_THRESHOLD}) ...",
      flush=True)
sym_to_ppi: dict[str, set[str]] = defaultdict(set)
n_total = 0
n_kept = 0
with gzip.open(STRING_LINKS, "rt") as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split(" ")
        if len(f) < 3: continue
        n_total += 1
        a, b, s = f[0], f[1], int(f[2])
        if s < STRING_THRESHOLD: continue
        n_kept += 1
        ea = a.split(".", 1)[1] if "." in a else a
        eb = b.split(".", 1)[1] if "." in b else b
        sa = ensp_to_sym.get(ea); sb = ensp_to_sym.get(eb)
        if not sa or not sb: continue
        sym_to_ppi[sa].add(sb)
        sym_to_ppi[sb].add(sa)
print(f"  {n_total:,} total edges; {n_kept:,} kept at ≥{STRING_THRESHOLD}")
print(f"  {len(sym_to_ppi):,} genes with ≥1 HC PPI partner")

# ---------------------------------------------------------
# 3. Build Variant→{genes, pathways, ppi} from GTEx eQTL cache
# ---------------------------------------------------------
print("\n[3/3] Combining with GTEx chr22 eQTL cache ...", flush=True)
eqtl = json.load(open(EQTL_PATH))
print(f"  {len(eqtl):,} variants in GTEx chr22 cache")

# GTEx cache has ENSG, but STRING is in symbols. We need a gene symbol
# per ENSG. Use preferred_name heuristic: compute ENSG→symbol from
# ENSP→symbol via first-pass STRING info approximation — not perfect,
# but reasonable for chr22 analysis since HBP only cares about graph
# connectivity not ID types.
# Simpler: treat the ENSG from GTEx as the gene node itself and build
# ENSG↔ENSG PPI by matching partner genes' ENSGs.
# Since STRING is ENSP, we'd need a full ENSG↔ENSP map. Without it,
# we use the ENSG string as the gene node and leave PPI per-variant
# empty. Then we add symbol-based PPI via a HGNC alias step.
# For pragmatism, we pull gene symbol for each GTEx ENSG using
# /mnt/data/GraphGWAS/data/annotations/gencode/ if available, else
# skip PPI coupling and keep genes/pathways only.
gencode = ANNO / "gencode"
ensg_to_sym: dict[str, str] = {}
# Search for a GTF or a symbol lookup
import glob
for fp in glob.glob(str(gencode / "*.gtf*")):
    print(f"  trying gencode file {fp}")
    opener = gzip.open if fp.endswith(".gz") else open
    with opener(fp, "rt") as fh:
        for line in fh:
            if line.startswith("#") or "\tgene\t" not in line: continue
            if "chr22" not in line and "\t22\t" not in line: continue
            attrs = line.rstrip().split("\t")[-1]
            ensg = None; sym = None
            for kv in attrs.split(";"):
                kv = kv.strip()
                if kv.startswith("gene_id "):
                    ensg = kv.split('"')[1].split(".")[0]
                elif kv.startswith("gene_name "):
                    sym = kv.split('"')[1]
            if ensg and sym:
                ensg_to_sym[ensg] = sym
    break  # first GTF file
print(f"  ENSG→symbol map: {len(ensg_to_sym):,} entries")

cache: dict[str, dict] = {}
for vid, info in eqtl.items():
    gene_ensg = info.get("gene", "").split(".")[0]
    if not gene_ensg:
        continue
    sym = ensg_to_sym.get(gene_ensg)
    # Try stripping version suffix even if map missing
    genes = [sym] if sym else [gene_ensg]
    # PPI: all high-confidence partners of this gene's symbol
    ppi_partners: list[str] = []
    if sym and sym in sym_to_ppi:
        # Cap per-variant fan-out to keep cache small
        ppi_partners = list(sym_to_ppi[sym])[:20]
    if not genes and not ppi_partners:
        continue
    cache[vid] = {
        "genes": genes,
        "pathways": [],  # Reactome not parsed on this host
        "ppi": ppi_partners,
    }

print(f"\n  Output cache: {len(cache):,} variants")
print(f"  Writing {OUT} ...", flush=True)
with OUT.open("w") as f:
    json.dump(cache, f)
print(f"  Size: {OUT.stat().st_size / 1e6:.1f} MB")
print(f"  Done.")
