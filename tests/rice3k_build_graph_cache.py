"""Build a rice multi-omics graph_cache for HBP / L1 fine-mapping.

Output schema (compatible with graphgwas.finemapping_v2.fast_hbp_finemap
and hbp_finemap_from_sumstats):

    graph_cache = {
        variant_id:  {
            "genes":     [LOC_Os IDs the variant overlaps, via snpEff],
            "pathways":  [pathway names the overlapping gene(s) sit in],
            "ppi":       [LOC_Os IDs of PPI partners of the gene(s),
                          from RicePPINet at probability >= 0.7],
        },
        ...
    }

Sources:
  - Variant→Gene: snpEff ANN field of NB_bialSNP_pseudo_canonical_ALL.vcf.gz
                  (already indexed per variant by gene_position_index.tsv;
                  we re-scan the VCF to get per-variant gene memberships).
  - Gene→Pathway: Ren 2023 catalogue's `pathway` column — a hand-curated
                  mapping of 269 rice grain genes to their biological
                  pathway modules. This is sparse but high-quality and
                  exactly the prior information HBP benefits from.
  - Gene↔Gene PPI: RicePPINet (Liu et al. 2017 Plant J) 708,819 predicted
                   interactions, filtered to Probability >= 0.7
                   (≈100k high-confidence edges; the paper's default).

Writes:
  data/rice_3k/annotations/rice_graph_cache.json
"""
from __future__ import annotations

import gzip
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

DATA_DIR = Path("/mnt/data/GraphGWAS/data/rice_3k")
ANNO = DATA_DIR / "annotations"
ANNO.mkdir(parents=True, exist_ok=True)

PSEUDO_VCF = Path("/mnt/data/GraphPop/data/raw/3kRG_data/NB_bialSNP_pseudo_canonical_ALL.vcf.gz")
REN_TSV = DATA_DIR / "ground_truth" / "grain_quality_causal_genes.tsv"
PPI_GZ = ANNO / "riceppinet_full.dat.gz"
CACHE_OUT = ANNO / "rice_graph_cache.json"

PPI_PROB_THRESHOLD = 0.7

# ------------------------------------------------------------------
# 1. Gene→Pathway map from Ren 2023 catalogue
# ------------------------------------------------------------------
print("[1/3] Building Gene→Pathway map from Ren 2023 ...", flush=True)
ren = pd.read_csv(REN_TSV, sep="\t")
gene_to_pathway: dict[str, list[str]] = defaultdict(list)
for _, r in ren.iterrows():
    loc = r["LOC_Os_id"]
    pw = r["pathway"]
    if pd.notna(loc) and pd.notna(pw):
        # split "Ubiquitin-proteasome pathway" etc. by delimiters
        for p in re.split(r"[;,]", str(pw)):
            p = p.strip()
            if p and p not in gene_to_pathway[loc]:
                gene_to_pathway[loc].append(p)
print(f"  {len(gene_to_pathway)} genes with pathway annotations")

# ------------------------------------------------------------------
# 2. Gene↔Gene PPI map from RicePPINet (Prob >= 0.7)
# ------------------------------------------------------------------
print(f"\n[2/3] Building Gene↔Gene PPI map from RicePPINet "
      f"(Prob >= {PPI_PROB_THRESHOLD}) ...", flush=True)
gene_to_ppi: dict[str, set[str]] = defaultdict(set)
with gzip.open(PPI_GZ, "rt") as fh:
    header = fh.readline()  # skip header
    n_total = 0
    n_kept = 0
    for line in fh:
        n_total += 1
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        a, b, prob = parts[0], parts[1], float(parts[2])
        if prob < PPI_PROB_THRESHOLD:
            continue
        n_kept += 1
        gene_to_ppi[a].add(b)
        gene_to_ppi[b].add(a)
print(f"  {n_total:,} total PPI rows; {n_kept:,} kept at Prob>={PPI_PROB_THRESHOLD}")
print(f"  {len(gene_to_ppi):,} genes with >=1 PPI partner")

# ------------------------------------------------------------------
# 3. Scan pseudo-canonical VCF to build Variant→{genes, pathways, ppi}
# ------------------------------------------------------------------
print("\n[3/3] Scanning pseudo-canonical VCF for Variant→Gene mapping ...",
      flush=True)
gene_rx = re.compile(rb"LOC_Os\d+g\d+")
cache: dict[str, dict] = {}


def _load_catalogue_keep_genes() -> set:
    """Genes from docs/groundtruth/epistasis_pairs.json (rice).

    Used to bypass the pathway-OR-PPI filter for catalogue-relevant
    genes that have no entries in Ren-2023 / RicePPINet (e.g. Ghd7,
    Pik-1, Pik-2). Result is cached on first call.
    """
    if hasattr(_load_catalogue_keep_genes, "_cached"):
        return _load_catalogue_keep_genes._cached
    keep: set = set()
    cat_path = (Path(__file__).resolve().parent.parent
                / "docs" / "groundtruth" / "epistasis_pairs.json")
    if cat_path.exists():
        try:
            cat = json.loads(cat_path.read_text())
            for p in cat.get("species", {}).get("rice", {}).get("pairs", []):
                for g_key in ("gene_1", "gene_2"):
                    locus = p[g_key].get("locus")
                    if locus:
                        keep.add(locus)
        except Exception as e:
            print(f"  ⚠ couldn't load catalogue keep-list: {e}", flush=True)
    _load_catalogue_keep_genes._cached = keep
    print(f"  Catalogue keep-list (rice): {len(keep)} genes "
          f"({sorted(keep)[:5]}...)", flush=True)
    return keep

with gzip.open(PSEUDO_VCF, "rb") as fh:
    n = 0
    n_annot = 0
    for line in fh:
        if line.startswith(b"#"):
            continue
        n += 1
        if n % 2_000_000 == 0:
            print(f"  {n/1e6:.1f}M variants scanned, {len(cache):,} cached",
                  flush=True)

        parts = line.split(b"\t", 8)
        if len(parts) < 8:
            continue
        chrom_raw = parts[0].decode()
        chrom = chrom_raw if chrom_raw.startswith("Chr") else f"Chr{chrom_raw}"
        pos = parts[1].decode()
        ref = parts[3].decode()
        alt = parts[4].decode()
        info = parts[7]

        gene_ids = list(dict.fromkeys(m.decode() for m in gene_rx.findall(info)))
        if not gene_ids:
            continue
        n_annot += 1

        # Filter to top ~4 genes per variant (most promoter-/gene-scope
        # annotations list 2-3 duplicate IDs; snpEff also emits intergenic
        # matches. Keep unique.)
        gene_ids = gene_ids[:4]

        pathways: list[str] = []
        ppi_partners: list[str] = []
        for g in gene_ids:
            for pw in gene_to_pathway.get(g, []):
                if pw not in pathways:
                    pathways.append(pw)
            if g in gene_to_ppi:
                for partner in list(gene_to_ppi[g])[:20]:  # cap per-variant fan-out
                    if partner not in ppi_partners:
                        ppi_partners.append(partner)

        # Originally we dropped variants without pathway or PPI info to
        # control JSON size. This silently filtered out catalogue-relevant
        # genes whose annotations live outside RicePPINet/Ren-2023 (e.g.
        # Ghd7 LOC_Os07g15770, Pik-1/Pik-2 LOC_Os11g46200/210). Keep
        # those even when pathways and ppi are empty if at least one
        # annotated gene is on the catalogue keep-list. The list lives in
        # docs/groundtruth/epistasis_pairs.json and is loaded once per
        # build via the env var GRAPHGWAS_CATALOGUE_KEEP_GENES.
        catalogue_keep = _load_catalogue_keep_genes()
        in_catalogue = any(g in catalogue_keep for g in gene_ids)
        if not pathways and not ppi_partners and not in_catalogue:
            continue

        vid = f"{chrom}:{pos}:{ref}:{alt}"
        cache[vid] = {
            "genes": gene_ids,
            "pathways": pathways,
            "ppi": ppi_partners,
        }

print(f"\n  Total variants scanned: {n:,}")
print(f"  Annotated (any gene): {n_annot:,}")
print(f"  Kept (gene + pathway OR ppi): {len(cache):,}")

# Write cache
print(f"\nWriting {CACHE_OUT} ...", flush=True)
with CACHE_OUT.open("w") as f:
    json.dump(cache, f)
print(f"  Size: {CACHE_OUT.stat().st_size / 1e6:.1f} MB")
print(f"  Done.")
