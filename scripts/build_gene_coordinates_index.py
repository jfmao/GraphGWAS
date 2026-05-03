"""Pre-scan all species chromosome caches and emit a gene-systematic →
(chrom, pos_min, pos_max, n_variants_in_cache) lookup table.

Used by tests/multi_pair_catalogue_benchmark.py to discover which
chromosome each catalogue gene lives on without re-loading the per-chromosome
caches each time. One-time cost: ~30 sec, peak memory ~600 MB (one cache
loaded at a time and freed).

Output: docs/groundtruth/gene_coordinates_index.json
"""
from __future__ import annotations

import gc
import json
import sys
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
INDEX_OUT = REPO_ROOT / "docs" / "groundtruth" / "gene_coordinates_index.json"


SPECIES_CACHE_FILES = {
    "yeast": [REPO_ROOT / "data" / "yeast" / "yeast_graph_cache_v2.json"],
    "arabidopsis": [
        REPO_ROOT / "data" / "arabidopsis" / f"arabidopsis_graph_cache_v3_chr{c}.json"
        for c in (1, 2, 3, 4, 5)
    ],
    "rice": [
        REPO_ROOT / "data" / "rice_3k" / "annotations" / f"rice_graph_cache_v2_Chr{c}.json"
        for c in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
    ],
    "human": [
        REPO_ROOT / "data" / "annotations" / f"human_graph_cache_v2_chr{c}.json"
        for c in [str(i) for i in range(1, 23)] + ["X"]
    ],
}


def index_one_cache(cache_path: Path) -> dict:
    """Return {gene_systematic: {chrom: {pos_min, pos_max, n_variants}}}.

    A gene may appear on multiple chromosomes if cache annotations cross-link
    (rare). We track per-chrom counts so the caller can pick the dominant one.
    """
    print(f"  scanning {cache_path.name} ({cache_path.stat().st_size / 1e6:.0f} MB)...",
          end="", flush=True)
    cache = json.loads(cache_path.read_text())
    n_var_total = len(cache)
    gene_chrom_pos: dict = defaultdict(lambda: defaultdict(list))
    for vid, entry in cache.items():
        # vid format: "chrom:pos:ref:alt"
        parts = vid.split(":")
        if len(parts) < 2:
            continue
        chrom = parts[0]
        try:
            pos = int(parts[1])
        except ValueError:
            continue
        for gene in entry.get("genes", []):
            gene_chrom_pos[gene][chrom].append(pos)
    print(f" {n_var_total:,} variants, {len(gene_chrom_pos):,} genes annotated")
    out: dict = {}
    for gene, by_chrom in gene_chrom_pos.items():
        out[gene] = {}
        for chrom, positions in by_chrom.items():
            out[gene][chrom] = {
                "pos_min": int(min(positions)),
                "pos_max": int(max(positions)),
                "n_variants": int(len(positions)),
            }
    del cache, gene_chrom_pos
    gc.collect()
    return out


def merge(into: dict, src: dict) -> None:
    for gene, by_chrom in src.items():
        if gene not in into:
            into[gene] = {}
        for chrom, stats in by_chrom.items():
            if chrom in into[gene]:
                # Same gene on same chrom in two caches — extend min/max
                into[gene][chrom]["pos_min"] = min(into[gene][chrom]["pos_min"],
                                                     stats["pos_min"])
                into[gene][chrom]["pos_max"] = max(into[gene][chrom]["pos_max"],
                                                     stats["pos_max"])
                into[gene][chrom]["n_variants"] += stats["n_variants"]
            else:
                into[gene][chrom] = stats


def main() -> None:
    index: dict = {"schema_version": "1.0", "species": {}}
    for sp, cache_files in SPECIES_CACHE_FILES.items():
        print(f"\n=== {sp.upper()} ({len(cache_files)} cache files) ===")
        sp_index: dict = {}
        for cf in cache_files:
            if not cf.exists():
                print(f"  ⚠ MISSING: {cf}")
                continue
            sub = index_one_cache(cf)
            merge(sp_index, sub)
            del sub
            gc.collect()
        index["species"][sp] = sp_index
        print(f"  total {sp} genes: {len(sp_index):,}")

    INDEX_OUT.parent.mkdir(parents=True, exist_ok=True)
    INDEX_OUT.write_text(json.dumps(index, indent=1))
    print(f"\n  Wrote {INDEX_OUT}  "
          f"({INDEX_OUT.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
