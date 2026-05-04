"""Populate the `pathways` field of human_graph_cache_v2_chr*.json from
Reactome (Homo sapiens subset).

The shipped human caches have `pathways: []` for all variants — they
were built before the Reactome integration step. This means M2's
`same_pathway` motif fires zero matches on human, severely handicapping
M2's recall (3/70 catalogue rank-1 in §Y.8). Patching adds Reactome
pathway names per gene; expected to substantially improve M2.

Source: data/reactome/ReactomePathways.gmt (Homo-sapiens-only pathways
filtered by R-HSA- prefix on col 2). Auto-downloads from Reactome if
missing.

Usage: python scripts/patch_human_cache_with_reactome.py [--dry-run]
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.request
import zipfile
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
REACTOME_DIR = REPO / "data" / "reactome"
GMT_PATH = REACTOME_DIR / "ReactomePathways.gmt"
GMT_ZIP = REACTOME_DIR / "ReactomePathways.gmt.zip"
HGNC_MAP_OUT = REACTOME_DIR / "hgnc_to_reactome_human.json"
CACHE_DIR = REPO / "data" / "annotations"


def ensure_gmt() -> None:
    """Download and extract ReactomePathways.gmt if missing."""
    if GMT_PATH.exists() and GMT_PATH.stat().st_size > 1000:
        return
    REACTOME_DIR.mkdir(parents=True, exist_ok=True)
    if not GMT_ZIP.exists() or GMT_ZIP.stat().st_size < 1000:
        url = "https://reactome.org/download/current/ReactomePathways.gmt.zip"
        print(f"Downloading {url}...")
        # Reactome requires a User-Agent header
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req) as r, open(GMT_ZIP, "wb") as f:
            f.write(r.read())
        print(f"  → {GMT_ZIP.stat().st_size/1e6:.2f} MB")
    print(f"Extracting {GMT_ZIP.name}...")
    with zipfile.ZipFile(GMT_ZIP) as z:
        z.extractall(REACTOME_DIR)
    print(f"  → {GMT_PATH.stat().st_size/1e6:.2f} MB")


def build_hgnc_pathway_map() -> dict[str, list[str]]:
    """Parse Reactome GMT, return {hgnc_symbol: [pathway_name, ...]} for human only."""
    print(f"Parsing {GMT_PATH.name} for Homo sapiens (R-HSA- prefix)...")
    out: dict[str, set] = defaultdict(set)
    n_pathways_human = 0
    n_pathways_total = 0
    with open(GMT_PATH) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            pathway_name, pathway_id = parts[0], parts[1]
            n_pathways_total += 1
            if not pathway_id.startswith("R-HSA-"):
                continue
            n_pathways_human += 1
            for gene in parts[2:]:
                gene = gene.strip()
                if gene:
                    out[gene].add(pathway_name)
    out_list = {g: sorted(p) for g, p in out.items()}
    print(f"  {n_pathways_human:,} human pathways (of {n_pathways_total:,} total)")
    print(f"  {len(out_list):,} unique HGNC symbols mapped")
    return out_list


def patch_caches(hgnc_to_pathways: dict[str, list[str]],
                  dry_run: bool = False) -> dict:
    """For each human_graph_cache_v2_chr*.json, populate `pathways`."""
    chrom_files = sorted(CACHE_DIR.glob("human_graph_cache_v2_chr*.json"))
    summary = []
    for cache_path in chrom_files:
        chrom = cache_path.stem.replace("human_graph_cache_v2_", "")
        t0 = time.time()
        cache = json.loads(cache_path.read_text())
        n_total = len(cache)
        n_modified = 0
        n_with_pathways_before = 0
        n_with_pathways_after = 0
        max_pathways_per_var = 0
        for vid, entry in cache.items():
            if entry.get("pathways"):
                n_with_pathways_before += 1
            existing = set(entry.get("pathways") or [])
            new_paths: set = set()
            for g in entry.get("genes", []) or []:
                if g in hgnc_to_pathways:
                    new_paths.update(hgnc_to_pathways[g])
            # Cap per-variant fan-out at 10 (avoid bloat for hub genes;
            # M2 motif filtering already deduplicates per-pathway anyway)
            new_paths_capped = sorted(new_paths)[:10]
            merged = sorted(existing | set(new_paths_capped))[:10]
            if merged and merged != entry.get("pathways"):
                entry["pathways"] = merged
                n_modified += 1
            if entry.get("pathways"):
                n_with_pathways_after += 1
                max_pathways_per_var = max(max_pathways_per_var,
                                            len(entry["pathways"]))
        if not dry_run:
            cache_path.write_text(json.dumps(cache))
        summary.append({
            "chrom": chrom,
            "n_total": n_total,
            "n_modified": n_modified,
            "n_with_pathways_before": n_with_pathways_before,
            "n_with_pathways_after": n_with_pathways_after,
            "max_pathways_per_var": max_pathways_per_var,
            "elapsed_s": round(time.time() - t0, 1),
        })
        print(f"  {cache_path.name}: {n_total:>7,} variants  "
              f"with_paths {n_with_pathways_before:>6,}→{n_with_pathways_after:>6,}  "
              f"max_paths/var={max_pathways_per_var:>2}  "
              f"in {summary[-1]['elapsed_s']}s")
    return summary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true",
                    help="don't write patched caches; only report counts")
    args = ap.parse_args()

    ensure_gmt()
    hgnc_map = build_hgnc_pathway_map()
    HGNC_MAP_OUT.write_text(json.dumps(hgnc_map, indent=1))
    print(f"Wrote {HGNC_MAP_OUT} ({HGNC_MAP_OUT.stat().st_size/1e6:.2f} MB)")
    print()
    print(f"Patching {len(list(CACHE_DIR.glob('human_graph_cache_v2_chr*.json')))} "
          f"human cache files{' (dry run)' if args.dry_run else ''}...")
    summary = patch_caches(hgnc_map, dry_run=args.dry_run)
    print()
    total_modified = sum(s["n_modified"] for s in summary)
    total_with_paths = sum(s["n_with_pathways_after"] for s in summary)
    total_var = sum(s["n_total"] for s in summary)
    print(f"Total variants modified: {total_modified:,} "
          f"({100*total_modified/total_var:.1f}% of all human variants)")
    print(f"Total variants with pathways after: {total_with_paths:,} "
          f"({100*total_with_paths/total_var:.1f}%)")


if __name__ == "__main__":
    main()
