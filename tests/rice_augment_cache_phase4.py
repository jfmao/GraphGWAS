"""Phase 4 — augment rice graph_cache with per-variant prior_score.

v1 cache: Ren-2023 pathway labels + RicePPINet PPI (already heterogeneous).
v2 cache: + prior_score derived from snpEff impact class (HIGH=1.0,
          MODERATE=0.7, LOW=0.4, MODIFIER=0.1) by re-scanning the pseudo-
          canonical VCF; HIGH-impact = stop_gained, frameshift, missense, etc.

This adds the most heterogeneous per-variant prior we have available for rice
without external downloads. Future Phase 4.B will add Cao 2020 / Tu 2020
ChIP-seq peak overlap when we have time to download.
"""
from __future__ import annotations

import bisect
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path

DATA = Path("/mnt/data/GraphGWAS/data/rice_3k")
ANN_DIR = DATA / "annotations"
PSEUDO_VCF = Path("/mnt/data/GraphPop/data/raw/3kRG_data/NB_bialSNP_pseudo_canonical_ALL.vcf.gz")
IMPACT_SCORE = {b"HIGH": 1.0, b"MODERATE": 0.7, b"LOW": 0.4, b"MODIFIER": 0.1}


def main():
    """Stream VCF, extract per-variant max snpEff impact, augment v1 caches."""
    print("Streaming pseudo-canonical VCF for snpEff impact class ...", flush=True)
    impact_re = re.compile(rb"\|(HIGH|MODERATE|LOW|MODIFIER)\|")
    var_score: dict[str, float] = {}
    n = 0; n_with = 0
    with gzip.open(PSEUDO_VCF, "rb") as fh:
        for line in fh:
            if line.startswith(b"#"): continue
            n += 1
            if n % 2_000_000 == 0:
                print(f"  {n/1e6:.1f}M scanned, {n_with:,} with impact", flush=True)
            parts = line.split(b"\t", 8)
            if len(parts) < 8: continue
            chrom_raw = parts[0].decode()
            chrom = chrom_raw if chrom_raw.startswith("Chr") else f"Chr{chrom_raw}"
            pos = parts[1].decode()
            ref = parts[3].decode()
            alt = parts[4].decode()
            info = parts[7]
            # Find max impact across ANN annotations
            max_score = 0.0
            for m in impact_re.findall(info):
                s = IMPACT_SCORE.get(m, 0.0)
                if s > max_score:
                    max_score = s
            if max_score > 0:
                vid = f"{chrom}:{pos}:{ref}:{alt}"
                var_score[vid] = max_score
                n_with += 1
    print(f"\nTotal variants with snpEff impact: {len(var_score):,}")

    # Augment per-chromosome cache files in place
    print("\nAugmenting per-chrom v1 caches → v2 ...", flush=True)
    for chrom_n in range(1, 13):
        chrom = f"Chr{chrom_n}"
        in_p = ANN_DIR / f"rice_graph_cache_{chrom}.json"
        if not in_p.exists():
            print(f"  {chrom}: missing v1, skip"); continue
        cache = json.load(in_p.open())
        n_aug = 0
        for vid, info in cache.items():
            sc = var_score.get(vid)
            if sc is not None and sc > 0:
                info["prior_score"] = float(sc)
                n_aug += 1
        out = ANN_DIR / f"rice_graph_cache_v2_{chrom}.json"
        with out.open("w") as fh:
            json.dump(cache, fh)
        print(f"  {chrom}: {n_aug:,}/{len(cache):,} variants got prior_score "
              f"({100*n_aug/len(cache):.0f}%); wrote {out.name} "
              f"({out.stat().st_size/1e6:.0f} MB)", flush=True)


if __name__ == "__main__":
    main()
