"""Phase 1.A — augment human graph_cache with local-file functional annotations.

Adds:
  - regulatory: ENCODE cCRE element ID + class (pELS, dELS, PLS, CTCF-bound, etc.)
                from /data/annotations/encode/encodeCcreCombined.bed
  - prior_score: 1.0 if variant overlaps any cCRE; 0.5 if PLS; 1.5 if PLS+CTCF
                 (per-variant heterogeneous score driving the heterogeneity test)

The cCRE catalog has 1.06M elements with 9 classes (pELS, dELS, PLS, CTCF-bound,
etc.); per-variant cCRE classes will be highly heterogeneous, exactly the kind
of annotation our heterogeneity-not-density principle predicts will help HBP.

Reads existing data/annotations/human_graph_cache_chr<N>.json (genome-wide GTEx + STRING),
augments in-place to data/annotations/human_graph_cache_v2_chr<N>.json.
"""
from __future__ import annotations

import bisect
import json
from collections import defaultdict
from pathlib import Path

ANN = Path("/mnt/data/GraphGWAS/data/annotations")
CCRE_BED = ANN / "encode" / "encodeCcreCombined.bed"
CACHE_GLOB = "human_graph_cache_chr*.json"
OUT_PREFIX = "human_graph_cache_v2_chr"

# cCRE class → prior_score (heterogeneous: PLS = strong promoter, pELS = mid, etc.)
CLASS_SCORE = {
    "PLS,CTCF-bound": 1.5,
    "PLS": 1.0,
    "pELS,CTCF-bound": 1.2,
    "pELS": 0.8,
    "dELS,CTCF-bound": 1.0,
    "dELS": 0.6,
    "CTCF-only,CTCF-bound": 0.7,
    "DNase-H3K4me3,CTCF-bound": 0.9,
    "DNase-H3K4me3": 0.5,
}


def load_ccre_bed() -> dict[str, list[tuple[int,int,str,str,float]]]:
    """Returns {chrom: sorted list of (start, end, ccre_id, ccre_class, prior_score)}."""
    print(f"Loading {CCRE_BED} ...", flush=True)
    out: dict[str, list] = defaultdict(list)
    n = 0
    with CCRE_BED.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11: continue
            chrom = f[0]
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                continue
            ccre_id = f[3]
            cls = f[10]   # 11th column = class string
            score = CLASS_SCORE.get(cls, 0.4)  # default mild for unknown
            out[chrom].append((start, end, ccre_id, cls, score))
            n += 1
    for c in out:
        out[c].sort()
    print(f"  {n:,} cCRE elements across {len(out)} chromosomes", flush=True)
    return out


def augment_chrom_cache(chrom: str, cache: dict, ccres: list) -> dict:
    """Add `regulatory` + `prior_score` fields to each variant in the cache."""
    starts = [c[0] for c in ccres]
    n_aug = 0
    n_with_prior = 0
    for vid, info in cache.items():
        try:
            pos = int(vid.split(":", 2)[1])
        except (IndexError, ValueError):
            continue
        # Find any cCRE element overlapping pos
        # Use bisect to find candidate range
        idx = bisect.bisect_right(starts, pos) - 1
        regs = []
        max_score = 0.0
        for j in range(max(0, idx-3), min(len(ccres), idx+2)):
            s, e, cid, cls, sc = ccres[j]
            if s <= pos <= e:
                regs.append(f"cCRE:{cid}:{cls}")
                if sc > max_score:
                    max_score = sc
        if regs:
            info["regulatory"] = regs
            info["prior_score"] = float(max_score)
            n_with_prior += 1
        n_aug += 1
    return {"n_total": len(cache), "n_with_prior": n_with_prior}


def main():
    ccres = load_ccre_bed()
    caches = sorted(ANN.glob(CACHE_GLOB))
    print(f"\n{len(caches)} per-chrom cache files to augment\n")
    summary = []
    for cf in caches:
        chrom = cf.stem.replace("human_graph_cache_", "")
        # Skip the v2 outputs if they exist (don't double-augment)
        if "_v2_" in cf.name: continue
        if not chrom.startswith("chr"): continue
        cache = json.load(cf.open())
        ccre_list = ccres.get(chrom, [])
        if not ccre_list:
            print(f"  [{chrom}] no cCRE; skipping")
            continue
        stats = augment_chrom_cache(chrom, cache, ccre_list)
        out = ANN / f"{OUT_PREFIX}{chrom}.json".replace("v2_chrchr","v2_chr")
        out = ANN / f"human_graph_cache_v2_{chrom}.json"
        with out.open("w") as fh:
            json.dump(cache, fh)
        size_mb = out.stat().st_size / 1e6
        n_with = stats["n_with_prior"]
        n_total = stats["n_total"]
        pct = 100.0 * n_with / n_total if n_total else 0
        print(f"  [{chrom}] {n_with:,}/{n_total:,} variants got cCRE overlap "
              f"({pct:.0f}%); wrote {out.name} ({size_mb:.0f} MB)", flush=True)
        summary.append({"chrom": chrom, "n_total": n_total,
                        "n_with_prior": n_with, "frac": pct})
    # Save summary
    import pandas as pd
    pd.DataFrame(summary).to_csv(ANN / "human_graph_cache_v2_summary.tsv",
                                  sep="\t", index=False)


if __name__ == "__main__":
    main()
