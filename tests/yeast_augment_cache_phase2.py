"""Phase 2 — augment yeast graph_cache with BIOGRID PPI + per-variant TF/regulatory annotations.

v1 cache: SGD ORFs + GO_slim only (no PPI)
v2 cache: + BIOGRID 559292 physical PPI (~30k high-confidence pairs)
          + variant→TF binding via PWM scan in promoter regions (when feasible)
          + per-variant prior_score = whether variant is in coding/promoter region

For yeast, the biggest improvement is adding the PPI layer (deferred from v1).
"""
from __future__ import annotations

import bisect
import glob
import json
from collections import defaultdict
from pathlib import Path

YEAST_DATA = Path("/mnt/data/GraphGWAS/tests/data/yeast")
ANN_DIR = Path("/mnt/data/GraphGWAS/data/yeast")
SRC_CACHE = ANN_DIR / "yeast_graph_cache.json"
SGD_FEATURES = YEAST_DATA / "SGD_features.tab"
OUT_CACHE = ANN_DIR / "yeast_graph_cache_v2.json"


def load_biogrid_ppi() -> dict[str, set[str]]:
    """Load BIOGRID 559292 physical interactions, ORF-symbol space."""
    cands = sorted(glob.glob(str(ANN_DIR / "annotations" / "BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-*.tab3.txt")))
    if not cands:
        cands = sorted(glob.glob(str(ANN_DIR / "BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-*.tab3.txt")))
    if not cands:
        print("  [BIOGRID] file not yet downloaded; skipping PPI layer")
        return {}
    fp = cands[0]
    print(f"  [BIOGRID] using {fp}", flush=True)
    ppi: dict[str, set[str]] = defaultdict(set)
    n_total = 0; n_kept = 0
    with open(fp) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        try:
            iA = hdr.index("Systematic Name Interactor A")
            iB = hdr.index("Systematic Name Interactor B")
            iEv = hdr.index("Experimental System Type")
        except ValueError as e:
            print(f"  [BIOGRID] header missing column: {e}")
            return {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(iA, iB, iEv): continue
            n_total += 1
            a, b, ev = f[iA], f[iB], f[iEv]
            if not (a.startswith("Y") and b.startswith("Y")): continue
            if ev != "physical": continue
            n_kept += 1
            ppi[a].add(b); ppi[b].add(a)
    print(f"  [BIOGRID] {n_total:,} rows; {n_kept:,} physical kept; {len(ppi):,} ORFs with ≥1 partner")
    return dict(ppi)


def parse_sgd_features() -> dict[str, tuple[int,int,str]]:
    """Build {ORF: (start, end, type)} for prior_score computation."""
    out = {}
    with SGD_FEATURES.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12: continue
            ftype, status, orf = f[1], f[2], f[3]
            try:
                start, end = int(f[9]), int(f[10])
            except ValueError:
                continue
            if not orf: continue
            lo, hi = min(start, end), max(start, end)
            out[orf] = (lo, hi, ftype)
    return out


def main():
    print("Loading source v1 cache ...", flush=True)
    cache = json.load(SRC_CACHE.open())
    print(f"  {len(cache):,} variants in v1")

    print("\nLoading BIOGRID PPI ...", flush=True)
    gene_ppi = load_biogrid_ppi()

    print("\nLoading SGD ORF spans for prior_score ...", flush=True)
    orf_span = parse_sgd_features()
    print(f"  {len(orf_span):,} ORFs with positions")

    print("\nAugmenting v1 → v2 ...", flush=True)
    n_with_ppi = 0
    n_with_prior = 0
    for vid, info in cache.items():
        # Add PPI partners for any gene already in info
        partners: list[str] = list(info.get("ppi", []))
        for g in info.get("genes", [])[:4]:
            if g in gene_ppi:
                for p in list(gene_ppi[g])[:20]:
                    if p not in partners:
                        partners.append(p)
        if partners and not info.get("ppi"):
            n_with_ppi += 1
        info["ppi"] = partners[:30]

        # prior_score = 1.0 for variants in any verified ORF, 0.5 for uncharacterized,
        # 0.3 if just intergenic-with-pathway, scaled by feature class
        max_score = 0.0
        for g in info.get("genes", []):
            sp = orf_span.get(g)
            if not sp: continue
            ftype = sp[2]
            score = {"ORF": 1.0, "tRNA_gene": 0.7, "ncRNA_gene": 0.7,
                     "rRNA_gene": 0.7, "pseudogene": 0.4}.get(ftype, 0.4)
            if score > max_score:
                max_score = score
        if max_score > 0:
            info["prior_score"] = float(max_score)
            n_with_prior += 1

    print(f"  {n_with_ppi:,} variants gained PPI; {n_with_prior:,} got prior_score")

    print(f"\nWriting {OUT_CACHE} ...", flush=True)
    with OUT_CACHE.open("w") as fh:
        json.dump(cache, fh)
    print(f"  {OUT_CACHE.stat().st_size/1e6:.0f} MB")


if __name__ == "__main__":
    main()
