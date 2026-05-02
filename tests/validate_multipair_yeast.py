"""§Y.4 multi-pair extension — sweep all Tier-A 1011-segregating yeast pairs.

The original §Y.4 benchmark tested only the BCY1×TPK1 anchor.  This script
sweeps the five Tier-A panel-segregating yeast pairs from the ground-truth
catalogue (``docs/groundtruth/epistasis_pairs.json``):

  1. BCY1 × TPK1   — anchor (cAMP-PKA reg:cat)
  2. WHI5 × CLN3   — Goldstein 2025-confirmed segregating; G1/S commitment
  3. MKT1 × GPA1   — natural hub-QTL pleiotropic with multiple loci
  4. MIP1 × SAL1   — mitonuclear epistasis (mt-DNA polymerase × ADP/ATP carrier)
  5. TKL1 × NQM1   — pentose-phosphate paralog redundancy (CellMap-validated)

For each pair, simulate Scenario A (positive-interaction quantitative DGP,
β=3, R² ≈ 10%) and run all 5 methods (M1-M5).  Produce a 5-pair × 5-method
ground-truth recovery rank grid.

Canonical-PPI injection ON by default (else most Tier-A pairs miss because
the cache PPI is incomplete; see §Y.5).

Usage:
    python tests/validate_multipair_yeast.py [--scenario A] [--no-inject-ppi]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src" / "python"))
sys.path.insert(0, str(REPO_ROOT / "tests"))

from graphgwas.canonical_ppi import (  # noqa: E402
    YEAST_CANONICAL_EDGES,
    inject_canonical_ppi_edges,
)
from graphgwas.epistasis_v2 import (  # noqa: E402
    InteractionResult,
    dark_matter_epistasis_from_data,
    differential_subgraph_from_data,
    ld_pruned_cooccurrence_from_data,
    motif_filtered_epistasis_from_data,
)
from validate_m2_yeast import (  # noqa: E402
    YEAST_CACHE, YEAST_VCF,
    _bipartite_from_dict, _run_m5_then_test, _subsample_with_truth_priority,
    load_yeast_chroms, simulate_with_truth,
)


def determine_required_chroms(cache: dict, pairs: list[dict]) -> list[str]:
    """Find every chromosome that hosts any Tier-A pair gene.

    Walks the cache to find each gene's chromosome (every variant in that
    gene shares a chromosome).  Returns sorted unique chromosome names
    (``chromosome9``, ``chromosome10``, etc).
    """
    needed_genes = set()
    for p in pairs:
        needed_genes.add(p["gene_1"]["systematic"])
        needed_genes.add(p["gene_2"]["systematic"])
    chr_for_gene: dict[str, str] = {}
    for vid, entry in cache.items():
        for g in entry.get("genes", []):
            if g in needed_genes and g not in chr_for_gene:
                chr_for_gene[g] = vid.split(":", 1)[0]
                if len(chr_for_gene) == len(needed_genes):
                    break
        if len(chr_for_gene) == len(needed_genes):
            break
    return sorted(set(chr_for_gene.values()))


RESULTS_DIR = REPO_ROOT / "results" / "paper2_epistasis"
GROUNDTRUTH_PATH = REPO_ROOT / "docs" / "groundtruth" / "epistasis_pairs.json"


def find_recovery_rank(results: list[InteractionResult],
                       gene_1_variants: set[str],
                       gene_2_variants: set[str]) -> int | None:
    for r, ir in enumerate(results, 1):
        v1, v2 = ir.variant_1, ir.variant_2
        if (v1 in gene_1_variants and v2 in gene_2_variants) or \
           (v1 in gene_2_variants and v2 in gene_1_variants):
            return r
    return None


def best_variant_for_gene(cache: dict,
                          gene_systematic: str,
                          variant_id_to_col: dict,
                          dosages: np.ndarray) -> tuple[str, int] | tuple[None, None]:
    """Return (variant_id, dosage_col) for the highest-MAF variant in the gene."""
    candidates = [
        (vid, idx) for vid, idx in variant_id_to_col.items()
        if gene_systematic in cache.get(vid, {}).get("genes", [])
    ]
    if not candidates:
        return None, None
    af = np.array([np.nanmean(dosages[:, idx]) / 2.0 for _, idx in candidates])
    maf = np.minimum(af, 1.0 - af)
    k = int(np.argmax(maf))
    return candidates[k][0], candidates[k][1]


def run_one_pair(pair: dict, dosages, variant_ids, sample_ids, cache,
                 motifs_list: list[str], beta: float,
                 rng: np.random.Generator) -> dict:
    """Run all 5 methods on one pair under Scenario A.  Return per-method ranks."""
    g1_sys = pair["gene_1"]["systematic"]
    g1_name = pair["gene_1"]["standard"]
    g2_sys = pair["gene_2"]["systematic"]
    g2_name = pair["gene_2"]["standard"]
    print(f"\n{'=' * 70}")
    print(f"  PAIR: {g1_name} ({g1_sys}) × {g2_name} ({g2_sys})  [{pair['id']}]")
    print(f"{'=' * 70}")

    variant_id_to_col = {v: i for i, v in enumerate(variant_ids)}
    g1v, g1_idx = best_variant_for_gene(cache, g1_sys, variant_id_to_col, dosages)
    g2v, g2_idx = best_variant_for_gene(cache, g2_sys, variant_id_to_col, dosages)
    if g1v is None or g2v is None:
        print(f"  ✗ Skipping {pair['id']} — no representative variant for "
              f"{g1_sys if g1v is None else g2_sys}")
        return {"pair_id": pair["id"], "g1": g1_name, "g2": g2_name,
                "skipped": True, "reason": "no representative variant"}

    af_g1 = np.nanmean(dosages[:, g1_idx]) / 2.0
    af_g2 = np.nanmean(dosages[:, g2_idx]) / 2.0
    print(f"  {g1_name}: {g1v}  MAF={min(af_g1, 1-af_g1):.3f}")
    print(f"  {g2_name}: {g2v}  MAF={min(af_g2, 1-af_g2):.3f}")

    # Recovery sets
    gene1_variants = {vid for vid in variant_ids
                       if g1_sys in cache.get(vid, {}).get("genes", [])}
    gene2_variants = {vid for vid in variant_ids
                       if g2_sys in cache.get(vid, {}).get("genes", [])}
    print(f"  Recovery sets: {g1_name}={len(gene1_variants)} variants, "
          f"{g2_name}={len(gene2_variants)} variants")

    # Phenotype simulation (Scenario A)
    truth = {"gene_1_idx": g1_idx, "gene_2_idx": g2_idx,
             "gene_1_variant": g1v, "gene_2_variant": g2v,
             "gene_1_systematic": g1_sys, "gene_2_systematic": g2_sys,
             "gene_1_symbol": g1_name, "gene_2_symbol": g2_name,
             "shared_pathways": []}
    phenotype = simulate_with_truth(dosages, truth, beta=beta, n_nuisance=0, rng=rng)

    out: dict = {
        "pair_id": pair["id"], "g1": g1_name, "g2": g2_name,
        "g1_systematic": g1_sys, "g2_systematic": g2_sys,
        "g1_variant": g1v, "g1_maf": float(min(af_g1, 1 - af_g1)),
        "g2_variant": g2v, "g2_maf": float(min(af_g2, 1 - af_g2)),
        "n_g1_variants": len(gene1_variants),
        "n_g2_variants": len(gene2_variants),
        "methods": {},
    }

    # ----- M1 -----
    print(f"\n  [M1]")
    t0 = time.time()
    af_vec = np.nanmean(dosages, axis=0) / 2.0
    m1_variants = [
        {"variantId": variant_ids[k], "pos": int(variant_ids[k].split(":")[1]),
         "af_total": float(af_vec[k])}
        for k in range(len(variant_ids))
    ]
    m1_dosage_list = [dosages[:, k] for k in range(dosages.shape[1])]
    m1_res = ld_pruned_cooccurrence_from_data(
        variants=m1_variants, dosage_list=m1_dosage_list, phenotype=phenotype,
        r2_prune=0.5, min_cocarriers=1, min_distance_bp=10_000,
        max_variants=2000, verbose=False,
    )
    rank = find_recovery_rank(m1_res, gene1_variants, gene2_variants)
    out["methods"]["M1"] = {"rank": rank, "n_results": len(m1_res),
                             "runtime_sec": time.time() - t0}
    print(f"    M1 rank: {rank} ({len(m1_res)} pairs in {out['methods']['M1']['runtime_sec']:.1f}s)")

    # ----- M2 -----
    print(f"  [M2 motifs={'+'.join(motifs_list)}]")
    t0 = time.time()
    m2_res = motif_filtered_epistasis_from_data(
        dosages=dosages, variant_ids=variant_ids, phenotype=phenotype,
        graph_cache=cache, motifs=motifs_list, mac_min=10, correction="BH",
        max_pairs_per_entity=500, max_pairs_total=200_000, verbose=False,
    )
    rank = find_recovery_rank(m2_res, gene1_variants, gene2_variants)
    out["methods"]["M2"] = {"rank": rank, "n_results": len(m2_res),
                             "runtime_sec": time.time() - t0}
    print(f"    M2 rank: {rank} ({len(m2_res)} pairs in {out['methods']['M2']['runtime_sec']:.1f}s)")

    # ----- M3 -----
    print(f"  [M3]")
    t0 = time.time()
    q75 = np.percentile(phenotype, 75)
    q25 = np.percentile(phenotype, 25)
    case_mask = phenotype >= q75
    ctrl_mask = phenotype <= q25
    sub_idx = _subsample_with_truth_priority(
        variant_ids, gene1_variants, gene2_variants, n_target=1500, rng=rng,
    )
    sub_dos = dosages[:, sub_idx]
    sub_vids = [variant_ids[k] for k in sub_idx]
    m3_res = differential_subgraph_from_data(
        dosages=sub_dos, variant_ids=sub_vids,
        case_mask=case_mask, control_mask=ctrl_mask,
        r2_prune=0.5, min_distance_bp=10_000, min_cocarriers=2, mac_min=10,
        max_variants=500, verbose=False,
    )
    rank = find_recovery_rank(m3_res, gene1_variants, gene2_variants)
    out["methods"]["M3"] = {"rank": rank, "n_results": len(m3_res),
                             "runtime_sec": time.time() - t0}
    print(f"    M3 rank: {rank} ({len(m3_res)} pairs in {out['methods']['M3']['runtime_sec']:.1f}s)")

    # ----- M4 -----
    print(f"  [M4]")
    t0 = time.time()
    m4_res = dark_matter_epistasis_from_data(
        dosages=sub_dos, variant_ids=sub_vids,
        case_mask=case_mask, control_mask=ctrl_mask,
        maf_min=0.05, depletion_threshold=0.7, min_expected=0.5,
        mac_min=10, max_variants=500, verbose=False,
    )
    rank = find_recovery_rank(m4_res, gene1_variants, gene2_variants)
    out["methods"]["M4"] = {"rank": rank, "n_results": len(m4_res),
                             "runtime_sec": time.time() - t0}
    print(f"    M4 rank: {rank} ({len(m4_res)} pairs in {out['methods']['M4']['runtime_sec']:.1f}s)")

    # ----- M5 -----
    print(f"  [M5]")
    t0 = time.time()
    m5_res = _run_m5_then_test(
        dosages=dosages, variant_ids=variant_ids, phenotype=phenotype,
        cache=cache, truth=truth, alpha=0.15, n_seeds_extra="local_chroms",
        top_k=2000, mac_min=10, rng=rng,
    )
    rank = find_recovery_rank(m5_res, gene1_variants, gene2_variants)
    out["methods"]["M5"] = {"rank": rank, "n_results": len(m5_res),
                             "runtime_sec": time.time() - t0}
    print(f"    M5 rank: {rank} ({len(m5_res)} pairs in {out['methods']['M5']['runtime_sec']:.1f}s)")

    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--beta", type=float, default=3.0)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--motifs", default="same_gene,same_pathway,protein_interaction")
    p.add_argument("--no-inject-ppi", action="store_true",
                   help="Skip canonical PPI injection (most Tier-A pairs "
                        "will fail without it).")
    p.add_argument("--limit", type=int, default=None,
                   help="Run only the first N pairs (for testing)")
    args = p.parse_args()
    motifs_list = [m.strip() for m in args.motifs.split(",") if m.strip()]

    print("=== §Y.4 multi-pair benchmark on yeast Scenario A ===\n")

    # Load catalog
    catalog = json.loads(GROUNDTRUTH_PATH.read_text())
    yeast_pairs = [
        p for p in catalog["species"]["yeast"]["pairs"]
        if p.get("panel_natural_variation", False)
    ]
    if args.limit:
        yeast_pairs = yeast_pairs[:args.limit]
    print(f"Pairs to test: {len(yeast_pairs)}")
    for p_ in yeast_pairs:
        print(f"  - {p_['id']} ({p_['gene_1']['standard']}×{p_['gene_2']['standard']})")
    print()

    # Load cache first (needed to determine required chromosomes)
    print("Loading cache...")
    cache = json.loads(YEAST_CACHE.read_text())

    # Determine which chromosomes to load based on the pairs' gene locations
    chroms = determine_required_chroms(cache, yeast_pairs)
    print(f"Required chromosomes: {chroms}")

    # Load genotypes for the union of needed chromosomes
    print("Loading genotypes...")
    dosages, variant_ids, sample_ids = load_yeast_chroms(
        YEAST_VCF, chroms, maf_min=0.05,
    )
    if not args.no_inject_ppi:
        cache, diag = inject_canonical_ppi_edges(
            cache, edges=YEAST_CANONICAL_EDGES, inplace=True, verbose=False,
        )
        print(f"Canonical PPI injected: {diag['n_edges_applied']}/{diag['n_edges_supplied']} "
              f"edges, {diag['n_variants_modified']} variants modified")

    # Run each pair
    rng = np.random.default_rng(args.seed)
    results = {}
    for pair in yeast_pairs:
        results[pair["id"]] = run_one_pair(
            pair, dosages, variant_ids, sample_ids, cache,
            motifs_list, args.beta, rng,
        )

    # Print final grid
    print(f"\n\n{'=' * 70}")
    print(f"  FINAL: PAIR × METHOD RECOVERY RANK GRID (Scenario A, β={args.beta})")
    print(f"{'=' * 70}\n")
    methods = ["M1", "M2", "M3", "M4", "M5"]
    print(f"  {'Pair':<32} " + " ".join(f"{m:>8}" for m in methods))
    print(f"  {'-' * 32} " + " ".join("-" * 8 for _ in methods))
    for pair_id, data in results.items():
        if data.get("skipped"):
            cells = [f"{'SKIP':>8}"] * len(methods)
        else:
            cells = []
            for m in methods:
                r = data["methods"].get(m, {}).get("rank")
                cells.append(f"{(str(r) if r else 'NF'):>8}")
        label = f"{data.get('g1','?')}×{data.get('g2','?')}"
        print(f"  {label:<32} " + " ".join(cells))
    print()

    # Save
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out_path = RESULTS_DIR / "y4_multipair_yeast.json"
    payload = {
        "scenario": "A",
        "params": {"beta": args.beta, "motifs": motifs_list,
                   "inject_canonical_ppi": not args.no_inject_ppi,
                   "seed": args.seed},
        "n_samples": int(dosages.shape[0]),
        "n_variants_loaded": int(dosages.shape[1]),
        "results": results,
    }
    out_path.write_text(json.dumps(payload, indent=2, default=str))
    print(f"  Wrote {out_path}")


if __name__ == "__main__":
    main()
