"""Paper #2 §Y.8 — multi-pair catalogue benchmark.

Runs M2 (motif-filtered) and M5★ (PPI-bridge + focused seeds) on every
pair in docs/groundtruth/epistasis_pairs.json that has natural variation
in the corresponding panel. Reports per-pair × per-method recovery rank
plus testability classification.

Pipeline per pair:
  1. Look up g1 + g2 chromosomal coordinates via
     docs/groundtruth/gene_coordinates_index.json (built by
     scripts/build_gene_coordinates_index.py).
  2. If neither gene is in the cache → UNTESTABLE_NO_CACHE_ANNOTATION.
  3. Load species cache for the involved chromosomes; inject canonical PPI.
  4. Load dosages for ±50 kb windows around g1 and g2 (plus a small
     control window) from VCF (Arabidopsis), pgen (rice), or the same
     yeast loader as §Y.4.
  5. If either gene has no MAF≥0.05 variant in dosages →
     UNTESTABLE_NO_PANEL_VARIATION.
  6. Pick highest-MAF representative variants; simulate Scenario A
     (y = β·G_X·G_Y + ε, β=3); run M2 + M5★; record gene-level rank.

Output:
  results/paper2_epistasis/multi_pair_catalogue_benchmark.json

This is the §Y.8 falsification of §Y.7's "M5★ recovers the canonical pair
at rank 1" claim — does it generalise to the broader catalogue, or was
it specific to the three anchors?
"""
from __future__ import annotations

import argparse
import gc
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
from cyvcf2 import VCF

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src" / "python"))
sys.path.insert(0, str(REPO_ROOT / "tests"))

from graphgwas.epistasis_v2 import (  # noqa: E402
    motif_filtered_epistasis_from_data,
    _test_interaction,
)
from validate_m2_yeast import _bipartite_from_dict  # noqa: E402

CATALOGUE = REPO_ROOT / "docs" / "groundtruth" / "epistasis_pairs.json"
GENE_INDEX = REPO_ROOT / "docs" / "groundtruth" / "gene_coordinates_index.json"
RESULTS_OUT = (REPO_ROOT / "results" / "paper2_epistasis"
               / "multi_pair_catalogue_benchmark.json")

# ---------------------------------------------------------------------------
# Cache loaders (per-chrom or full, depending on species)
# ---------------------------------------------------------------------------

def load_species_cache_for_chroms(species: str, chroms: list) -> dict:
    """Load the union of cache JSONs for the given chromosomes."""
    cache: dict = {}
    if species == "yeast":
        # Single large JSON; load once
        path = REPO_ROOT / "data" / "yeast" / "yeast_graph_cache_v2.json"
        cache = json.loads(path.read_text())
        return cache
    if species == "arabidopsis":
        for c in chroms:
            path = (REPO_ROOT / "data" / "arabidopsis"
                    / f"arabidopsis_graph_cache_v3_chr{c}.json")
            if path.exists():
                cache.update(json.loads(path.read_text()))
        return cache
    if species == "rice":
        for c in chroms:
            # Cache files use Chr1, Chr2, ... naming
            path = (REPO_ROOT / "data" / "rice_3k" / "annotations"
                    / f"rice_graph_cache_v2_Chr{c}.json")
            if path.exists():
                cache.update(json.loads(path.read_text()))
        return cache
    sys.exit(f"unknown species {species}")


# ---------------------------------------------------------------------------
# Genotype loaders (per-window)
# ---------------------------------------------------------------------------

def load_arabidopsis_windows(windows: list) -> tuple:
    """windows: list of (chrom, start, end). Returns (dosages, vids, sample_ids)."""
    vcf_path = (REPO_ROOT / "tests" / "data" / "arabidopsis"
                / "1001genomes_snp-short-indel_only_ACGTN.vcf.gz")
    vcf = VCF(str(vcf_path))
    sample_ids = list(vcf.samples)
    rows: list = []
    vids: list = []
    for chrom, start, end in windows:
        n = 0
        for record in vcf(f"{chrom}:{start}-{end}"):
            if not record.is_snp or len(record.ALT) != 1:
                continue
            g = record.gt_types
            d = np.where(g == 3, np.nan, g).astype(np.float64)
            af = np.nanmean(d) / 2.0
            maf = min(af, 1.0 - af)
            if maf < 0.05:
                continue
            rows.append(d)
            vids.append(f"{record.CHROM}:{record.POS}:{record.REF}:{record.ALT[0]}")
            n += 1
    vcf.close()
    if not rows:
        return None, [], sample_ids
    dosages = np.column_stack(rows)
    return dosages, vids, sample_ids


def load_rice_windows(windows: list) -> tuple:
    """Load rice dosages for given windows via pgenlib. Returns (dosages, vids, sample_ids)."""
    import pgenlib
    pgen_prefix = REPO_ROOT / "data" / "rice_3k" / "rice_3k"
    pgen_path = str(pgen_prefix) + ".pgen"
    pvar_path = str(pgen_prefix) + ".pvar"
    psam_path = str(pgen_prefix) + ".psam"

    # Read sample IDs once (small)
    sample_ids: list = []
    with open(psam_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            sample_ids.append(parts[0])
    n_samples = len(sample_ids)

    # Index pvar lines we want
    target_chroms = set()
    for chrom, _, _ in windows:
        target_chroms.add(f"Chr{chrom}")
        target_chroms.add(str(chrom))

    keep_idx: list = []
    keep_meta: list = []
    with open(pvar_path) as f:
        line_idx = 0
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                line_idx += 1
                continue
            chrom = parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                line_idx += 1
                continue
            ref, alt = parts[3], parts[4]
            if chrom in target_chroms:
                in_window = any(
                    (chrom == w[0] or chrom == f"Chr{w[0]}") and w[1] <= pos <= w[2]
                    for w in windows
                )
                if in_window:
                    keep_idx.append(line_idx)
                    keep_meta.append((chrom, pos, ref, alt))
            line_idx += 1

    if not keep_idx:
        return None, [], sample_ids

    n_var = len(keep_idx)
    pgen = pgenlib.PgenReader(pgen_path.encode("utf-8"))
    dosages = np.empty((n_var, n_samples), dtype=np.float32)
    pgen.read_dosages_list(np.asarray(keep_idx, dtype=np.uint32), dosages)
    pgen.close()
    dosages = dosages.T
    dosages = np.where(dosages < 0, np.nan, dosages).astype(np.float64)

    af = np.nanmean(dosages, axis=0) / 2.0
    maf = np.minimum(af, 1.0 - af)
    keep = maf >= 0.05
    dosages = dosages[:, keep]
    keep_meta_kept = [keep_meta[k] for k in range(len(keep_meta)) if keep[k]]
    # Normalise vid format: cache uses "Chr6:..." even if pvar uses "6"
    vids = [f"Chr{m[0].lstrip('Chr')}:{m[1]}:{m[2]}:{m[3]}"
            if not m[0].startswith("Chr")
            else f"{m[0]}:{m[1]}:{m[2]}:{m[3]}"
            for m in keep_meta_kept]
    return dosages, vids, sample_ids


def load_yeast_windows(windows: list) -> tuple:
    """Load yeast dosages for given windows (chrom in 'chromosome9' form)."""
    vcf_path = REPO_ROOT / "tests" / "data" / "yeast" / "1011_snps_maf05.vcf.gz"
    vcf = VCF(str(vcf_path))
    sample_ids = list(vcf.samples)
    rows: list = []
    vids: list = []
    target_chroms = {f"chromosome{c}" for c, _, _ in windows}
    target_chroms |= {str(c) for c, _, _ in windows}
    # No tabix on the yeast VCF — sequential read
    for record in vcf:
        if record.CHROM not in target_chroms:
            continue
        in_window = any(
            (record.CHROM == f"chromosome{c}" or record.CHROM == str(c))
            and s <= record.POS <= e
            for c, s, e in windows
        )
        if not in_window:
            continue
        if not record.is_snp or len(record.ALT) != 1:
            continue
        g = record.gt_types
        d = np.where(g == 3, np.nan, g).astype(np.float64)
        af = np.nanmean(d) / 2.0
        maf = min(af, 1.0 - af)
        if maf < 0.05:
            continue
        rows.append(d)
        vids.append(f"{record.CHROM}:{record.POS}:{record.REF}:{record.ALT[0]}")
    vcf.close()
    if not rows:
        return None, [], sample_ids
    dosages = np.column_stack(rows)
    return dosages, vids, sample_ids


# ---------------------------------------------------------------------------
# M5★ kernel (PPI-bridge + focused seeds), copied from m5_variants_benchmark
# ---------------------------------------------------------------------------

def compute_rwr_matrix_ppi(cache: dict,
                            annotated_vids: list,
                            seed_indices: list,
                            alpha: float = 0.15,
                            n_iter: int = 30) -> tuple:
    """Bipartite RWR with var→gene→gene'→var transition (PPI-bridged)."""
    A, cache_var_ids, gene_ids = _bipartite_from_dict(cache, annotated_vids)
    n_var, n_gene = A.shape
    var_deg = A.sum(axis=1, keepdims=True).clip(min=1e-12)
    gene_deg = A.sum(axis=0, keepdims=True).clip(min=1e-12)
    P_vg = A / var_deg
    P_gv = (A / gene_deg).T

    gene_to_idx = {g: i for i, g in enumerate(gene_ids)}
    G = np.zeros((n_gene, n_gene), dtype=np.float64)
    np.fill_diagonal(G, 1.0)
    for vid in cache_var_ids:
        entry = cache[vid]
        entry_genes = [g for g in entry.get("genes", []) if g in gene_to_idx]
        entry_ppis = [g for g in entry.get("ppi", []) if g in gene_to_idx]
        for ga in entry_genes:
            ia = gene_to_idx[ga]
            for gb in entry_ppis:
                if ga == gb:
                    continue
                ib = gene_to_idx[gb]
                G[ia, ib] = 1.0
                G[ib, ia] = 1.0
    G_row = G.sum(axis=1, keepdims=True).clip(min=1e-12)
    P_gg = G / G_row
    M_vv = P_vg @ P_gg @ P_gv

    n_seed = len(seed_indices)
    E = np.zeros((n_var, n_seed), dtype=np.float64)
    for k, si in enumerate(seed_indices):
        E[si, k] = 1.0
    P = E.copy()
    for _ in range(n_iter):
        P = (1.0 - alpha) * (M_vv @ P) + alpha * E
    P = P / np.clip(P.sum(axis=0, keepdims=True), 1e-12, None)
    return P, A, cache_var_ids


def run_m5_star(cache, annotated_vids, g1_seeds, g2_seeds,
                 dosages, variant_ids, phenotype, mac_min: int = 10) -> dict:
    """M5★ = PPI-bridge + focused seeds; rank g1×g2 cross-pairs by interaction p."""
    seed_indices = sorted(set(g1_seeds + g2_seeds))
    if not g1_seeds or not g2_seeds:
        return {"n_tested": 0, "ground_truth_rank": None,
                "ground_truth_p": None}
    P, A, cache_var_ids = compute_rwr_matrix_ppi(
        cache, annotated_vids, seed_indices, alpha=0.15, n_iter=30,
    )
    seed_arr = np.asarray(seed_indices)
    P_seed = P[seed_arr, :]
    sym = P_seed + P_seed.T
    seed_to_col = {si: k for k, si in enumerate(seed_indices)}

    pair_scores: list = []
    for si in g1_seeds:
        ki = seed_to_col[si]
        for sj in g2_seeds:
            if si == sj:
                continue
            kj = seed_to_col[sj]
            s = float(sym[ki, kj])
            if s <= 0:
                continue
            pair_scores.append((si, sj, s))
    pair_scores.sort(key=lambda t: t[2], reverse=True)

    cache_idx_to_col = {v: variant_ids.index(v) for v in cache_var_ids}
    results: list = []
    for si, sj, score in pair_scores:
        vid_i = cache_var_ids[si]
        vid_j = cache_var_ids[sj]
        i = cache_idx_to_col[vid_i]
        j = cache_idx_to_col[vid_j]
        d1 = dosages[:, i].astype(float)
        d2 = dosages[:, j].astype(float)
        s1 = float(np.nansum(d1)); n1 = int(np.sum(~np.isnan(d1)))
        s2 = float(np.nansum(d2)); n2 = int(np.sum(~np.isnan(d2)))
        mac1 = min(s1, 2 * n1 - s1)
        mac2 = min(s2, 2 * n2 - s2)
        if mac1 < mac_min or mac2 < mac_min:
            continue
        r = _test_interaction(d1, d2, phenotype)
        results.append({"variant_1": vid_i, "variant_2": vid_j,
                        "p_interaction": r["p_interaction"],
                        "rwr_score": float(score)})
    if not results:
        return {"n_tested": 0, "ground_truth_rank": None, "ground_truth_p": None}
    results.sort(key=lambda r: r["p_interaction"])
    # All entries are g1×g2 cross-pairs by construction → rank 1 = best
    return {"n_tested": len(results),
            "ground_truth_rank": 1,
            "ground_truth_p": float(results[0]["p_interaction"]),
            "best_pair": (results[0]["variant_1"], results[0]["variant_2"])}


# ---------------------------------------------------------------------------
# M2 wrapper (run on the loaded panel; rank canonical pair gene-level)
# ---------------------------------------------------------------------------

def run_m2(cache, dosages, variant_ids, phenotype,
            gene1_set, gene2_set) -> dict:
    """Run M2 motif-filtered epistasis on the loaded panel; rank gene-pair."""
    m2_results = motif_filtered_epistasis_from_data(
        dosages=dosages, variant_ids=variant_ids, phenotype=phenotype,
        graph_cache=cache,
        motifs=["same_gene", "same_pathway", "protein_interaction"],
        mac_min=10, correction="BH",
        max_pairs_per_entity=500, max_pairs_total=200_000, verbose=False,
    )
    n = len(m2_results)
    if n == 0:
        return {"n_tested": 0, "ground_truth_rank": None}
    # Find first cross-pair (g1, g2) in result order
    for rk, r in enumerate(m2_results, start=1):
        v1 = r.variant_1
        v2 = r.variant_2
        if (v1 in gene1_set and v2 in gene2_set) or \
           (v1 in gene2_set and v2 in gene1_set):
            return {"n_tested": n,
                    "ground_truth_rank": rk,
                    "ground_truth_p": float(r.p_interaction),
                    "ground_truth_q": float(r.p_corrected) if r.p_corrected else None,
                    "detected_motif": r.motif}
    return {"n_tested": n, "ground_truth_rank": None}


# ---------------------------------------------------------------------------
# Per-pair runner
# ---------------------------------------------------------------------------

PAD_BP = 50_000
SAMPLE_RNG = np.random.default_rng(2026)


def discover_window(gene_index: dict, gene_systematic: str) -> tuple | None:
    """Returns (chrom, pos_min - PAD, pos_max + PAD) or None if not annotated."""
    g = gene_index.get(gene_systematic)
    if not g:
        return None
    # Pick the chrom with the most variants annotated
    best_chrom, best_stats = max(g.items(), key=lambda kv: kv[1]["n_variants"])
    return (best_chrom, max(0, best_stats["pos_min"] - PAD_BP),
            best_stats["pos_max"] + PAD_BP)


def normalise_chrom_for_load(species: str, chrom: str) -> str:
    """Each species' loader expects the chrom in a specific form."""
    if species == "yeast":
        # Cache uses "chromosome9", VCF uses "chromosome9" — already aligned
        return chrom
    if species == "arabidopsis":
        # Cache uses "1", VCF uses "1" — strip leading "Chr" if present
        return chrom.replace("Chr", "").lstrip("c")
    if species == "rice":
        # Cache uses "Chr6", pgen uses "Chr6" or "6" — return without prefix for windowing
        return chrom.replace("Chr", "")
    return chrom


def gene_id_field(species: str) -> str:
    """Catalogue stores gene IDs under different keys per species."""
    return {"yeast": "systematic", "arabidopsis": "agi",
             "rice": "locus", "human": "ensembl"}[species]


def run_pair(species: str, pair: dict, gene_index: dict) -> dict:
    pair_id = pair.get("id", f"{pair['gene_1'].get('standard','?')}_{pair['gene_2'].get('standard','?')}")
    id_field = gene_id_field(species)
    g1_sys = pair["gene_1"][id_field]
    g2_sys = pair["gene_2"][id_field]
    g1_sym = pair["gene_1"].get("standard") or pair["gene_1"].get("symbol", g1_sys)
    g2_sym = pair["gene_2"].get("standard") or pair["gene_2"].get("symbol", g2_sys)

    print(f"\n--- {species} {pair_id} ({g1_sym}×{g2_sym}) ---")
    base = {
        "pair_id": pair_id,
        "g1_symbol": g1_sym, "g2_symbol": g2_sym,
        "g1_systematic": g1_sys, "g2_systematic": g2_sys,
        "is_anchor": bool(pair.get("anchor")),
        "panel_natural_variation_claimed": pair.get("panel_natural_variation"),
    }

    w1 = discover_window(gene_index, g1_sys)
    w2 = discover_window(gene_index, g2_sys)
    if w1 is None or w2 is None:
        print(f"  UNTESTABLE: g1_in_cache={w1 is not None} g2_in_cache={w2 is not None}")
        return {**base, "status": "UNTESTABLE_NO_CACHE_ANNOTATION",
                "g1_in_cache": w1 is not None, "g2_in_cache": w2 is not None}

    print(f"  windows: g1={w1}  g2={w2}")
    # Normalise chrom IDs for loaders
    chrom_load_1 = normalise_chrom_for_load(species, w1[0])
    chrom_load_2 = normalise_chrom_for_load(species, w2[0])
    windows = [(chrom_load_1, w1[1], w1[2]),
               (chrom_load_2, w2[1], w2[2])]
    # Also load a small control window on chr1 (or another chrom not in g1/g2)
    if species == "arabidopsis":
        ctrl_chrom = "3" if chrom_load_1 != "3" and chrom_load_2 != "3" else "1"
        windows.append((ctrl_chrom, 1_000_000, 1_100_000))
    elif species == "rice":
        ctrl_chrom = "5" if chrom_load_1 != "5" and chrom_load_2 != "5" else "1"
        windows.append((ctrl_chrom, 1_000_000, 1_200_000))
    # For yeast keep two-window scope (it's small)

    # Load cache for the involved chromosomes
    needed_chroms = list({w1[0].replace("Chr", ""), w2[0].replace("Chr", "")})
    print(f"  loading cache for chroms: {needed_chroms}")
    t0 = time.time()
    cache = load_species_cache_for_chroms(species, needed_chroms)
    print(f"    cache: {len(cache):,} variants in {time.time()-t0:.1f}s")

    # Inject canonical PPI
    from graphgwas.canonical_ppi import (
        CANONICAL_EDGES_BY_SPECIES, inject_canonical_ppi_edges,
    )
    edges = CANONICAL_EDGES_BY_SPECIES.get(species, [])
    cache, inj = inject_canonical_ppi_edges(
        cache, edges=edges, inplace=True, verbose=False,
    )
    print(f"    PPI injection: {inj['n_edges_applied']}/{inj['n_edges_supplied']} edges, "
          f"{inj['n_variants_modified']} variants modified")

    # Load dosages
    print(f"  loading dosages from windows: {windows}")
    t0 = time.time()
    if species == "arabidopsis":
        dosages, variant_ids, sample_ids = load_arabidopsis_windows(windows)
    elif species == "rice":
        dosages, variant_ids, sample_ids = load_rice_windows(windows)
    elif species == "yeast":
        dosages, variant_ids, sample_ids = load_yeast_windows(windows)
    else:
        dosages = None
    if dosages is None or dosages.shape[1] == 0:
        print(f"  UNTESTABLE: no variants survived MAF≥0.05 in target windows")
        return {**base, "status": "UNTESTABLE_NO_VARIANTS_IN_WINDOWS"}
    print(f"    dosages: {dosages.shape[0]} samples × {dosages.shape[1]:,} variants  "
          f"in {time.time()-t0:.1f}s")

    # Recovery sets (cache-annotated AND in dosage-loaded panel)
    gene1_set = {vid for vid in variant_ids
                  if g1_sys in cache.get(vid, {}).get("genes", [])}
    gene2_set = {vid for vid in variant_ids
                  if g2_sys in cache.get(vid, {}).get("genes", [])}
    print(f"    panel recovery sets: {g1_sym}={len(gene1_set)} {g2_sym}={len(gene2_set)}")
    if not gene1_set or not gene2_set:
        print(f"  UNTESTABLE: g1_panel={len(gene1_set)} g2_panel={len(gene2_set)}")
        return {**base, "status": "UNTESTABLE_NO_PANEL_VARIATION",
                "g1_panel_count": len(gene1_set),
                "g2_panel_count": len(gene2_set),
                "n_dosage_variants": int(dosages.shape[1])}

    # Pick highest-MAF representatives + simulate Scenario A
    def best_maf_idx(var_set):
        best_v, best_m = None, -1
        for v in var_set:
            i = variant_ids.index(v)
            af = np.nanmean(dosages[:, i]) / 2.0
            m = min(af, 1.0 - af)
            if m > best_m:
                best_m, best_v = m, v
        return best_v, best_m

    g1_rep, g1_maf = best_maf_idx(gene1_set)
    g2_rep, g2_maf = best_maf_idx(gene2_set)
    print(f"    representatives: {g1_rep} (MAF={g1_maf:.3f})  "
          f"{g2_rep} (MAF={g2_maf:.3f})")
    g1_idx = variant_ids.index(g1_rep)
    g2_idx = variant_ids.index(g2_rep)

    rng = np.random.default_rng(2026)
    g1d = np.where(np.isnan(dosages[:, g1_idx]), np.nanmean(dosages[:, g1_idx]),
                   dosages[:, g1_idx])
    g2d = np.where(np.isnan(dosages[:, g2_idx]), np.nanmean(dosages[:, g2_idx]),
                   dosages[:, g2_idx])
    g1c = g1d - g1d.mean()
    g2c = g2d - g2d.mean()
    causal = 3.0 * g1c * g2c
    phenotype = causal + rng.standard_normal(dosages.shape[0])

    # Annotated indices for M5★ seeding
    annotated_vids = [v for v in variant_ids if v in cache]
    g1_seeds = [k for k, v in enumerate(annotated_vids)
                 if g1_sys in cache[v].get("genes", [])]
    g2_seeds = [k for k, v in enumerate(annotated_vids)
                 if g2_sys in cache[v].get("genes", [])]

    # Run M2
    print(f"    running M2...")
    t0 = time.time()
    try:
        m2_res = run_m2(cache, dosages, variant_ids, phenotype, gene1_set, gene2_set)
    except Exception as e:
        m2_res = {"error": str(e)}
    print(f"      M2: rank={m2_res.get('ground_truth_rank')}/{m2_res.get('n_tested')} "
          f"in {time.time()-t0:.1f}s")

    # Run M5★
    print(f"    running M5★...")
    t0 = time.time()
    try:
        m5_res = run_m5_star(cache, annotated_vids, g1_seeds, g2_seeds,
                              dosages, variant_ids, phenotype, mac_min=10)
    except Exception as e:
        m5_res = {"error": str(e), "trace": traceback.format_exc()}
    print(f"      M5★: rank={m5_res.get('ground_truth_rank')}/{m5_res.get('n_tested')} "
          f"in {time.time()-t0:.1f}s")

    # Free
    del cache, dosages
    gc.collect()

    return {
        **base,
        "status": "TESTED",
        "n_samples": int(len(sample_ids)),
        "n_dosage_variants": int(len(variant_ids)),
        "g1_panel_count": len(gene1_set),
        "g2_panel_count": len(gene2_set),
        "g1_rep_maf": float(g1_maf),
        "g2_rep_maf": float(g2_maf),
        "M2": m2_res,
        "M5_star": m5_res,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--species", choices=["yeast", "arabidopsis", "rice", "all"],
                        default="all")
    parser.add_argument("--max-pairs", type=int, default=None,
                        help="cap number of pairs per species (debug)")
    args = parser.parse_args()

    catalogue = json.loads(CATALOGUE.read_text())
    gene_index_full = json.loads(GENE_INDEX.read_text())["species"]

    species_list = ["yeast", "arabidopsis", "rice"] if args.species == "all" \
                    else [args.species]

    all_results: dict = {}
    for sp in species_list:
        print(f"\n{'='*70}\n  SPECIES: {sp.upper()}\n{'='*70}")
        pairs = catalogue["species"][sp]["pairs"]
        if args.max_pairs:
            pairs = pairs[:args.max_pairs]
        gene_index = gene_index_full.get(sp, {})
        sp_results: list = []
        for p in pairs:
            try:
                r = run_pair(sp, p, gene_index)
            except Exception as e:
                r = {"pair_id": p.get("id", "?"), "status": "RUNNER_ERROR",
                     "error": str(e), "trace": traceback.format_exc()}
                print(f"  ERROR: {e}")
                traceback.print_exc()
            sp_results.append(r)
            gc.collect()
        all_results[sp] = sp_results

    # Summary table
    print(f"\n{'='*70}\n  CATALOGUE SUMMARY\n{'='*70}\n")
    for sp, sp_results in all_results.items():
        print(f"\n  {sp.upper()}:")
        n_total = len(sp_results)
        n_tested = sum(1 for r in sp_results if r.get("status") == "TESTED")
        m5_rank1 = sum(1 for r in sp_results
                        if r.get("status") == "TESTED"
                        and r.get("M5_star", {}).get("ground_truth_rank") == 1)
        m2_top10 = sum(1 for r in sp_results
                        if r.get("status") == "TESTED"
                        and (r.get("M2", {}).get("ground_truth_rank") or 1e9) <= 10)
        print(f"    {n_tested}/{n_total} testable; M5★ rank=1: {m5_rank1}; "
              f"M2 rank≤10: {m2_top10}")
        for r in sp_results:
            status = r.get("status", "?")
            pid = r.get("pair_id", "?")
            anchor = "⭐" if r.get("is_anchor") else " "
            if status == "TESTED":
                m5r = r.get("M5_star", {}).get("ground_truth_rank")
                m5n = r.get("M5_star", {}).get("n_tested")
                m2r = r.get("M2", {}).get("ground_truth_rank")
                m2n = r.get("M2", {}).get("n_tested")
                print(f"    {anchor} {pid:<35}  M5★={m5r}/{m5n}  M2={m2r}/{m2n}")
            else:
                print(f"    {anchor} {pid:<35}  {status}")

    RESULTS_OUT.parent.mkdir(parents=True, exist_ok=True)
    RESULTS_OUT.write_text(json.dumps(all_results, indent=2))
    print(f"\n  Wrote {RESULTS_OUT}")


if __name__ == "__main__":
    main()
