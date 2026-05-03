"""Paper #2 §Y.6 — cross-species replication of M2 epistasis recovery.

Tests whether the M2 motif-filtered epistasis result on yeast (§Y.4 —
recovers BCY1 × TPK1 at q < 10⁻⁴) generalises to other species:

  - Arabidopsis (1001 Genomes): FT (AT1G65480, chr1) × FLC (AT5G10140, chr5)
    canonical flowering-time epistatic pair.
  - Rice (3,000 Rice Genomes): Hd1 (LOC_Os06g16370, chr6) × Hd3a
    (LOC_Os06g06320, chr6) canonical heading-date epistatic pair.

Same recipe as ``validate_m2_yeast.py``:
  1. Load genotypes restricted to the chromosomes of the ground-truth pair.
  2. Match against the species' graph cache (gene/pathway/PPI annotations).
  3. Pick highest-MAF representative variants from each gene.
  4. Simulate ``y = beta * G_X * G_Y + eps`` (Scenario A — positive interaction).
  5. Run all 5 methods (M1, M2, M3, M4, M5) on the same data.
  6. Use gene-level recovery matching: any (gene_X-var × gene_Y-var) pair
     in a method's result list counts as recovery.

Output:
  results/paper2_epistasis/cross_species_validation.json
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from cyvcf2 import VCF

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src" / "python"))
sys.path.insert(0, str(REPO_ROOT / "tests"))

from graphgwas.epistasis_v2 import (  # noqa: E402
    InteractionResult,
    _test_interaction,
    dark_matter_epistasis_from_data,
    differential_subgraph_from_data,
    ld_pruned_cooccurrence_from_data,
    motif_filtered_epistasis_from_data,
)

# Reuse helpers from the yeast script
from validate_m2_yeast import (  # noqa: E402
    _bipartite_from_dict,
    _print_comparison_table,
    _run_m5_then_test,
    _subsample_with_truth_priority,
    _summarise_method,
)

RESULTS_DIR = REPO_ROOT / "results" / "paper2_epistasis"


# ===========================================================================
# Per-species configuration
# ===========================================================================

SPECIES_CONFIG = {
    "arabidopsis": {
        "vcf_path": REPO_ROOT / "tests" / "data" / "arabidopsis" /
                    "1001genomes_snp-short-indel_only_ACGTN.vcf.gz",
        # Cache: one JSON per chromosome
        "cache_paths": {
            "1": REPO_ROOT / "data" / "arabidopsis" /
                 "arabidopsis_graph_cache_v3_chr1.json",
            "5": REPO_ROOT / "data" / "arabidopsis" /
                 "arabidopsis_graph_cache_v3_chr5.json",
        },
        # Variant ID format in cache + genotype: "1:24331373:A:C"
        "vid_format": "{chrom}:{pos}:{ref}:{alt}",
        "ground_truth": {
            "gene_1_systematic": "AT1G65480",
            "gene_1_symbol": "FT",
            "gene_1_chrom": "1",
            "gene_1_window": (24_300_000, 24_400_000),  # FT is around 24.34Mb
            "gene_2_systematic": "AT5G10140",
            "gene_2_symbol": "FLC",
            "gene_2_chrom": "5",
            "gene_2_window": (3_150_000, 3_250_000),    # FLC ~3.18Mb
            "biological_role": (
                "FT and FLC are the canonical flowering-time epistatic pair "
                "in Arabidopsis. FLC represses FT directly via chromatin "
                "binding; vernalisation suppresses FLC, releasing FT to "
                "promote flowering. FT × FLC interactions appear in 100+ "
                "Arabidopsis QTL studies."
            ),
            "literature": "Michaels & Amasino 1999 (Plant Cell), "
                          "Searle et al. 2006 (Genes Dev)",
        },
    },
    "rice": {
        # Rice uses PLINK 2 pgen format — loaded via pgenlib
        "pgen_prefix": REPO_ROOT / "data" / "rice_3k" / "rice_3k",
        "cache_paths": {
            "Chr6": REPO_ROOT / "data" / "rice_3k" / "annotations" /
                    "rice_graph_cache_v2_Chr6.json",
        },
        # Variant ID format: "Chr6:126124:A:C"
        "vid_format": "Chr{chrom_num}:{pos}:{ref}:{alt}",
        "ground_truth": {
            "gene_1_systematic": "LOC_Os06g16370",
            "gene_1_symbol": "Hd1",
            "gene_1_chrom": "6",
            "gene_1_window": (9_300_000, 9_400_000),  # Hd1 ~9.34Mb
            "gene_2_systematic": "LOC_Os06g06320",
            "gene_2_symbol": "Hd3a",
            "gene_2_chrom": "6",
            "gene_2_window": (2_800_000, 2_900_000),  # Hd3a ~2.83Mb
            "biological_role": (
                "Hd1 and Hd3a are the canonical heading-date epistatic pair "
                "in rice. Hd3a is the rice ortholog of FT; Hd1 is a "
                "photoperiod-sensitive transcription factor that promotes "
                "Hd3a under short days but represses it under long days. "
                "Their epistasis controls flowering-time variation across "
                "the 3K Rice Genomes."
            ),
            "literature": "Yano et al. 2000 (Plant Cell); "
                          "Kojima et al. 2002 (Plant Cell Physiol); "
                          "Tamaki et al. 2007 (Science)",
        },
    },
}


# ===========================================================================
# Genotype loaders (per species)
# ===========================================================================

def load_arabidopsis(config: dict) -> tuple[np.ndarray, list, list]:
    """Load Arabidopsis genotypes around FT (chr1) + FLC (chr5)."""
    vcf_path = config["vcf_path"]
    if not vcf_path.exists():
        sys.exit(f"VCF missing: {vcf_path}")

    gt = config["ground_truth"]
    windows = [
        (gt["gene_1_chrom"], *gt["gene_1_window"]),  # FT region
        (gt["gene_2_chrom"], *gt["gene_2_window"]),  # FLC region
    ]
    # Add a third "control" window for nuisance variation — picks chr3 1.0-1.2Mb
    # (was 1-2Mb but caused OOM in M4 with 24k variants × 1135 samples)
    windows.append(("3", 1_000_000, 1_200_000))

    print(f"[1/6] Loading Arabidopsis VCF {vcf_path.name}")
    vcf = VCF(str(vcf_path))
    sample_ids = list(vcf.samples)
    rows: list = []
    vids: list = []
    for chrom, start, end in windows:
        n_chr = 0
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
            n_chr += 1
        print(f"      chr{chrom}:{start}-{end} → {n_chr:,} bi-allelic SNPs (MAF≥0.05)")
    vcf.close()
    if not rows:
        sys.exit("No Arabidopsis variants loaded")
    dosages = np.column_stack(rows)
    print(f"      total: {dosages.shape[0]} samples × {dosages.shape[1]:,} variants")
    return dosages, vids, sample_ids


def load_rice(config: dict) -> tuple[np.ndarray, list, list]:
    """Load rice chr6 genotypes via pgenlib."""
    import pgenlib

    pgen_prefix = config["pgen_prefix"]
    pgen_path = str(pgen_prefix) + ".pgen"
    pvar_path = str(pgen_prefix) + ".pvar"
    psam_path = str(pgen_prefix) + ".psam"
    if not Path(pgen_path).exists():
        sys.exit(f"pgen missing: {pgen_path}")

    print(f"[1/6] Loading rice pgen {Path(pgen_path).name}")
    # 1. Read .psam to get sample IDs
    sample_ids: list = []
    with open(psam_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            sample_ids.append(parts[0])  # IID is column 0 (or 1 if FID present)
    n_samples = len(sample_ids)
    print(f"      n_samples = {n_samples}")

    # 2. Read .pvar to find chr6 variants in narrow windows around the two
    # target genes plus one nuisance window.  The full chr6:1-12 Mb window
    # yields 183K MAF-filtered variants, which causes M1's per-variant
    # dosage_list copy to allocate 4.4 GB and accumulates to OOM during M3.
    # Mirror the Arabidopsis loader: ±50 kb around each gene + a 200 kb
    # control region — produces ~5-15 K variants, comparable to Arabidopsis.
    gt = config["ground_truth"]
    g1_start, g1_end = gt["gene_1_window"]
    g2_start, g2_end = gt["gene_2_window"]
    # Pad ±50 kb beyond the listed window
    PAD = 50_000
    rice_windows = [
        (gt["gene_1_chrom"], g1_start - PAD, g1_end + PAD),
        (gt["gene_2_chrom"], g2_start - PAD, g2_end + PAD),
        # Control region: chr6:5.5-5.7 Mb (between Hd3a and Hd1, neither gene)
        ("6", 5_500_000, 5_700_000),
    ]
    target_chroms = {"Chr6", "6", str(gt["gene_1_chrom"]), str(gt["gene_2_chrom"])}
    print(f"      scanning .pvar for {len(rice_windows)} chr6 windows: "
          f"{[(w[0], f'{w[1]/1e6:.2f}-{w[2]/1e6:.2f}Mb') for w in rice_windows]}")
    chr6_idx: list = []
    chr6_meta: list = []  # (pos, ref, alt)
    with open(pvar_path) as f:
        line_idx = 0
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                line_idx += 1
                continue
            chrom, pos, _id, ref, alt = parts[0], int(parts[1]), parts[2], parts[3], parts[4]
            if chrom in target_chroms:
                in_any_window = any(
                    s <= pos <= e for (_c, s, e) in rice_windows
                )
                if in_any_window:
                    chr6_idx.append(line_idx)
                    chr6_meta.append((chrom, pos, ref, alt))
            line_idx += 1
    print(f"      chr6 variants in target windows: {len(chr6_idx):,}")

    # 3. Read genotypes for chr6 from .pgen
    n_var = len(chr6_idx)
    if n_var == 0:
        sys.exit("No chr6 variants in rice pvar")

    pgen = pgenlib.PgenReader(pgen_path.encode("utf-8"))
    dosages = np.empty((n_var, n_samples), dtype=np.float32)
    pgen.read_dosages_list(np.asarray(chr6_idx, dtype=np.uint32), dosages)
    pgen.close()
    dosages = dosages.T  # → (n_samples, n_var)
    # pgenlib uses -9 for missing → convert to NaN
    dosages = np.where(dosages < 0, np.nan, dosages).astype(np.float64)

    # 4. MAF filter
    af = np.nanmean(dosages, axis=0) / 2.0
    maf = np.minimum(af, 1.0 - af)
    keep = maf >= 0.05
    dosages = dosages[:, keep]
    chr6_meta_kept = [chr6_meta[k] for k in range(len(chr6_meta)) if keep[k]]
    # Normalize chromosome prefix to match the cache's "Chr6:..." vid format
    # (pvar may use "6" while the cache always uses "Chr6").
    vids = [f"Chr6:{m[1]}:{m[2]}:{m[3]}" for m in chr6_meta_kept]
    print(f"      after MAF≥0.05: {dosages.shape[1]:,} variants")
    return dosages, vids, sample_ids


# ===========================================================================
# Cache-aggregator (single dict from per-chromosome JSONs)
# ===========================================================================

def load_cache_for_species(config: dict) -> dict:
    cache: dict = {}
    for chrom, path in config["cache_paths"].items():
        if not path.exists():
            print(f"      ⚠ cache missing for chr{chrom}: {path}")
            continue
        sub = json.loads(path.read_text())
        cache.update(sub)
        print(f"      chr{chrom} cache: {len(sub):,} variants")
    return cache


# ===========================================================================
# Pick representative variants
# ===========================================================================

def pick_representatives_cross_species(cache: dict,
                                        variant_id_to_col: dict,
                                        dosages: np.ndarray,
                                        gt: dict) -> dict:
    print(f"\n[2/6] Picking representative variants for "
          f"{gt['gene_1_symbol']} × {gt['gene_2_symbol']}")

    def find_best(gene_systematic: str):
        candidates = [
            (vid, idx) for vid, idx in variant_id_to_col.items()
            if gene_systematic in cache.get(vid, {}).get("genes", [])
        ]
        if not candidates:
            return None, None
        af = np.array([np.nanmean(dosages[:, idx]) / 2.0 for _, idx in candidates])
        maf = np.minimum(af, 1.0 - af)
        return candidates[int(np.argmax(maf))][0], candidates[int(np.argmax(maf))][1]

    g1_vid, g1_idx = find_best(gt["gene_1_systematic"])
    g2_vid, g2_idx = find_best(gt["gene_2_systematic"])
    if g1_vid is None or g2_vid is None:
        sys.exit(f"Could not find representative variants for one of the "
                 f"target genes (g1={gt['gene_1_systematic']} / "
                 f"g2={gt['gene_2_systematic']})")

    g1_paths = set(cache[g1_vid].get("pathways", []))
    g2_paths = set(cache[g2_vid].get("pathways", []))
    g1_ppi = set(cache[g1_vid].get("ppi", []))
    g2_ppi = set(cache[g2_vid].get("ppi", []))
    shared_paths = sorted(g1_paths & g2_paths)
    shared_ppi = sorted(g1_ppi & g2_ppi)
    g1_in_g2_ppi = gt["gene_1_systematic"] in g2_ppi
    g2_in_g1_ppi = gt["gene_2_systematic"] in g1_ppi

    def maf_of(idx):
        af = np.nanmean(dosages[:, idx]) / 2.0
        return min(af, 1.0 - af)

    print(f"      {gt['gene_1_symbol']} ({gt['gene_1_systematic']}): "
          f"{g1_vid}  MAF={maf_of(g1_idx):.3f}")
    print(f"      {gt['gene_2_symbol']} ({gt['gene_2_systematic']}): "
          f"{g2_vid}  MAF={maf_of(g2_idx):.3f}")
    print(f"      shared pathways: {shared_paths}")
    print(f"      shared PPI partners: {shared_ppi}")
    print(f"      direct PPI link "
          f"({gt['gene_1_symbol']} ↔ {gt['gene_2_symbol']}): "
          f"{g1_in_g2_ppi or g2_in_g1_ppi}")

    return {
        "gene_1_variant": g1_vid, "gene_1_idx": g1_idx,
        "gene_2_variant": g2_vid, "gene_2_idx": g2_idx,
        "shared_pathways": shared_paths,
        "shared_ppi": shared_ppi,
        "direct_ppi_link": bool(g1_in_g2_ppi or g2_in_g1_ppi),
    }


# ===========================================================================
# Phenotype simulator (Scenario A)
# ===========================================================================

def simulate_scenario_a(dosages: np.ndarray,
                         truth: dict,
                         beta: float,
                         rng: np.random.Generator) -> np.ndarray:
    print(f"\n[3/6] Simulating Scenario A (positive interaction): beta={beta}")
    n_samples = dosages.shape[0]
    g1 = dosages[:, truth["gene_1_idx"]]
    g2 = dosages[:, truth["gene_2_idx"]]
    g1f = np.where(np.isnan(g1), np.nanmean(g1), g1)
    g2f = np.where(np.isnan(g2), np.nanmean(g2), g2)
    g1c = g1f - g1f.mean()
    g2c = g2f - g2f.mean()
    causal = beta * g1c * g2c
    eps = rng.standard_normal(n_samples)
    y = causal + eps
    var_causal = causal.var()
    var_total = y.var()
    print(f"      Var decomposition:  causal={var_causal:.3f}  noise≈1.0  "
          f"total={var_total:.3f}")
    print(f"      Causal R^2 (interaction):  {var_causal / var_total:.3f}")
    return y


# ===========================================================================
# Run all 5 methods on one species
# ===========================================================================

def run_species(species: str, beta: float, seed: int,
                 motifs_list: list[str], inject_canonical_ppi: bool,
                 skip_heavy: bool = False) -> dict:
    print(f"\n{'=' * 70}\n  SPECIES: {species.upper()}\n{'=' * 70}")
    config = SPECIES_CONFIG[species]
    rng = np.random.default_rng(seed)

    # 1. Load genotypes
    if species == "arabidopsis":
        dosages, variant_ids, sample_ids = load_arabidopsis(config)
    elif species == "rice":
        dosages, variant_ids, sample_ids = load_rice(config)
    else:
        sys.exit(f"unknown species: {species}")

    # 2. Cache
    print(f"      Loading {species} cache...")
    cache = load_cache_for_species(config)
    n_in_cache = sum(1 for v in variant_ids if v in cache)
    print(f"      {n_in_cache:,}/{len(variant_ids):,} variants in cache")

    # 2b. Optional canonical-PPI injection
    inject_diag = None
    if inject_canonical_ppi:
        from graphgwas.canonical_ppi import (
            CANONICAL_EDGES_BY_SPECIES,
            inject_canonical_ppi_edges,
        )
        edges = CANONICAL_EDGES_BY_SPECIES.get(species, [])
        cache, inject_diag = inject_canonical_ppi_edges(
            cache, edges=edges, inplace=True, verbose=False,
        )
        print(f"      ✓ injected canonical PPI: "
              f"{inject_diag['n_edges_applied']}/{inject_diag['n_edges_supplied']} "
              f"edges, {inject_diag['n_variants_modified']} variants modified")
        if inject_diag["edges_with_zero_variants"]:
            print(f"      ⚠ edges dropped (no variant in cache): "
                  f"{inject_diag['edges_with_zero_variants']}")

    variant_id_to_col = {v: i for i, v in enumerate(variant_ids)}

    # 3. Truth + simulate
    truth = pick_representatives_cross_species(
        cache, variant_id_to_col, dosages, config["ground_truth"],
    )
    phenotype = simulate_scenario_a(dosages, truth, beta, rng)

    # 4. Recovery sets (gene-level)
    gt = config["ground_truth"]
    g1v = truth["gene_1_variant"]
    g2v = truth["gene_2_variant"]
    gene1_variants = {
        vid for vid in variant_ids
        if gt["gene_1_systematic"] in cache.get(vid, {}).get("genes", [])
    }
    gene2_variants = {
        vid for vid in variant_ids
        if gt["gene_2_systematic"] in cache.get(vid, {}).get("genes", [])
    }
    print(f"\n      Recovery sets: {gt['gene_1_symbol']} = "
          f"{len(gene1_variants)} variants, {gt['gene_2_symbol']} = "
          f"{len(gene2_variants)} variants (cache-annotated + dosage-loaded)")
    if not gene1_variants or not gene2_variants:
        sys.exit(f"Empty recovery set — increase the load-window for one of "
                 f"the genes")

    # Case/control masks (25/75 percentile binarisation)
    q75 = np.percentile(phenotype, 75)
    q25 = np.percentile(phenotype, 25)
    case_mask = phenotype >= q75
    ctrl_mask = phenotype <= q25

    methods_out: list = []

    # ----- M1 -----
    print(f"\n[M1/{species}] LD-pruned + interaction test")
    t0 = time.time()
    af_vec = np.nanmean(dosages, axis=0) / 2.0
    m1_variants = [
        {"variantId": variant_ids[k],
         # Position parsing: chr-prefix-stripped second field
         "pos": int(variant_ids[k].split(":")[1]),
         "af_total": float(af_vec[k])}
        for k in range(len(variant_ids))
    ]
    m1_dosage_list = [dosages[:, k] for k in range(dosages.shape[1])]
    m1_results = ld_pruned_cooccurrence_from_data(
        variants=m1_variants, dosage_list=m1_dosage_list, phenotype=phenotype,
        r2_prune=0.5, min_cocarriers=1, min_distance_bp=10_000,
        max_variants=2000, verbose=True,
    )
    m1_runtime = time.time() - t0
    methods_out.append(_summarise_method(
        "M1", "LD-pruned co-occurrence (Bonferroni)",
        m1_results, m1_runtime, g1v, g2v,
        gene1_variants=gene1_variants, gene2_variants=gene2_variants,
    ))

    # ----- M2 -----
    print(f"\n[M2/{species}] Motif-filtered ({'+'.join(motifs_list)})")
    t0 = time.time()
    m2_results = motif_filtered_epistasis_from_data(
        dosages=dosages, variant_ids=variant_ids, phenotype=phenotype,
        graph_cache=cache, motifs=motifs_list,
        mac_min=10, correction="BH",
        max_pairs_per_entity=500, max_pairs_total=200_000, verbose=True,
    )
    m2_runtime = time.time() - t0
    methods_out.append(_summarise_method(
        "M2", "Motif-filtered (BH-FDR)",
        m2_results, m2_runtime, g1v, g2v,
        gene1_variants=gene1_variants, gene2_variants=gene2_variants,
    ))

    # ----- M3 -----
    print(f"\n[M3/{species}] Differential subgraph (case/control)")
    t0 = time.time()
    sub_idx = _subsample_with_truth_priority(
        variant_ids, gene1_variants, gene2_variants, n_target=1500, rng=rng,
    )
    sub_dosages = dosages[:, sub_idx]
    sub_vids = [variant_ids[k] for k in sub_idx]
    print(f"      subsampled to {len(sub_idx)} variants (truth genes forced in)")
    m3_results = differential_subgraph_from_data(
        dosages=sub_dosages, variant_ids=sub_vids,
        case_mask=case_mask, control_mask=ctrl_mask,
        r2_prune=0.5, min_distance_bp=10_000,
        min_cocarriers=2, mac_min=10, max_variants=500, verbose=True,
    )
    m3_runtime = time.time() - t0
    methods_out.append(_summarise_method(
        "M3", "Differential subgraph (Fisher's exact)",
        m3_results, m3_runtime, g1v, g2v,
        gene1_variants=gene1_variants, gene2_variants=gene2_variants,
    ))

    # ----- M4 -----
    if skip_heavy:
        print(f"\n[M4/{species}] Skipped (--skip-heavy)")
    else:
        print(f"\n[M4/{species}] Dark matter / synthetic incompatibility")
        t0 = time.time()
        m4_results = dark_matter_epistasis_from_data(
            dosages=sub_dosages, variant_ids=sub_vids,
            case_mask=case_mask, control_mask=ctrl_mask,
            maf_min=0.05, depletion_threshold=0.7, min_expected=0.5,
            mac_min=10, max_variants=500, verbose=True,
        )
        m4_runtime = time.time() - t0
        methods_out.append(_summarise_method(
            "M4", "Dark matter (Poisson)",
            m4_results, m4_runtime, g1v, g2v,
            gene1_variants=gene1_variants, gene2_variants=gene2_variants,
        ))

    # ----- M5 -----
    print(f"\n[M5/{species}] RWR + interaction test")
    t0 = time.time()
    truth_for_m5 = {
        "gene_1_systematic": gt["gene_1_systematic"],
        "gene_2_systematic": gt["gene_2_systematic"],
    }
    # _run_m5_then_test from the yeast script hardcodes BCY1/TPK1 systematic IDs;
    # we mirror its logic inline here for cross-species generality.
    m5_results = _run_m5_then_test_general(
        dosages=dosages, variant_ids=variant_ids, phenotype=phenotype,
        cache=cache, gt=gt,
        alpha=0.15, top_k=2000, mac_min=10, rng=rng,
    )
    m5_runtime = time.time() - t0
    methods_out.append(_summarise_method(
        "M5", "RWR + interaction test (BH-FDR)",
        m5_results, m5_runtime, g1v, g2v,
        gene1_variants=gene1_variants, gene2_variants=gene2_variants,
    ))

    # Per-species comparison
    print(f"\n[{species} summary]")
    _print_comparison_table(methods_out)

    return {
        "species": species,
        "ground_truth": config["ground_truth"],
        "n_samples": int(dosages.shape[0]),
        "n_variants_loaded": int(dosages.shape[1]),
        "n_variants_in_cache": int(n_in_cache),
        "shared_pathways_g1_g2": truth["shared_pathways"],
        "shared_ppi_g1_g2": truth["shared_ppi"],
        "direct_ppi_link": truth["direct_ppi_link"],
        "methods": methods_out,
    }


# ===========================================================================
# M5 wrapper, generalised across species (no BCY1/TPK1 hardcoding)
# ===========================================================================

def _run_m5_then_test_general(dosages: np.ndarray,
                                variant_ids: list,
                                phenotype: np.ndarray,
                                cache: dict,
                                gt: dict,
                                alpha: float,
                                top_k: int,
                                mac_min: int,
                                rng: np.random.Generator) -> list:
    """M5 RWR + interaction test, parameterised by species ground-truth."""
    annotated_vids = [v for v in variant_ids if v in cache]
    print(f"      bipartite graph: {len(annotated_vids):,} annotated variants")
    A, cache_var_ids, _ = _bipartite_from_dict(cache, annotated_vids)

    cache_idx_to_col = {
        cache_var_ids[k]: variant_ids.index(cache_var_ids[k])
        for k in range(len(cache_var_ids))
    }
    g1_seeds = [k for k, v in enumerate(cache_var_ids)
                if gt["gene_1_systematic"] in cache[v].get("genes", [])]
    g2_seeds = [k for k, v in enumerate(cache_var_ids)
                if gt["gene_2_systematic"] in cache[v].get("genes", [])]
    other_pool = [k for k in range(len(cache_var_ids))
                  if k not in set(g1_seeds) | set(g2_seeds)]
    extras = list(other_pool)  # all annotated variants
    seed_indices = sorted(set(g1_seeds + g2_seeds + extras))
    print(f"      seeds: {len(g1_seeds)} {gt['gene_1_symbol']} + "
          f"{len(g2_seeds)} {gt['gene_2_symbol']} + {len(extras)} = "
          f"{len(seed_indices)}")

    # Batched RWR
    print(f"      computing batched RWR...")
    n_var, n_gene = A.shape
    var_deg = A.sum(axis=1, keepdims=True).clip(min=1e-12)
    gene_deg = A.sum(axis=0, keepdims=True).clip(min=1e-12)
    P_vg = A / var_deg
    P_gv = (A / gene_deg).T
    M_vv = P_vg @ P_gv
    n_seed = len(seed_indices)
    E = np.zeros((n_var, n_seed), dtype=np.float64)
    for k, si in enumerate(seed_indices):
        E[si, k] = 1.0
    P = E.copy()
    for _ in range(30):
        P = (1.0 - alpha) * (M_vv @ P) + alpha * E
    P = P / np.clip(P.sum(axis=0, keepdims=True), 1e-12, None)

    # Pair scoring
    seed_arr = np.asarray(seed_indices)
    P_seed = P[seed_arr, :]
    mutual = P_seed * P_seed.T
    iu = np.triu_indices(n_seed, k=1)
    pair_idx_pairs = list(zip(iu[0], iu[1]))
    scores = mutual[iu]
    order = np.argsort(scores)[::-1]
    pair_scores: list = []
    for o in order:
        s = float(scores[o])
        if s <= 0:
            break
        ki, kj = pair_idx_pairs[o]
        pair_scores.append((seed_indices[ki], seed_indices[kj], s))
        if len(pair_scores) >= top_k:
            break
    print(f"      top {len(pair_scores):,} pairs by mutual-RWR (out of {(scores>0).sum():,} non-zero)")

    # Test interactions on top-K
    af_vec = np.nanmean(dosages, axis=0) / 2.0
    maf_vec = np.minimum(af_vec, 1.0 - af_vec)
    results: list = []
    for i_cache, j_cache, rwr_score in pair_scores:
        vid_i = cache_var_ids[i_cache]
        vid_j = cache_var_ids[j_cache]
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
        results.append(InteractionResult(
            variant_1=vid_i, variant_2=vid_j,
            motif="rwr", shared_entity=f"rwr_score={rwr_score:.4g}",
            beta_marginal_1=r["beta_1"], beta_marginal_2=r["beta_2"],
            beta_interaction=r["beta_interaction"],
            se_interaction=r["se_interaction"],
            p_interaction=r["p_interaction"],
            p_corrected=None, n_samples=int(r["n"]),
            maf_1=float(maf_vec[i]), maf_2=float(maf_vec[j]),
        ))
    if not results:
        return []
    pvals = np.array([r.p_interaction for r in results])
    n_tests = len(pvals)
    order = np.argsort(pvals)
    ranks = np.empty(n_tests, dtype=int)
    ranks[order] = np.arange(1, n_tests + 1)
    corrected = np.minimum(1.0, pvals * n_tests / ranks)
    for k in range(n_tests - 2, -1, -1):
        corrected[order[k]] = min(corrected[order[k]], corrected[order[k + 1]])
    for k, r in enumerate(results):
        r.p_corrected = float(corrected[k])
    results.sort(key=lambda r: r.p_interaction)
    return results


# ===========================================================================
# Main
# ===========================================================================

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--species", choices=["arabidopsis", "rice", "all"],
                   default="all", help="species to run")
    p.add_argument("--beta", type=float, default=3.0)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--motifs", default="same_gene,same_pathway,protein_interaction",
                   help="Comma-separated M2 motifs (default all 3)")
    p.add_argument("--no-inject-ppi", action="store_true",
                   help="Skip canonical-PPI injection")
    p.add_argument("--skip-heavy", action="store_true",
                   help="Skip M4 (avoid OOM at 3K samples × 500 vars).  "
                        "M2 + M3 still run; useful for rice 3K RG.")
    args = p.parse_args()
    motifs_list = [m.strip() for m in args.motifs.split(",") if m.strip()]

    if args.species == "all":
        species_list = ["arabidopsis", "rice"]
    else:
        species_list = [args.species]

    all_results: dict = {}
    for sp in species_list:
        all_results[sp] = run_species(
            sp, args.beta, args.seed,
            motifs_list=motifs_list,
            inject_canonical_ppi=not args.no_inject_ppi,
            skip_heavy=args.skip_heavy,
        )

    # Final cross-species grid
    print(f"\n{'=' * 70}\n  CROSS-SPECIES SUMMARY\n{'=' * 70}\n")
    methods = ["M1", "M2", "M3", "M4", "M5"]
    print(f"  {'Species':<14}  " + "  ".join(f"{m:>10}" for m in methods))
    print(f"  {'-' * 14}  " + "  ".join("-" * 10 for _ in methods))
    for sp, sp_data in all_results.items():
        m_by_name = {m["method"]: m for m in sp_data["methods"]}
        cells = []
        for m in methods:
            entry = m_by_name.get(m)
            if entry is None:
                cells.append("—"); continue
            r = entry["ground_truth_rank"]
            cells.append(f"{r:,}" if r is not None else "NF")
        print(f"  {sp:<14}  " + "  ".join(f"{c:>10}" for c in cells))
    print()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out = RESULTS_DIR / "cross_species_validation.json"
    out.write_text(json.dumps(all_results, indent=2))
    print(f"  Wrote {out}")


if __name__ == "__main__":
    main()
