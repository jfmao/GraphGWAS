"""Paper #2 §Y.4 — yeast real-data validation of M2 motif-filtered epistasis.

Tests whether M2 motif-filtering recovers a literature-validated epistatic
pair on real yeast genotypes.

Ground-truth pair: BCY1 (YIL033C) -- TPK1 (YJL164C), the regulatory and
catalytic subunits of yeast Protein Kinase A.  These two genes interact
through cAMP signalling (Phillips 2008 \\cite{phillips2008epistasis};
literature reviewed in Bloom 2015 yeast QTL studies).

Setup:
- Load yeast 1011 Genomes MAF-≥0.05 SNP genotypes from
  tests/data/yeast/1011_snps_maf05.vcf.gz (1,011 samples, ~69K SNPs)
- Restrict to chr9, chr10, chr12 (the chromosomes harbouring BCY1, TPK1,
  HSP104) for tractability
- Match against data/yeast/yeast_graph_cache_v2.json
  (gene/pathway/PPI annotations from SGD + BIOGRID, with prior_score)
- Pick representative BCY1 and TPK1 variants (highest-MAF SNP in each gene)
- Simulate a phenotype y = beta * G_BCY1 * G_TPK1 + Z_nuisance * theta + eps
  where the BCY1 x TPK1 product is the ONLY causal interaction (everything
  else is null nuisance)
- Run motif_filtered_epistasis_from_data with motifs=["same_gene","same_pathway"]
- Report the rank of (BCY1*, TPK1*) in the result list

Expectation: BCY1 and TPK1 share 6 pathways in the graph cache including
generic terms ("chromatin", "cytoplasm", "nucleus", "response to stress").
The same_pathway motif should enumerate the (BCY1*, TPK1*) pair and the
strong simulated interaction should put it at rank 1.

Honest-disclosure mode: if the pair fails to surface at rank 1, the script
reports its actual rank and the top-3 pairs ahead of it for diagnostic.

Usage:
    python tests/validate_m2_yeast.py [--beta 0.5] [--seed 2026]

Output:
    results/paper2_epistasis/yeast_validation.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from cyvcf2 import VCF

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src" / "python"))

from graphgwas.epistasis_v2 import motif_filtered_epistasis_from_data  # noqa: E402

YEAST_VCF = REPO_ROOT / "tests" / "data" / "yeast" / "1011_snps_maf05.vcf.gz"
YEAST_CACHE = REPO_ROOT / "data" / "yeast" / "yeast_graph_cache_v2.json"
RESULTS_DIR = REPO_ROOT / "results" / "paper2_epistasis"

# Literature-validated yeast epistatic pair
GROUND_TRUTH = {
    "gene_1_systematic": "YIL033C",
    "gene_1_symbol": "BCY1",
    "gene_1_chrom": "chromosome9",
    "gene_2_systematic": "YJL164C",
    "gene_2_symbol": "TPK1",
    "gene_2_chrom": "chromosome10",
    "biological_role": (
        "BCY1 is the regulatory subunit of yeast Protein Kinase A; "
        "TPK1 is one of three catalytic subunits.  Their interaction "
        "is the canonical cAMP-PKA signalling switch."
    ),
    "literature": "Phillips 2008 (Nat Rev Genet); Toda et al. 1987 (Cell)",
}

# Limit VCF reading to chromosomes containing the targets to keep runtime
# under a minute on a 1011×69K-SNP panel.
CHROMS_TO_LOAD = ["chromosome9", "chromosome10", "chromosome12"]


# ===========================================================================
# VCF loading
# ===========================================================================

def load_yeast_chroms(vcf_path: Path,
                      chroms: list,
                      maf_min: float = 0.05) -> tuple:
    """Load yeast SNP dosages for the given chromosomes.

    Returns:
        (dosages, variant_ids, sample_ids)
        dosages: shape (n_samples, n_variants), float64, NaN for missing
        variant_ids: list of "chromosome{N}:{pos}:{REF}:{ALT}" strings
        sample_ids: list of length n_samples
    """
    print(f"[1/6] Loading VCF from {vcf_path}")
    print(f"      chromosomes: {chroms}")

    vcf = VCF(str(vcf_path))
    sample_ids = list(vcf.samples)
    n_samples = len(sample_ids)

    rows: list = []
    vids: list = []
    n_skipped_non_biallelic = 0
    n_skipped_low_maf = 0
    n_skipped_other_chrom = 0
    chrom_set = set(chroms)
    per_chr_counts: dict = {c: 0 for c in chroms}

    # Sequential read (no tabix index required) — filter on CHROM inside.
    for record in vcf:
        if record.CHROM not in chrom_set:
            n_skipped_other_chrom += 1
            continue
        if not record.is_snp:
            continue
        if len(record.ALT) != 1:
            n_skipped_non_biallelic += 1
            continue
        g = record.gt_types
        d = np.where(g == 3, np.nan, g).astype(np.float64)
        af = np.nanmean(d) / 2.0
        maf = min(af, 1.0 - af)
        if maf < maf_min:
            n_skipped_low_maf += 1
            continue
        rows.append(d)
        vids.append(f"{record.CHROM}:{record.POS}:{record.REF}:{record.ALT[0]}")
        per_chr_counts[record.CHROM] += 1
    for c in chroms:
        print(f"      {c}: {per_chr_counts[c]:,} biallelic SNPs after MAF≥{maf_min}")
    vcf.close()

    if not rows:
        sys.exit("No variants loaded")
    dosages = np.column_stack(rows)
    print(f"      total: {dosages.shape[0]} samples × {dosages.shape[1]:,} variants")
    print(f"      skipped: {n_skipped_non_biallelic:,} non-biallelic, "
          f"{n_skipped_low_maf:,} low-MAF")
    return dosages, vids, sample_ids


# ===========================================================================
# Ground-truth variant selection
# ===========================================================================

def pick_representative_variants(cache: dict,
                                 variant_id_to_col: dict,
                                 dosages: np.ndarray,
                                 ground_truth: dict) -> dict:
    """For each ground-truth gene, pick the variant with highest MAF.

    Highest MAF maximises power for the simulated interaction.

    Returns:
        dict with gene_1_variant, gene_1_idx, gene_2_variant, gene_2_idx,
        plus shared_pathways for QC.
    """
    print(f"\n[2/6] Picking representative variants for ground-truth pair")

    g1_gene = ground_truth["gene_1_systematic"]
    g2_gene = ground_truth["gene_2_systematic"]

    def find_best(gene: str):
        candidates = [
            (vid, idx) for vid, idx in variant_id_to_col.items()
            if gene in cache.get(vid, {}).get("genes", [])
        ]
        if not candidates:
            return None, None
        # Highest MAF first
        af = np.array([np.nanmean(dosages[:, idx]) / 2.0 for _, idx in candidates])
        maf = np.minimum(af, 1.0 - af)
        best_k = int(np.argmax(maf))
        return candidates[best_k][0], candidates[best_k][1]

    g1_vid, g1_idx = find_best(g1_gene)
    g2_vid, g2_idx = find_best(g2_gene)
    if g1_vid is None or g2_vid is None:
        sys.exit(f"Could not find variants for {g1_gene} or {g2_gene} "
                 f"in cache + dosage data")

    # Shared pathways (sanity QC for same_pathway motif)
    g1_paths = set(cache[g1_vid]["pathways"])
    g2_paths = set(cache[g2_vid]["pathways"])
    shared = sorted(g1_paths & g2_paths)
    print(f"      {ground_truth['gene_1_symbol']} ({g1_gene}): "
          f"{g1_vid}  MAF={min(np.nanmean(dosages[:, g1_idx])/2.0, 1.0 - np.nanmean(dosages[:, g1_idx])/2.0):.3f}")
    print(f"      {ground_truth['gene_2_symbol']} ({g2_gene}): "
          f"{g2_vid}  MAF={min(np.nanmean(dosages[:, g2_idx])/2.0, 1.0 - np.nanmean(dosages[:, g2_idx])/2.0):.3f}")
    print(f"      shared pathways: {len(shared)} ({shared[:3]}{'...' if len(shared)>3 else ''})")
    return {
        "gene_1_variant": g1_vid, "gene_1_idx": g1_idx,
        "gene_2_variant": g2_vid, "gene_2_idx": g2_idx,
        "shared_pathways": shared,
    }


# ===========================================================================
# Phenotype simulation (single causal pair)
# ===========================================================================

def simulate_with_truth(dosages: np.ndarray,
                        truth: dict,
                        beta: float,
                        n_nuisance: int,
                        rng: np.random.Generator) -> np.ndarray:
    """Generate y = beta * g_BCY1 * g_TPK1 + Z_nuisance * theta + eps.

    Args:
        dosages: shape (n_samples, n_variants).
        truth: output of pick_representative_variants().
        beta: causal-interaction effect size.
        n_nuisance: number of random nuisance interaction features.
        rng: numpy Generator.
    """
    print(f"\n[3/6] Simulating phenotype: beta={beta}, n_nuisance={n_nuisance}")
    n_samples, n_variants = dosages.shape

    # Causal interaction
    g1 = dosages[:, truth["gene_1_idx"]]
    g2 = dosages[:, truth["gene_2_idx"]]
    g1_filled = np.where(np.isnan(g1), np.nanmean(g1), g1)
    g2_filled = np.where(np.isnan(g2), np.nanmean(g2), g2)
    g1c = g1_filled - g1_filled.mean()
    g2c = g2_filled - g2_filled.mean()
    causal = beta * g1c * g2c

    # Optional nuisance Z (random non-causal pair products)
    if n_nuisance > 0:
        forbidden = {truth["gene_1_idx"], truth["gene_2_idx"]}
        candidates = np.array([i for i in range(n_variants) if i not in forbidden])
        a_idx = rng.choice(candidates, size=n_nuisance, replace=False)
        b_idx = rng.choice(candidates, size=n_nuisance, replace=False)
        Z = dosages[:, a_idx] * dosages[:, b_idx]
        Z = np.where(np.isnan(Z), np.nanmean(Z, axis=0, keepdims=True), Z)
        Z = Z - Z.mean(axis=0, keepdims=True)
        Z_std = Z.std(axis=0, keepdims=True).clip(min=1e-8)
        Z = Z / Z_std
        # Each theta ~ N(0, 0.01²) → expected Var(Z·θ) = n_nuisance × 1e-4
        theta = rng.standard_normal(n_nuisance) * 0.01
        nuisance = Z @ theta
    else:
        nuisance = np.zeros(n_samples, dtype=np.float64)

    # Noise
    eps = rng.standard_normal(n_samples) * 1.0

    y = causal + nuisance + eps
    var_causal = causal.var()
    var_nuisance = nuisance.var()
    var_total = y.var()
    print(f"      Var decomposition:  causal={var_causal:.3f}  "
          f"nuisance={var_nuisance:.3f}  noise≈1.0  total={var_total:.3f}")
    print(f"      Causal R^2 (interaction):  {var_causal / var_total:.3f}")
    return y


# ===========================================================================
# Main
# ===========================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--beta", type=float, default=3.0,
                        help="Causal interaction effect size (default 3.0). "
                             "BCY1 and TPK1 representatives have MAF≈0.06; "
                             "Var(G1*G2)≈0.014 ⇒ beta=3 ⇒ causal R²≈11%%.")
    parser.add_argument("--n-nuisance", type=int, default=0,
                        help="Number of nuisance interaction features in DGP "
                             "(default 0 = pure null background; only the "
                             "BCY1×TPK1 interaction is non-null).")
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--max-pairs-per-entity", type=int, default=500,
                        help="Per-pathway/gene cap on enumerated pairs")
    args = parser.parse_args()

    print("=== Yeast §Y.4 validation: M2 motif-filtered epistasis ===")
    print(f"VCF        : {YEAST_VCF}")
    print(f"Graph cache: {YEAST_CACHE}")
    print(f"Beta       : {args.beta}")
    print(f"Seed       : {args.seed}")
    print()

    if not YEAST_VCF.exists():
        sys.exit(f"VCF missing: {YEAST_VCF}")
    if not YEAST_CACHE.exists():
        sys.exit(f"Graph cache missing: {YEAST_CACHE}")

    rng = np.random.default_rng(args.seed)

    # 1. Load genotypes
    dosages, variant_ids, sample_ids = load_yeast_chroms(
        YEAST_VCF, CHROMS_TO_LOAD, maf_min=0.05,
    )

    # Match against cache
    print(f"\n      Loading graph cache...")
    cache = json.loads(YEAST_CACHE.read_text())
    print(f"      cache has {len(cache):,} variants total")

    variant_id_to_col = {v: i for i, v in enumerate(variant_ids)}
    n_in_cache = sum(1 for v in variant_ids if v in cache)
    print(f"      {n_in_cache:,}/{len(variant_ids):,} loaded variants are in the cache")
    if n_in_cache < 100:
        sys.exit("Too few cache hits — cache/VCF format mismatch?")

    # 2. Pick representative variants for ground-truth pair
    truth = pick_representative_variants(cache, variant_id_to_col, dosages, GROUND_TRUTH)

    # 3. Simulate phenotype
    phenotype = simulate_with_truth(dosages, truth, args.beta, args.n_nuisance, rng)

    # 4. Run M2
    print(f"\n[4/6] Running M2 motif-filtered epistasis (no-Neo4j path)")
    results = motif_filtered_epistasis_from_data(
        dosages=dosages,
        variant_ids=variant_ids,
        phenotype=phenotype,
        graph_cache=cache,
        motifs=["same_gene", "same_pathway"],
        mac_min=10,
        correction="BH",
        max_pairs_per_entity=args.max_pairs_per_entity,
        max_pairs_total=200_000,
        verbose=True,
    )

    # 5. Find rank of ground-truth pair
    print(f"\n[5/6] Ranking ground-truth pair {GROUND_TRUTH['gene_1_symbol']} × "
          f"{GROUND_TRUTH['gene_2_symbol']}")

    g1v = truth["gene_1_variant"]
    g2v = truth["gene_2_variant"]
    truth_rank = None
    truth_record = None
    for k, r in enumerate(results, start=1):
        if {r.variant_1, r.variant_2} == {g1v, g2v}:
            truth_rank = k
            truth_record = r
            break

    # 6. Report + write JSON
    print(f"\n[6/6] Result")
    out: dict = {
        "ground_truth": GROUND_TRUTH,
        "params": {
            "beta": args.beta,
            "n_nuisance": args.n_nuisance,
            "seed": args.seed,
            "max_pairs_per_entity": args.max_pairs_per_entity,
            "chroms_loaded": CHROMS_TO_LOAD,
        },
        "n_samples": int(dosages.shape[0]),
        "n_variants_loaded": int(dosages.shape[1]),
        "n_variants_in_cache": int(n_in_cache),
        "shared_pathways_g1_g2": truth["shared_pathways"],
        "n_pairs_tested": int(len(results)),
        "n_pairs_significant_q05": sum(
            1 for r in results
            if r.p_corrected is not None and r.p_corrected < 0.05
        ),
    }

    if truth_rank is not None:
        print(f"      ✓ ground-truth pair found at RANK {truth_rank} of "
              f"{len(results):,} tested pairs")
        print(f"        p_interaction: {truth_record.p_interaction:.2e}")
        print(f"        BH-FDR q-val : {truth_record.p_corrected:.2e}")
        print(f"        β_interaction: {truth_record.beta_interaction:+.4f}")
        print(f"        motif        : {truth_record.motif}")
        print(f"        shared_entity: {truth_record.shared_entity}")

        out["ground_truth_recovered"] = True
        out["ground_truth_rank"] = int(truth_rank)
        out["ground_truth_p_interaction"] = float(truth_record.p_interaction)
        out["ground_truth_q_value"] = float(truth_record.p_corrected)
        out["ground_truth_beta"] = float(truth_record.beta_interaction)
        out["ground_truth_motif"] = truth_record.motif
        out["ground_truth_shared_entity"] = truth_record.shared_entity
    else:
        print(f"      ✗ ground-truth pair NOT found in M2's tested pool")
        print(f"        Possible reasons:")
        print(f"        - The motif enumeration's per-entity cap excluded it "
              f"(try --max-pairs-per-entity {args.max_pairs_per_entity * 4})")
        print(f"        - MAC filter excluded it (BCY1 or TPK1 representative "
              f"is too rare)")
        print(f"        - cache lacks a shared pathway (we found "
              f"{len(truth['shared_pathways'])})")
        out["ground_truth_recovered"] = False
        out["ground_truth_rank"] = None

    # Top-5 pairs in any case (diagnostic / paper figure source)
    out["top_5_pairs"] = [
        {
            "rank": k,
            "variant_1": r.variant_1,
            "variant_2": r.variant_2,
            "motif": r.motif,
            "shared_entity": r.shared_entity,
            "p_interaction": float(r.p_interaction),
            "q_value": float(r.p_corrected) if r.p_corrected else None,
            "beta_interaction": float(r.beta_interaction),
        }
        for k, r in enumerate(results[:5], start=1)
    ]

    print()
    print(f"      Top-5 pairs:")
    for entry in out["top_5_pairs"]:
        marker = " ★" if {entry["variant_1"], entry["variant_2"]} == {g1v, g2v} else ""
        print(f"        #{entry['rank']:<2}  q={entry['q_value']:.2e}  "
              f"{entry['variant_1']} × {entry['variant_2']}  "
              f"[{entry['motif']}]{marker}")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out_path = RESULTS_DIR / "yeast_validation.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\n      Wrote {out_path}")


if __name__ == "__main__":
    main()
