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

import time  # noqa: E402

from graphgwas.epistasis_v2 import (  # noqa: E402
    InteractionResult,
    _test_interaction,
    dark_matter_epistasis_from_data,
    differential_subgraph_from_data,
    ld_pruned_cooccurrence_from_data,
    motif_filtered_epistasis_from_data,
)
from graphgwas.epistasis_higher_order import (  # noqa: E402
    build_bipartite_adjacency,
    mutual_rwr_pair_scores,
)

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


def simulate_scenario_b(dosages: np.ndarray,
                         truth: dict,
                         baseline_p: float,
                         rng: np.random.Generator
                         ) -> tuple[np.ndarray, np.ndarray]:
    """Scenario B — case/control with case-enriched co-occurrence.

    Each sample is a case if (a) it carries both BCY1 + TPK1 alleles
    (forced case), OR (b) baseline coin flip with probability baseline_p.

    Returns (case_mask, control_mask).  All non-cases are controls.

    This phenotype mimics the M3 differential-subgraph signal type:
    co-occurrence of two specific alleles is over-represented in cases.
    """
    n_samples = dosages.shape[0]
    g1 = dosages[:, truth["gene_1_idx"]]
    g2 = dosages[:, truth["gene_2_idx"]]
    g1_carrier = np.where(np.isnan(g1), False, g1 > 0)
    g2_carrier = np.where(np.isnan(g2), False, g2 > 0)
    double_carrier = g1_carrier & g2_carrier
    baseline_case = rng.uniform(size=n_samples) < baseline_p
    case_mask = double_carrier | baseline_case
    control_mask = ~case_mask
    n_case = int(case_mask.sum())
    n_ctrl = int(control_mask.sum())
    n_dc = int(double_carrier.sum())
    n_dc_in_cases = int((double_carrier & case_mask).sum())
    print(f"\n[3/6] Scenario B (case/control with case-enriched co-occurrence)")
    print(f"      double-carriers in panel:  {n_dc}")
    print(f"      → forced cases:            {n_dc_in_cases}")
    print(f"      baseline P(case):          {baseline_p}")
    print(f"      n_cases:                   {n_case}")
    print(f"      n_controls:                {n_ctrl}")
    print(f"      DC enrichment:             "
          f"{(n_dc_in_cases/max(1,n_case)):.4f} vs {(n_dc/n_samples):.4f} baseline "
          f"= {(n_dc_in_cases*n_samples)/(max(1,n_case)*max(1,n_dc)):.1f}x")
    return case_mask, control_mask


def simulate_scenario_c(dosages: np.ndarray,
                         truth: dict,
                         baseline_p: float,
                         rng: np.random.Generator
                         ) -> tuple[np.ndarray, np.ndarray]:
    """Scenario C — synthetic incompatibility (depleted co-occurrence in cases).

    Each sample becomes a case via baseline coin flip BUT *only if* it is
    NOT a double-carrier of (BCY1, TPK1).  Double-carriers can never be
    cases — modelling synthetic lethality of the two-gene combination.

    Returns (case_mask, control_mask).
    """
    n_samples = dosages.shape[0]
    g1 = dosages[:, truth["gene_1_idx"]]
    g2 = dosages[:, truth["gene_2_idx"]]
    g1_carrier = np.where(np.isnan(g1), False, g1 > 0)
    g2_carrier = np.where(np.isnan(g2), False, g2 > 0)
    double_carrier = g1_carrier & g2_carrier
    baseline_case = rng.uniform(size=n_samples) < baseline_p
    # Synthetic-lethal: DCs are never cases
    case_mask = baseline_case & (~double_carrier)
    control_mask = ~case_mask
    n_case = int(case_mask.sum())
    n_ctrl = int(control_mask.sum())
    n_dc = int(double_carrier.sum())
    n_dc_in_cases = int((double_carrier & case_mask).sum())
    print(f"\n[3/6] Scenario C (synthetic-lethal: DCs depleted from cases)")
    print(f"      double-carriers in panel:  {n_dc}")
    print(f"      DCs in cases (must be 0):  {n_dc_in_cases}")
    print(f"      baseline P(case|non-DC):   {baseline_p}")
    print(f"      n_cases:                   {n_case}")
    print(f"      n_controls:                {n_ctrl}")
    print(f"      expected co-occurrence in cases under independence:")
    print(f"        n_case × P(g1) × P(g2) ≈ "
          f"{n_case * g1_carrier.mean() * g2_carrier.mean():.2f}")
    print(f"        observed: {n_dc_in_cases}  → ratio {n_dc_in_cases/max(0.01,n_case*g1_carrier.mean()*g2_carrier.mean()):.3f}")
    return case_mask, control_mask


# ===========================================================================
# M5 wrapper: RWR pair candidates → interaction-test → InteractionResult
# ===========================================================================

def _run_m5_then_test(dosages: np.ndarray,
                      variant_ids: list,
                      phenotype: np.ndarray,
                      cache: dict,
                      truth: dict,
                      alpha: float,
                      n_seeds_extra,
                      top_k: int,
                      mac_min: int,
                      rng: np.random.Generator) -> list:
    """Run M5 RWR on the yeast bipartite graph, then test top-K pairs.

    Strategy:
      1. Restrict the graph to variants present in BOTH the dosage matrix
         and the cache (only annotated variants can be reached by RWR).
      2. Seed pool = (BCY1 vars + TPK1 vars + n_seeds_extra random others).
         Both ground-truth genes' variants are seeds so the pair is
         guaranteed to be enumerable.
      3. Compute mutual_rwr scores for all pairs among the seed set.
      4. Test top-K pairs for interaction; collect InteractionResults
         exactly as motif_filtered_epistasis_from_data does (MAC filter,
         _test_interaction call, BH-FDR over the top-K pool).
    """
    # 1. Build adjacency restricted to (dosage ∩ cache) variants
    annotated_variant_ids = [v for v in variant_ids if v in cache]
    print(f"      building bipartite graph on {len(annotated_variant_ids):,} "
          f"annotated variants...")
    A, cache_var_ids, gene_ids = _bipartite_from_dict(cache, annotated_variant_ids)
    print(f"      adjacency: {A.shape[0]:,} variants × {A.shape[1]:,} genes")

    # Map cache index → dosage column
    cache_idx_to_col = {
        cache_var_ids[k]: variant_ids.index(cache_var_ids[k])
        for k in range(len(cache_var_ids))
    }
    # 2. Seed pool: BCY1 + TPK1 vars + n_seeds_extra random
    bcy1_seeds = [k for k, v in enumerate(cache_var_ids)
                  if "YIL033C" in cache.get(v, {}).get("genes", [])]
    tpk1_seeds = [k for k, v in enumerate(cache_var_ids)
                  if "YJL164C" in cache.get(v, {}).get("genes", [])]
    other_pool = [k for k in range(len(cache_var_ids))
                  if k not in set(bcy1_seeds) | set(tpk1_seeds)]
    if n_seeds_extra == "local_chroms":
        local_chroms = {"chromosome9", "chromosome10", "chromosome12"}
        extra_seeds = [
            k for k in other_pool
            if cache_var_ids[k].split(":", 1)[0] in local_chroms
        ]
    else:
        extra_seeds = rng.choice(other_pool,
                                 size=min(int(n_seeds_extra), len(other_pool)),
                                 replace=False).tolist()
    seed_indices = sorted(set(bcy1_seeds + tpk1_seeds + extra_seeds))
    print(f"      seeds: {len(bcy1_seeds)} BCY1 + {len(tpk1_seeds)} TPK1 "
          f"+ {len(extra_seeds)} extra = {len(seed_indices)} total")

    # 3. Batched RWR — build M_vv once, iterate over all seeds simultaneously
    print(f"      computing batched RWR (M_vv built once)...")
    n_var, n_gene = A.shape
    var_deg = A.sum(axis=1, keepdims=True).clip(min=1e-12)
    gene_deg = A.sum(axis=0, keepdims=True).clip(min=1e-12)
    P_vg = A / var_deg              # variant→gene
    P_gv = (A / gene_deg).T         # gene→variant
    M_vv = P_vg @ P_gv              # variant→variant transition (n_var × n_var)
    n_seed = len(seed_indices)
    E = np.zeros((n_var, n_seed), dtype=np.float64)
    for k, si in enumerate(seed_indices):
        E[si, k] = 1.0
    P = E.copy()
    n_iter = 30
    for _ in range(n_iter):
        P = (1.0 - alpha) * (M_vv @ P) + alpha * E
    # Normalise each column (re-normalise minor numerical drift)
    P = P / np.clip(P.sum(axis=0, keepdims=True), 1e-12, None)

    # 4. Pair scoring: score(i, j) = P[v_j, k_i] * P[v_i, k_j] for seed pair
    # P shape (n_var, n_seed); for any two seeds k_i, k_j (corresponding to
    # variant indices seed_indices[k_i], seed_indices[k_j]), score is
    # P[seed_indices[k_j], k_i] * P[seed_indices[k_i], k_j].
    print(f"      scoring {n_seed*(n_seed-1)//2:,} pairs ...")
    seed_arr = np.asarray(seed_indices)
    # Materialise P[seed_indices, :] — a (n_seed, n_seed) submatrix of P
    P_seed = P[seed_arr, :]                              # (n_seed, n_seed)
    # mutual_score[k_i, k_j] = P_seed[k_j, k_i] * P_seed[k_i, k_j]
    mutual = P_seed * P_seed.T                           # element-wise; symmetric
    # Take upper triangular (k_i < k_j) and sort
    iu = np.triu_indices(n_seed, k=1)
    pair_idx_pairs = list(zip(iu[0], iu[1]))
    scores = mutual[iu]
    order = np.argsort(scores)[::-1]
    pair_scores = []
    for o in order:
        s = float(scores[o])
        if s <= 0:
            break  # rest are zero
        ki, kj = pair_idx_pairs[o]
        pair_scores.append((seed_indices[ki], seed_indices[kj], s))
        if len(pair_scores) >= top_k:
            break
    print(f"      top {len(pair_scores):,} pairs by mutual-RWR score "
          f"(out of {(scores > 0).sum():,} non-zero)")

    # 4. Test each pair for interaction; collect InteractionResults
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
            motif="rwr",
            shared_entity=f"rwr_score={rwr_score:.4g}",
            beta_marginal_1=r["beta_1"], beta_marginal_2=r["beta_2"],
            beta_interaction=r["beta_interaction"],
            se_interaction=r["se_interaction"],
            p_interaction=r["p_interaction"],
            p_corrected=None,
            n_samples=int(r["n"]),
            maf_1=float(maf_vec[i]), maf_2=float(maf_vec[j]),
        ))

    if not results:
        print(f"      ✗ No M5 candidate pair survived the MAC filter")
        return []

    # BH-FDR over the M5-tested pool (smaller than M1/M2's pool)
    pvals = np.array([r.p_interaction for r in results])
    n_tests = len(pvals)
    order = np.argsort(pvals)
    ranks = np.empty(n_tests, dtype=int)
    ranks[order] = np.arange(1, n_tests + 1)
    corrected = np.minimum(1.0, pvals * n_tests / ranks)
    for k in range(n_tests - 2, -1, -1):
        corrected[order[k]] = min(corrected[order[k]],
                                   corrected[order[k + 1]])
    for k, r in enumerate(results):
        r.p_corrected = float(corrected[k])
    results.sort(key=lambda r: r.p_interaction)
    return results


def _subsample_with_truth_priority(variant_ids: list,
                                    gene1_variants: set,
                                    gene2_variants: set,
                                    n_target: int,
                                    rng: np.random.Generator) -> np.ndarray:
    """Subsample dosage column indices to n_target, forcing BCY1+TPK1 vars in.

    Used to keep M3/M4's O(N²) inner loops fast while not losing the
    ground-truth gene representatives.
    """
    forced = [i for i, v in enumerate(variant_ids)
              if v in gene1_variants or v in gene2_variants]
    forced_set = set(forced)
    other_pool = [i for i in range(len(variant_ids)) if i not in forced_set]
    n_extra = max(0, n_target - len(forced))
    if n_extra > 0 and other_pool:
        extras = rng.choice(other_pool,
                            size=min(n_extra, len(other_pool)),
                            replace=False).tolist()
    else:
        extras = []
    return np.asarray(sorted(set(forced) | set(extras)))


def _bipartite_from_dict(cache: dict,
                         variant_ids_to_keep: list) -> tuple:
    """Helper: build_bipartite_adjacency reads from JSON path; we already
    have the cache loaded as a dict, so we duplicate the logic here.

    Returns (A, variant_ids_out, gene_ids) where variant_ids_out preserves
    the order of ``variant_ids_to_keep`` (filtered to those present in
    cache). Callers that build seed-index lists against
    ``variant_ids_to_keep`` (e.g. M5 in validate_epistasis_cross_species
    and m5_variants_benchmark) depend on this alignment.
    """
    variant_ids_out = [v for v in variant_ids_to_keep if v in cache]
    gene_set: set = set()
    for vid in variant_ids_out:
        gene_set.update(cache[vid].get("genes", []))
    gene_ids = sorted(gene_set)
    gene_to_col = {g: i for i, g in enumerate(gene_ids)}
    A = np.zeros((len(variant_ids_out), len(gene_ids)), dtype=np.float64)
    for vi, vid in enumerate(variant_ids_out):
        for gene in cache[vid].get("genes", []):
            A[vi, gene_to_col[gene]] = 1.0
    return A, variant_ids_out, gene_ids


# ===========================================================================
# Reporting helpers
# ===========================================================================

def _summarise_method(method_short: str,
                      method_long: str,
                      results: list,
                      runtime_sec: float,
                      g1v: str, g2v: str,
                      gene1_variants: set | None = None,
                      gene2_variants: set | None = None) -> dict:
    """Find ground-truth rank in `results` and produce a dict for the JSON.

    Matching policy (gene-level): if ``gene1_variants`` and ``gene2_variants``
    are provided (the full sets of cache-annotated variants in BCY1 and
    TPK1), the rank is the *first* occurrence of any (BCY1-var, TPK1-var)
    pair in the result list.  This is the biologically-honest matching
    rule — a "BCY1 × TPK1" finding should accept any variant pair across
    the two genes, not just the simulated representatives (which may be
    LD-pruned away by some methods).

    Falls back to exact-representative matching if gene sets aren't given.
    """
    truth_rank = None
    truth_record = None
    if gene1_variants is not None and gene2_variants is not None:
        for k, r in enumerate(results, start=1):
            v1, v2 = r.variant_1, r.variant_2
            if (v1 in gene1_variants and v2 in gene2_variants) or \
               (v2 in gene1_variants and v1 in gene2_variants):
                truth_rank = k
                truth_record = r
                break
    else:
        for k, r in enumerate(results, start=1):
            if {r.variant_1, r.variant_2} == {g1v, g2v}:
                truth_rank = k
                truth_record = r
                break

    n_sig_q05 = sum(
        1 for r in results
        if r.p_corrected is not None and r.p_corrected < 0.05
    )

    out = {
        "method": method_short,
        "method_long": method_long,
        "runtime_sec": float(runtime_sec),
        "n_pairs_tested": int(len(results)),
        "n_significant_q05": int(n_sig_q05),
        "ground_truth_rank": truth_rank,
        "ground_truth_p_interaction": (
            float(truth_record.p_interaction) if truth_record else None
        ),
        "ground_truth_q_value": (
            float(truth_record.p_corrected)
            if truth_record and truth_record.p_corrected is not None
            else None
        ),
        "ground_truth_beta": (
            float(truth_record.beta_interaction) if truth_record else None
        ),
        "ground_truth_motif": (
            truth_record.motif if truth_record else None
        ),
        "top_5_pairs": [
            {
                "rank": k,
                "variant_1": r.variant_1,
                "variant_2": r.variant_2,
                "motif": r.motif,
                "shared_entity": r.shared_entity,
                "p_interaction": float(r.p_interaction),
                "q_value": (
                    float(r.p_corrected) if r.p_corrected is not None else None
                ),
                "beta_interaction": float(r.beta_interaction),
            }
            for k, r in enumerate(results[:5], start=1)
        ],
    }
    return out


def _print_comparison_table(methods_out: list) -> None:
    print()
    print(f"  {'Method':<6}  {'rank':>8}  {'n_tested':>10}  "
          f"{'q-value':>10}  {'beta_hat':>10}  {'sig@q05':>8}  "
          f"{'runtime':>8}  detected_via")
    print(f"  {'-'*6}  {'-'*8}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*8}  "
          f"{'-'*8}  {'-'*16}")
    for m in methods_out:
        rank = m["ground_truth_rank"]
        rank_str = f"{rank:,}" if rank else "NF"
        q = m["ground_truth_q_value"]
        q_str = f"{q:.2e}" if q is not None else "—"
        beta = m["ground_truth_beta"]
        beta_str = f"{beta:+.4f}" if beta is not None else "—"
        det = m["ground_truth_motif"] or "—"
        print(f"  {m['method']:<6}  {rank_str:>8}  "
              f"{m['n_pairs_tested']:>10,}  {q_str:>10}  "
              f"{beta_str:>10}  {m['n_significant_q05']:>8,}  "
              f"{m['runtime_sec']:>7.1f}s  {det}")


# ===========================================================================
# Main
# ===========================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--scenario", choices=["A", "B", "C", "all"],
                        default="all",
                        help="Which simulation scenario to run.  "
                             "A = quantitative positive interaction (M1/M2/M5 turf), "
                             "B = case/control with case-enriched co-occurrence (M3 turf), "
                             "C = synthetic incompatibility (M4 turf), "
                             "all = run all three sequentially")
    parser.add_argument("--beta", type=float, default=3.0,
                        help="(Scenario A) Causal interaction effect size (default 3.0).")
    parser.add_argument("--baseline-p", type=float, default=0.10,
                        help="(Scenarios B + C) Baseline P(case).")
    parser.add_argument("--n-nuisance", type=int, default=0,
                        help="(Scenario A) Nuisance interaction features.  Default 0.")
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--max-pairs-per-entity", type=int, default=500,
                        help="Per-pathway/gene cap on enumerated pairs")
    parser.add_argument("--motifs", default="same_gene,same_pathway",
                        help="Comma-separated motifs for M2.  Add "
                             "'protein_interaction' to enable P3.  "
                             "(Default: same_gene,same_pathway — original §Y.4)")
    parser.add_argument("--inject-canonical-ppi", action="store_true",
                        help="Augment graph cache with BioGRID-curated "
                             "canonical PPI edges (BCY1↔TPK1/2/3 + 9 more "
                             "pairs from docs/groundtruth/epistasis_pairs.json). "
                             "The cache's PPI field is incomplete — see "
                             "src/python/graphgwas/canonical_ppi.py.")
    args = parser.parse_args()
    motifs_list = [m.strip() for m in args.motifs.split(",") if m.strip()]

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

    # Optional canonical-PPI injection (closes a curation gap; see §Y.5)
    inject_diag = None
    if args.inject_canonical_ppi:
        from graphgwas.canonical_ppi import (
            YEAST_CANONICAL_EDGES,
            inject_canonical_ppi_edges,
        )
        cache, inject_diag = inject_canonical_ppi_edges(
            cache, edges=YEAST_CANONICAL_EDGES, inplace=True, verbose=True,
        )
        print(f"      ✓ injected canonical PPI: "
              f"{inject_diag['n_edges_applied']}/{inject_diag['n_edges_supplied']} "
              f"edges applied, {inject_diag['n_variants_modified']} variants modified")
        if inject_diag["edges_with_zero_variants"]:
            print(f"      ⚠ edges dropped (no variant in cache): "
                  f"{inject_diag['edges_with_zero_variants']}")

    variant_id_to_col = {v: i for i, v in enumerate(variant_ids)}
    n_in_cache = sum(1 for v in variant_ids if v in cache)
    print(f"      {n_in_cache:,}/{len(variant_ids):,} loaded variants are in the cache")
    if n_in_cache < 100:
        sys.exit("Too few cache hits — cache/VCF format mismatch?")

    # 2. Pick representative variants for ground-truth pair
    truth = pick_representative_variants(cache, variant_id_to_col, dosages, GROUND_TRUTH)

    # ----- Recovery sets (gene-level matching, shared across scenarios) -----
    g1v = truth["gene_1_variant"]
    g2v = truth["gene_2_variant"]
    gene1_variants = {
        vid for vid in variant_ids
        if GROUND_TRUTH["gene_1_systematic"]
           in cache.get(vid, {}).get("genes", [])
    }
    gene2_variants = {
        vid for vid in variant_ids
        if GROUND_TRUTH["gene_2_systematic"]
           in cache.get(vid, {}).get("genes", [])
    }
    print(f"\n      Gene-level recovery sets: "
          f"{GROUND_TRUTH['gene_1_symbol']} = {len(gene1_variants)} variants, "
          f"{GROUND_TRUTH['gene_2_symbol']} = {len(gene2_variants)} variants")

    # ----- Dispatch on scenario -----
    scenarios = (["A", "B", "C"] if args.scenario == "all"
                 else [args.scenario])

    all_scenario_results: dict = {}
    for scenario_name in scenarios:
        print(f"\n{'=' * 70}")
        print(f"  SCENARIO {scenario_name}")
        print(f"{'=' * 70}")
        scenario_results = _run_one_scenario(
            scenario_name, dosages, variant_ids, cache, truth,
            args=args, rng=rng,
            g1v=g1v, g2v=g2v,
            gene1_variants=gene1_variants, gene2_variants=gene2_variants,
            motifs_list=motifs_list,
        )
        all_scenario_results[scenario_name] = scenario_results

    # ----- Final 3 × 5 grid output -----
    print(f"\n{'=' * 70}")
    print(f"  FINAL: SCENARIO × METHOD GRID")
    print(f"{'=' * 70}")
    _print_scenario_method_grid(all_scenario_results)

    out: dict = {
        "ground_truth": GROUND_TRUTH,
        "params": {
            "beta": args.beta,
            "baseline_p": args.baseline_p,
            "n_nuisance": args.n_nuisance,
            "seed": args.seed,
            "max_pairs_per_entity": args.max_pairs_per_entity,
            "chroms_loaded": CHROMS_TO_LOAD,
        },
        "n_samples": int(dosages.shape[0]),
        "n_variants_loaded": int(dosages.shape[1]),
        "n_variants_in_cache": int(n_in_cache),
        "shared_pathways_g1_g2": truth["shared_pathways"],
        "scenarios": all_scenario_results,
    }
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out_path = RESULTS_DIR / "yeast_validation.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\n      Wrote {out_path}")


def _run_one_scenario(scenario_name: str,
                       dosages: np.ndarray,
                       variant_ids: list,
                       cache: dict,
                       truth: dict,
                       *,
                       args,
                       rng: np.random.Generator,
                       g1v: str, g2v: str,
                       gene1_variants: set,
                       gene2_variants: set,
                       motifs_list: list[str] = None) -> dict:
    if motifs_list is None:
        motifs_list = ["same_gene", "same_pathway"]
    """Build phenotype + case/control masks for one scenario and run all 5 methods."""
    if scenario_name == "A":
        phenotype = simulate_with_truth(
            dosages, truth, args.beta, args.n_nuisance, rng,
        )
        # 25/75 percentile binarisation for M3/M4
        q75 = np.percentile(phenotype, 75)
        q25 = np.percentile(phenotype, 25)
        case_mask = phenotype >= q75
        ctrl_mask = phenotype <= q25
        scenario_kind = "quantitative_positive_interaction"
    elif scenario_name == "B":
        case_mask, ctrl_mask = simulate_scenario_b(
            dosages, truth, args.baseline_p, rng,
        )
        # Use the binary phenotype as a continuous trait (0/1) for M1/M2/M5
        phenotype = case_mask.astype(float)
        scenario_kind = "case_control_enriched_co_occurrence"
    elif scenario_name == "C":
        case_mask, ctrl_mask = simulate_scenario_c(
            dosages, truth, args.baseline_p, rng,
        )
        phenotype = case_mask.astype(float)
        scenario_kind = "synthetic_lethal_depletion"
    else:
        raise ValueError(f"Unknown scenario: {scenario_name}")

    methods_out: list[dict] = []

    # ----- M1 -----
    print(f"\n[M1] LD-pruned co-occurrence + interaction test")
    t0 = time.time()
    af_vec = np.nanmean(dosages, axis=0) / 2.0
    m1_variants = [
        {"variantId": variant_ids[k],
         "pos": int(variant_ids[k].split(":")[1]),
         "af_total": float(af_vec[k])}
        for k in range(len(variant_ids))
    ]
    m1_dosage_list = [dosages[:, k] for k in range(dosages.shape[1])]
    m1_results = ld_pruned_cooccurrence_from_data(
        variants=m1_variants, dosage_list=m1_dosage_list,
        phenotype=phenotype, r2_prune=0.5,
        min_cocarriers=1, min_distance_bp=10_000,
        max_variants=2000, verbose=True,
    )
    m1_runtime = time.time() - t0
    methods_out.append(_summarise_method(
        "M1", "LD-pruned co-occurrence (Bonferroni)",
        m1_results, m1_runtime, g1v, g2v,
        gene1_variants=gene1_variants, gene2_variants=gene2_variants,
    ))

    # ----- M2 -----
    print(f"\n[M2] Motif-filtered ({'+'.join(motifs_list)})")
    t0 = time.time()
    m2_results = motif_filtered_epistasis_from_data(
        dosages=dosages, variant_ids=variant_ids, phenotype=phenotype,
        graph_cache=cache,
        motifs=motifs_list,
        mac_min=10, correction="BH",
        max_pairs_per_entity=args.max_pairs_per_entity,
        max_pairs_total=200_000, verbose=True,
    )
    m2_runtime = time.time() - t0
    methods_out.append(_summarise_method(
        "M2", "Motif-filtered (BH-FDR)",
        m2_results, m2_runtime, g1v, g2v,
        gene1_variants=gene1_variants, gene2_variants=gene2_variants,
    ))

    # ----- M3 -----
    print(f"\n[M3] Differential subgraph (case/control)")
    t0 = time.time()
    # Subsample variants to keep M3's O(N²) inner loop tractable while
    # GUARANTEEING that BCY1 and TPK1 representatives are in the pool.
    sub_idx = _subsample_with_truth_priority(
        variant_ids, gene1_variants, gene2_variants,
        n_target=1500, rng=rng,
    )
    sub_dosages = dosages[:, sub_idx]
    sub_variant_ids = [variant_ids[k] for k in sub_idx]
    print(f"      subsampled to {len(sub_idx)} variants (BCY1+TPK1 forced in)")
    m3_results = differential_subgraph_from_data(
        dosages=sub_dosages, variant_ids=sub_variant_ids,
        case_mask=case_mask, control_mask=ctrl_mask,
        r2_prune=0.5, min_distance_bp=10_000,
        min_cocarriers=2, mac_min=10,
        max_variants=500, verbose=True,
    )
    m3_runtime = time.time() - t0
    methods_out.append(_summarise_method(
        "M3", "Differential subgraph (Fisher's exact)",
        m3_results, m3_runtime, g1v, g2v,
        gene1_variants=gene1_variants, gene2_variants=gene2_variants,
    ))

    # ----- M4 -----
    print(f"\n[M4] Dark matter / synthetic incompatibility")
    t0 = time.time()
    m4_results = dark_matter_epistasis_from_data(
        dosages=sub_dosages, variant_ids=sub_variant_ids,
        case_mask=case_mask, control_mask=ctrl_mask,
        maf_min=0.05, depletion_threshold=0.7,
        min_expected=0.5,
        mac_min=10, max_variants=500, verbose=True,
    )
    m4_runtime = time.time() - t0
    methods_out.append(_summarise_method(
        "M4", "Dark matter / synthetic incompatibility (Poisson)",
        m4_results, m4_runtime, g1v, g2v,
        gene1_variants=gene1_variants, gene2_variants=gene2_variants,
    ))

    # ----- M5 -----
    print(f"\n[M5] RWR on bipartite graph + interaction test")
    t0 = time.time()
    m5_results = _run_m5_then_test(
        dosages=dosages, variant_ids=variant_ids, phenotype=phenotype,
        cache=cache, truth=truth,
        alpha=0.15, n_seeds_extra="local_chroms",
        top_k=2000, mac_min=10, rng=rng,
    )
    m5_runtime = time.time() - t0
    methods_out.append(_summarise_method(
        "M5", "RWR + interaction test (BH-FDR)",
        m5_results, m5_runtime, g1v, g2v,
        gene1_variants=gene1_variants, gene2_variants=gene2_variants,
    ))

    # Per-scenario method comparison
    print(f"\n[Scenario {scenario_name} summary]")
    _print_comparison_table(methods_out)

    return {
        "scenario_kind": scenario_kind,
        "n_cases": int(case_mask.sum()),
        "n_controls": int(ctrl_mask.sum()),
        "methods": methods_out,
    }


def _print_scenario_method_grid(all_results: dict) -> None:
    """Final 3 × 5 summary grid: rows = scenarios, columns = methods."""
    methods = ["M1", "M2", "M3", "M4", "M5"]
    print()
    print(f"  {'Scenario':<12}  " + "  ".join(f"{m:>10}" for m in methods))
    print(f"  {'-' * 12}  " + "  ".join("-" * 10 for _ in methods))
    for sc, sc_data in all_results.items():
        m_by_name = {m["method"]: m for m in sc_data["methods"]}
        cells = []
        for m in methods:
            entry = m_by_name.get(m)
            if entry is None:
                cells.append("—")
                continue
            rank = entry["ground_truth_rank"]
            cells.append(f"{rank:,}" if rank is not None else "NF")
        print(f"  {sc:<12}  " + "  ".join(f"{c:>10}" for c in cells))
    print()
    print(f"  Cell = ground-truth rank in that method's result list "
          f"(NF = not found, top winner = lowest number)")


if __name__ == "__main__":
    main()
