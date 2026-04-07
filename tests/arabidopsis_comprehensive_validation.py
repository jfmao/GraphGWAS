"""Comprehensive Arabidopsis GWAS validation: 536 traits + eQTL.

Validates GraphGWAS on Arabidopsis 1001 Genomes (1,135 accessions):
  1. GWAS on 536 organismal traits (vectorized, with population correction)
  2. Heritability estimation across traits
  3. Comparison to AraGWAS ground truth (275 Bonferroni hits)
  4. eQTL scan (100 random genes from 24K expression matrix)
  5. Multi-trait genetic correlation (ionomics: 18 elements)
  6. Pleiotropy analysis

Usage: python tests/arabidopsis_comprehensive_validation.py
"""

import csv
import glob
import gzip
import json
import os
import sys
import time

import numpy as np
from scipy import stats as sp_stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "python"))

DATA_DIR = os.path.join(os.path.dirname(__file__), "data", "arabidopsis")
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results", "arabidopsis_validation")
os.makedirs(RESULTS_DIR, exist_ok=True)


# ===================================================================
# Data loading
# ===================================================================

def load_genotypes_vcf(chrom=None, max_variants=200000):
    """Load genotypes from VCF using cyvcf2. Returns dosage matrix."""
    from cyvcf2 import VCF

    vcf_path = os.path.join(DATA_DIR,
                            "1001genomes_snp-short-indel_only_ACGTN.vcf.gz")
    vcf = VCF(vcf_path)
    sample_ids = vcf.samples

    print(f"  Loading genotypes from VCF ({len(sample_ids)} accessions)...")
    if chrom:
        print(f"    Chromosome: {chrom}")

    t0 = time.time()
    positions = []
    chroms = []
    dosages = []

    iterator = vcf(chrom) if chrom else vcf
    n = 0
    for v in iterator:
        if len(v.ALT) != 1:  # biallelic only
            continue
        gt = v.gt_types  # 0=HOM_REF, 1=HET, 2=UNKNOWN, 3=HOM_ALT
        # Convert to dosage: 0→0, 1→1, 3→2, 2→NaN
        dosage = np.where(gt == 0, 0.0,
                 np.where(gt == 1, 1.0,
                 np.where(gt == 3, 2.0, np.nan)))

        # MAF filter
        af = np.nanmean(dosage) / 2
        maf = min(af, 1 - af)
        if maf < 0.05:  # stricter MAF for speed
            continue

        dosages.append(dosage)
        positions.append(v.POS)
        chroms.append(v.CHROM)
        n += 1

        if n >= max_variants:
            break
        if n % 50000 == 0:
            print(f"    ...{n:,} variants loaded ({time.time()-t0:.0f}s)")

    vcf.close()
    geno = np.array(dosages, dtype=np.float64)
    elapsed = time.time() - t0
    print(f"    {geno.shape[0]:,} variants × {geno.shape[1]} accessions "
          f"({elapsed:.0f}s)")

    return geno, sample_ids, positions, chroms


def load_phenotype_csv(path):
    """Load one AraPheno phenotype CSV. Returns {accession_id: value}."""
    data = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            acc_id = row.get("accession_id", "")
            val = row.get("phenotype_value", "")
            if acc_id and val:
                try:
                    data[acc_id] = float(val)
                except ValueError:
                    pass
    return data


def load_all_phenotypes(min_samples=50):
    """Load all AraPheno phenotype files."""
    pheno_dir = os.path.join(DATA_DIR, "phenotypes")
    catalog_path = os.path.join(DATA_DIR, "phenotype_catalog.csv")

    # Load catalog for names
    catalog = {}
    if os.path.exists(catalog_path):
        with open(catalog_path) as f:
            for row in csv.DictReader(f):
                catalog[row.get("phenotype_id", "")] = row

    phenotypes = {}
    for fname in sorted(os.listdir(pheno_dir)):
        if not fname.endswith(".csv"):
            continue
        pid = fname.replace("pheno_", "").replace(".csv", "")
        data = load_phenotype_csv(os.path.join(pheno_dir, fname))
        if len(data) >= min_samples:
            name = catalog.get(pid, {}).get("name", pid)
            phenotypes[name] = {"data": data, "pid": pid, "n": len(data)}

    return phenotypes


def load_aragwas_ground_truth():
    """Load AraGWAS Bonferroni-significant associations."""
    path = os.path.join(DATA_DIR, "aragwas_bonf_associations.csv")
    if not os.path.exists(path):
        return {}
    hits = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            pheno = row.get("phenotype", "")
            if pheno not in hits:
                hits[pheno] = []
            hits[pheno].append({
                "chr": row.get("snp_chr", ""),
                "pos": int(row.get("snp_pos", 0)),
                "pvalue": float(row.get("pvalue", 0)),
            })
    return hits


def align_data(geno, sample_ids, pheno_dict):
    """Align genotype matrix to phenotype dict by accession ID."""
    # AraPheno uses numeric accession IDs; VCF uses the same
    shared = [s for s in sample_ids if s in pheno_dict]
    if len(shared) < 30:
        return None, None
    id_to_idx = {s: i for i, s in enumerate(sample_ids)}
    cols = np.array([id_to_idx[s] for s in shared])
    pheno_vec = np.array([pheno_dict[s] for s in shared])
    return geno[:, cols], pheno_vec


# ===================================================================
# Main
# ===================================================================

def main():
    t_start = time.time()
    print("=" * 70)
    print("  GraphGWAS Arabidopsis Comprehensive Validation")
    print("  536 traits + eQTL × 1,135 accessions")
    print("=" * 70)

    # --- Load genotypes (chr4 first — contains FLC/FRI) ---
    # Use chr4 for focused validation, then expand
    print("\n--- Loading genotypes ---")
    geno, sample_ids, positions, chroms_list = load_genotypes_vcf(
        chrom="4", max_variants=200000
    )

    from graphgwas.sv import fast_gwas
    from graphgwas.popstruct import compute_grm, corrected_gwas, lambda_gc

    # --- Load phenotypes ---
    print("\n--- Loading phenotypes ---")
    phenotypes = load_all_phenotypes(min_samples=50)
    print(f"  {len(phenotypes)} traits loaded")

    ground_truth = load_aragwas_ground_truth()
    print(f"  AraGWAS ground truth: {sum(len(v) for v in ground_truth.values())} "
          f"hits across {len(ground_truth)} traits")

    # --- Compute GRM for population correction ---
    print("\n--- Computing GRM ---")
    # Use a random subset for GRM (for speed)
    rng = np.random.default_rng(42)
    grm_idx = rng.choice(geno.shape[0], min(20000, geno.shape[0]), replace=False)
    grm = compute_grm(geno[grm_idx], verbose=True)

    # =================================================================
    # Phase 1: GWAS on all traits (with PCA correction)
    # =================================================================
    print(f"\n{'='*70}")
    print(f"  Phase 1: GWAS on {len(phenotypes)} traits (chr4, PCA-corrected)")
    print(f"{'='*70}")

    gwas_results = {}
    t1 = time.time()

    for ti, (name, info) in enumerate(sorted(phenotypes.items())):
        geno_a, pheno_a = align_data(geno, sample_ids, info["data"])
        if geno_a is None:
            continue

        # Quick PCA-corrected GWAS
        # Align GRM to this trait's samples
        shared = [s for s in sample_ids if s in info["data"]]
        id_to_idx = {s: i for i, s in enumerate(sample_ids)}
        cols = np.array([id_to_idx[s] for s in shared])

        # Use corrected_gwas with auto method
        result = corrected_gwas(geno_a, pheno_a, method="pca",
                                n_pcs=5, verbose=False)

        gwas_results[name] = {
            "n_samples": geno_a.shape[1],
            "lambda_gc": result["lambda_gc"],
            "n_bonf": result["n_bonferroni"],
            "n_sug": result["n_suggestive"],
            "min_p": result["min_p"],
            "pvals": result["pvals"],
        }

        if (ti + 1) % 100 == 0 or (ti + 1) == len(phenotypes):
            elapsed = time.time() - t1
            print(f"  ...{ti+1}/{len(phenotypes)} traits ({elapsed:.0f}s)")

    elapsed = time.time() - t1
    print(f"\n  GWAS complete: {len(gwas_results)} traits in {elapsed:.0f}s")
    print(f"  Traits with Bonf hits: "
          f"{sum(1 for r in gwas_results.values() if r['n_bonf'] > 0)}")
    print(f"  Total Bonf hits: "
          f"{sum(r['n_bonf'] for r in gwas_results.values()):,}")

    # Lambda GC distribution
    lgcs = [r["lambda_gc"] for r in gwas_results.values()]
    print(f"  Lambda GC: mean={np.mean(lgcs):.3f}, "
          f"median={np.median(lgcs):.3f}, "
          f"range=[{min(lgcs):.3f}, {max(lgcs):.3f}]")

    # Top traits
    top_traits = sorted(gwas_results.items(), key=lambda x: x[1]["min_p"])[:15]
    print(f"\n  Top 15 traits by min p-value:")
    for name, r in top_traits:
        print(f"    {name:<35s} min_p={r['min_p']:.2e}  "
              f"bonf={r['n_bonf']}  λ={r['lambda_gc']:.2f}  "
              f"n={r['n_samples']}")

    # =================================================================
    # Phase 2: Heritability
    # =================================================================
    print(f"\n{'='*70}")
    print(f"  Phase 2: Heritability estimation")
    print(f"{'='*70}")

    h2_results = {}
    for name, info in sorted(phenotypes.items()):
        geno_a, pheno_a = align_data(geno, sample_ids, info["data"])
        if geno_a is None or geno_a.shape[1] < 50:
            continue

        n_sam = geno_a.shape[1]
        train = np.arange(0, n_sam, 2)
        test = np.arange(1, n_sam, 2)

        betas_t, pvals_t = fast_gwas(geno_a[:, train], pheno_a[train])
        top_idx = np.argsort(pvals_t)[:30]

        G_test = geno_a[:, test].copy()
        rm = np.nanmean(G_test, axis=1)
        for i in range(G_test.shape[0]):
            mask = np.isnan(G_test[i])
            if np.any(mask):
                G_test[i, mask] = rm[i] if not np.isnan(rm[i]) else 0

        prs = np.zeros(len(test))
        for idx in top_idx:
            prs += betas_t[idx] * G_test[idx]

        if np.std(prs) > 0 and np.std(pheno_a[test]) > 0:
            r2 = np.corrcoef(prs, pheno_a[test])[0, 1] ** 2
        else:
            r2 = 0.0
        h2_results[name] = float(r2)

    sorted_h2 = sorted(h2_results.items(), key=lambda x: x[1], reverse=True)
    print(f"  Top 15 heritable traits (PRS R²):")
    for name, h2 in sorted_h2[:15]:
        print(f"    {name:<35s} h²={h2:.4f}")

    h2_vals = [h for _, h in sorted_h2]
    print(f"\n  h² range: [{min(h2_vals):.4f}, {max(h2_vals):.4f}]")
    print(f"  h² mean: {np.mean(h2_vals):.4f}")

    # =================================================================
    # Phase 3: AraGWAS ground truth comparison
    # =================================================================
    if ground_truth:
        print(f"\n{'='*70}")
        print(f"  Phase 3: AraGWAS ground truth comparison")
        print(f"{'='*70}")

        overlap = 0
        tested = 0
        for pheno_name, hits in ground_truth.items():
            chr4_hits = [h for h in hits if h["chr"] == "4"]
            if not chr4_hits:
                continue

            # Find matching trait in our results
            matched_trait = None
            for name in gwas_results:
                if pheno_name.lower().strip() in name.lower():
                    matched_trait = name
                    break

            if matched_trait is None:
                continue

            our_pvals = gwas_results[matched_trait]["pvals"]
            our_sig_idx = set(np.where(our_pvals < 1e-4)[0])

            for hit in chr4_hits:
                tested += 1
                for oi in our_sig_idx:
                    if oi < len(positions):
                        if abs(positions[oi] - hit["pos"]) < 50000:
                            overlap += 1
                            break

        print(f"  Chr4 ground truth hits tested: {tested}")
        print(f"  Overlapping (±50kb): {overlap}")
        if tested > 0:
            print(f"  Overlap rate: {overlap/tested:.1%}")

    # =================================================================
    # Phase 4: eQTL scan (100 genes)
    # =================================================================
    print(f"\n{'='*70}")
    print(f"  Phase 4: eQTL scan (100 genes)")
    print(f"{'='*70}")

    expr_path = os.path.join(DATA_DIR, "expression",
                             "ath1001_tx_norm_counts.tsv.gz")
    if os.path.exists(expr_path):
        print("  Loading expression matrix...")
        with gzip.open(expr_path, "rt") as f:
            header = f.readline().strip().split("\t")
            expr_accessions = [h.replace("X", "") for h in header[1:]]

            # Read 100 random genes
            all_lines = f.readlines()

        rng = np.random.default_rng(42)
        gene_idx = rng.choice(len(all_lines), min(100, len(all_lines)),
                              replace=False)

        eqtl_results = {}
        t4 = time.time()
        for i, gi in enumerate(gene_idx):
            parts = all_lines[gi].strip().split("\t")
            gene_id = parts[0]
            expr_vals = {}
            for j, acc in enumerate(expr_accessions):
                try:
                    expr_vals[acc] = float(parts[j + 1])
                except (ValueError, IndexError):
                    pass

            geno_a, pheno_a = align_data(geno, sample_ids, expr_vals)
            if geno_a is None:
                continue

            _, pvals = fast_gwas(geno_a, pheno_a)
            n_sig = int(np.sum(pvals < 0.05 / geno.shape[0]))
            eqtl_results[gene_id] = {"n_sig": n_sig,
                                     "min_p": float(np.nanmin(pvals))}

            if (i + 1) % 25 == 0:
                print(f"  ...{i+1}/100 genes ({time.time()-t4:.0f}s)")

        n_with_eqtl = sum(1 for r in eqtl_results.values() if r["n_sig"] > 0)
        print(f"\n  {len(eqtl_results)} genes scanned, "
              f"{n_with_eqtl} with significant eQTL on chr4")

        top_eqtl = sorted(eqtl_results.items(),
                          key=lambda x: x[1]["min_p"])[:10]
        print(f"  Top 10:")
        for gene, r in top_eqtl:
            print(f"    {gene:<15s} min_p={r['min_p']:.2e}  hits={r['n_sig']}")

    # =================================================================
    # Phase 5: Ionomics genetic correlation
    # =================================================================
    print(f"\n{'='*70}")
    print(f"  Phase 5: Ionomics genetic correlation")
    print(f"{'='*70}")

    ion_traits = [n for n in gwas_results
                  if any(e in n for e in ["Li7", "B11", "Na23", "Mg25",
                         "P31", "S34", "K39", "Ca43", "Mn55", "Fe57",
                         "Co59", "Ni60", "Cu65", "Zn66", "As75",
                         "Se82", "Rb85", "Mo98", "Cd111"])]

    if len(ion_traits) >= 5:
        T = len(ion_traits)
        z_mat = np.zeros((geno.shape[0], T))
        for ti, name in enumerate(ion_traits):
            pvals = gwas_results[name]["pvals"]
            betas = np.zeros_like(pvals)  # we don't store betas, use sign from pvals
            z = sp_stats.norm.isf(pvals / 2)
            z[~np.isfinite(z)] = 0
            z_mat[:, ti] = z

        rg = np.corrcoef(z_mat.T)
        print(f"  Ionomics genetic correlation ({T} elements):")
        labels = [n[:8] for n in ion_traits[:10]]
        header = "          " + "  ".join(f"{l:>8}" for l in labels)
        print(f"  {header}")
        for i in range(min(10, T)):
            row = f"  {labels[i]:<8}" + "  ".join(
                f"{rg[i,j]:8.3f}" for j in range(min(10, T)))
            print(row)

        off = rg[np.triu_indices(T, k=1)]
        print(f"\n  Off-diagonal: mean={np.mean(off):.3f}, "
              f"range=[{off.min():.3f}, {off.max():.3f}]")

    # =================================================================
    # Phase 6: Pleiotropy
    # =================================================================
    print(f"\n{'='*70}")
    print(f"  Phase 6: Pleiotropy")
    print(f"{'='*70}")

    n_var = geno.shape[0]
    pleio = np.zeros(n_var, dtype=int)
    for r in gwas_results.values():
        pleio += (r["pvals"] < 1e-4).astype(int)

    print(f"  Variants sig for ≥2 traits:  {int(np.sum(pleio >= 2)):,}")
    print(f"  Variants sig for ≥5 traits:  {int(np.sum(pleio >= 5)):,}")
    print(f"  Variants sig for ≥10 traits: {int(np.sum(pleio >= 10)):,}")
    print(f"  Max traits per variant:      {int(pleio.max())}")

    top_pleio = np.argsort(-pleio)[:10]
    print(f"\n  Top 10 pleiotropic variants (chr4):")
    for idx in top_pleio:
        if idx < len(positions):
            print(f"    chr4:{positions[idx]}  → {pleio[idx]} traits")

    # =================================================================
    # Summary
    # =================================================================
    elapsed = time.time() - t_start
    print(f"\n{'='*70}")
    print(f"  SUMMARY")
    print(f"{'='*70}")
    print(f"  Time: {elapsed:.0f}s ({elapsed/60:.1f} min)")
    print(f"  Traits analysed: {len(gwas_results)}")
    print(f"  Bonf-sig traits: "
          f"{sum(1 for r in gwas_results.values() if r['n_bonf'] > 0)}")
    print(f"  Mean λ_GC: {np.mean(lgcs):.3f}")
    print(f"  Mean h²: {np.mean(h2_vals):.4f}")
    print(f"  Pleiotropic variants (≥2): {int(np.sum(pleio >= 2)):,}")

    # Save
    summary = {
        "n_traits": len(gwas_results),
        "n_bonf_traits": sum(1 for r in gwas_results.values() if r["n_bonf"] > 0),
        "mean_lambda_gc": float(np.mean(lgcs)),
        "mean_h2": float(np.mean(h2_vals)),
        "n_pleio2": int(np.sum(pleio >= 2)),
        "elapsed": elapsed,
    }
    with open(os.path.join(RESULTS_DIR, "validation_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n  Saved to {RESULTS_DIR}/")


if __name__ == "__main__":
    main()
