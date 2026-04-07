"""Generate simulation data with known ground truth using msprime.

Creates a VCF with known causal variants and a phenotype CSV,
suitable for validating the entire GraphGWAS pipeline without Neo4j.

Ground truth:
- 1000 samples, ~50K variants on a 1Mb chromosome
- 5 causal variants with known effect sizes
- Binary phenotype (case/control) derived from liability model
- Quantitative phenotype from additive genetic model

Output:
    tests/data/sim/simulation.vcf
    tests/data/sim/phenotypes.csv
    tests/data/sim/ground_truth.json

Usage:
    python tests/generate_simulation.py
"""

import json
import os
from pathlib import Path

import msprime
import numpy as np


def generate_simulation(output_dir: str = "tests/data/sim",
                        n_samples: int = 1000,
                        seq_length: int = 1_000_000,
                        n_causal: int = 5,
                        heritability: float = 0.5,
                        prevalence: float = 0.3,
                        seed: int = 42):
    """Generate simulation with known ground truth.

    Uses msprime to simulate realistic LD structure, then plants
    causal variants with known effect sizes.

    Args:
        output_dir: directory for output files.
        n_samples: number of diploid individuals.
        seq_length: chromosome length in bp.
        n_causal: number of causal variants.
        heritability: narrow-sense h² for the phenotype.
        prevalence: disease prevalence for binary trait.
        seed: random seed.
    """
    os.makedirs(output_dir, exist_ok=True)
    rng = np.random.default_rng(seed)

    print(f"=== Generating Simulation Data ===")
    print(f"  Samples: {n_samples}, Seq length: {seq_length:,}")
    print(f"  Causal variants: {n_causal}, h²: {heritability}")

    # --- 1. Simulate with msprime ---
    print("\n1. Running msprime coalescent simulation...")
    ts = msprime.sim_ancestry(
        samples=n_samples,
        sequence_length=seq_length,
        recombination_rate=1e-8,
        population_size=10_000,
        random_seed=seed,
    )
    ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=seed + 1)

    n_variants = ts.num_sites
    print(f"   Generated {n_variants} variants for {n_samples} samples")

    # --- 2. Extract genotype matrix ---
    print("2. Extracting genotype matrix...")
    # Get variant positions and genotypes
    positions = []
    genotypes = []  # list of arrays, each (n_samples,) diploid dosage

    for var in ts.variants():
        pos = int(var.site.position)
        # Dosage: count alt alleles per diploid individual
        gt = var.genotypes  # haploid (2*n_samples,)
        dosage = gt[0::2] + gt[1::2]  # diploid dosage (n_samples,)
        af = np.mean(dosage) / 2
        # Filter: only biallelic, MAF > 0.01
        if len(var.alleles) == 2 and 0.01 < af < 0.99:
            positions.append(pos)
            genotypes.append(dosage.astype(np.float64))

    genotypes = np.array(genotypes)  # (n_variants_filtered, n_samples)
    positions = np.array(positions)
    n_var_filtered = len(positions)
    print(f"   Filtered to {n_var_filtered} biallelic variants (MAF > 0.01)")

    # --- 3. Select causal variants ---
    print(f"3. Selecting {n_causal} causal variants...")
    # Pick causal variants spread across the chromosome
    causal_indices = np.linspace(0, n_var_filtered - 1, n_causal + 2,
                                  dtype=int)[1:-1]
    # Adjust to exactly n_causal
    causal_indices = causal_indices[:n_causal]

    # Effect sizes: drawn from normal, scaled by MAF
    causal_afs = np.array([np.mean(genotypes[ci]) / 2 for ci in causal_indices])
    raw_effects = rng.normal(0, 1, n_causal)
    # Scale effects so that genetic variance = heritability × total variance
    # Var(genetic) = Σ 2*p*(1-p)*β² for additive model
    var_factors = 2 * causal_afs * (1 - causal_afs)
    # We want h² = Var(G) / (Var(G) + Var(E)), so Var(G) = h²/(1-h²) * Var(E)
    # With Var(E) = 1: Var(G) = h²/(1-h²)
    target_var_g = heritability / (1 - heritability)
    current_var_g = np.sum(var_factors * raw_effects ** 2)
    if current_var_g > 0:
        scale = np.sqrt(target_var_g / current_var_g)
        effects = raw_effects * scale
    else:
        effects = raw_effects

    causal_info = []
    for i, (ci, eff) in enumerate(zip(causal_indices, effects)):
        info = {
            "index": int(ci),
            "position": int(positions[ci]),
            "af": float(causal_afs[i]),
            "effect_size": float(eff),
            "variant_id": f"sim_chr1_{positions[ci]}",
        }
        causal_info.append(info)
        print(f"   Causal {i+1}: pos={positions[ci]:,}, AF={causal_afs[i]:.3f}, "
              f"β={eff:.4f}")

    # --- 4. Generate phenotypes ---
    print("4. Generating phenotypes...")
    # Genetic value: Σ β_j × dosage_ij
    genetic_value = np.zeros(n_samples)
    for ci, eff in zip(causal_indices, effects):
        genetic_value += eff * genotypes[ci]

    # Environmental noise
    env_noise = rng.normal(0, 1, n_samples)

    # Quantitative phenotype
    quant_pheno = genetic_value + env_noise
    realized_h2 = np.var(genetic_value) / np.var(quant_pheno)
    print(f"   Quantitative: realized h²={realized_h2:.4f} (target={heritability})")

    # Binary phenotype via liability threshold
    liability = genetic_value + env_noise
    threshold = np.quantile(liability, 1 - prevalence)
    case_status = (liability >= threshold).astype(int)
    n_cases = int(np.sum(case_status))
    n_controls = n_samples - n_cases
    print(f"   Binary: {n_cases} cases, {n_controls} controls "
          f"(prevalence={n_cases/n_samples:.2f})")

    # --- 5. Write VCF ---
    vcf_path = os.path.join(output_dir, "simulation.vcf")
    print(f"5. Writing VCF to {vcf_path}...")

    sample_names = [f"SIM_{i:04d}" for i in range(n_samples)]

    with open(vcf_path, "w") as f:
        # Header
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##contig=<ID=chr1,length={seq_length}>\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join(sample_names) + "\n")

        # Variant records
        for vi in range(n_var_filtered):
            pos = positions[vi]
            vid = f"sim_chr1_{pos}"
            dosage = genotypes[vi].astype(int)

            # Convert dosage to GT strings
            gt_strings = []
            for d in dosage:
                if d == 0:
                    gt_strings.append("0/0")
                elif d == 1:
                    gt_strings.append("0/1")
                else:
                    gt_strings.append("1/1")

            is_causal = vi in causal_indices
            info = "CAUSAL" if is_causal else "."

            f.write(f"chr1\t{pos}\t{vid}\tA\tT\t.\tPASS\t{info}\tGT\t"
                    + "\t".join(gt_strings) + "\n")

    print(f"   Written {n_var_filtered} variants × {n_samples} samples")

    # --- 6. Write phenotype CSV ---
    pheno_path = os.path.join(output_dir, "phenotypes.csv")
    print(f"6. Writing phenotypes to {pheno_path}...")

    with open(pheno_path, "w") as f:
        f.write("sample_id,case_control,quantitative_trait,genetic_value\n")
        for i in range(n_samples):
            f.write(f"{sample_names[i]},{case_status[i]},"
                    f"{quant_pheno[i]:.6f},{genetic_value[i]:.6f}\n")

    # --- 7. Write ground truth ---
    truth_path = os.path.join(output_dir, "ground_truth.json")
    print(f"7. Writing ground truth to {truth_path}...")

    ground_truth = {
        "n_samples": n_samples,
        "n_variants": n_var_filtered,
        "seq_length": seq_length,
        "heritability_target": heritability,
        "heritability_realized": float(realized_h2),
        "prevalence": float(n_cases / n_samples),
        "n_cases": n_cases,
        "n_controls": n_controls,
        "causal_variants": causal_info,
        "seed": seed,
    }
    with open(truth_path, "w") as f:
        json.dump(ground_truth, f, indent=2)

    print(f"\n=== Simulation Complete ===")
    print(f"  VCF: {vcf_path} ({os.path.getsize(vcf_path) / 1e6:.1f} MB)")
    print(f"  Phenotypes: {pheno_path}")
    print(f"  Ground truth: {truth_path}")

    return ground_truth


if __name__ == "__main__":
    generate_simulation()
