"""Generate test phenotype data for GraphGWAS Phase 1 validation.

Creates two phenotype files:
1. phenotypes_alzheimer_sim.csv — simulated Alzheimer's-like phenotype
   with rs429358 (APOE) as causal variant
2. phenotypes_null.csv — random case/control (null control)

Requires Neo4j connection to read sample IDs and genotypes.
"""

import csv
import sys
import os

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "python"))

from graphgwas.db import GraphGWASConnection
from graphgwas.config import NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD, N_SAMPLES
from graphgwas.genotype import unpack_genotypes


def generate_alzheimer_sim(conn, output_path):
    """Generate simulated Alzheimer's phenotype based on rs429358 (APOE).

    Case probability:
    - TT (hom_ref, GT=0): P(case) = 0.05
    - TC (het, GT=1):      P(case) = 0.15 (OR ~3.4)
    - CC (hom_alt, GT=2):  P(case) = 0.45 (OR ~15)
    """
    # Get all sample IDs and packed indices
    result = conn.execute_read(
        """
        MATCH (s:Sample)
        RETURN s.sampleId AS sid, s.packed_index AS idx, s.population AS pop, s.sex AS sex
        ORDER BY s.packed_index
        """
    )
    samples = [dict(r) for r in result]

    # Find rs429358 or a variant near APOE locus
    result = conn.execute_read(
        """
        MATCH (v:Variant)
        WHERE v.chr = 'chr19' AND v.pos >= 45410000 AND v.pos <= 45420000
        RETURN v.variantId AS vid, v.pos AS pos, v.gt_packed AS gtp, v.af_total AS af
        ORDER BY v.af_total DESC
        LIMIT 5
        """
    )
    apoe_variants = [dict(r) for r in result]

    if not apoe_variants:
        print("WARNING: No APOE-region variants found. Using random phenotype.")
        return generate_null(conn, output_path)

    # Use the variant with highest AF in the APOE region
    apoe = apoe_variants[0]
    print(f"Using APOE variant: {apoe['vid']} (AF={apoe['af']:.4f})")

    gt = unpack_genotypes(apoe["gtp"], N_SAMPLES)

    # Assign case probability based on genotype
    np.random.seed(42)
    case_prob = {0: 0.05, 1: 0.15, 2: 0.45, 3: 0.0}  # 3=missing → not assigned

    rows = []
    for s in samples:
        idx = s["idx"]
        g = int(gt[idx])
        p = case_prob.get(g, 0.0)

        # Add population-correlated noise (confounding)
        pop = s.get("pop", "")
        if pop in ("YRI", "LWK", "GWD", "MSL", "ESN", "ACB", "ASW"):
            p += 0.02  # slight African ancestry effect
        elif pop in ("CHB", "JPT", "CHS", "CDX", "KHV"):
            p -= 0.01  # slight East Asian effect

        p = max(0.01, min(0.99, p))
        is_case = np.random.random() < p
        age = int(np.random.uniform(40, 85))
        sex = s.get("sex", 0)

        rows.append({
            "sample_id": s["sid"],
            "case_control": "case" if is_case else "control",
            "age": age,
            "sex": sex,
            "population": pop,
        })

    n_cases = sum(1 for r in rows if r["case_control"] == "case")
    n_controls = len(rows) - n_cases
    print(f"Generated: {n_cases} cases, {n_controls} controls")

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample_id", "case_control", "age", "sex", "population"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"Written to {output_path}")


def generate_null(conn, output_path):
    """Generate random case/control phenotype (no genetic effect)."""
    result = conn.execute_read(
        """
        MATCH (s:Sample)
        RETURN s.sampleId AS sid, s.population AS pop, s.sex AS sex
        ORDER BY s.sampleId
        """
    )
    samples = [dict(r) for r in result]

    np.random.seed(123)
    rows = []
    for s in samples:
        is_case = np.random.random() < 0.15  # ~15% case rate
        rows.append({
            "sample_id": s["sid"],
            "case_control": "case" if is_case else "control",
            "age": int(np.random.uniform(20, 80)),
            "sex": s.get("sex", 0),
            "population": s.get("pop", ""),
        })

    n_cases = sum(1 for r in rows if r["case_control"] == "case")
    print(f"Null phenotype: {n_cases} cases, {len(rows) - n_cases} controls")

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample_id", "case_control", "age", "sex", "population"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"Written to {output_path}")


if __name__ == "__main__":
    data_dir = os.path.join(os.path.dirname(__file__), "data")
    os.makedirs(data_dir, exist_ok=True)

    with GraphGWASConnection(NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD) as conn:
        generate_alzheimer_sim(conn, os.path.join(data_dir, "phenotypes_alzheimer_sim.csv"))
        generate_null(conn, os.path.join(data_dir, "phenotypes_null.csv"))
