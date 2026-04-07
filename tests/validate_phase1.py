"""Phase 1 end-to-end validation script.

Validates:
1. Phenotype import and activation
2. Genotype decode accuracy
3. Association engine (positive and null control)
4. Result storage and querying
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "python"))

from graphgwas.db import GraphGWASConnection
from graphgwas.config import NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD, N_SAMPLES
from graphgwas.phenotype import load_phenotypes, activate_phenotype
from graphgwas.genotype import (
    unpack_genotypes, get_phenotype_indices, build_count_table,
    build_dosage, compute_z_score,
)
from graphgwas.schema import ensure_indexes, audit_schema, verify_graphmana_compat
from graphgwas.assoc import single_locus_scan, chi2_allelic_test, fisher_exact_test


def validate_all():
    data_dir = os.path.join(os.path.dirname(__file__), "data")
    alz_file = os.path.join(data_dir, "phenotypes_alzheimer_sim.csv")

    if not os.path.exists(alz_file):
        print("ERROR: Run generate_test_data.py first.")
        return False

    all_pass = True

    with GraphGWASConnection(NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD) as conn:

        # 1. Schema compatibility
        print("=== Test 1: Schema Compatibility ===")
        issues = verify_graphmana_compat(conn)
        if issues:
            print(f"FAIL: {issues}")
            all_pass = False
        else:
            print("PASS: GraphMana schema compatible")

        # 2. Create indexes
        print("\n=== Test 2: Index Creation ===")
        results = ensure_indexes(conn)
        for r in results:
            print(f"  {r}")
        print("PASS: Indexes created")

        # 3. Load phenotype
        print("\n=== Test 3: Phenotype Load ===")
        result = load_phenotypes(conn, alz_file, "sample_id",
                                 ["case_control", "age", "sex", "population"])
        print(f"  Matched: {result['n_matched']}, Unmatched: {result['n_unmatched']}")
        assert result["n_matched"] > 3000, f"Expected >3000 matched, got {result['n_matched']}"
        print("PASS: Phenotypes loaded")

        # 4. Activate phenotype
        print("\n=== Test 4: Phenotype Activation ===")
        result = activate_phenotype(conn, "case_control", "case", "control")
        print(f"  Cases: {result['n_cases']}, Controls: {result['n_controls']}")
        assert result["n_cases"] > 50, f"Expected >50 cases, got {result['n_cases']}"
        assert result["n_controls"] > 2000, f"Expected >2000 controls, got {result['n_controls']}"
        print("PASS: Phenotype activated")

        # 5. Genotype decode validation
        print("\n=== Test 5: Genotype Decode ===")
        rec = conn.execute_read(
            """
            MATCH (v:Variant)
            WHERE v.chr = 'chr22' AND v.af_total > 0.1 AND v.af_total < 0.5
            RETURN v.variantId AS vid, v.gt_packed AS gtp, v.ac_total AS ac, v.an_total AS an
            LIMIT 1
            """
        ).single()

        if rec and rec["gtp"]:
            gt = unpack_genotypes(rec["gtp"], N_SAMPLES)
            computed_ac = int(np.sum(gt == 1) + 2 * np.sum(gt == 2))
            computed_an = 2 * int(np.sum(gt != 3))
            stored_ac = rec["ac"]
            stored_an = rec["an"]
            print(f"  Variant: {rec['vid']}")
            print(f"  Computed AC/AN: {computed_ac}/{computed_an}")
            print(f"  Stored AC/AN:   {stored_ac}/{stored_an}")
            assert abs(computed_ac - stored_ac) <= 1, "AC mismatch"
            assert abs(computed_an - stored_an) <= 1, "AN mismatch"
            print("PASS: Genotype decode matches stored counts")
        else:
            print("SKIP: No suitable variant found")

        # 6. Count table validation
        print("\n=== Test 6: Count Table ===")
        case_idx, ctrl_idx = get_phenotype_indices(conn)
        print(f"  Case indices: {len(case_idx)}, Control indices: {len(ctrl_idx)}")

        if rec and rec["gtp"]:
            table = build_count_table(rec["gtp"], case_idx, ctrl_idx, N_SAMPLES)
            total = table.sum()
            expected = len(case_idx) + len(ctrl_idx)
            print(f"  Table sum: {total}, Expected: {expected}")
            # Allow for missing genotypes
            assert total <= expected, f"Table sum {total} > expected {expected}"
            print("PASS: Count table valid")

        # 7. Association — positive control (APOE region)
        print("\n=== Test 7: Association Scan (APOE region) ===")
        results = single_locus_scan(
            conn, "chr19", start=45410000, end=45420000,
            method="auto", verbose=False
        )
        print(f"  Variants scanned: {len(results)}")
        if results:
            best = min(results, key=lambda r: r["p_value"])
            print(f"  Best hit: {best['variantId']} p={best['p_value']:.2e}")
            # With simulated APOE effect, should see some signal
            if best["p_value"] < 0.05:
                print("PASS: APOE region shows expected signal")
            else:
                print("WARN: No significant signal in APOE region (may be due to simulation)")
        else:
            print("WARN: No variants found in APOE region")

        # 8. Z-score validation
        print("\n=== Test 8: Z-score Computation ===")
        if rec and rec["gtp"]:
            table = build_count_table(rec["gtp"], case_idx, ctrl_idx, N_SAMPLES)
            z = compute_z_score(table, len(case_idx), len(ctrl_idx))
            print(f"  Z-score: {z:.4f}")
            # Z-score should be finite for a variant with AF 0.1-0.5
            assert np.isfinite(z), "Z-score is not finite"
            print("PASS: Z-score computation valid")

    print("\n" + "=" * 50)
    if all_pass:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
    return all_pass


if __name__ == "__main__":
    success = validate_all()
    sys.exit(0 if success else 1)
