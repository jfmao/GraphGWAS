"""Pre-association QC with stored metrics.

Runs a full-genome scan computing:
- HWE p-value in controls (stored on Variant nodes)
- MAC in cases/controls (stored on Variant nodes)
- Per-sample missingness and heterozygosity (stored on Sample nodes)
"""

from __future__ import annotations


import numpy as np

from .config import (
    HWE_P_THRESHOLD, MAC_MIN,
    SAMPLE_MISSINGNESS_MAX, HET_RATE_SD, N_SAMPLES,
)
from .db import GraphGWASConnection
from .genotype import unpack_genotypes, get_phenotype_indices


def hwe_exact_test(n_het: int, n_hom_minor: int, n_total: int) -> float:
    """Wigginton et al. (2005) exact test for Hardy-Weinberg equilibrium.

    Returns mid-p value. Small p-value = departure from HWE.
    """
    if n_total == 0:
        return 1.0

    n_hom_major = n_total - n_het - n_hom_minor
    if n_hom_major < 0:
        return 1.0

    # Rare homozygote count
    n_min = min(n_hom_minor, n_hom_major)
    n_max = max(n_hom_minor, n_hom_major)
    n_het_obs = n_het

    # Total allele count
    n_alleles = 2 * n_total
    n_minor_alleles = 2 * n_min + n_het

    # Compute probability distribution of het count under HWE
    # Using the recursive method from Wigginton et al.
    het_probs = np.zeros(n_minor_alleles + 1)

    # Start from mid value
    mid = int(n_minor_alleles * (n_alleles - n_minor_alleles) / n_alleles)
    if mid % 2 != n_minor_alleles % 2:
        mid += 1
    mid = min(mid, n_minor_alleles)

    het_probs[mid] = 1.0
    prob_sum = 1.0

    # Recurse downward
    curr_het = mid
    while curr_het > 1:
        curr_hom_minor = (n_minor_alleles - curr_het) // 2
        curr_hom_major = n_total - curr_het - curr_hom_minor
        het_probs[curr_het - 2] = (
            het_probs[curr_het]
            * curr_het
            * (curr_het - 1)
            / (4.0 * (curr_hom_minor + 1) * (curr_hom_major + 1))
        )
        prob_sum += het_probs[curr_het - 2]
        curr_het -= 2

    # Recurse upward
    curr_het = mid
    while curr_het < n_minor_alleles - 1:
        curr_hom_minor = (n_minor_alleles - curr_het) // 2
        curr_hom_major = n_total - curr_het - curr_hom_minor
        het_probs[curr_het + 2] = (
            het_probs[curr_het]
            * 4.0
            * curr_hom_minor
            * curr_hom_major
            / ((curr_het + 1) * (curr_het + 2))
        )
        prob_sum += het_probs[curr_het + 2]
        curr_het += 2

    # Normalize
    het_probs /= prob_sum

    # Mid-p value: sum of probs for het counts <= observed
    p_value = 0.0
    for i in range(0, n_minor_alleles + 1, 2 if n_minor_alleles % 2 == 0 else 1):
        if i % 2 == n_het_obs % 2:  # same parity
            if het_probs[i] <= het_probs[n_het_obs] + 1e-12:
                p_value += het_probs[i]

    return min(p_value, 1.0)


def run_variant_qc(conn: GraphGWASConnection, chromosomes: list[str] | None = None,
                   batch_size: int = 5000, verbose: bool = True) -> dict:
    """Run variant QC: HWE in controls, MAC in case/ctrl.

    Stores hwe_pvalue_controls, mac_cases, mac_controls on Variant nodes.
    Also accumulates per-sample missingness and het counts.

    Returns summary dict.
    """
    case_idx, ctrl_idx = get_phenotype_indices(conn)
    n_case = len(case_idx)
    n_ctrl = len(ctrl_idx)
    all_idx = np.concatenate([case_idx, ctrl_idx])
    n_all = len(all_idx)

    if verbose:
        print(f"QC: {n_case} cases, {n_ctrl} controls, {n_all} total samples")

    # Per-sample accumulators
    sample_missing = np.zeros(N_SAMPLES, dtype=np.int64)
    sample_het = np.zeros(N_SAMPLES, dtype=np.int64)
    sample_total = np.zeros(N_SAMPLES, dtype=np.int64)

    # Get chromosome list
    if chromosomes is None:
        result = conn.execute_read(
            "MATCH (v:Variant) RETURN DISTINCT v.chr AS chr ORDER BY chr"
        )
        chromosomes = [r["chr"] for r in result]

    n_variants_total = 0
    n_hwe_fail = 0
    n_mac_fail = 0

    for chrom in chromosomes:
        if verbose:
            print(f"  Scanning {chrom}...", end="", flush=True)

        # Use batched iteration to avoid loading entire chromosome
        qc_batch_size = 5000
        qc_offset = 0
        updates = []
        chrom_count = 0
        more_data = True

        while more_data:
            result = conn.execute_read(
                """
                MATCH (v:Variant)
                WHERE v.chr = $chr
                RETURN v.variantId AS vid, v.gt_packed AS gtp, v.call_rate AS cr
                ORDER BY v.pos
                SKIP $skip LIMIT $limit
                """,
                {"chr": chrom, "skip": qc_offset, "limit": qc_batch_size},
            )
            records = list(result)
            if not records:
                break
            if len(records) < qc_batch_size:
                more_data = False
            qc_offset += len(records)

            for rec in records:
                gt_packed = rec["gtp"]
                if gt_packed is None:
                    continue

                all_gt = unpack_genotypes(gt_packed, N_SAMPLES)
                chrom_count += 1

                # Accumulate per-sample stats (vectorized — no Python loop)
                gt_subset = all_gt[all_idx]
                sample_total[all_idx] += 1
                sample_missing[all_idx] += (gt_subset == 3).astype(np.int64)
                sample_het[all_idx] += (gt_subset == 1).astype(np.int64)

                # Control genotype counts for HWE
                ctrl_gt = all_gt[ctrl_idx]
                ctrl_het = int(np.sum(ctrl_gt == 1))
                ctrl_hom_alt = int(np.sum(ctrl_gt == 2))
                ctrl_non_missing = int(np.sum(ctrl_gt != 3))

                hwe_p = hwe_exact_test(ctrl_het, ctrl_hom_alt, ctrl_non_missing)

                # MAC in cases and controls
                case_gt = all_gt[case_idx]
                case_ac = int(np.sum(case_gt == 1) + 2 * np.sum(case_gt == 2))
                case_an = 2 * int(np.sum(case_gt != 3))
                ctrl_ac = int(np.sum(ctrl_gt == 1) + 2 * np.sum(ctrl_gt == 2))
                ctrl_an = 2 * ctrl_non_missing

                mac_case = min(case_ac, case_an - case_ac) if case_an > 0 else 0
                mac_ctrl = min(ctrl_ac, ctrl_an - ctrl_ac) if ctrl_an > 0 else 0

                if hwe_p < HWE_P_THRESHOLD:
                    n_hwe_fail += 1
                if mac_case < MAC_MIN:
                    n_mac_fail += 1

                updates.append({
                    "vid": rec["vid"],
                    "hwe_p": hwe_p,
                    "mac_case": mac_case,
                    "mac_ctrl": mac_ctrl,
                })

                # Batch write
                if len(updates) >= batch_size:
                    _write_variant_qc_batch(conn, updates)
                    updates = []

        if updates:
            _write_variant_qc_batch(conn, updates)

        n_variants_total += chrom_count
        if verbose:
            print(f" {chrom_count} variants")

    # Write per-sample QC
    _write_sample_qc(conn, all_idx, sample_missing, sample_het, sample_total)

    return {
        "n_variants_scanned": n_variants_total,
        "n_hwe_fail": n_hwe_fail,
        "n_mac_fail": n_mac_fail,
        "n_cases": n_case,
        "n_controls": n_ctrl,
    }


def _write_variant_qc_batch(conn: GraphGWASConnection, updates: list[dict]):
    """Batch-write HWE and MAC to Variant nodes."""
    conn.execute_write(
        """
        UNWIND $batch AS row
        MATCH (v:Variant {variantId: row.vid})
        SET v.hwe_pvalue_controls = row.hwe_p,
            v.mac_cases = row.mac_case,
            v.mac_controls = row.mac_ctrl
        """,
        {"batch": updates},
    )


def _write_sample_qc(conn: GraphGWASConnection, all_idx: np.ndarray,
                      missing: np.ndarray, het: np.ndarray, total: np.ndarray):
    """Write per-sample QC metrics and apply exclusion flags."""
    updates = []
    het_rates = []

    for idx in all_idx:
        t = int(total[idx])
        if t == 0:
            continue
        miss_rate = float(missing[idx]) / t
        het_rate = float(het[idx]) / t
        het_rates.append(het_rate)
        updates.append({
            "idx": int(idx),
            "miss": miss_rate,
            "het": het_rate,
        })

    # Compute het_rate thresholds
    het_arr = np.array(het_rates)
    het_mean = np.mean(het_arr) if len(het_arr) > 0 else 0
    het_std = np.std(het_arr) if len(het_arr) > 0 else 1
    het_lo = het_mean - HET_RATE_SD * het_std
    het_hi = het_mean + HET_RATE_SD * het_std

    # Add exclusion flags
    for u in updates:
        excluded = False
        reason = []
        if u["miss"] > SAMPLE_MISSINGNESS_MAX:
            excluded = True
            reason.append(f"missingness={u['miss']:.3f}")
        if u["het"] < het_lo or u["het"] > het_hi:
            excluded = True
            reason.append(f"het_rate={u['het']:.4f}")
        u["excluded"] = excluded
        u["reason"] = "; ".join(reason) if reason else None

    # Batch write
    for i in range(0, len(updates), 1000):
        batch = updates[i : i + 1000]
        conn.execute_write(
            """
            UNWIND $batch AS row
            MATCH (s:Sample {packed_index: row.idx})
            SET s.genotype_missingness = row.miss,
                s.het_rate = row.het,
                s.is_excluded_qc = row.excluded,
                s.exclusion_reason = row.reason
            """,
            {"batch": batch},
        )
