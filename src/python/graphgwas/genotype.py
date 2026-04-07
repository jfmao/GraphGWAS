"""Genotype access layer — decode gt_packed for GWAS computation.

Provides 7 output formats mapped to every statistical method:
- decode_genotypes: raw int8 {0,1,2,3}
- build_dosage: float64 {0,1,2,NaN}
- build_count_table: int[2][3] case/ctrl × HomRef/Het/HomAlt
- build_carrier_set: bool[N]
- compute_z_score: standardized z-score
- compute_ld_r2: pairwise LD r²
- variant_iterator: stream variants from Neo4j
"""

from __future__ import annotations

import numpy as np

from .config import GT_HOM_REF, GT_HET, GT_HOM_ALT, GT_MISSING
from .db import GraphGWASConnection


# ---------------------------------------------------------------------------
# Core decode — replicates GraphMana's genotype_packer.unpack_genotypes()
# ---------------------------------------------------------------------------

def unpack_genotypes(gt_packed: bytes, n_samples: int) -> np.ndarray:
    """Unpack gt_packed bytes to per-sample genotype codes.

    Returns int8 array: 0=HomRef, 1=Het, 2=HomAlt, 3=Missing.
    """
    arr = np.frombuffer(gt_packed, dtype=np.uint8)
    bits = np.empty(len(arr) * 4, dtype=np.uint8)
    bits[0::4] = arr & 0x03
    bits[1::4] = (arr >> 2) & 0x03
    bits[2::4] = (arr >> 4) & 0x03
    bits[3::4] = (arr >> 6) & 0x03
    return bits[:n_samples].astype(np.int8)


# ---------------------------------------------------------------------------
# Phenotype-stratified sample indices (built ONCE per GWAS session)
# ---------------------------------------------------------------------------

def get_phenotype_indices(conn: GraphGWASConnection) -> tuple[np.ndarray, np.ndarray]:
    """Return (case_packed_indices, control_packed_indices) for active phenotype.

    These arrays map into the gt_packed byte array via PackedGenotypeReader logic.
    Built once, reused for every variant in the scan.
    """
    result = conn.execute_read(
        """
        MATCH (s:Sample)
        WHERE s.is_case = true AND (s.is_excluded_qc IS NULL OR s.is_excluded_qc = false)
        RETURN s.packed_index AS idx
        ORDER BY s.packed_index
        """
    )
    case_idx = np.array([r["idx"] for r in result], dtype=np.int32)

    result = conn.execute_read(
        """
        MATCH (s:Sample)
        WHERE s.is_control = true AND (s.is_excluded_qc IS NULL OR s.is_excluded_qc = false)
        RETURN s.packed_index AS idx
        ORDER BY s.packed_index
        """
    )
    ctrl_idx = np.array([r["idx"] for r in result], dtype=np.int32)

    return case_idx, ctrl_idx


def get_all_indices(conn: GraphGWASConnection) -> np.ndarray:
    """Return packed_indices for all non-excluded samples."""
    result = conn.execute_read(
        """
        MATCH (s:Sample)
        WHERE s.is_excluded_qc IS NULL OR s.is_excluded_qc = false
        RETURN s.packed_index AS idx
        ORDER BY s.packed_index
        """
    )
    return np.array([r["idx"] for r in result], dtype=np.int32)


# ---------------------------------------------------------------------------
# Output format 1: Raw genotype array (for GNN export, full analysis)
# ---------------------------------------------------------------------------

def decode_genotypes(gt_packed: bytes, indices: np.ndarray, n_samples: int) -> np.ndarray:
    """Decode gt_packed for specific sample indices.

    Args:
        gt_packed: packed genotype bytes from Variant node.
        indices: array of packed_index values to extract.
        n_samples: total samples encoded in gt_packed.

    Returns:
        int8 array of length len(indices): {0,1,2,3}.
    """
    all_gt = unpack_genotypes(gt_packed, n_samples)
    return all_gt[indices]


# ---------------------------------------------------------------------------
# Output format 2: Dosage vector (for regression, burden, GRM)
# ---------------------------------------------------------------------------

def build_dosage(gt_packed: bytes, indices: np.ndarray, n_samples: int) -> np.ndarray:
    """Build float64 dosage vector: 0.0, 1.0, 2.0, NaN (missing).

    This is the primary input to logistic/linear regression.
    """
    gt = decode_genotypes(gt_packed, indices, n_samples).astype(np.float64)
    gt[gt == GT_MISSING] = np.nan
    return gt


# ---------------------------------------------------------------------------
# Output format 3: 2×3 count table (for chi2, Fisher, HWE)
# ---------------------------------------------------------------------------

def build_count_table(gt_packed: bytes, case_idx: np.ndarray,
                      ctrl_idx: np.ndarray, n_samples: int) -> np.ndarray:
    """Build 2×3 genotype count table: case/ctrl × HomRef/Het/HomAlt.

    Returns:
        int array shape (2, 3): [[case_hom_ref, case_het, case_hom_alt],
                                  [ctrl_hom_ref, ctrl_het, ctrl_hom_alt]]
    """
    all_gt = unpack_genotypes(gt_packed, n_samples)
    case_gt = all_gt[case_idx]
    ctrl_gt = all_gt[ctrl_idx]

    table = np.zeros((2, 3), dtype=np.int64)
    # Vectorized with np.bincount (C-level, ~5-10x faster than Python loop)
    valid_case = case_gt[case_gt < 3].astype(np.intp)
    valid_ctrl = ctrl_gt[ctrl_gt < 3].astype(np.intp)
    if len(valid_case) > 0:
        table[0] = np.bincount(valid_case, minlength=3)[:3]
    if len(valid_ctrl) > 0:
        table[1] = np.bincount(valid_ctrl, minlength=3)[:3]
    return table


def count_table_to_allelic(table: np.ndarray) -> np.ndarray:
    """Convert 2×3 genotype table to 2×2 allelic table.

    Returns:
        int array shape (2, 2): [[case_ref, case_alt],
                                  [ctrl_ref, ctrl_alt]]
    """
    # REF alleles = 2*HomRef + Het, ALT alleles = 2*HomAlt + Het
    allelic = np.zeros((2, 2), dtype=np.int64)
    allelic[0, 0] = 2 * table[0, 0] + table[0, 1]  # case REF
    allelic[0, 1] = 2 * table[0, 2] + table[0, 1]  # case ALT
    allelic[1, 0] = 2 * table[1, 0] + table[1, 1]  # ctrl REF
    allelic[1, 1] = 2 * table[1, 2] + table[1, 1]  # ctrl ALT
    return allelic


# ---------------------------------------------------------------------------
# Output format 4: Carrier set (for epistasis, co-occurrence)
# ---------------------------------------------------------------------------

def build_carrier_set(gt_packed: bytes, indices: np.ndarray,
                      n_samples: int) -> np.ndarray:
    """Bool array: True where genotype is HET or HOM_ALT (carrier)."""
    gt = decode_genotypes(gt_packed, indices, n_samples)
    return (gt == GT_HET) | (gt == GT_HOM_ALT)


# ---------------------------------------------------------------------------
# Output format 5: Standardized z-score (for MPAT, Manhattan)
# ---------------------------------------------------------------------------

def compute_z_score(table: np.ndarray, n_case: int, n_ctrl: int) -> float:
    """Compute standardized z-score from 2×3 count table.

    z = (p_case - p_ctrl) / SE, where SE = sqrt(p_hat*(1-p_hat)*(1/AN_case + 1/AN_ctrl))
    """
    # Allele counts from genotype table
    ac_case = table[0, 1] + 2 * table[0, 2]
    an_case = 2 * (table[0, 0] + table[0, 1] + table[0, 2])
    ac_ctrl = table[1, 1] + 2 * table[1, 2]
    an_ctrl = 2 * (table[1, 0] + table[1, 1] + table[1, 2])

    if an_case == 0 or an_ctrl == 0:
        return 0.0

    p_case = ac_case / an_case
    p_ctrl = ac_ctrl / an_ctrl
    p_hat = (ac_case + ac_ctrl) / (an_case + an_ctrl)

    if p_hat == 0.0 or p_hat == 1.0:
        return 0.0

    se = np.sqrt(p_hat * (1 - p_hat) * (1 / an_case + 1 / an_ctrl))
    if se == 0:
        return 0.0

    return (p_case - p_ctrl) / se


# ---------------------------------------------------------------------------
# Output format 6: Pairwise LD r² (for MPAT null, gcSKAT)
# ---------------------------------------------------------------------------

def compute_ld_r2(gt_packed_v1: bytes, gt_packed_v2: bytes,
                  indices: np.ndarray, n_samples: int) -> float:
    """Compute pairwise LD r² between two variants.

    Uses Pearson correlation of dosages (0,1,2) for non-missing samples.
    """
    d1 = build_dosage(gt_packed_v1, indices, n_samples)
    d2 = build_dosage(gt_packed_v2, indices, n_samples)

    # Remove samples where either is missing
    valid = ~(np.isnan(d1) | np.isnan(d2))
    if np.sum(valid) < 10:
        return 0.0

    d1 = d1[valid]
    d2 = d2[valid]

    corr = np.corrcoef(d1, d2)[0, 1]
    if np.isnan(corr):
        return 0.0
    return corr ** 2


# ---------------------------------------------------------------------------
# Output format 7: Variant iterator (stream from Neo4j)
# ---------------------------------------------------------------------------

def variant_iterator(conn: GraphGWASConnection, chr: str,
                     start: int | None = None, end: int | None = None,
                     batch_size: int = 10000):
    """Yield variant dicts from Neo4j for a genomic region.

    Each dict contains: variantId, chr, pos, ref, alt, af_total, call_rate, gt_packed.

    Uses cursor-based pagination (WHERE pos > last_pos) instead of SKIP/LIMIT.
    SKIP/LIMIT degrades to O(N*offset) on large chromosomes because Neo4j must
    scan past all skipped rows. Cursor-based pagination is O(batch_size) per query
    regardless of position in the chromosome.
    """
    fields = ("v.variantId AS variantId, v.chr AS chr, v.pos AS pos, "
              "v.ref AS ref, v.alt AS alt, v.af_total AS af_total, "
              "v.call_rate AS call_rate, v.gt_packed AS gt_packed")

    cursor_pos = (start - 1) if start is not None else -1
    end_pos = end

    while True:
        if end_pos is not None:
            query = f"""
                MATCH (v:Variant)
                WHERE v.chr = $chr AND v.pos > $cursor AND v.pos <= $end
                RETURN {fields}
                ORDER BY v.pos
                LIMIT $limit
            """
            params = {"chr": chr, "cursor": cursor_pos, "end": end_pos, "limit": batch_size}
        else:
            query = f"""
                MATCH (v:Variant)
                WHERE v.chr = $chr AND v.pos > $cursor
                RETURN {fields}
                ORDER BY v.pos
                LIMIT $limit
            """
            params = {"chr": chr, "cursor": cursor_pos, "limit": batch_size}

        result = conn.execute_read(query, params)
        records = list(result)
        if not records:
            break
        for rec in records:
            yield dict(rec)
        # Advance cursor to the last position in this batch
        cursor_pos = records[-1]["pos"]
        if len(records) < batch_size:
            break


# ---------------------------------------------------------------------------
# PCA import and covariate matrix (Phase 2)
# ---------------------------------------------------------------------------

def import_pca(conn: GraphGWASConnection, eigenvec_path: str,
               n_components: int = 10) -> int:
    """Import PCA coordinates from PLINK2 eigenvec file into Sample nodes.

    Sets properties pca_1, pca_2, ..., pca_k on Sample nodes.
    Returns number of samples updated.
    """
    import csv

    rows = []
    with open(eigenvec_path, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)  # #FID IID PC1 PC2 ...
        # Find PC columns
        pc_cols = [i for i, h in enumerate(header) if h.startswith("PC")]
        id_col = 1  # IID column in PLINK2 eigenvec

        for row in reader:
            sample_id = row[id_col]
            pcs = {}
            for j, col_idx in enumerate(pc_cols[:n_components]):
                pcs[f"pca_{j + 1}"] = float(row[col_idx])
            rows.append({"sample_id": sample_id, **pcs})

    # Build dynamic SET clause
    set_parts = [f"s.pca_{i+1} = row.pca_{i+1}" for i in range(min(n_components, len(pc_cols)))]
    set_clause = ", ".join(set_parts)

    n_updated = 0
    batch_size = 500
    for i in range(0, len(rows), batch_size):
        batch = rows[i: i + batch_size]
        result = conn.execute_write(
            f"""
            UNWIND $batch AS row
            MATCH (s:Sample {{sampleId: row.sample_id}})
            SET {set_clause}
            RETURN count(s) AS n
            """,
            {"batch": batch},
        )
        rec = result.single()
        n_updated += rec["n"] if rec else 0

    return n_updated


def get_covariate_matrix(conn: GraphGWASConnection,
                         covariate_names: list[str],
                         sample_indices: np.ndarray) -> np.ndarray:
    """Build covariate matrix aligned to sample index order.

    Args:
        conn: database connection.
        covariate_names: property names on Sample nodes (e.g., ['pca_1', 'pca_2', 'sex']).
        sample_indices: packed_index array (case_idx + ctrl_idx order).

    Returns:
        (N, k) float64 array. NaN for missing values.
    """
    # Build RETURN clause for requested covariates
    return_parts = ["s.packed_index AS idx"]
    for name in covariate_names:
        safe = name.replace(" ", "_").replace("-", "_")
        return_parts.append(f"s.{safe} AS {safe}")
    return_clause = ", ".join(return_parts)

    result = conn.execute_read(
        f"""
        MATCH (s:Sample)
        WHERE s.packed_index IN $indices
        RETURN {return_clause}
        """,
        {"indices": [int(i) for i in sample_indices]},
    )

    # Build lookup: packed_index → covariate values
    cov_lookup = {}
    for rec in result:
        idx = rec["idx"]
        vals = []
        for name in covariate_names:
            safe = name.replace(" ", "_").replace("-", "_")
            v = rec[safe]
            vals.append(float(v) if v is not None else np.nan)
        cov_lookup[idx] = vals

    # Align to sample_indices order
    matrix = np.full((len(sample_indices), len(covariate_names)), np.nan)
    for i, idx in enumerate(sample_indices):
        if int(idx) in cov_lookup:
            matrix[i] = cov_lookup[int(idx)]

    return matrix


def get_phenotype_values(conn: GraphGWASConnection,
                         sample_indices: np.ndarray) -> np.ndarray:
    """Get gwas_value for samples, aligned to index order.

    For binary traits: 1.0 (case), 0.0 (control).
    For quantitative traits: continuous value.
    """
    result = conn.execute_read(
        """
        MATCH (s:Sample)
        WHERE s.packed_index IN $indices
        RETURN s.packed_index AS idx, s.gwas_value AS val
        """,
        {"indices": [int(i) for i in sample_indices]},
    )

    val_lookup = {}
    for rec in result:
        v = rec["val"]
        val_lookup[rec["idx"]] = float(v) if v is not None else np.nan

    values = np.full(len(sample_indices), np.nan)
    for i, idx in enumerate(sample_indices):
        if int(idx) in val_lookup:
            values[i] = val_lookup[int(idx)]

    return values
