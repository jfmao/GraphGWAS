"""Phenotype import and activation for GWAS analysis.

Loads phenotype data from CSV into Sample.phenotypes MAP property,
then activates a specific trait by setting indexed is_case/is_control flags.
"""

import csv
import sys

import click

from .db import GraphGWASConnection


def load_phenotypes(conn: GraphGWASConnection, csv_path: str,
                    sample_id_col: str, trait_cols: list[str] | None = None) -> dict:
    """Load phenotype CSV into Sample.phenotypes MAP.

    Args:
        conn: active database connection.
        csv_path: path to phenotype CSV file.
        sample_id_col: column name for sample ID.
        trait_cols: columns to load (None = all except sample_id_col).

    Returns:
        dict with n_matched, n_unmatched, n_total, trait_summary.
    """
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        if sample_id_col not in reader.fieldnames:
            raise ValueError(f"Column '{sample_id_col}' not found. Available: {reader.fieldnames}")

        if trait_cols is None:
            trait_cols = [c for c in reader.fieldnames if c != sample_id_col]

        rows = []
        for row in reader:
            phenotypes = {}
            for col in trait_cols:
                val = row.get(col, "")
                if val == "":
                    continue
                # Try numeric conversion
                try:
                    phenotypes[col] = float(val)
                    if phenotypes[col] == int(phenotypes[col]):
                        phenotypes[col] = int(phenotypes[col])
                except ValueError:
                    phenotypes[col] = val
            rows.append({"sample_id": row[sample_id_col], "phenotypes": phenotypes})

    # Batch update — store each trait as a separate property (pheno_<trait>)
    # Neo4j property values cannot be nested MAPs, so we prefix trait names.
    n_matched = 0
    n_unmatched = 0
    batch_size = 500

    # Build SET clause dynamically for each trait
    set_parts = []
    for col in trait_cols:
        safe_col = col.replace(" ", "_").replace("-", "_")
        set_parts.append(f"s.pheno_{safe_col} = row.pheno_{safe_col}")
    set_clause = ", ".join(set_parts)

    for i in range(0, len(rows), batch_size):
        batch = []
        for r in rows[i : i + batch_size]:
            entry = {"sample_id": r["sample_id"]}
            for col in trait_cols:
                safe_col = col.replace(" ", "_").replace("-", "_")
                entry[f"pheno_{safe_col}"] = r["phenotypes"].get(col)
            batch.append(entry)

        result = conn.execute_write(
            f"""
            UNWIND $batch AS row
            OPTIONAL MATCH (s:Sample {{sampleId: row.sample_id}})
            WITH row, s
            WHERE s IS NOT NULL
            SET {set_clause}, s.phenotypes_loaded = true
            RETURN count(s) AS matched
            """,
            {"batch": batch},
        )
        rec = result.single()
        matched = rec["matched"] if rec else 0
        n_matched += matched
        n_unmatched += len(batch) - matched

    # Build trait summary
    trait_summary = {}
    for col in trait_cols:
        values = [r["phenotypes"].get(col) for r in rows if col in r["phenotypes"]]
        numeric = [v for v in values if isinstance(v, (int, float))]
        if numeric:
            trait_summary[col] = {
                "type": "numeric",
                "n": len(numeric),
                "min": min(numeric),
                "max": max(numeric),
                "mean": sum(numeric) / len(numeric),
            }
        else:
            counts = {}
            for v in values:
                counts[str(v)] = counts.get(str(v), 0) + 1
            trait_summary[col] = {"type": "categorical", "n": len(values), "levels": counts}

    return {
        "n_total": len(rows),
        "n_matched": n_matched,
        "n_unmatched": n_unmatched,
        "traits_loaded": trait_cols,
        "trait_summary": trait_summary,
    }


def activate_phenotype(conn: GraphGWASConnection, trait: str,
                       case_value: str | None = None,
                       control_value: str | None = None) -> dict:
    """Activate a trait for GWAS by setting is_case/is_control/gwas_value.

    For binary traits: provide case_value and control_value.
    For quantitative traits: omit case_value/control_value (gwas_value = numeric phenotype).

    Returns:
        dict with n_cases, n_controls, n_missing, trait, trait_type.
    """
    is_binary = case_value is not None

    # Trait is stored as property pheno_<trait> on Sample nodes
    safe_trait = trait.replace(" ", "_").replace("-", "_")
    prop_name = f"pheno_{safe_trait}"

    if is_binary:
        # Binary trait: set is_case/is_control from pheno_<trait> property
        result = conn.execute_write(
            f"""
            MATCH (s:Sample)
            WHERE s.{prop_name} IS NOT NULL
            WITH s, toString(s.{prop_name}) AS val
            SET s.gwas_phenotype = $trait,
                s.is_case = CASE WHEN val = $case_val THEN true ELSE false END,
                s.is_control = CASE WHEN val = $ctrl_val THEN true ELSE false END,
                s.gwas_value = CASE
                    WHEN val = $case_val THEN 1.0
                    WHEN val = $ctrl_val THEN 0.0
                    ELSE null
                END
            RETURN
                sum(CASE WHEN s.is_case = true THEN 1 ELSE 0 END) AS n_cases,
                sum(CASE WHEN s.is_control = true THEN 1 ELSE 0 END) AS n_controls,
                sum(CASE WHEN s.gwas_value IS NULL THEN 1 ELSE 0 END) AS n_missing
            """,
            {"trait": trait, "case_val": case_value, "ctrl_val": control_value},
        )
        rec = result.single()
        return {
            "trait": trait,
            "trait_type": "binary",
            "n_cases": rec["n_cases"],
            "n_controls": rec["n_controls"],
            "n_missing": rec["n_missing"],
        }
    else:
        # Quantitative trait: gwas_value = numeric, no case/control
        result = conn.execute_write(
            f"""
            MATCH (s:Sample)
            WHERE s.{prop_name} IS NOT NULL
            WITH s, s.{prop_name} AS val
            SET s.gwas_phenotype = $trait,
                s.is_case = null,
                s.is_control = null,
                s.gwas_value = CASE WHEN val IS NOT NULL THEN toFloat(val) ELSE null END
            RETURN
                count(CASE WHEN s.gwas_value IS NOT NULL THEN 1 END) AS n_with_value,
                count(CASE WHEN s.gwas_value IS NULL THEN 1 END) AS n_missing
            """,
            {"trait": trait},
        )
        rec = result.single()
        return {
            "trait": trait,
            "trait_type": "quantitative",
            "n_with_value": rec["n_with_value"],
            "n_missing": rec["n_missing"],
        }


def clear_phenotype(conn: GraphGWASConnection):
    """Clear all phenotype activation flags (not the phenotypes MAP itself)."""
    conn.execute_write(
        """
        MATCH (s:Sample)
        REMOVE s.is_case, s.is_control, s.gwas_phenotype, s.gwas_value
        """
    )
