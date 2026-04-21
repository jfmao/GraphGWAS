"""Import PLINK2/regenie GWAS summary statistics into Neo4j.

Bridges the hybrid architecture: external tools run fast GWAS on BGEN/PGEN
genotypes, producing summary statistics, and this module persists those results
as graph-queryable AssociationResult nodes linked to Variant/GWASStudy.

Usage:
    from graphgwas.summary_import import import_plink2, import_regenie
    import_plink2("height.glm.linear", run_id="ukb_height_2026",
                  phenotype_key="height", method="plink2_linear")

Or CLI:
    python -m graphgwas.summary_import plink2 --file height.glm.linear \\
      --run-id ukb_height_2026 --phenotype height
"""
from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable

from . import config
from .db import GraphGWASConnection
from .results import create_gwas_study, link_study, store_results


def _open_maybe_gz(path: Path):
    if str(path).endswith(".gz"):
        import gzip
        return gzip.open(path, "rt")
    return open(path, "r")


def _build_variant_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Produce chr:pos:ref:alt variantId matching the GraphMana import convention."""
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    return f"{chrom}:{pos}:{ref}:{alt}"


def _safe_float(x: str) -> float | None:
    if x in ("", "NA", "nan", "NaN"):
        return None
    try:
        return float(x)
    except ValueError:
        return None


def _safe_int(x: str) -> int | None:
    v = _safe_float(x)
    return int(v) if v is not None else None


def parse_plink2_glm(path: Path, p_threshold: float = 1.0) -> Iterable[dict]:
    """Yield normalized association result rows from a PLINK2 .glm.{linear,logistic} file.

    Filters rows by P < p_threshold. Distinguishes logistic vs linear by header columns
    (logistic reports OR, linear reports BETA).
    """
    with _open_maybe_gz(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        # PLINK2 sometimes prefixes header with # — normalize
        fieldnames = [fn.lstrip("#") for fn in (reader.fieldnames or [])]
        reader.fieldnames = fieldnames
        is_logistic = "OR" in fieldnames
        for row in reader:
            test = row.get("TEST", "ADD")
            if test != "ADD":
                continue
            p = _safe_float(row.get("P", ""))
            if p is None or p >= p_threshold:
                continue
            chrom = row["CHROM"]
            pos = _safe_int(row["POS"])
            ref = row.get("REF", "")
            alt = row.get("ALT", "")
            a1 = row.get("A1", alt)
            # PLINK2 may swap REF/ALT after BGEN round-trip (see PROVISIONAL_REF?).
            # Prefer the original ID column if it looks like chr:pos:ref:alt.
            orig_id = row.get("ID", "")
            if orig_id and orig_id.count(":") >= 3:
                variant_id_str = orig_id if orig_id.startswith("chr") else f"chr{orig_id}"
            else:
                variant_id_str = _build_variant_id(chrom, pos, ref, alt)
            # PLINK2 defaults A1 to ALT. BETA refers to A1 effect.
            flip = a1 != alt
            if is_logistic:
                or_val = _safe_float(row.get("OR", ""))
                beta = math.log(or_val) if (or_val and or_val > 0) else None
            else:
                or_val = None
                beta = _safe_float(row.get("BETA", ""))
            if beta is not None and flip:
                beta = -beta
                or_val = None if or_val is None else 1.0 / or_val
            se = _safe_float(row.get("SE", row.get("LOG(OR)_SE", "")))
            n = _safe_int(row.get("OBS_CT", "0")) or 0
            log10p = -math.log10(p) if p > 0 else 999.0
            yield dict(
                variantId=variant_id_str,
                method="plink2_logistic" if is_logistic else "plink2_linear",
                beta=beta,
                se=se,
                p_value=p,
                p_value_log10=log10p,
                odds_ratio=or_val,
                n_cases=n if is_logistic else 0,
                n_controls=0,
                af=_safe_float(row.get("A1_FREQ", "")),
            )


def parse_regenie(path: Path, p_threshold: float = 1.0) -> Iterable[dict]:
    """Yield normalized association result rows from a regenie .regenie file.

    regenie logs LOG10P (= -log10(P)). ALLELE1 is the effect allele.
    """
    with _open_maybe_gz(path) as fh:
        # regenie uses space-delimited with header starting on first non-# line
        header = None
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            header = line.rstrip("\n").split()
            break
        if not header:
            return
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split()
            chrom = parts[idx["CHROM"]]
            pos = int(parts[idx["GENPOS"]])
            a0 = parts[idx["ALLELE0"]]
            a1 = parts[idx["ALLELE1"]]
            log10p = _safe_float(parts[idx["LOG10P"]])
            if log10p is None:
                continue
            p = 10.0 ** (-log10p)
            if p >= p_threshold:
                continue
            yield dict(
                variantId=_build_variant_id(chrom, pos, a0, a1),
                method="regenie",
                beta=_safe_float(parts[idx["BETA"]]),
                se=_safe_float(parts[idx["SE"]]),
                p_value=p,
                p_value_log10=log10p,
                odds_ratio=None,
                n_cases=_safe_int(parts[idx["N"]]) or 0,
                n_controls=0,
                af=_safe_float(parts[idx.get("A1FREQ", -1)]) if "A1FREQ" in idx else None,
            )


def import_plink2(
    path: Path,
    run_id: str,
    phenotype_key: str,
    p_threshold: float = 1.0,
    batch_size: int = 5000,
    notes: str = "",
) -> int:
    """Import a PLINK2 .glm output into Neo4j as AssociationResult + GWASStudy."""
    results = list(parse_plink2_glm(path, p_threshold=p_threshold))
    if not results:
        return 0
    method = results[0]["method"]
    with GraphGWASConnection(config.NEO4J_URI, config.NEO4J_USER, config.NEO4J_PASSWORD) as conn:
        create_gwas_study(
            conn, run_id=run_id, phenotype_key=phenotype_key, method=method,
            n_cases=0, n_controls=0, n_variants_tested=len(results),
            region="genome-wide", notes=notes or f"imported from {path.name}",
        )
        n = store_results(conn, results, run_id=run_id,
                          phenotype_key=phenotype_key, batch_size=batch_size)
        link_study(conn, run_id=run_id)
    return n


def import_regenie(
    path: Path,
    run_id: str,
    phenotype_key: str,
    p_threshold: float = 1.0,
    batch_size: int = 5000,
    notes: str = "",
) -> int:
    """Import a regenie .regenie output into Neo4j."""
    results = list(parse_regenie(path, p_threshold=p_threshold))
    if not results:
        return 0
    with GraphGWASConnection(config.NEO4J_URI, config.NEO4J_USER, config.NEO4J_PASSWORD) as conn:
        create_gwas_study(
            conn, run_id=run_id, phenotype_key=phenotype_key, method="regenie",
            n_cases=0, n_controls=0, n_variants_tested=len(results),
            region="genome-wide", notes=notes or f"imported from {path.name}",
        )
        n = store_results(conn, results, run_id=run_id,
                          phenotype_key=phenotype_key, batch_size=batch_size)
        link_study(conn, run_id=run_id)
    return n


def _cli() -> None:
    import argparse
    ap = argparse.ArgumentParser(description="Import PLINK2/regenie summary stats to Neo4j")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p2 = sub.add_parser("plink2", help="Import PLINK2 .glm output")
    p2.add_argument("--file", required=True, type=Path)
    p2.add_argument("--run-id", required=True)
    p2.add_argument("--phenotype", required=True)
    p2.add_argument("--p-threshold", type=float, default=1.0,
                    help="Only keep results with P < threshold (default 1.0, i.e. all)")
    p2.add_argument("--batch-size", type=int, default=5000)
    p2.add_argument("--notes", default="")

    rg = sub.add_parser("regenie", help="Import regenie output")
    rg.add_argument("--file", required=True, type=Path)
    rg.add_argument("--run-id", required=True)
    rg.add_argument("--phenotype", required=True)
    rg.add_argument("--p-threshold", type=float, default=1.0)
    rg.add_argument("--batch-size", type=int, default=5000)
    rg.add_argument("--notes", default="")

    args = ap.parse_args()
    if args.cmd == "plink2":
        n = import_plink2(args.file, args.run_id, args.phenotype,
                          p_threshold=args.p_threshold, batch_size=args.batch_size,
                          notes=args.notes)
        print(f"Imported {n:,} PLINK2 results (run_id={args.run_id})")
    elif args.cmd == "regenie":
        n = import_regenie(args.file, args.run_id, args.phenotype,
                           p_threshold=args.p_threshold, batch_size=args.batch_size,
                           notes=args.notes)
        print(f"Imported {n:,} regenie results (run_id={args.run_id})")


if __name__ == "__main__":
    _cli()
