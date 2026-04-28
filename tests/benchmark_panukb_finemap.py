"""End-to-end Pan-UKB cross-ancestry fine-mapping (sumstats + 1KG LD fallback).

Demonstrates the full sumstats-only pipeline by combining:
  1. Real Pan-UKB summary statistics (streamed via tabix, no Hail required)
  2. Ancestry-matched LD computed locally from 1KG Phase 3 BGEN

This is a *fallback* configuration for environments without Hail. When Hail
is installed and Pan-UKB LD BlockMatrices are accessible, replace the
_compute_1kg_ld_matched() call with graphgwas.panukb.fetch_ld_slice() for
in-sample Pan-UKB LD. The downstream fine-mapping code is identical.

Example usage:
    python tests/benchmark_panukb_finemap.py --locus FTO --out results/panukb/fto
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd
from pyliftover import LiftOver

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src" / "python"))

from graphgwas.bgen_reader import BgenReader  # noqa: E402
from graphgwas.finemapping_v2 import (  # noqa: E402
    hbp_finemap_from_sumstats,
    l1_finemap_from_sumstats,
)
from graphgwas.panukb import (  # noqa: E402
    ANCESTRY_N,
    LocusSumstats,
    fetch_sumstats_locus,
)


# ===================================================================
# Canonical loci for cross-ancestry demonstration (GRCh37)
# ===================================================================

LOCI = {
    "FTO":   {"trait": "BMI",           "phenocode": "21001", "trait_type": "continuous",
              "modifier": "irnt", "chr": "16", "center": 53_820_527, "window": 100_000},
    "APOA5": {"trait": "Triglycerides", "phenocode": "30870", "trait_type": "biomarkers",
              "modifier": "irnt", "chr": "11", "center": 116_662_407, "window": 100_000},
    "LDLR":  {"trait": "LDL",           "phenocode": "30780", "trait_type": "biomarkers",
              "modifier": "irnt", "chr": "19", "center":  11_200_038, "window": 100_000},
    "HMGA2": {"trait": "Height",        "phenocode": "50",    "trait_type": "continuous",
              "modifier": "irnt", "chr": "12", "center":  66_360_071, "window": 100_000},
}

# Map Pan-UKB ancestry → 1KG superpopulation. CSA ↔ SAS is the best approximation:
# Pan-UKB CSA pools Indian/Pakistani/Bangladeshi (UKB field 21000), 1KG SAS pools
# the same population groups (GIH, PJL, ITU, STU, BEB).
ANC_TO_1KG_SUPERPOP = {
    "EUR": "EUR", "CSA": "SAS", "AFR": "AFR", "EAS": "EAS",
}

# Pan-UKB sumstats are GRCh37; our local 1KG BGEN (NYGC 30× high-coverage
# release, used by plink2 --export bgen-1.2) is GRCh38. The LiftOver step
# maps Pan-UKB variants to the BGEN build so sumstats × LD intersect.
_LIFTOVER = None  # lazy singleton


def _lift_sumstats_to_b38(sumstats: LocusSumstats) -> LocusSumstats:
    """Return a copy of sumstats with positions lifted GRCh37 → GRCh38.

    Variants that fail liftover or land on a different chromosome are
    dropped. The variant_id is rebuilt from the new position; allele
    orientation is preserved.
    """
    global _LIFTOVER
    if _LIFTOVER is None:
        _LIFTOVER = LiftOver("hg19", "hg38")

    df = sumstats.variants.copy()
    new_pos = np.full(len(df), -1, dtype=np.int64)
    for i, row in enumerate(df.itertuples(index=False)):
        chr_tag = f"chr{row.chr}"
        res = _LIFTOVER.convert_coordinate(chr_tag, int(row.pos))
        if res and res[0][0] == chr_tag:
            new_pos[i] = int(res[0][1])
    df["pos"] = new_pos
    df = df[df["pos"] >= 0].copy()
    df["variant_id"] = (df["chr"].astype(str) + ":" + df["pos"].astype(str)
                        + ":" + df["ref"] + ":" + df["alt"])
    # End window too — use the max lifted pos as the new end (approximate)
    return LocusSumstats(
        chr=sumstats.chr, start=int(df["pos"].min()) if len(df) else sumstats.start,
        end=int(df["pos"].max()) if len(df) else sumstats.end,
        ancestry=sumstats.ancestry, variants=df.reset_index(drop=True),
        n_samples=sumstats.n_samples,
    )


# ===================================================================
# Ancestry-matched 1KG LD
# ===================================================================

def _compute_1kg_ld_matched(
    chr: str, start: int, end: int,
    ancestry: str,
    bgen_dir: Path,
    popmap_path: Path,
    min_af: float = 0.01,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Compute r² LD matrix from 1KG BGEN, restricted to samples matching
    the requested Pan-UKB ancestry.

    Returns (variant_df, R_sq) where:
      variant_df columns: chr, pos, a1, a2, af, variant_id
      R_sq: n × n squared-correlation matrix (float32)
    """
    superpop = ANC_TO_1KG_SUPERPOP[ancestry]
    popmap = pd.read_csv(popmap_path, sep="\t")
    kg_samples = popmap[popmap.superpopulation == superpop].iloc[:, 0].tolist()

    reader = BgenReader(bgen_dir)
    all_samples = reader.samples(chr)
    sample_mask = np.isin(all_samples, kg_samples)
    if sample_mask.sum() < 30:
        raise RuntimeError(
            f"Too few 1KG {superpop} samples ({sample_mask.sum()}) "
            f"match BGEN samples — check popmap/BGEN alignment.",
        )

    variants, dosage = reader.load_locus(chr, start, end)
    if len(variants) == 0:
        return variants, np.zeros((0, 0), dtype=np.float32)
    D = dosage[sample_mask, :]  # (n_anc_samples, n_var)

    af = D.mean(axis=0) / 2.0
    keep = (af >= min_af) & (af <= 1 - min_af)
    variants = variants.loc[keep].reset_index(drop=True)
    D = D[:, keep]
    if len(variants) == 0:
        return variants, np.zeros((0, 0), dtype=np.float32)

    # Standardize
    mu = D.mean(axis=0, keepdims=True)
    std = D.std(axis=0, keepdims=True)
    std = np.where(std < 1e-10, 1.0, std)
    Z = (D - mu) / std
    # r² = corr²
    R = (Z.T @ Z) / Z.shape[0]
    R_sq = np.clip(R * R, 0, 1).astype(np.float32)
    np.fill_diagonal(R_sq, 1.0)

    variants["af"] = af[keep]
    # 1KG BGEN plink2 convention: a1 = ALT, a2 = REF
    variants["variant_id"] = (variants["chr"].astype(str).str.replace("chr", "")
                              + ":" + variants["pos"].astype(str)
                              + ":" + variants["a2"] + ":" + variants["a1"])
    return variants, R_sq


# ===================================================================
# Sumstats × LD intersection
# ===================================================================

def _intersect_sumstats_ld(
    sumstats: LocusSumstats,
    ld_variants: pd.DataFrame,
    R_sq: np.ndarray,
) -> tuple[list[dict], np.ndarray, np.ndarray]:
    """Align a Pan-UKB LocusSumstats with the 1KG LD variant set.

    Matches by chr:pos:ref:alt; tries the allele swap chr:pos:alt:ref
    (with z sign flip) if direct match fails. Returns:
        variants: list of dicts with variantId/chr/pos/ref/alt/af_total
        z: aligned z-score array
        R_sq_aligned: square LD sub-matrix for the intersected variants
    """
    ss = sumstats.variants.set_index("variant_id")
    out_variants: list[dict] = []
    z_list: list[float] = []
    ld_idx_list: list[int] = []
    for i, row in ld_variants.iterrows():
        vid = row["variant_id"]
        flip = 1.0
        rec = None
        if vid in ss.index:
            rec = ss.loc[vid]
        else:
            parts = vid.split(":", 3)
            if len(parts) != 4:
                continue
            c, p, a, b = parts
            alt_vid = f"{c}:{p}:{b}:{a}"
            if alt_vid in ss.index:
                rec = ss.loc[alt_vid]
                flip = -1.0
        if rec is None:
            continue
        # Handle multiple rows in sumstats for the same variant_id (should be rare)
        if isinstance(rec, pd.DataFrame):
            rec = rec.iloc[0]
        # Skip if z is non-finite (can happen near the window edges)
        if not np.isfinite(float(rec["z"])):
            continue
        out_variants.append({
            "variantId": vid,
            "chr": str(row["chr"]).replace("chr", ""),
            "pos": int(row["pos"]),
            "ref": str(row["a2"]),
            "alt": str(row["a1"]),
            "af_total": float(row["af"]),
        })
        z_list.append(float(rec["z"]) * flip)
        ld_idx_list.append(i)

    if not out_variants:
        return [], np.zeros(0), np.zeros((0, 0))

    idx = np.asarray(ld_idx_list, dtype=int)
    R_sub = R_sq[np.ix_(idx, idx)]
    return out_variants, np.asarray(z_list, dtype=np.float64), R_sub.astype(np.float64)


# ===================================================================
# Runner
# ===================================================================

def run_locus(
    locus_name: str,
    ancestries: list[str] = ["EUR", "CSA", "AFR", "EAS"],
    bgen_dir: Path = Path("/mnt/data/GraphGWAS/tests/data/human/1kGP_bgen"),
    popmap_path: Path = Path("/mnt/data/GraphGWAS/tests/data/human/1kg_popmap.tsv"),
    out_dir: Path | None = None,
    verbose: bool = True,
) -> dict:
    """Run Pan-UKB cross-ancestry fine-mapping for one locus."""
    locus = LOCI[locus_name]
    start = locus["center"] - locus["window"]
    end = locus["center"] + locus["window"]

    # Liftover center to GRCh38 for the 1KG BGEN query window
    global _LIFTOVER
    if _LIFTOVER is None:
        _LIFTOVER = LiftOver("hg19", "hg38")
    b38 = _LIFTOVER.convert_coordinate(f"chr{locus['chr']}", locus["center"])
    if not b38:
        raise RuntimeError(f"liftOver failed for {locus_name} center")
    ld_center_b38 = b38[0][1]
    ld_start_b38 = ld_center_b38 - locus["window"]
    ld_end_b38 = ld_center_b38 + locus["window"]

    # 1) Pan-UKB sumstats (GRCh37) → lifted to GRCh38 to match local BGEN
    t0 = time.time()
    sumstats_b37 = fetch_sumstats_locus(
        phenocode=locus["phenocode"], chr=locus["chr"],
        start=start, end=end,
        ancestries=ancestries, trait_type=locus["trait_type"],
        modifier=locus["modifier"],
    )
    sumstats = {anc: _lift_sumstats_to_b38(s) for anc, s in sumstats_b37.items()}
    t_sumstats = time.time() - t0
    if verbose:
        print(f"[{locus_name}] fetched sumstats for "
              f"{list(sumstats)} in {t_sumstats:.1f}s (with GRCh37→38 liftover)")

    results: dict[str, dict] = {}
    for anc in ancestries:
        if anc not in sumstats:
            if verbose:
                print(f"  {anc}: no sumstats returned, skipping")
            continue

        # 2) Ancestry-matched 1KG LD (GRCh38 coordinates)
        t1 = time.time()
        try:
            ld_vars, R_sq = _compute_1kg_ld_matched(
                chr=locus["chr"], start=ld_start_b38, end=ld_end_b38,
                ancestry=anc, bgen_dir=bgen_dir, popmap_path=popmap_path,
            )
        except Exception as e:
            print(f"  {anc}: LD computation failed: {e}")
            continue
        t_ld = time.time() - t1
        if len(ld_vars) == 0:
            if verbose:
                print(f"  {anc}: no LD variants, skipping")
            continue

        # 3) Intersect sumstats with LD
        variants, z, R_aligned = _intersect_sumstats_ld(sumstats[anc], ld_vars, R_sq)
        n_int = len(variants)
        if n_int < 5:
            if verbose:
                print(f"  {anc}: only {n_int} intersecting variants, skipping")
            continue

        # 4) Sumstats-only fine-mapping (HBP and L1, no graph cache here;
        # this runs the algorithms with z + R only, no functional prior)
        t2 = time.time()
        hbp_out = hbp_finemap_from_sumstats(
            variants, z, R_aligned, graph_cache={}, chr_name=f"chr{locus['chr']}",
        )
        t_hbp = time.time() - t2
        t3 = time.time()
        l1_out = l1_finemap_from_sumstats(
            variants, z, R_aligned, z_func=None, alpha=1.0,
            chr_name=f"chr{locus['chr']}",
        )
        t_l1 = time.time() - t3

        top_hbp = hbp_out[0] if hbp_out else None
        top_l1 = l1_out[0] if l1_out else None
        n_cs_hbp = sum(1 for c in hbp_out if c.in_credible_set)
        n_cs_l1 = sum(1 for c in l1_out if c.in_credible_set)
        lead_sumstats_vid = sumstats[anc].variants.loc[
            sumstats[anc].variants["log10p"].idxmax(), "variant_id"
        ]

        results[anc] = {
            "N_sumstats": int(sumstats[anc].n_samples),
            "n_sumstats_variants": int(len(sumstats[anc].variants)),
            "n_ld_variants": int(len(ld_vars)),
            "n_intersect": int(n_int),
            "lead_sumstats_variant": lead_sumstats_vid,
            "top_hbp_variant": top_hbp.variant_id if top_hbp else None,
            "top_hbp_pip": float(top_hbp.pip) if top_hbp else None,
            "n_cs_hbp": int(n_cs_hbp),
            "top_l1_variant": top_l1.variant_id if top_l1 else None,
            "top_l1_pip": float(top_l1.pip) if top_l1 else None,
            "n_cs_l1": int(n_cs_l1),
            "t_ld_s": t_ld,
            "t_hbp_s": t_hbp,
            "t_l1_s": t_l1,
        }
        if verbose:
            r = results[anc]
            print(f"  {anc} N={r['N_sumstats']:>7,} intersect={n_int:>5}  "
                  f"lead={r['lead_sumstats_variant']:<28}  "
                  f"HBP_top={r['top_hbp_variant']:<28} PIP={r['top_hbp_pip']:.3f} "
                  f"CS={n_cs_hbp}  "
                  f"(LD {t_ld:.1f}s HBP {t_hbp*1000:.0f}ms L1 {t_l1*1000:.0f}ms)")

    record = {"locus": locus_name, **locus, "ancestries": results,
              "t_sumstats_s": t_sumstats}
    if out_dir:
        out_dir.mkdir(parents=True, exist_ok=True)
        with open(out_dir / f"{locus_name}.json", "w") as f:
            json.dump(record, f, indent=2)
    return record


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--locus", default="FTO", choices=sorted(LOCI))
    ap.add_argument("--all", action="store_true", help="Run all loci")
    ap.add_argument("--out", type=Path,
                    default=Path("/mnt/data/GraphGWAS/results/panukb"))
    args = ap.parse_args()

    if args.all:
        summaries = {name: run_locus(name, out_dir=args.out) for name in LOCI}
        with open(args.out / "summary.json", "w") as f:
            json.dump(summaries, f, indent=2)
        print(f"\nWrote {args.out}/summary.json with {len(summaries)} loci.")
    else:
        rec = run_locus(args.locus, out_dir=args.out)
        print(f"\nWrote {args.out}/{args.locus}.json")


if __name__ == "__main__":
    main()
