"""Summary statistics import and analysis for biobank-scale GWAS.

Mode 2 of GraphGWAS: instead of storing individual genotypes and running
GWAS internally, import pre-computed summary statistics from external
GWAS (Pan-UKBB, GWAS Catalog, FinnGen, etc.) and perform graph-native
analyses on them.

This enables GraphGWAS to work with 400K+ sample GWAS results across
thousands of traits without requiring individual-level genotype data.

Supported formats:
    - Pan-UKBB (per-phenotype TSV.bgz files)
    - Generic (chr, pos, ref, alt, beta, se, pval)
    - GWAS Catalog (standardized format)

Components:
    S1. parse_panukbb          — parse Pan-UKBB summary stats file
    S2. parse_generic          — parse generic GWAS summary stats
    S3. filter_significant     — extract significant hits
    S4. sumstats_to_zscores    — convert to z-scores for MPAT/LDSC
    S5. sumstats_gene_test     — gene-level test from summary stats (MAGMA-like)
    S6. sumstats_genetic_corr  — cross-trait genetic correlation from z-scores
    S7. sumstats_mr            — MR from two summary stats files
    S8. sumstats_pleiotropy    — cross-trait pleiotropy from multiple GWAS
    S9. sumstats_report        — integrated summary stats analysis report
"""

from __future__ import annotations

import gzip
import os
import time

import numpy as np
from scipy import stats as sp_stats


# ===================================================================
# S1. Pan-UKBB Parser
# ===================================================================

def parse_panukbb(path: str,
                  ancestry: str = "EUR",
                  p_threshold: float | None = None,
                  max_variants: int | None = None,
                  exclude_low_confidence: bool = True,
                  verbose: bool = True) -> dict:
    """Parse a Pan-UKBB per-phenotype summary statistics file.

    Pan-UKBB format: chr, pos, ref, alt, af_cases_{anc}, af_controls_{anc},
    beta_{anc}, se_{anc}, neglog10_pval_{anc}, low_confidence_{anc}

    Args:
        path: path to .tsv.bgz file.
        ancestry: ancestry group (EUR, AFR, EAS, SAS, AMR).
        p_threshold: if set, only return variants with p < threshold.
        max_variants: maximum variants to load (None = all).
        exclude_low_confidence: skip low-confidence variants.
        verbose: print progress.

    Returns:
        dict with arrays: chr, pos, ref, alt, beta, se, pval, neglog10p,
        af_cases, af_controls, n_variants.
    """
    if verbose:
        print(f"  Parsing Pan-UKBB: {os.path.basename(path)}")
        print(f"    Ancestry: {ancestry}, p_threshold: {p_threshold}")
        t0 = time.time()

    beta_col = f"beta_{ancestry}"
    se_col = f"se_{ancestry}"
    pval_col = f"neglog10_pval_{ancestry}"
    af_case_col = f"af_cases_{ancestry}"
    af_ctrl_col = f"af_controls_{ancestry}"
    lc_col = f"low_confidence_{ancestry}"

    chroms, positions, refs, alts = [], [], [], []
    betas, ses, neglog10ps = [], [], []
    af_cases, af_ctrls = [], []

    opener = gzip.open if path.endswith(".gz") or path.endswith(".bgz") else open

    with opener(path, "rt") as f:
        header = f.readline().strip().split("\t")
        col_idx = {name: i for i, name in enumerate(header)}

        # Verify columns exist
        for required in [beta_col, se_col, pval_col]:
            if required not in col_idx:
                raise ValueError(f"Column '{required}' not found. "
                                 f"Available: {header}")

        beta_i = col_idx[beta_col]
        se_i = col_idx[se_col]
        pval_i = col_idx[pval_col]
        af_case_i = col_idx.get(af_case_col)
        af_ctrl_i = col_idx.get(af_ctrl_col)
        lc_i = col_idx.get(lc_col)
        chr_i = col_idx["chr"]
        pos_i = col_idx["pos"]
        ref_i = col_idx["ref"]
        alt_i = col_idx["alt"]

        n_read = 0
        n_skipped_na = 0
        n_skipped_lc = 0
        n_skipped_p = 0

        for line in f:
            parts = line.strip().split("\t")

            # Skip NA values
            beta_str = parts[beta_i]
            if beta_str == "NA" or beta_str == "":
                n_skipped_na += 1
                continue

            try:
                beta = float(beta_str)
                se = float(parts[se_i])
                neglog10p = float(parts[pval_i])
            except (ValueError, IndexError):
                n_skipped_na += 1
                continue

            # Skip low confidence
            if exclude_low_confidence and lc_i is not None:
                if parts[lc_i].lower() == "true":
                    n_skipped_lc += 1
                    continue

            # P-value threshold filter
            if p_threshold is not None:
                p = 10 ** (-neglog10p)
                if p > p_threshold:
                    n_skipped_p += 1
                    continue

            chroms.append(parts[chr_i])
            positions.append(int(parts[pos_i]))
            refs.append(parts[ref_i])
            alts.append(parts[alt_i])
            betas.append(beta)
            ses.append(se)
            neglog10ps.append(neglog10p)

            if af_case_i is not None:
                try:
                    af_cases.append(float(parts[af_case_i]))
                except ValueError:
                    af_cases.append(np.nan)
            if af_ctrl_i is not None:
                try:
                    af_ctrls.append(float(parts[af_ctrl_i]))
                except ValueError:
                    af_ctrls.append(np.nan)

            n_read += 1
            if max_variants and n_read >= max_variants:
                break

            if verbose and n_read % 1000000 == 0:
                print(f"    ...{n_read:,} variants loaded", flush=True)

    # Convert to numpy
    result = {
        "chr": np.array(chroms),
        "pos": np.array(positions, dtype=np.int64),
        "ref": np.array(refs),
        "alt": np.array(alts),
        "beta": np.array(betas, dtype=np.float64),
        "se": np.array(ses, dtype=np.float64),
        "neglog10p": np.array(neglog10ps, dtype=np.float64),
        "pval": 10.0 ** (-np.array(neglog10ps, dtype=np.float64)),
        "n_variants": n_read,
        "format": "panukbb",
        "ancestry": ancestry,
        "source": path,
    }

    if af_cases:
        result["af_cases"] = np.array(af_cases, dtype=np.float64)
    if af_ctrls:
        result["af_controls"] = np.array(af_ctrls, dtype=np.float64)

    if verbose:
        elapsed = time.time() - t0
        print(f"    Loaded: {n_read:,} variants ({elapsed:.0f}s)")
        print(f"    Skipped: {n_skipped_na:,} NA, {n_skipped_lc:,} low-confidence"
              f"{f', {n_skipped_p:,} above threshold' if p_threshold else ''}")

    return result


# ===================================================================
# S2. Generic Parser
# ===================================================================

def parse_generic(path: str,
                  chr_col: str = "chr", pos_col: str = "pos",
                  ref_col: str = "ref", alt_col: str = "alt",
                  beta_col: str = "beta", se_col: str = "se",
                  pval_col: str = "pval",
                  delimiter: str = "\t",
                  p_threshold: float | None = None,
                  verbose: bool = True) -> dict:
    """Parse generic GWAS summary statistics file.

    Flexible column mapping for any format.

    Args:
        path: path to summary stats file (.tsv, .csv, .gz).
        *_col: column names for each field.
        delimiter: column delimiter.
        p_threshold: optional p-value filter.
        verbose: print progress.

    Returns:
        dict with arrays (same structure as parse_panukbb).
    """
    opener = gzip.open if path.endswith(".gz") else open

    chroms, positions, refs, alts = [], [], [], []
    betas, ses, pvals = [], [], []

    with opener(path, "rt") as f:
        header = f.readline().strip().split(delimiter)
        col_idx = {name.strip(): i for i, name in enumerate(header)}

        chr_i = col_idx[chr_col]
        pos_i = col_idx[pos_col]
        beta_i = col_idx[beta_col]
        se_i = col_idx[se_col]
        pval_i = col_idx[pval_col]
        ref_i = col_idx.get(ref_col)
        alt_i = col_idx.get(alt_col)

        n = 0
        for line in f:
            parts = line.strip().split(delimiter)
            try:
                b = float(parts[beta_i])
                s = float(parts[se_i])
                p = float(parts[pval_i])
            except (ValueError, IndexError):
                continue

            if p_threshold and p > p_threshold:
                continue

            chroms.append(parts[chr_i])
            positions.append(int(parts[pos_i]))
            refs.append(parts[ref_i] if ref_i is not None else ".")
            alts.append(parts[alt_i] if alt_i is not None else ".")
            betas.append(b)
            ses.append(s)
            pvals.append(p)
            n += 1

    neglog10p = -np.log10(np.maximum(np.array(pvals), 1e-300))

    if verbose:
        print(f"  Generic parser: {n:,} variants from {os.path.basename(path)}")

    return {
        "chr": np.array(chroms),
        "pos": np.array(positions, dtype=np.int64),
        "ref": np.array(refs),
        "alt": np.array(alts),
        "beta": np.array(betas, dtype=np.float64),
        "se": np.array(ses, dtype=np.float64),
        "pval": np.array(pvals, dtype=np.float64),
        "neglog10p": neglog10p,
        "n_variants": n,
        "format": "generic",
        "source": path,
    }


# ===================================================================
# S3. Filter Significant
# ===================================================================

def filter_significant(sumstats: dict, p_threshold: float = 5e-8) -> dict:
    """Extract significant variants from parsed summary stats.

    Args:
        sumstats: output of parse_panukbb or parse_generic.
        p_threshold: significance threshold.

    Returns:
        Filtered sumstats dict with only significant variants.
    """
    mask = sumstats["pval"] < p_threshold
    n_sig = int(np.sum(mask))

    filtered = {}
    for key, val in sumstats.items():
        if isinstance(val, np.ndarray) and len(val) == sumstats["n_variants"]:
            filtered[key] = val[mask]
        else:
            filtered[key] = val

    filtered["n_variants"] = n_sig
    filtered["p_threshold"] = p_threshold

    return filtered


# ===================================================================
# S4. Z-scores
# ===================================================================

def sumstats_to_zscores(sumstats: dict) -> np.ndarray:
    """Convert summary stats to z-scores.

    z = beta / se, or sign(beta) × sqrt(chi2) from p-value.

    Args:
        sumstats: parsed summary stats dict.

    Returns:
        (N,) array of z-scores.
    """
    se = sumstats["se"]
    beta = sumstats["beta"]

    # Prefer beta/se when SE is available and nonzero
    valid_se = (se > 0) & np.isfinite(se)
    z = np.zeros(len(beta))
    z[valid_se] = beta[valid_se] / se[valid_se]

    # Fallback: from p-value
    invalid = ~valid_se & np.isfinite(sumstats["pval"])
    if np.any(invalid):
        chi2 = sp_stats.chi2.isf(sumstats["pval"][invalid], 1)
        z[invalid] = np.sign(beta[invalid]) * np.sqrt(np.maximum(chi2, 0))

    return z


# ===================================================================
# S5. Gene-Level Test (MAGMA-like)
# ===================================================================

def sumstats_gene_test(sumstats: dict,
                       gene_annotations: dict | str,
                       window: int = 10000,
                       verbose: bool = True) -> list[dict]:
    """Gene-level association test from summary statistics.

    Aggregates per-variant z-scores within gene boundaries (+window).
    Uses the mean chi² approach (simplified MAGMA without LD correction).

    For LD-corrected gene tests, use with the 1KG LD reference.

    Args:
        sumstats: parsed summary stats.
        gene_annotations: dict of {gene: (chr, start, end)} or path to
                          GNExT gene annotation TSV.
        window: bp window around gene boundaries.
        verbose: print progress.

    Returns:
        List of {gene, chr, start, end, n_snps, mean_chi2, p_value}.
    """
    # Load gene annotations if path
    if isinstance(gene_annotations, str):
        gene_annotations = _load_gene_annotations(gene_annotations)

    z = sumstats_to_zscores(sumstats)
    chi2 = z ** 2
    variant_chr = sumstats["chr"]
    variant_pos = sumstats["pos"]

    results = []
    for gene, (gchr, gstart, gend) in gene_annotations.items():
        # Find variants in gene ± window
        mask = ((variant_chr == str(gchr)) &
                (variant_pos >= gstart - window) &
                (variant_pos <= gend + window))
        n_snps = int(np.sum(mask))

        if n_snps == 0:
            continue

        gene_chi2 = chi2[mask]
        mean_chi2 = float(np.mean(gene_chi2))

        # P-value: under null, mean(chi²) × n_snps ~ chi²(n_snps)
        # Simplified: use the Brown method approximation
        # For now: mean chi² vs chi²(1) expectation
        # Gene-level stat = sum(chi²) / sqrt(n_snps) (accounts for LD approx)
        stat = np.sum(gene_chi2) / np.sqrt(n_snps)
        gene_p = float(sp_stats.norm.sf(stat / np.sqrt(2)))  # approximate

        results.append({
            "gene": gene,
            "chr": str(gchr),
            "start": gstart,
            "end": gend,
            "n_snps": n_snps,
            "mean_chi2": mean_chi2,
            "stat": float(stat),
            "p_value": gene_p,
        })

    results.sort(key=lambda r: r["p_value"])

    if verbose:
        n_sig = sum(1 for r in results if r["p_value"] < 0.05 / max(len(results), 1))
        print(f"  Gene-level test: {len(results)} genes, "
              f"{n_sig} Bonferroni-significant")
        for r in results[:10]:
            print(f"    {r['gene']:<15s} stat={r['stat']:.2f} "
                  f"p={r['p_value']:.2e} ({r['n_snps']} SNPs)")

    return results


def _load_gene_annotations(path: str) -> dict:
    """Load gene annotations from GNExT-format TSV.

    Format: ENSG_ID<tab>chr<tab>start<tab>end<tab>strand<tab>symbol
    """
    genes = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                symbol = parts[5]
                chr_val = parts[1]
                start = int(parts[2])
                end = int(parts[3])
                genes[symbol] = (chr_val, start, end)
    return genes


# ===================================================================
# S6. Cross-Trait Genetic Correlation
# ===================================================================

def sumstats_genetic_corr(sumstats_list: list[dict],
                          trait_names: list[str] | None = None,
                          verbose: bool = True) -> dict:
    """Cross-trait genetic correlation from multiple GWAS summary stats.

    Uses z-score correlation at shared variants as a proxy for genetic
    correlation (simplified cross-trait LDSC).

    Args:
        sumstats_list: list of parsed summary stats dicts.
        trait_names: labels for each trait.
        verbose: print progress.

    Returns:
        dict with correlation_matrix, trait_names.
    """
    T = len(sumstats_list)
    if trait_names is None:
        trait_names = [f"trait_{i}" for i in range(T)]

    # Build variant index from first sumstats (assume shared variant order)
    # For simplicity: use position-based matching
    # Create variant key: chr:pos
    variant_keys = {}
    for i, ss in enumerate(sumstats_list):
        for j in range(ss["n_variants"]):
            key = f"{ss['chr'][j]}:{ss['pos'][j]}"
            if key not in variant_keys:
                variant_keys[key] = {}
            z = ss["beta"][j] / ss["se"][j] if ss["se"][j] > 0 else 0
            variant_keys[key][i] = z

    # Build z-score matrix for shared variants
    shared_keys = [k for k, v in variant_keys.items() if len(v) == T]

    if len(shared_keys) < 100:
        if verbose:
            print(f"  Too few shared variants ({len(shared_keys)})")
        return {"error": "Too few shared variants",
                "n_shared": len(shared_keys)}

    z_mat = np.zeros((len(shared_keys), T))
    for i, key in enumerate(shared_keys):
        for j in range(T):
            z_mat[i, j] = variant_keys[key].get(j, 0)

    # Genetic correlation ≈ correlation of z-scores
    rg = np.corrcoef(z_mat.T)

    if verbose:
        print(f"  Genetic correlation ({T} traits, "
              f"{len(shared_keys):,} shared variants):")
        for i in range(min(T, 8)):
            row = "    " + "  ".join(f"{rg[i,j]:6.3f}" for j in range(min(T, 8)))
            print(f"  {trait_names[i][:12]:<12s} {row}")

    return {
        "correlation_matrix": rg.tolist(),
        "trait_names": trait_names,
        "n_shared_variants": len(shared_keys),
        "n_traits": T,
    }


# ===================================================================
# S7. MR from Summary Stats
# ===================================================================

def sumstats_mr(exposure_ss: dict, outcome_ss: dict,
                p_threshold: float = 5e-8,
                verbose: bool = True) -> dict:
    """Two-sample MR from summary statistics.

    Selects instruments from exposure, looks up effects in outcome,
    runs IVW, Egger, and weighted median.

    Args:
        exposure_ss: parsed summary stats for exposure trait.
        outcome_ss: parsed summary stats for outcome trait.
        p_threshold: instrument selection threshold.
        verbose: print progress.

    Returns:
        dict with IVW, Egger, weighted median estimates.
    """
    from .mr import ivw_estimate, egger_estimate, weighted_median_estimate

    # Index outcome by position
    outcome_idx = {}
    for i in range(outcome_ss["n_variants"]):
        key = f"{outcome_ss['chr'][i]}:{outcome_ss['pos'][i]}"
        outcome_idx[key] = i

    # Select instruments from exposure
    instruments = []
    for i in range(exposure_ss["n_variants"]):
        if exposure_ss["pval"][i] > p_threshold:
            continue
        if exposure_ss["se"][i] <= 0:
            continue

        key = f"{exposure_ss['chr'][i]}:{exposure_ss['pos'][i]}"
        if key not in outcome_idx:
            continue

        oi = outcome_idx[key]
        if outcome_ss["se"][oi] <= 0:
            continue

        instruments.append({
            "variant": f"{exposure_ss['chr'][i]}:{exposure_ss['pos'][i]}",
            "beta_exposure": float(exposure_ss["beta"][i]),
            "se_exposure": float(exposure_ss["se"][i]),
            "p_exposure": float(exposure_ss["pval"][i]),
            "beta_outcome": float(outcome_ss["beta"][oi]),
            "se_outcome": float(outcome_ss["se"][oi]),
            "p_outcome": float(outcome_ss["pval"][oi]),
        })

    if verbose:
        print(f"  Two-sample MR: {len(instruments)} instruments "
              f"(p < {p_threshold:.0e})")

    if len(instruments) < 3:
        return {"error": "Too few instruments", "n_instruments": len(instruments)}

    results = {
        "n_instruments": len(instruments),
        "ivw": ivw_estimate(instruments, verbose=verbose),
        "egger": egger_estimate(instruments, verbose=verbose),
        "weighted_median": weighted_median_estimate(instruments, verbose=verbose),
    }

    return results


# ===================================================================
# S8. Cross-Trait Pleiotropy
# ===================================================================

def sumstats_pleiotropy(sumstats_list: list[dict],
                        trait_names: list[str],
                        p_threshold: float = 5e-8,
                        verbose: bool = True) -> dict:
    """Identify pleiotropic loci from multiple GWAS summary stats.

    For each variant, count how many traits it's significant for.

    Args:
        sumstats_list: list of parsed summary stats.
        trait_names: labels.
        p_threshold: significance threshold.
        verbose: print progress.

    Returns:
        dict with pleiotropic variants and counts.
    """
    # Count significant traits per variant (by position)
    variant_traits = {}

    for ti, ss in enumerate(sumstats_list):
        for i in range(ss["n_variants"]):
            if ss["pval"][i] >= p_threshold:
                continue
            key = f"{ss['chr'][i]}:{ss['pos'][i]}"
            if key not in variant_traits:
                variant_traits[key] = []
            variant_traits[key].append(trait_names[ti])

    # Sort by pleiotropy count
    pleiotropic = []
    for key, traits in variant_traits.items():
        if len(traits) >= 2:
            chr_val, pos = key.split(":")
            pleiotropic.append({
                "variant": key,
                "chr": chr_val,
                "pos": int(pos),
                "n_traits": len(traits),
                "traits": traits,
            })

    pleiotropic.sort(key=lambda v: v["n_traits"], reverse=True)

    n_any = len(variant_traits)
    n_pleio = len(pleiotropic)

    if verbose:
        print(f"  Pleiotropy: {n_any:,} variants sig for any trait, "
              f"{n_pleio:,} for ≥2 traits")
        for v in pleiotropic[:10]:
            traits_str = ", ".join(v["traits"][:5])
            print(f"    {v['variant']}: {v['n_traits']} traits — {traits_str}")

    return {
        "pleiotropic_variants": pleiotropic[:1000],  # cap for memory
        "n_significant_any": n_any,
        "n_pleiotropic": n_pleio,
        "max_traits": pleiotropic[0]["n_traits"] if pleiotropic else 0,
    }


# ===================================================================
# S9. Summary Stats Report
# ===================================================================

def sumstats_report(sumstats: dict,
                    gene_annotations: dict | str | None = None,
                    verbose: bool = True) -> dict:
    """Single-trait summary statistics analysis report.

    Args:
        sumstats: parsed summary stats.
        gene_annotations: for gene-level test (optional).
        verbose: print progress.

    Returns:
        Comprehensive report dict.
    """
    if verbose:
        print("=" * 60)
        print(f"  Summary Statistics Report")
        print(f"  Source: {os.path.basename(sumstats.get('source', '?'))}")
        print("=" * 60)

    n = sumstats["n_variants"]
    pvals = sumstats["pval"]
    betas = sumstats["beta"]

    # Basic stats
    n_gw_sig = int(np.sum(pvals < 5e-8))
    n_suggestive = int(np.sum(pvals < 1e-5))
    min_p = float(np.min(pvals)) if n > 0 else 1.0

    # Lambda GC
    from .popstruct import lambda_gc
    lgc = lambda_gc(pvals)

    # Manhattan peaks (clump by 500kb windows)
    loci = []
    sig_idx = np.where(pvals < 5e-8)[0]
    if len(sig_idx) > 0:
        sig_sorted = sig_idx[np.argsort(pvals[sig_idx])]
        used = set()
        for idx in sig_sorted:
            chr_val = sumstats["chr"][idx]
            pos = sumstats["pos"][idx]
            # Check if within 500kb of existing locus
            is_new = True
            for u_chr, u_pos in used:
                if chr_val == u_chr and abs(pos - u_pos) < 500000:
                    is_new = False
                    break
            if is_new:
                used.add((chr_val, pos))
                loci.append({
                    "chr": str(chr_val),
                    "pos": int(pos),
                    "pval": float(pvals[idx]),
                    "beta": float(betas[idx]),
                    "neglog10p": float(sumstats["neglog10p"][idx]),
                })

    report = {
        "n_variants": n,
        "n_genome_wide_significant": n_gw_sig,
        "n_suggestive": n_suggestive,
        "min_pvalue": min_p,
        "lambda_gc": lgc,
        "n_independent_loci": len(loci),
        "top_loci": loci[:50],
    }

    if verbose:
        print(f"\n  Variants: {n:,}")
        print(f"  Genome-wide significant: {n_gw_sig:,}")
        print(f"  Suggestive: {n_suggestive:,}")
        print(f"  Lambda GC: {lgc:.3f}")
        print(f"  Independent loci: {len(loci)}")
        print(f"\n  Top loci:")
        for locus in loci[:10]:
            print(f"    chr{locus['chr']}:{locus['pos']} "
                  f"p={locus['pval']:.2e} β={locus['beta']:.4f}")

    # Gene-level test (if annotations available)
    if gene_annotations:
        if verbose:
            print(f"\n  Gene-level test:")
        gene_results = sumstats_gene_test(sumstats, gene_annotations,
                                          verbose=verbose)
        report["gene_results"] = gene_results[:100]

    return report
