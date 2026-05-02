"""Common interface to external epistasis-detection tools (paper #2 §Y.3).

Three baselines from the literature, all with the same call signature:

    run_boost(...)   -> pandas.DataFrame  # PLINK2's --epistasis (Wan 2010 algorithm)
    run_mapit(...)   -> pandas.DataFrame  # marginal epistasis test (Crawford 2017)
    run_mdr(...)     -> pandas.DataFrame  # multifactor dimensionality reduction

All return a DataFrame with columns:
    variant_1, variant_2, score, p_value, method

allowing direct comparison against GraphGWAS M2's output schema.

Installation (run once before using each tool):

    BOOST:   PLINK2 already supports the BOOST algorithm via --epistasis.
             We wrap PLINK2 instead of compiling the standalone BOOST C tool.
             Install:  mamba install -c bioconda plink2

    MAPIT:   R package; uses the official Crawford-lab MAPIT C++ via Rcpp.
             Install:  conda install -c conda-forge r-base
                       R -e 'install.packages("mapit", repos="http://cran.r-project.org")'
             Alternatively, the Python-only port mapit-py (less battle-tested).

    MDR:     scikit-mdr is the pure-Python re-implementation; pip-installable
             and easiest of the three.
             Install:  pip install scikit-mdr

Each wrapper checks for its dependency on first call and raises a clear
RuntimeError with the install command if missing.
"""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


# ===========================================================================
# Common output schema
# ===========================================================================

BASELINE_COLUMNS = ["variant_1", "variant_2", "score", "p_value", "method"]


def _empty_result(method: str) -> pd.DataFrame:
    return pd.DataFrame(columns=BASELINE_COLUMNS).astype({
        "variant_1": "string",
        "variant_2": "string",
        "score": "float64",
        "p_value": "float64",
        "method": "string",
    })


# ===========================================================================
# BOOST  (via PLINK2 --epistasis)
# ===========================================================================

def _check_plink2() -> str:
    p = shutil.which("plink2")
    if p is None:
        raise RuntimeError(
            "plink2 not on PATH.  Install:\n"
            "  mamba install -c bioconda plink2"
        )
    return p


def run_boost(
    bgen: Path,
    sample: Path,
    pheno: Path,
    *,
    out_dir: Path,
    p_threshold: float = 5e-8,
    extra: list[str] | None = None,
) -> pd.DataFrame:
    """BOOST-equivalent test via PLINK2's --epistasis.

    Wan 2010's BOOST algorithm: rapid 2x3x3 chi-squared screen followed
    by full LRT.  PLINK2's --epistasis uses the same boolean encoding.

    Args:
        bgen: BGEN file (whole-chromosome).
        sample: PLINK2 .sample file matching the BGEN.
        pheno: Phenotype file in PLINK2 format
            (header "FID IID PHENO_1"; values; one phenotype per call).
        out_dir: directory for plink2 output files.
        p_threshold: report only pairs at p < this.
        extra: extra plink2 flags.

    Returns:
        DataFrame in BASELINE_COLUMNS schema, sorted by p_value ascending.
    """
    plink = _check_plink2()
    out_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = out_dir / "boost_epistasis"
    cmd = [
        plink, "--bgen", str(bgen), "--sample", str(sample),
        "--pheno", str(pheno),
        "--epistasis", "set-by-set",  # BOOST-style mode
        "--epistasis-threshold", str(p_threshold),
        "--out", str(out_prefix),
    ]
    if extra:
        cmd.extend(extra)
    rc = subprocess.call(cmd)
    if rc != 0:
        raise RuntimeError(f"plink2 --epistasis failed (rc={rc})")

    epi_file = out_prefix.with_suffix(".epi.qt")
    if not epi_file.exists():
        # Try linear/case-control flavours
        for ext in (".epi.cc", ".epi.qt"):
            cand = out_prefix.with_suffix(ext)
            if cand.exists():
                epi_file = cand
                break
        else:
            return _empty_result("BOOST")

    df = pd.read_csv(epi_file, sep=r"\s+")
    # plink2 columns: CHR1, ID1, CHR2, ID2, OR/BETA, STAT, P
    out = pd.DataFrame({
        "variant_1": df["ID1"].astype(str),
        "variant_2": df["ID2"].astype(str),
        "score": df.get("STAT", df.get("CHISQ", np.nan)),
        "p_value": df["P"],
        "method": "BOOST/plink2",
    })
    return out.sort_values("p_value").reset_index(drop=True)


# ===========================================================================
# MAPIT  (via R)
# ===========================================================================

def _check_R(packages: list[str]) -> str:
    r = shutil.which("R") or shutil.which("Rscript")
    if r is None:
        raise RuntimeError(
            "R not on PATH.  Install:\n"
            "  conda install -c conda-forge r-base\n"
            f"  R -e 'install.packages(c({', '.join(repr(p) for p in packages)}))'"
        )
    return r


def run_mapit(
    genotype_matrix: np.ndarray,
    phenotype: np.ndarray,
    variant_ids: list[str],
    *,
    p_threshold: float = 0.05,
) -> pd.DataFrame:
    """Marginal epistasis test (Crawford et al. 2017, PLOS Genet 13:e1006869).

    For each variant g_k, MAPIT tests whether g_k participates in ANY
    pairwise interaction across the genome — without enumerating pairs.
    Output is per-variant, not per-pair.  We map this to the common
    schema by emitting (variant_k, "*", score, p_value) rows.

    Args:
        genotype_matrix: shape (n_samples, n_variants).
        phenotype: shape (n_samples,).
        variant_ids: list of n_variants identifiers.
        p_threshold: report only variants at p < this (per-variant, not pairwise).

    Returns:
        DataFrame in BASELINE_COLUMNS schema.

    [SCAFFOLD]: Currently calls Rscript via subprocess with an inline
    R script that loads the official 'mapit' R package.  Full
    implementation needs (a) MAPIT installation verified, (b) data
    serialisation via .rds or feather, (c) result parsing.

    TODO(paper2-Y3-mapit):
    - Add proper Rscript template at src/python/graphgwas/mapit_template.R
    - Wire genotype/phenotype via feather (faster than .rds)
    - Parse output, emit BASELINE_COLUMNS DataFrame
    """
    raise NotImplementedError(
        "MAPIT R wrapper not yet implemented.  See TODO in "
        "graphgwas.epistasis_baselines.run_mapit. "
        "Install plan: conda install -c conda-forge r-base; "
        "R -e 'install.packages(\"mapit\", repos=\"http://cran.r-project.org\")'"
    )


# ===========================================================================
# MDR  (via scikit-mdr)
# ===========================================================================

def run_mdr(
    genotype_matrix: np.ndarray,
    phenotype: np.ndarray,
    variant_ids: list[str],
    *,
    n_pairs_test: int = 1000,
    cv_folds: int = 10,
    seed: int = 0,
) -> pd.DataFrame:
    """Multifactor Dimensionality Reduction (scikit-mdr port).

    For each pair of variants in a tested set, MDR builds a 9-cell
    contingency-table classifier and scores its cross-validated balanced
    accuracy.  We rank pairs by accuracy.

    Args:
        genotype_matrix: shape (n_samples, n_variants), integer 0/1/2 dosages.
        phenotype: shape (n_samples,), binary case/control encoding 0/1.
        variant_ids: list of n_variants identifiers.
        n_pairs_test: number of random pairs to evaluate (full enumeration
            of millions of pairs is impractical; sample n_pairs_test).
        cv_folds: cross-validation folds for accuracy estimate.
        seed: RNG seed for pair sampling + CV split.

    Returns:
        DataFrame in BASELINE_COLUMNS schema, sorted by score (cv accuracy)
        descending.

    Install: pip install scikit-mdr
    """
    try:
        from mdr import MDR  # noqa: F401
    except ImportError:
        raise RuntimeError(
            "scikit-mdr not installed.  Install:\n"
            "  pip install scikit-mdr"
        )
    from mdr import MDR
    from sklearn.model_selection import cross_val_score

    rng = np.random.default_rng(seed)
    n_var = genotype_matrix.shape[1]
    if n_var < 2:
        return _empty_result("MDR")
    # Sample random pairs
    pair_indices = set()
    while len(pair_indices) < n_pairs_test:
        i, j = rng.integers(0, n_var, size=2)
        if i != j:
            pair_indices.add((min(i, j), max(i, j)))
    pair_indices = list(pair_indices)

    rows = []
    for (i, j) in pair_indices:
        X_pair = genotype_matrix[:, [i, j]].astype(np.int8)
        clf = MDR()
        try:
            scores = cross_val_score(clf, X_pair, phenotype, cv=cv_folds,
                                     scoring="balanced_accuracy")
            cv_acc = float(scores.mean())
        except Exception:  # noqa: BLE001
            continue
        rows.append({
            "variant_1": variant_ids[i],
            "variant_2": variant_ids[j],
            "score": cv_acc,
            "p_value": np.nan,  # MDR doesn't emit p-values; permutation needed
            "method": "MDR",
        })
    df = pd.DataFrame(rows, columns=BASELINE_COLUMNS)
    return df.sort_values("score", ascending=False).reset_index(drop=True)


# ===========================================================================
# CLI / smoke test
# ===========================================================================

if __name__ == "__main__":
    print("epistasis_baselines.py — wrappers for BOOST/MAPIT/MDR")
    print()
    print("Tool availability:")
    for name, check in [
        ("plink2  (BOOST)", lambda: shutil.which("plink2") is not None),
        ("R       (MAPIT)", lambda: shutil.which("R") is not None
                                     or shutil.which("Rscript") is not None),
        ("scikit-mdr (MDR)", lambda: __import__("importlib").util.find_spec("mdr") is not None),
    ]:
        try:
            avail = check()
        except Exception:  # noqa: BLE001
            avail = False
        print(f"  {name}: {'✓ available' if avail else '✗ not installed'}")
    print()
    print("See module docstring for install commands.")
