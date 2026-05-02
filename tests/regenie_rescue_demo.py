"""REGENIE rescue demonstration on 1KG chr22 — paper #2 §Y.2.

Reproduces the headline rescue claim: GraphGWAS M2 motif-pair detection,
fed back into REGENIE as covariates, collapses the spurious-significance
regime that omitted interactions create (Yelmen et al. 2026).

Pipeline:

    1. Simulate phenotypes under Yelmen's strict no-path null DGP:
           y = G_t @ alpha + Z @ theta_u + eps,  alpha = 0
       where target SNPs are excluded from Z's interaction features.

    2. Run REGENIE step 1 + step 2 on the simulated phenotype, scanning
       genome-wide.  Count SNPs at p < 5e-8 (these are spurious — the
       DGP has no marginal effect).

    3. Run GraphGWAS M2 motif-pair detection (P1 same-gene); add the
       detected pairs as interaction covariates to REGENIE step 2.

    4. Recount spurious significant SNPs.  Headline number:
       X -> Y % reduction in false-positive count.

PREREQUISITES (user must install before Step 2):

    mamba install -n graphmana -c bioconda regenie

    or

    conda install -c bioconda regenie

REGENIE 4.1.2 is available on bioconda; Yelmen et al. used v3.4 — both
are compatible with this script.  After install:

    regenie --version

should report.  Then run this script.

USAGE:
    python tests/regenie_rescue_demo.py --simulate-only   # Stage 1 only
    python tests/regenie_rescue_demo.py --regenie-only    # Stages 1 + 2
    python tests/regenie_rescue_demo.py                   # Full rescue demo

OUTPUT:
    paper/epistasis_v1/figures/fig3_regenie_rescue.{png,pdf,json}
    results/paper2_epistasis/regenie_*.{regenie,log}      # raw REGENIE output
    results/paper2_epistasis/spurious_counts.json         # before/after
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src" / "python"))

from graphgwas.bgen_reader import BgenReader  # noqa: E402
from graphgwas.epistasis_v2 import motif_filtered_epistasis_from_data  # noqa: E402

BGEN_DIR = REPO_ROOT / "tests" / "data" / "human" / "1kGP_bgen"
CACHE_PATH = REPO_ROOT / "data" / "annotations" / "human_graph_cache_v2_chr22.json"
RESULTS_DIR = REPO_ROOT / "results" / "paper2_epistasis"
FIG_DIR = REPO_ROOT / "paper" / "epistasis_v1" / "figures"


# ===========================================================================
# Stage 1: Phenotype simulator (strict no-path null DGP, Yelmen Eq. 1)
# ===========================================================================

def simulate_phenotype(
    dosages: np.ndarray,
    *,
    n_target: int = 1,
    n_interaction_features: int = 100,
    lambda_var: float = 0.1,
    rng: np.random.Generator,
) -> dict:
    """Simulate a strict-no-path-null phenotype.

    DGP (Yelmen et al. 2026 Eq. 1, simplified for demonstration):

        y = G_t @ alpha + Z @ theta_u + eps      (alpha = 0)

    where Z is built from random two-way SNP×SNP products and the target
    SNP G_t is forbidden from the interaction features (strict no-path).

    Args:
        dosages: dosage matrix, shape (n_samples, n_variants).
        n_target: number of target SNPs to simulate (each will become
            a separate phenotype column for batched REGENIE).
        n_interaction_features: width of Z, number of interaction-product columns.
        lambda_var: target variance fraction from the interaction term;
            controls signal strength of the bias.
        rng: numpy Generator.

    Returns:
        dict with keys:
            y                 : (n_samples, n_target) phenotype matrix
            target_idx        : (n_target,) target SNP column indices
            interaction_pairs : (n_interaction_features, 2) interaction SNP indices
            theta_u           : (n_interaction_features,) interaction coefficients
            lambda_var        : echoed
    """
    n_samples, n_variants = dosages.shape

    # Pick interaction SNPs first (m_int * 2 unique columns)
    int_pool_size = n_interaction_features * 2
    int_pool = rng.choice(n_variants, size=int_pool_size, replace=False)
    a_idx = int_pool[:n_interaction_features]
    b_idx = int_pool[n_interaction_features:]

    # Target SNPs: chosen DISJOINT from interaction pool (strict no-path)
    forbidden = set(int_pool.tolist())
    candidate_targets = np.array(
        [i for i in range(n_variants) if i not in forbidden]
    )
    target_idx = rng.choice(candidate_targets, size=n_target, replace=False)

    # Build Z = m_int columns of centred dosage products
    Z = dosages[:, a_idx] * dosages[:, b_idx]
    Z = Z - Z.mean(axis=0, keepdims=True)
    # Standardise per-column for predictable variance accounting
    Z = Z / np.std(Z, axis=0, keepdims=True).clip(min=1e-8)

    # theta_u from N(0, 1); scale so Var(Z @ theta_u) ≈ lambda_var / (1 - lambda_var)
    # (because lambda_var = Var(u) / (1 + Var(u)) ⇒ Var(u) = lambda_var/(1-lambda_var))
    target_var_u = lambda_var / (1.0 - lambda_var)
    theta_u = rng.standard_normal(n_interaction_features)
    u = Z @ theta_u
    # Rescale theta_u so empirical var(u) matches target_var_u
    u = u * np.sqrt(target_var_u / u.var())
    theta_u = theta_u * np.sqrt(target_var_u / (Z @ theta_u).var())

    # Noise: var = 1
    eps = rng.standard_normal((n_samples, n_target))

    # alpha = 0 by construction (strict null)
    # Phenotype: y = u (broadcast across n_target columns) + eps
    y = u[:, None] + eps

    return {
        "y": y.astype(np.float32),
        "target_idx": target_idx,
        "interaction_pairs": np.column_stack([a_idx, b_idx]),
        "theta_u": theta_u,
        "lambda_var": float(lambda_var),
        "n_samples": int(n_samples),
        "n_variants": int(n_variants),
    }


# ===========================================================================
# Stage 2: REGENIE wrapper (subprocess-based)
# ===========================================================================

def regenie_available() -> bool:
    return shutil.which("regenie") is not None


def write_regenie_pheno(y: np.ndarray, sample_ids: list[str], out: Path) -> Path:
    """Write a REGENIE-compatible .pheno file.

    Format: header line "FID IID PHENO_1 PHENO_2 ..."; rows are
    "0 <sample_id> <values...>".  REGENIE expects FID and IID columns;
    1KG samples don't have FIDs, so we use 0 as a placeholder.
    """
    out.parent.mkdir(parents=True, exist_ok=True)
    n_samples, n_pheno = y.shape
    assert len(sample_ids) == n_samples
    header = "FID IID " + " ".join(f"Y{i+1}" for i in range(n_pheno))
    lines = [header]
    for sid, row in zip(sample_ids, y):
        lines.append(f"0 {sid} " + " ".join(f"{v:.6f}" for v in row))
    out.write_text("\n".join(lines) + "\n")
    return out


def run_regenie_step1(
    bgen: Path, sample: Path, pheno: Path, out_prefix: Path, *,
    extra: list[str] | None = None,
) -> int:
    """Run REGENIE step 1 (whole-genome ridge regression).

    Output: <out_prefix>_pred.list and <out_prefix>_*.loco files.
    Returns subprocess exit code.
    """
    if not regenie_available():
        raise RuntimeError(
            "REGENIE not on PATH. Install via:\n"
            "  mamba install -n graphmana -c bioconda regenie\n"
            "Then 'conda activate graphmana' before re-running."
        )
    cmd = [
        "regenie", "--step", "1",
        "--bgen", str(bgen),
        "--sample", str(sample),
        "--phenoFile", str(pheno),
        "--bsize", "1000",
        "--out", str(out_prefix),
    ]
    if extra:
        cmd.extend(extra)
    print(f"  + {' '.join(cmd)}")
    return subprocess.call(cmd)


def run_regenie_step2(
    bgen: Path, sample: Path, pheno: Path,
    pred_list: Path, out_prefix: Path, *,
    covar: Path | None = None,
    extra: list[str] | None = None,
) -> int:
    """Run REGENIE step 2 (per-SNP association).

    Output: <out_prefix>_*.regenie summary statistics.
    """
    if not regenie_available():
        raise RuntimeError("REGENIE not on PATH (see step1 error message).")
    cmd = [
        "regenie", "--step", "2",
        "--bgen", str(bgen),
        "--sample", str(sample),
        "--phenoFile", str(pheno),
        "--pred", str(pred_list),
        "--bsize", "400",
        "--out", str(out_prefix),
    ]
    if covar is not None:
        cmd.extend(["--covarFile", str(covar)])
    if extra:
        cmd.extend(extra)
    print(f"  + {' '.join(cmd)}")
    return subprocess.call(cmd)


# ===========================================================================
# Stage 3: GraphGWAS M2 motif detection (no-Neo4j path)
# ===========================================================================

def detect_motif_interactions_no_neo4j(
    dosages: np.ndarray,
    variant_ids: list[str],
    phenotype: np.ndarray,
    graph_cache_path: Path,
    *,
    motifs: list[str] | None = None,
    fdr_threshold: float = 0.05,
    mac_min: int = 10,
    max_pairs_per_entity: int = 500,
    max_pairs_total: int = 50_000,
    verbose: bool = True,
) -> list:
    """Detect M2 (motif-filtered) interactions without Neo4j.

    Thin wrapper around
    :func:`graphgwas.epistasis_v2.motif_filtered_epistasis_from_data`
    that returns only the BH-FDR-significant pairs.

    Args:
        dosages: shape ``(n_samples, n_variants)``, NaN-tolerant.
        variant_ids: length ``n_variants`` in ``chrN:pos:REF:ALT`` format
            (matching the keys of ``graph_cache_path``).  When loading
            from BGEN, swap ``a1, a2`` → ``REF, ALT`` first
            (plink2 puts ALT in slot 0).
        phenotype: shape ``(n_samples,)``.
        graph_cache_path: JSON cache.
        motifs: default ``["same_gene", "same_pathway"]``.
        fdr_threshold: BH-FDR threshold for significance.
        mac_min, max_pairs_per_entity, max_pairs_total: passed through.
        verbose: print progress.

    Returns:
        List of significant ``InteractionResult`` records, sorted by
        ``p_interaction`` ascending.
    """
    if motifs is None:
        motifs = ["same_gene", "same_pathway"]
    results = motif_filtered_epistasis_from_data(
        dosages=dosages,
        variant_ids=variant_ids,
        phenotype=phenotype,
        graph_cache=graph_cache_path,
        motifs=motifs,
        mac_min=mac_min,
        correction="BH",
        max_pairs_per_entity=max_pairs_per_entity,
        max_pairs_total=max_pairs_total,
        verbose=verbose,
    )
    return [r for r in results
            if r.p_corrected is not None and r.p_corrected < fdr_threshold]


def write_regenie_covar(
    pair_indices: list[tuple[int, int]],
    dosages: np.ndarray,
    sample_ids: list[str],
    out: Path,
) -> Path:
    """Write a REGENIE-compatible covariate file built from motif-pair products.

    Each row: ``FID IID PAIR_001 PAIR_002 ... PAIR_K`` where each PAIR_k
    column is the centred dosage product of the (i, j) pair for that sample.

    Args:
        pair_indices: list of (i, j) variant-column indices.
        dosages: shape (n_samples, n_variants).
        sample_ids: length n_samples.
        out: destination path.

    Returns:
        ``out`` (unchanged), for convenience.
    """
    out.parent.mkdir(parents=True, exist_ok=True)
    n_samples = dosages.shape[0]
    if len(sample_ids) != n_samples:
        raise ValueError(
            f"sample_ids length ({len(sample_ids)}) != n_samples ({n_samples})"
        )
    cols = []
    for i, j in pair_indices:
        prod = dosages[:, i] * dosages[:, j]
        cols.append(prod - np.nanmean(prod))
    if cols:
        Z = np.column_stack(cols)
    else:
        Z = np.empty((n_samples, 0), dtype=np.float64)
    header = "FID IID " + " ".join(f"PAIR_{k+1:03d}" for k in range(Z.shape[1]))
    lines = [header]
    for sid, row in zip(sample_ids, Z):
        lines.append(f"0 {sid} " + " ".join(f"{v:.6f}" for v in row))
    out.write_text("\n".join(lines) + "\n")
    return out


# ===========================================================================
# Stage 4: Rescue counter
# ===========================================================================

def count_spurious_hits(regenie_output: Path, threshold: float = 5e-8) -> int:
    """Count SNPs with REGENIE p < threshold (= 'spurious significance').

    REGENIE writes one row per (variant, phenotype).  We count rows where
    the LOG10P column exceeds -log10(threshold).
    """
    import math
    cutoff_log10 = -math.log10(threshold)
    n_spurious = 0
    with regenie_output.open() as f:
        header = f.readline().split()
        try:
            log10p_col = header.index("LOG10P")
        except ValueError:
            raise RuntimeError(f"LOG10P column not in {regenie_output}: {header}")
        for line in f:
            parts = line.split()
            if len(parts) <= log10p_col:
                continue
            if parts[log10p_col] in ("NA", "."):
                continue
            if float(parts[log10p_col]) >= cutoff_log10:
                n_spurious += 1
    return n_spurious


# ===========================================================================
# Driver
# ===========================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--simulate-only", action="store_true",
                        help="Stage 1 only — write pheno file, exit before REGENIE.")
    parser.add_argument("--regenie-only", action="store_true",
                        help="Stages 1 and 2 only — skip M2 detection + rescue.")
    parser.add_argument("--m2-only", action="store_true",
                        help="Stages 1 and 3 only — simulate phenotype + run M2 motif "
                             "detection + write covariate file; skip REGENIE.")
    parser.add_argument("--fdr-threshold", type=float, default=0.05,
                        help="BH-FDR threshold for M2 significance (Stage 3).")
    parser.add_argument("--lambda-var", type=float, default=0.1,
                        help="Variance fraction from interactions (Yelmen lambda)")
    parser.add_argument("--n-target", type=int, default=1,
                        help="Number of target SNPs / phenotypes to simulate")
    parser.add_argument("--m-int", type=int, default=100,
                        help="Number of interaction-feature columns in Z")
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()

    print("=== REGENIE rescue demo ===")
    print(f"BGEN dir       : {BGEN_DIR}")
    print(f"Results dir    : {RESULTS_DIR}")
    print(f"REGENIE on PATH: {regenie_available()}")
    print(f"lambda_var     : {args.lambda_var}")
    print(f"m_int          : {args.m_int}")
    print(f"n_target       : {args.n_target}")
    print()

    rng = np.random.default_rng(args.seed)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # ----- Stage 1: simulate phenotype -----
    print("[1/4] Loading chr22 dosages + simulating phenotype...")
    with BgenReader(BGEN_DIR) as reader:
        df, dos = reader.load_locus("22", 30_000_000, 35_000_000, format="dosage")
        sample_ids = reader.samples("22").tolist()
    af = np.nanmean(dos, axis=0) / 2.0
    keep = np.minimum(af, 1.0 - af) >= 0.05
    dos = dos[:, keep].astype(np.float64)
    df = df[keep].reset_index(drop=True)
    # Build variant IDs in 'chrN:pos:REF:ALT' format (matching graph_cache).
    # BGEN reader's `a1, a2` are (ALT, REF) per plink2 convention — swap.
    variant_ids = [
        f"{r.chr}:{int(r.pos)}:{r.a2}:{r.a1}"
        for r in df.itertuples()
    ]
    print(f"  dosages: n_samples={dos.shape[0]}, n_variants={dos.shape[1]}")

    sim = simulate_phenotype(
        dos,
        n_target=args.n_target,
        n_interaction_features=args.m_int,
        lambda_var=args.lambda_var,
        rng=rng,
    )

    pheno_path = RESULTS_DIR / "sim_pheno.txt"
    write_regenie_pheno(sim["y"], sample_ids, pheno_path)
    print(f"  ✓ wrote {pheno_path}")
    print(f"    target SNP indices (chr22): {sim['target_idx'][:5].tolist()}...")
    print(f"    interaction pairs sampled : {len(sim['interaction_pairs'])}")
    print(f"    achieved Var(u)/Var(y)    : "
          f"{sim['theta_u'].dot(sim['theta_u']) / (1.0 + sim['theta_u'].dot(sim['theta_u'])):.3f}")

    if args.simulate_only:
        return

    # ----- Stage 2: REGENIE -----
    if args.m2_only:
        print("\n[2/4] Skipping REGENIE (--m2-only)")
    else:
        print("\n[2/4] Running REGENIE step 1 + step 2...")
        if not regenie_available():
            print("  ⚠ regenie not installed — skipping Stages 2+4.")
            print("    To enable later:  mamba install -n graphmana -c bioconda regenie")
            print("  Continuing with Stage 3 (M2 motif detection)...")
            args.m2_only = True  # graceful fallthrough to Stage 3

    n_baseline = None
    if not args.m2_only:
        bgen_path = BGEN_DIR / "chr22.bgen"
        sample_path = BGEN_DIR / "chr22.sample"
        step1_prefix = RESULTS_DIR / "regenie_step1"
        step2_prefix = RESULTS_DIR / "regenie_step2_baseline"

        rc = run_regenie_step1(bgen_path, sample_path, pheno_path, step1_prefix)
        if rc != 0:
            sys.exit(f"REGENIE step1 failed with exit code {rc}")
        rc = run_regenie_step2(
            bgen_path, sample_path, pheno_path,
            step1_prefix.with_name(step1_prefix.name + "_pred.list"),
            step2_prefix,
        )
        if rc != 0:
            sys.exit(f"REGENIE step2 failed with exit code {rc}")

        baseline_path = step2_prefix.with_suffix(".regenie")
        n_baseline = count_spurious_hits(baseline_path)
        print(f"  ✓ baseline spurious significant SNPs: {n_baseline:,}")

    if args.regenie_only:
        return

    # ----- Stage 3: M2 motif detection -----
    print("\n[3/4] Running M2 motif-pair detection (no-Neo4j path)...")
    print(f"  graph cache: {CACHE_PATH.name}")
    if not CACHE_PATH.exists():
        sys.exit(f"  ✗ graph cache missing: {CACHE_PATH}")

    sig_results = detect_motif_interactions_no_neo4j(
        dosages=dos,
        variant_ids=variant_ids,
        phenotype=sim["y"][:, 0],
        graph_cache_path=CACHE_PATH,
        motifs=["same_gene", "same_pathway"],
        fdr_threshold=args.fdr_threshold,
        mac_min=10,
        max_pairs_per_entity=500,
        max_pairs_total=50_000,
        verbose=True,
    )
    print(f"  ✓ {len(sig_results):,} pairs at BH-FDR < {args.fdr_threshold}")
    if sig_results:
        best = sig_results[0]
        print(f"    best: {best.variant_1} × {best.variant_2}  "
              f"motif={best.motif}  q={best.p_corrected:.2e}")

    # Map InteractionResult variant IDs back to dosage column indices
    vid_to_col = {v: i for i, v in enumerate(variant_ids)}
    pair_indices = [
        (vid_to_col[r.variant_1], vid_to_col[r.variant_2])
        for r in sig_results
    ]

    # Emit REGENIE covariate file
    covar_path = RESULTS_DIR / "motif_covariates.txt"
    write_regenie_covar(pair_indices, dos, sample_ids, covar_path)
    print(f"  ✓ wrote covariate file: {covar_path}")
    print(f"    columns: FID IID + {len(pair_indices)} pair-product columns")

    # Persist a JSON sidecar with the M2 detection result
    sidecar = RESULTS_DIR / "stage3_m2_results.json"
    sidecar.write_text(json.dumps({
        "n_significant": len(sig_results),
        "fdr_threshold": args.fdr_threshold,
        "motifs": ["same_gene", "same_pathway"],
        "lambda_var": args.lambda_var,
        "n_samples": int(dos.shape[0]),
        "n_variants": int(dos.shape[1]),
        "top_pairs": [
            {
                "variant_1": r.variant_1, "variant_2": r.variant_2,
                "motif": r.motif, "shared_entity": r.shared_entity,
                "p_interaction": float(r.p_interaction),
                "p_corrected": float(r.p_corrected) if r.p_corrected else None,
                "beta_interaction": float(r.beta_interaction),
            }
            for r in sig_results[:20]
        ],
    }, indent=2))
    print(f"  ✓ wrote {sidecar}")

    if args.m2_only:
        print("\n[done] --m2-only specified; stopping before Stage 4.")
        return

    # ----- Stage 4: rescue REGENIE step 2 with motif covariates + count -----
    print("\n[4/4] Re-running REGENIE step 2 with motif covariates...")
    if not regenie_available() or n_baseline is None:
        print("  ⚠ baseline not produced (REGENIE missing); cannot compute rescue Δ.")
        return

    rescue_prefix = RESULTS_DIR / "regenie_step2_rescue"
    rc = run_regenie_step2(
        BGEN_DIR / "chr22.bgen",
        BGEN_DIR / "chr22.sample",
        pheno_path,
        (RESULTS_DIR / "regenie_step1_pred.list"),
        rescue_prefix,
        covar=covar_path,
    )
    if rc != 0:
        sys.exit(f"REGENIE rescue step2 failed with exit code {rc}")

    rescue_path = rescue_prefix.with_suffix(".regenie")
    n_rescue = count_spurious_hits(rescue_path)
    delta_pct = 100 * (n_baseline - n_rescue) / max(1, n_baseline)
    print(f"  ✓ rescue spurious significant SNPs: {n_rescue:,}")
    print(f"  Δ (baseline → rescue): {n_baseline} → {n_rescue} "
          f"({delta_pct:+.1f}% reduction)")

    (RESULTS_DIR / "spurious_counts.json").write_text(json.dumps({
        "baseline": int(n_baseline),
        "rescue": int(n_rescue),
        "delta_pct_reduction": float(delta_pct),
        "n_motif_covariates": len(pair_indices),
        "fdr_threshold": args.fdr_threshold,
        "lambda_var": args.lambda_var,
    }, indent=2))


if __name__ == "__main__":
    main()
