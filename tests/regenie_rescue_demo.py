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

BGEN_DIR = REPO_ROOT / "tests" / "data" / "human" / "1kGP_bgen"
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
    sample_ids: list[str],
    phenotype: np.ndarray,
    graph_cache_path: Path,
    *,
    fdr_threshold: float = 0.05,
) -> list[tuple[int, int, float]]:
    """Detect M2 (motif-filtered) interactions without Neo4j.

    Uses graph_cache JSON file (gene-membership annotations) to enumerate
    same-gene pairs (motif P1), then fits y ~ G1 + G2 + G1*G2 and reports
    pairs with BH-FDR < fdr_threshold.

    [SCAFFOLD] full implementation deferred — for the rescue demo we
    currently use random-pair detection as a placeholder.  Wiring to
    epistasis_v2.motif_filtered_epistasis() requires building the
    Cypher-equivalent same-gene index from graph_cache, then calling
    the existing _test_interaction() helper.

    TODO(paper2-Y2):
    - Adapt epistasis_v2.motif_filtered_epistasis to accept graph_cache
      JSON instead of Neo4j connection.
    - Wire phenotype + dosages here.
    - Return list of (i, j, p_interaction) tuples.

    Returns:
        Empty list (placeholder).  Once implemented, returns
        [(i, j, p), ...] for each motif pair surviving FDR cutoff.
    """
    raise NotImplementedError(
        "M2 motif detection (no-Neo4j path) not yet implemented. "
        "See TODO in tests/regenie_rescue_demo.py:detect_motif_interactions_no_neo4j. "
        "Plan: lift epistasis_v2.motif_filtered_epistasis() into a "
        "graph_cache-driven variant, then re-enable this function."
    )


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
    print("\n[2/4] Running REGENIE step 1 + step 2...")
    if not regenie_available():
        print("  ⚠ regenie not installed — skipping. To enable:")
        print("    mamba install -n graphmana -c bioconda regenie")
        print("  Stages 3+4 also skipped.")
        return

    bgen_path = BGEN_DIR / "chr22.bgen"
    sample_path = BGEN_DIR / "chr22.sample"
    step1_prefix = RESULTS_DIR / "regenie_step1"
    step2_prefix = RESULTS_DIR / "regenie_step2_baseline"

    rc = run_regenie_step1(bgen_path, sample_path, pheno_path, step1_prefix)
    if rc != 0:
        sys.exit(f"REGENIE step1 failed with exit code {rc}")
    rc = run_regenie_step2(bgen_path, sample_path, pheno_path,
                           step1_prefix.with_name(step1_prefix.name + "_pred.list"),
                           step2_prefix)
    if rc != 0:
        sys.exit(f"REGENIE step2 failed with exit code {rc}")

    baseline_path = step2_prefix.with_suffix(".regenie")
    n_baseline = count_spurious_hits(baseline_path)
    print(f"  ✓ baseline spurious significant SNPs: {n_baseline:,}")

    if args.regenie_only:
        return

    # ----- Stage 3: M2 motif detection -----
    print("\n[3/4] Running M2 motif-pair detection...")
    print("  ⚠ M2 no-Neo4j wrapper is a TODO — see "
          "detect_motif_interactions_no_neo4j() docstring.")
    print("  This stage is intentionally a scaffold.  Once the wrapper "
          "exists, it will:")
    print("    - call detect_motif_interactions_no_neo4j(dos, sample_ids, "
          "sim['y'][:,0], CACHE_PATH)")
    print("    - return motif pairs + their interaction p-values")
    print("    - emit a covariate file: per-pair G1*G2 dosage product columns")
    print("  STOPPING here. Re-run after Stage 3 wrapper is implemented.")
    return

    # ----- Stage 4: rescue + comparison -----
    # (Not reached until Stage 3 is implemented.)


if __name__ == "__main__":
    main()
