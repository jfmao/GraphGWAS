"""Statistical association engine for GraphGWAS.

Methods:
- Chi-squared allelic test (fast screening, MAC >= 20)
- Fisher's exact test (rare variants, MAC < 5)
- Logistic regression (binary traits, with covariates)
- Firth penalized logistic regression (rare variants, 5 <= MAC < 20)
- Linear regression (quantitative traits, with covariates)
- Single-locus genome scan (auto-selects method per variant)
- Variant z-score computation (foundation for MPAT)
"""

from __future__ import annotations

import uuid
from datetime import datetime, timezone

import numpy as np
from scipy import stats as sp_stats
from scipy.linalg import solve_triangular, cho_factor, cho_solve

from . import config as _cfg
from .db import GraphGWASConnection
from .genotype import (
    build_count_table,
    build_dosage,
    compute_z_score,
    count_table_to_allelic,
    get_all_indices,
    get_covariate_matrix,
    get_phenotype_indices,
    get_phenotype_values,
    variant_iterator,
)


# ---------------------------------------------------------------------------
# A1.1 Chi-Squared Allelic Test
# ---------------------------------------------------------------------------

def chi2_allelic_test(count_table: np.ndarray) -> dict:
    """Chi-squared test on 2x2 allelic table derived from 2x3 genotype table."""
    allelic = count_table_to_allelic(count_table)
    a, b = allelic[0]  # case REF, case ALT
    c, d = allelic[1]  # ctrl REF, ctrl ALT

    if b == 0 or c == 0:
        odds_ratio = float("inf") if (a * d) > 0 else 0.0
        log_or = float("nan")
    else:
        odds_ratio = (a * d) / (b * c)
        log_or = np.log(odds_ratio)

    n = a + b + c + d
    if n == 0:
        return _null_result("chi2")

    chi2, p_value = sp_stats.chi2_contingency(allelic, correction=True)[:2]

    se_log_or = float("nan")
    ci_lower, ci_upper = float("nan"), float("nan")
    if a > 0 and b > 0 and c > 0 and d > 0:
        se_log_or = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
        ci_lower = np.exp(log_or - 1.96 * se_log_or)
        ci_upper = np.exp(log_or + 1.96 * se_log_or)

    return {
        "chi2": float(chi2), "p_value": float(p_value),
        "odds_ratio": float(odds_ratio),
        "beta": float(log_or) if not np.isnan(log_or) else 0.0,
        "se": float(se_log_or),
        "ci_lower": float(ci_lower), "ci_upper": float(ci_upper),
        "method": "chi2",
    }


# ---------------------------------------------------------------------------
# A1.2 Fisher's Exact Test
# ---------------------------------------------------------------------------

def fisher_exact_test(count_table: np.ndarray) -> dict:
    """Fisher's exact test on 2x2 allelic table."""
    allelic = count_table_to_allelic(count_table)
    odds_ratio, p_value = sp_stats.fisher_exact(allelic, alternative="two-sided")

    a, b = allelic[0]
    c, d = allelic[1]
    log_or = np.log(odds_ratio) if odds_ratio > 0 and not np.isinf(odds_ratio) else 0.0
    se = float("nan")
    if a > 0 and b > 0 and c > 0 and d > 0:
        se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)

    return {
        "odds_ratio": float(odds_ratio), "p_value": float(p_value),
        "beta": float(log_or), "se": float(se), "method": "fisher",
    }


# ---------------------------------------------------------------------------
# A1.3 Logistic Regression
# ---------------------------------------------------------------------------

def logistic_regression(dosage: np.ndarray, phenotype: np.ndarray,
                        covariates: np.ndarray | None = None) -> dict:
    """Standard logistic regression: logit(P(Y=1)) = beta0 + beta1*G + covariates."""
    valid = ~np.isnan(dosage)
    if covariates is not None:
        valid &= ~np.any(np.isnan(covariates), axis=1)
    if np.sum(valid) < 10:
        return _null_result("logistic")

    g = dosage[valid]
    y = phenotype[valid]
    if np.std(g) == 0 or np.std(y) == 0:
        return _null_result("logistic")

    if covariates is not None:
        X = np.column_stack([np.ones(len(g)), g, covariates[valid]])
    else:
        X = np.column_stack([np.ones(len(g)), g])

    try:
        import statsmodels.api as sm
        model = sm.Logit(y, X)
        result = model.fit(disp=0, method="newton", maxiter=25)
        beta = result.params[1]
        se = result.bse[1]
        p_value = result.pvalues[1]
    except Exception:
        return _score_test_logistic(g, y)

    odds_ratio = np.exp(beta)
    return {
        "beta": float(beta), "se": float(se), "p_value": float(p_value),
        "odds_ratio": float(odds_ratio),
        "ci_lower": float(np.exp(beta - 1.96 * se)),
        "ci_upper": float(np.exp(beta + 1.96 * se)),
        "method": "logistic",
    }


def _score_test_logistic(g: np.ndarray, y: np.ndarray) -> dict:
    """Score test fallback when full logistic regression fails."""
    p_bar = np.mean(y)
    if p_bar == 0 or p_bar == 1:
        return _null_result("logistic_score")
    U = np.sum((y - p_bar) * g)
    V = p_bar * (1 - p_bar) * np.sum((g - np.mean(g)) ** 2)
    if V == 0:
        return _null_result("logistic_score")
    chi2 = U ** 2 / V
    p_value = float(sp_stats.chi2.sf(chi2, 1))
    beta = U / V
    se = 1.0 / np.sqrt(V)
    return {
        "beta": float(beta), "se": float(se), "p_value": p_value,
        "odds_ratio": float(np.exp(beta)), "method": "logistic_score",
    }


# ---------------------------------------------------------------------------
# A1.5 Firth Penalized Logistic Regression (NEW — Phase 2)
# ---------------------------------------------------------------------------

def _sigmoid(x: np.ndarray) -> np.ndarray:
    """Numerically stable sigmoid."""
    return np.where(x >= 0, 1 / (1 + np.exp(-x)), np.exp(x) / (1 + np.exp(x)))


def firth_logistic_regression(dosage: np.ndarray, phenotype: np.ndarray,
                              covariates: np.ndarray | None = None,
                              max_iter: int = 25, tol: float = 1e-6) -> dict:
    """Firth penalized logistic regression.

    Penalized log-likelihood: L*(beta) = L(beta) + 0.5 * log|I(beta)|
    Uses IRLS with modified score incorporating hat matrix diagonal.
    P-value from penalized likelihood ratio test.
    """
    valid = ~np.isnan(dosage)
    if covariates is not None:
        valid &= ~np.any(np.isnan(covariates), axis=1)
    if np.sum(valid) < 10:
        return _null_result("firth")

    g = dosage[valid]
    y = phenotype[valid]
    n = len(y)

    if np.std(g) == 0:
        return _null_result("firth")

    # Design matrix: [intercept, genotype, covariates...]
    if covariates is not None:
        X = np.column_stack([np.ones(n), g, covariates[valid]])
    else:
        X = np.column_stack([np.ones(n), g])
    p = X.shape[1]

    # --- Fit full model ---
    beta = np.zeros(p)
    converged = False

    for iteration in range(max_iter):
        eta = X @ beta
        mu = _sigmoid(eta)
        mu = np.clip(mu, 1e-10, 1 - 1e-10)
        W = mu * (1 - mu)
        W_sqrt = np.sqrt(W)

        # Fisher information: X^T W X
        XtWX = (X * W[:, None]).T @ X

        try:
            L, low = cho_factor(XtWX)
        except np.linalg.LinAlgError:
            return _null_result("firth")

        # Hat diagonal: h_i = diag(W^{1/2} X (X^T W X)^{-1} X^T W^{1/2})
        # Efficient: Q = L^{-1} (X * sqrt(W))^T, h = colsum(Q^2)
        XW = (X * W_sqrt[:, None]).T  # (p, n)
        Q = solve_triangular(L, XW, lower=low)
        h = np.sum(Q ** 2, axis=0)  # hat diagonal, length n

        # Penalized score: U_j = sum_i x_ij * (y_i - mu_i + h_i * (0.5 - mu_i))
        working_residual = y - mu + h * (0.5 - mu)
        U = X.T @ working_residual

        # Newton step: delta = (X^T W X)^{-1} U
        delta = cho_solve((L, low), U)
        beta += delta

        if np.max(np.abs(delta)) < tol:
            converged = True
            break

    if not converged:
        return _null_result("firth")

    # Penalized log-likelihood at full model
    eta_full = X @ beta
    mu_full = _sigmoid(eta_full)
    mu_full = np.clip(mu_full, 1e-10, 1 - 1e-10)
    ll_full = (np.sum(y * np.log(mu_full) + (1 - y) * np.log(1 - mu_full))
               + 0.5 * np.log(np.linalg.det((X * (mu_full * (1 - mu_full))[:, None]).T @ X)))

    # --- Fit null model (beta_1 = 0) for LRT ---
    X_null = np.delete(X, 1, axis=1)  # remove genotype column
    beta_null = np.zeros(X_null.shape[1])

    for iteration in range(max_iter):
        eta_n = X_null @ beta_null
        mu_n = _sigmoid(eta_n)
        mu_n = np.clip(mu_n, 1e-10, 1 - 1e-10)
        W_n = mu_n * (1 - mu_n)

        XtWX_n = (X_null * W_n[:, None]).T @ X_null
        try:
            L_n, low_n = cho_factor(XtWX_n)
        except np.linalg.LinAlgError:
            break

        W_sqrt_n = np.sqrt(W_n)
        Q_n = solve_triangular(L_n, (X_null * W_sqrt_n[:, None]).T, lower=low_n)
        h_n = np.sum(Q_n ** 2, axis=0)

        U_n = X_null.T @ (y - mu_n + h_n * (0.5 - mu_n))
        delta_n = cho_solve((L_n, low_n), U_n)
        beta_null += delta_n

        if np.max(np.abs(delta_n)) < tol:
            break

    eta_null = X_null @ beta_null
    mu_null = _sigmoid(eta_null)
    mu_null = np.clip(mu_null, 1e-10, 1 - 1e-10)
    ll_null = (np.sum(y * np.log(mu_null) + (1 - y) * np.log(1 - mu_null))
               + 0.5 * np.log(np.linalg.det((X_null * (mu_null * (1 - mu_null))[:, None]).T @ X_null)))

    # Penalized LRT
    lr_stat = 2 * (ll_full - ll_null)
    lr_stat = max(lr_stat, 0.0)
    p_value = float(sp_stats.chi2.sf(lr_stat, 1))

    # Wald SE from Fisher information
    W_full = mu_full * (1 - mu_full)
    XtWX_full = (X * W_full[:, None]).T @ X
    try:
        var_beta = np.linalg.inv(XtWX_full)
        se = np.sqrt(var_beta[1, 1])
    except np.linalg.LinAlgError:
        se = float("nan")

    beta_g = beta[1]
    odds_ratio = np.exp(beta_g)

    return {
        "beta": float(beta_g), "se": float(se), "p_value": p_value,
        "odds_ratio": float(odds_ratio),
        "ci_lower": float(np.exp(beta_g - 1.96 * se)) if not np.isnan(se) else float("nan"),
        "ci_upper": float(np.exp(beta_g + 1.96 * se)) if not np.isnan(se) else float("nan"),
        "method": "firth", "converged": converged, "n_iter": iteration + 1,
    }


# ---------------------------------------------------------------------------
# A1.4 Linear Regression
# ---------------------------------------------------------------------------

def linear_regression(dosage: np.ndarray, phenotype: np.ndarray,
                      covariates: np.ndarray | None = None) -> dict:
    """Linear regression: Y = beta0 + beta1*G + covariates + epsilon."""
    valid = ~np.isnan(dosage) & ~np.isnan(phenotype)
    if covariates is not None:
        valid &= ~np.any(np.isnan(covariates), axis=1)
    if np.sum(valid) < 10:
        return _null_result("linear")

    g = dosage[valid]
    y = phenotype[valid]
    if np.std(g) == 0:
        return _null_result("linear")

    if covariates is not None:
        X = np.column_stack([np.ones(len(g)), g, covariates[valid]])
    else:
        X = np.column_stack([np.ones(len(g)), g])

    try:
        beta_hat, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta_hat
        n, p = X.shape
        mse = np.sum(resid ** 2) / (n - p)
        var_beta = mse * np.linalg.inv(X.T @ X)
        se = np.sqrt(np.diag(var_beta))
        t_stat = beta_hat[1] / se[1]
        p_value = float(2 * sp_stats.t.sf(abs(t_stat), n - p))
        ss_res = np.sum(resid ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    except Exception:
        return _null_result("linear")

    return {
        "beta": float(beta_hat[1]), "se": float(se[1]), "p_value": p_value,
        "r_squared": float(r_squared), "method": "linear",
    }


# ---------------------------------------------------------------------------
# Single-locus genome scan (Phase 2: Firth, covariates, streaming writes)
# ---------------------------------------------------------------------------

def single_locus_scan(conn: GraphGWASConnection, chr: str,
                      start: int | None = None, end: int | None = None,
                      method: str = "auto",
                      covariate_names: list[str] | None = None,
                      store_results_fn=None, run_id: str | None = None,
                      phenotype_key: str | None = None,
                      verbose: bool = True) -> list[dict]:
    """Run single-locus association for all variants in a region.

    Args:
        conn: database connection.
        chr: chromosome (e.g., "chr19").
        start, end: optional position range.
        method: 'auto', 'chi2', 'fisher', 'logistic', 'firth', 'linear'.
        covariate_names: property names for covariates (e.g., ['pca_1','pca_2','sex']).
        store_results_fn: if provided, callable(conn, batch, run_id, phenotype_key) for streaming writes.
        run_id: GWAS run identifier.
        phenotype_key: active phenotype name.
        verbose: print progress.

    Returns:
        list of result dicts.
    """
    if method == "linear":
        # Quantitative trait: use all samples with gwas_value
        all_idx = get_all_indices(conn)
        pheno = get_phenotype_values(conn, all_idx)
        # Filter to samples that have non-NaN phenotype values
        valid = ~np.isnan(pheno)
        all_idx = all_idx[valid]
        pheno = pheno[valid]
        n_case, n_ctrl = len(all_idx), 0  # no case/control distinction
        case_idx, ctrl_idx = all_idx, np.array([], dtype=np.int32)
    else:
        case_idx, ctrl_idx = get_phenotype_indices(conn)
        n_case, n_ctrl = len(case_idx), len(ctrl_idx)
        all_idx = np.concatenate([case_idx, ctrl_idx])
        pheno = np.concatenate([np.ones(n_case), np.zeros(n_ctrl)])

    # Build covariate matrix (once)
    covariates = None
    if covariate_names:
        covariates = get_covariate_matrix(conn, covariate_names, all_idx)

    if verbose:
        region = f"{chr}:{start}-{end}" if start else chr
        cov_str = f", covariates={covariate_names}" if covariate_names else ""
        print(f"Scanning {region} | {n_case} cases, {n_ctrl} controls | method={method}{cov_str}")

    results = []
    write_buffer = []
    n_scanned = 0
    WRITE_BATCH = 5000

    for variant in variant_iterator(conn, chr, start, end):
        gt_packed = variant["gt_packed"]
        if gt_packed is None:
            continue

        table = build_count_table(gt_packed, case_idx, ctrl_idx, _cfg.N_SAMPLES)
        n_scanned += 1

        # MAC computation
        ac_case = table[0, 1] + 2 * table[0, 2]
        ac_ctrl = table[1, 1] + 2 * table[1, 2]
        an_case = 2 * int(table[0].sum())
        an_ctrl = 2 * int(table[1].sum())
        mac_case = min(ac_case, an_case - ac_case) if an_case > 0 else 0
        mac_ctrl = min(ac_ctrl, an_ctrl - ac_ctrl) if an_ctrl > 0 else 0
        mac = min(mac_case, mac_ctrl)
        af = (ac_case + ac_ctrl) / (an_case + an_ctrl) if (an_case + an_ctrl) > 0 else 0

        # Method selection
        if method == "auto":
            if mac >= 20 and np.all(table >= 5):
                result = chi2_allelic_test(table)
            elif mac >= 5:
                # Firth for intermediate MAC
                dosage = build_dosage(gt_packed, all_idx, _cfg.N_SAMPLES)
                result = firth_logistic_regression(dosage, pheno, covariates)
            else:
                result = fisher_exact_test(table)
        elif method == "chi2":
            result = chi2_allelic_test(table)
        elif method == "fisher":
            result = fisher_exact_test(table)
        elif method == "logistic":
            dosage = build_dosage(gt_packed, all_idx, _cfg.N_SAMPLES)
            result = logistic_regression(dosage, pheno, covariates)
        elif method == "firth":
            dosage = build_dosage(gt_packed, all_idx, _cfg.N_SAMPLES)
            result = firth_logistic_regression(dosage, pheno, covariates)
        elif method == "linear":
            dosage = build_dosage(gt_packed, all_idx, _cfg.N_SAMPLES)
            result = linear_regression(dosage, pheno, covariates)
        else:
            raise ValueError(f"Unknown method: {method}")

        # Add variant metadata
        result["variantId"] = variant["variantId"]
        result["chr"] = variant["chr"]
        result["pos"] = variant["pos"]
        result["ref"] = variant["ref"]
        result["alt"] = variant["alt"]
        result["af"] = float(af)
        result["mac"] = int(mac)
        result["n_cases"] = n_case
        result["n_controls"] = n_ctrl
        result["p_value_log10"] = -np.log10(result["p_value"]) if result["p_value"] > 0 else 300.0

        results.append(result)

        # Streaming writes
        if store_results_fn and run_id:
            write_buffer.append(result)
            if len(write_buffer) >= WRITE_BATCH:
                store_results_fn(conn, write_buffer, run_id, phenotype_key)
                write_buffer = []

        if verbose and n_scanned % 10000 == 0:
            print(f"  ...{n_scanned} variants scanned", flush=True)

    # Flush remaining writes
    if store_results_fn and run_id and write_buffer:
        store_results_fn(conn, write_buffer, run_id, phenotype_key)

    if verbose:
        sig = sum(1 for r in results if r["p_value"] < 5e-8)
        print(f"Done: {n_scanned} variants, {sig} genome-wide significant")

    return results


# ---------------------------------------------------------------------------
# MPAT z-score computation
# ---------------------------------------------------------------------------

def compute_variant_z_scores(conn: GraphGWASConnection, chr: str,
                             start: int | None = None,
                             end: int | None = None) -> list[dict]:
    """Compute per-variant GWAS z-scores for MPAT."""
    case_idx, ctrl_idx = get_phenotype_indices(conn)
    results = []
    for variant in variant_iterator(conn, chr, start, end):
        gt_packed = variant["gt_packed"]
        if gt_packed is None:
            continue
        table = build_count_table(gt_packed, case_idx, ctrl_idx, _cfg.N_SAMPLES)
        z = compute_z_score(table, len(case_idx), len(ctrl_idx))
        ac_case = table[0, 1] + 2 * table[0, 2]
        an_case = 2 * int(table[0].sum())
        ac_ctrl = table[1, 1] + 2 * table[1, 2]
        an_ctrl = 2 * int(table[1].sum())
        results.append({
            "variantId": variant["variantId"], "chr": variant["chr"],
            "pos": variant["pos"], "z_score": float(z),
            "af_case": ac_case / an_case if an_case > 0 else 0,
            "af_ctrl": ac_ctrl / an_ctrl if an_ctrl > 0 else 0,
            "mac_case": min(ac_case, an_case - ac_case) if an_case > 0 else 0,
        })
    return results


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _null_result(method: str) -> dict:
    return {
        "beta": 0.0, "se": float("nan"), "p_value": 1.0,
        "odds_ratio": 1.0, "method": method,
    }


def generate_run_id() -> str:
    ts = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    return f"gwas_{ts}_{uuid.uuid4().hex[:8]}"
