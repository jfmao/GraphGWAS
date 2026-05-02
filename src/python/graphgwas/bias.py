"""GWAS test-statistic bias under omitted interaction terms.

Implements the closed-form bias-quantification framework from Yelmen
et~al.\ 2026 (bioRxiv 10.1101/2025.11.21.689603). When the true
data-generating process contains an interaction term that the
fitted single-marker LMM omits, the null distribution of the
$t$-statistic is shifted in mean and variance:

    t  ~  N(mu, 1 / sigma_res^2)         under H0: alpha = 0
    mu          = rho * sqrt(lambda * n) / sqrt(1 - lambda * rho^2)
    sigma_res^2 = (1 - lambda * rho^2) / (1 - lambda)

where:
    rho        = correlation between target SNP g and the realised
                 interaction signal u in col(Z), bounded above by
                 rho_max  = ||P_H g|| / ||g||  with P_H the orthogonal
                 projector onto col(Z).
    lambda     = Var(u) / (1 + Var(u)),  the fraction of phenotypic
                 variance from interactions, in [0, 1).
    n          = sample size.

The conservativeness ratio R(x) = p_true(|x|) / p_nominal(|x|)
quantifies how anti-conservative (R > 1) or conservative (R < 1) the
nominal p-value is at threshold |x|.  At the canonical GWAS
significance threshold |x| = 5.45 (p_nominal = 5e-8), Yelmen et~al.\
show R > 1 across virtually all realistic (n, lambda, rho).

This module's contribution beyond Yelmen et~al.\ is the substitution
of a graph-typed Z (built from M2 motif-pair products) for the random
Z they analysed.  See ``graphgwas.epistasis_v2.motif_interaction_matrix``.
"""
from __future__ import annotations

import numpy as np
from scipy import stats as sp_stats


# ---------------------------------------------------------------------------
# Preprocessing
# ---------------------------------------------------------------------------

def preprocess(g: np.ndarray,
               Z: np.ndarray,
               X: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Residualise g and Z on covariates X, then mean-center.

    Mirrors the simplification applied in Yelmen et~al.\ 2026 §"Practical
    computation of rho_max" when no GRM is structured-out: the operator
    T = M = I - X (X^T X)^{-1} X^T is the orthogonal projector that
    removes the covariate subspace.  Working in the residualised space
    means the noise covariance is (approximately) spherical and the
    closed-form rho_max applies directly.

    Args:
        g: target SNP dosage vector, shape (n,).
        Z: interaction-feature matrix, shape (n, m).
        X: covariate matrix, shape (n, q), optional.  If omitted, only
           mean-centring is applied.

    Returns:
        (g_resid, Z_resid) with the covariate subspace removed and
        mean zero per column.
    """
    g = np.asarray(g, dtype=float).ravel()
    Z = np.asarray(Z, dtype=float)
    if Z.ndim == 1:
        Z = Z.reshape(-1, 1)
    if X is not None:
        X = np.asarray(X, dtype=float)
        if X.ndim == 1:
            X = X.reshape(-1, 1)
        # Solve XtX b = Xt y for residualisation; numerically stable least-squares
        Xt_g, *_ = np.linalg.lstsq(X, g, rcond=None)
        g = g - X @ Xt_g
        Xt_Z, *_ = np.linalg.lstsq(X, Z, rcond=None)
        Z = Z - X @ Xt_Z
    g = g - g.mean()
    Z = Z - Z.mean(axis=0, keepdims=True)
    return g, Z


# ---------------------------------------------------------------------------
# rho_max  (Yelmen et al. 2026 Eq. 12)
# ---------------------------------------------------------------------------

def rho_max(g: np.ndarray,
            Z: np.ndarray,
            X: np.ndarray | None = None,
            preprocess_inputs: bool = True) -> float:
    """Maximum |correlation(g, u)| over u in col(Z).

    Yelmen et al. 2026 Eq. 12::

        rho_max  =  ||P_H g|| / ||g||

    where P_H = Z (Z^T Z)^+ Z^T is the orthogonal projector onto col(Z).
    Computed in practice by solving the least-squares problem
    b_hat = argmin_b ||Z b - g||^2 (avoiding explicit pseudo-inverse).

    Args:
        g: target SNP, shape (n,).
        Z: interaction-feature matrix, shape (n, m).
        X: optional covariates to residualise on.
        preprocess_inputs: if True (default), apply :func:`preprocess`
            before computing rho_max.  If False, assume g and Z are
            already in the preprocessed space.

    Returns:
        rho_max in [0, 1].
    """
    if preprocess_inputs:
        g, Z = preprocess(g, Z, X)
    g_norm = float(np.linalg.norm(g))
    if g_norm == 0.0:
        return 0.0
    b_hat, *_ = np.linalg.lstsq(Z, g, rcond=None)
    proj_norm = float(np.linalg.norm(Z @ b_hat))
    return min(1.0, proj_norm / g_norm)


# ---------------------------------------------------------------------------
# Bias quantities  (Yelmen et al. 2026 Eq. 10, 11)
# ---------------------------------------------------------------------------

def mu_shift(rho: float, lambda_var: float, n: int) -> float:
    """Mean shift of the t-statistic null under omitted interactions.

    Yelmen et al. 2026 Eq. 10::

        mu  =  rho * sqrt(lambda_var * n) / sqrt(1 - lambda_var * rho^2)

    Args:
        rho: target-SNP / interaction-subspace correlation.
        lambda_var: variance fraction from interactions, in [0, 1).
        n: sample size.
    """
    if not (0.0 <= lambda_var < 1.0):
        raise ValueError(f"lambda_var must lie in [0, 1); got {lambda_var}")
    if abs(rho) >= 1.0:
        raise ValueError(f"|rho| must be < 1; got {rho}")
    if n <= 0:
        raise ValueError(f"n must be positive; got {n}")
    denom = np.sqrt(1.0 - lambda_var * rho * rho)
    return float(rho * np.sqrt(lambda_var * n) / denom)


def sigma_res_sq(rho: float, lambda_var: float) -> float:
    """Residual-variance multiplier of the t-statistic null.

    Yelmen et al. 2026 Eq. 10::

        sigma_res^2  =  (1 - lambda_var * rho^2) / (1 - lambda_var)

    Returns:
        Variance scale factor.  Var(t) under the true null = 1 / sigma_res_sq().
    """
    if not (0.0 <= lambda_var < 1.0):
        raise ValueError(f"lambda_var must lie in [0, 1); got {lambda_var}")
    return float((1.0 - lambda_var * rho * rho) / (1.0 - lambda_var))


def conservativeness_ratio(x: float,
                           n: int,
                           lambda_var: float,
                           rho: float) -> float:
    """R(x) = p_true(|x|) / p_nominal(|x|).

    Yelmen et al. 2026 Eq. 11.  R > 1 = anti-conservative
    (more spurious significance than nominal p-value suggests);
    R < 1 = conservative.

    Args:
        x: t-statistic value (sign ignored, two-sided).
        n: sample size.
        lambda_var: variance fraction from interactions.
        rho: target-SNP / interaction-subspace correlation.

    Returns:
        R(x).  When p_nominal underflows numerically, returns +inf if
        p_true > 0 else 1.0.
    """
    x_abs = abs(float(x))
    mu = mu_shift(rho, lambda_var, n)
    sigma_res = float(np.sqrt(sigma_res_sq(rho, lambda_var)))
    p_nominal = 2.0 * sp_stats.norm.cdf(-x_abs)
    p_true = (sp_stats.norm.cdf(-sigma_res * (x_abs - mu)) +
              sp_stats.norm.cdf(-sigma_res * (x_abs + mu)))
    if p_nominal <= 0.0:
        return float("inf") if p_true > 0.0 else 1.0
    return float(p_true / p_nominal)


# ---------------------------------------------------------------------------
# Convenience entry point
# ---------------------------------------------------------------------------

def bias_diagnostic(g: np.ndarray,
                    Z: np.ndarray,
                    n: int,
                    lambda_var: float,
                    x: float = 5.45,
                    X: np.ndarray | None = None) -> dict:
    """Compute (rho_max, mu, sigma_res, R(x)) for one target SNP and Z matrix.

    Args:
        g: target SNP.
        Z: interaction-feature matrix.
        n: sample size for the LMM analysis.
        lambda_var: assumed variance fraction from interactions.
        x: t-statistic threshold (default 5.45 ≈ p = 5e-8).
        X: optional covariates.

    Returns:
        dict with keys ``rho_max``, ``mu``, ``sigma_res``, ``R``.
    """
    rho = rho_max(g, Z, X=X)
    return {
        "rho_max": rho,
        "mu": mu_shift(rho, lambda_var, n),
        "sigma_res": float(np.sqrt(sigma_res_sq(rho, lambda_var))),
        "R": conservativeness_ratio(x, n, lambda_var, rho),
    }
