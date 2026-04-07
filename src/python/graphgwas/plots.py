"""Manhattan and QQ plots + TSV export for GraphGWAS results."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .db import GraphGWASConnection


# Chromosome order and color scheme
CHR_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
CHR_COLORS = ["#1f77b4", "#aec7e8"]  # alternating blue shades


def manhattan_plot(conn: GraphGWASConnection, run_id: str, output_path: str,
                   title: str = "", figsize: tuple = (16, 6),
                   subsample_threshold: float = 0.05, subsample_frac: float = 0.01):
    """Generate Manhattan plot from stored AssociationResult nodes.

    Args:
        conn: database connection.
        run_id: GWAS run identifier.
        output_path: path for output PNG/PDF.
        title: plot title.
        subsample_threshold: p-values above this are subsampled for performance.
        subsample_frac: fraction to keep when subsampling.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Query results
    result = conn.execute_read(
        """
        MATCH (ar:AssociationResult {run_id: $run_id})-[:FOR_VARIANT]->(v:Variant)
        RETURN v.chr AS chr, v.pos AS pos, ar.p_value_log10 AS log10p, ar.p_value AS p
        ORDER BY v.chr, v.pos
        """,
        {"run_id": run_id},
    )
    data = result.data()
    if not data:
        print("No results found for run_id:", run_id)
        return

    # Organize by chromosome
    chr_data = {}
    for d in data:
        c = d["chr"]
        if c not in chr_data:
            chr_data[c] = {"pos": [], "log10p": [], "p": []}
        chr_data[c]["pos"].append(d["pos"])
        chr_data[c]["log10p"].append(d["log10p"])
        chr_data[c]["p"].append(d["p"])

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    x_offset = 0
    x_ticks = []
    x_labels = []

    for i, chrom in enumerate(CHR_ORDER):
        if chrom not in chr_data:
            continue

        pos = np.array(chr_data[chrom]["pos"])
        log10p = np.array(chr_data[chrom]["log10p"])
        p_vals = np.array(chr_data[chrom]["p"])

        x = pos - pos.min() + x_offset

        # Subsample non-significant points for performance
        if subsample_threshold > 0:
            sig_mask = p_vals < subsample_threshold
            nonsig_mask = ~sig_mask
            if np.sum(nonsig_mask) > 1000:
                np.random.seed(42)
                keep = np.random.random(np.sum(nonsig_mask)) < subsample_frac
                nonsig_x = x[nonsig_mask][keep]
                nonsig_y = log10p[nonsig_mask][keep]
            else:
                nonsig_x = x[nonsig_mask]
                nonsig_y = log10p[nonsig_mask]

            color = CHR_COLORS[i % 2]
            ax.scatter(nonsig_x, nonsig_y, c=color, s=1, alpha=0.5, rasterized=True)
            ax.scatter(x[sig_mask], log10p[sig_mask], c=color, s=4, alpha=0.8, rasterized=True)
        else:
            color = CHR_COLORS[i % 2]
            ax.scatter(x, log10p, c=color, s=1, alpha=0.5, rasterized=True)

        mid = x_offset + (pos.max() - pos.min()) / 2
        x_ticks.append(mid)
        x_labels.append(chrom.replace("chr", ""))
        x_offset += (pos.max() - pos.min()) + 1e6  # gap between chromosomes

    # Significance lines
    ax.axhline(y=-np.log10(5e-8), color="red", linestyle="--", linewidth=0.8, label="5×10⁻⁸")
    ax.axhline(y=-np.log10(1e-5), color="blue", linestyle="--", linewidth=0.5, alpha=0.5, label="1×10⁻⁵")

    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels, fontsize=8)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("-log₁₀(p)")
    ax.set_title(title or f"Manhattan Plot (run: {run_id})")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"Manhattan plot saved to {output_path}")


def qq_plot(conn: GraphGWASConnection, run_id: str, output_path: str,
            title: str = "", figsize: tuple = (7, 7)):
    """Generate QQ plot with genomic inflation factor lambda_GC.

    Args:
        conn: database connection.
        run_id: GWAS run identifier.
        output_path: path for output PNG/PDF.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    result = conn.execute_read(
        """
        MATCH (ar:AssociationResult {run_id: $run_id})
        RETURN ar.p_value AS p
        ORDER BY ar.p_value
        """,
        {"run_id": run_id},
    )
    p_values = np.array([r["p"] for r in result if r["p"] is not None and r["p"] > 0])

    if len(p_values) == 0:
        print("No p-values found for run_id:", run_id)
        return

    n = len(p_values)
    p_values.sort()

    # Expected p-values under null
    expected = (np.arange(1, n + 1) - 0.5) / n

    # Observed and expected -log10
    obs_log10 = -np.log10(p_values)
    exp_log10 = -np.log10(expected)

    # Genomic inflation factor
    chi2_obs = sp_stats_chi2_from_p(p_values)
    lambda_gc = np.median(chi2_obs) / 0.4549

    # 95% confidence band
    ci_upper = -np.log10(sp_stats_beta_ppf(0.025, np.arange(1, n+1), n - np.arange(1, n+1) + 1))
    ci_lower = -np.log10(sp_stats_beta_ppf(0.975, np.arange(1, n+1), n - np.arange(1, n+1) + 1))

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Confidence band
    ax.fill_between(exp_log10, ci_lower, ci_upper, color="lightgray", alpha=0.5, label="95% CI")

    # Identity line
    max_val = max(obs_log10.max(), exp_log10.max())
    ax.plot([0, max_val], [0, max_val], color="red", linewidth=0.8, linestyle="--")

    # Points
    ax.scatter(exp_log10, obs_log10, c="black", s=3, alpha=0.5, rasterized=True)

    ax.set_xlabel("Expected -log₁₀(p)")
    ax.set_ylabel("Observed -log₁₀(p)")
    ax.set_title(title or f"QQ Plot (run: {run_id})")
    ax.text(0.05, 0.95, f"λ_GC = {lambda_gc:.3f}", transform=ax.transAxes,
            fontsize=12, verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
    ax.legend(loc="lower right", fontsize=8)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"QQ plot saved to {output_path} (λ_GC = {lambda_gc:.3f})")


def sp_stats_chi2_from_p(p_values):
    """Convert p-values to chi-squared statistics (1 df)."""
    from scipy import stats as sp_stats
    return sp_stats.chi2.isf(p_values, 1)


def sp_stats_beta_ppf(q, a, b):
    """Beta distribution percent point function for QQ CI."""
    from scipy import stats as sp_stats
    return sp_stats.beta.ppf(q, a, b)


def results_to_tsv(conn: GraphGWASConnection, run_id: str, output_path: str):
    """Export association results to TSV (GNExT/PheWeb compatible format).

    Columns: variant, chr, pos, ref, alt, beta, se, p_value, af, n_case, n_control, method
    """
    result = conn.execute_read(
        """
        MATCH (ar:AssociationResult {run_id: $run_id})-[:FOR_VARIANT]->(v:Variant)
        RETURN ar.variant_id AS variant, v.chr AS chr, v.pos AS pos,
               v.ref AS ref, v.alt AS alt,
               ar.beta AS beta, ar.se AS se, ar.p_value AS p_value,
               ar.maf AS af, ar.n_cases AS n_case, ar.n_controls AS n_control,
               ar.method AS method
        ORDER BY v.chr, v.pos
        """,
        {"run_id": run_id},
    )
    data = result.data()

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t",
                                fieldnames=["variant", "chr", "pos", "ref", "alt",
                                            "beta", "se", "p_value", "af",
                                            "n_case", "n_control", "method"])
        writer.writeheader()
        for row in data:
            writer.writerow(row)

    print(f"Exported {len(data)} results to {output_path}")
