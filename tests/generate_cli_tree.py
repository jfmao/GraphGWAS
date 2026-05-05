"""Visualise the GraphGWAS CLI as a compact tree with command descriptions.

Produces Supplementary Figure S3 for the accompanying Nature Genetics
manuscript. Structure follows the GraphPop sfig_cli_tree.py template:
a vertical tree with colour-coded functional domains on the left and
a right-aligned description column.

Regenerate with:
    python tests/generate_cli_tree.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch


plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
})


# ===================================================================
# Command catalogue: organised by functional domain
# ===================================================================

DOMAINS = [
    {
        "name": "Fine-mapping\n(this paper)",
        "color": "#C62828", "bg": "#FFCDD2",
        "tier": "benchmarked",
        "functions": [
            ("Per-locus", [
                ("finemap", "GAFM / HBP / CLGF / GLEM fine-mapping on a locus"),
                ("finemap --source panukb", "stream Pan-UKB sumstats directly"),
            ]),
            ("Gene-level", [
                ("assoc mpat", "Multi-marker gene-level association test"),
            ]),
        ],
    },
    {
        "name": "Single-locus\n& QC",
        "color": "#1565C0", "bg": "#BBDEFB",
        "tier": "infrastructure",
        "functions": [
            ("GWAS", [
                ("assoc top-hits", "Top signals from a GWAS run"),
                ("phenotype activate", "Activate a trait for scanning"),
            ]),
            ("QC", [
                ("qc", "HWE, MAC, missingness, het-rate"),
            ]),
        ],
    },
    {
        "name": "Population\nstructure",
        "color": "#2E7D32", "bg": "#C8E6C9",
        "tier": "infrastructure",
        "functions": [
            ("Spectral PCs", [
                ("popstruct spectral-pcs", "Graph-spectral PCs from rare variants"),
                ("spectral correction", "Phenotype correction via spectral PCs"),
            ]),
        ],
    },
    {
        "name": "Multi-locus\n(paper #2)",
        "color": "#EF6C00", "bg": "#FFE0B2",
        "tier": "preview",
        "functions": [
            ("Epistasis", [
                ("epistasis scan", "LPCE LD-pruned co-occurrence (further methods forthcoming)"),
            ]),
            ("Pathway flow", [
                ("flow architecture", "Max-flow disease architecture"),
            ]),
        ],
    },
    {
        "name": "Heritability\n(platform)",
        "color": "#6A1B9A", "bg": "#E1BEE7",
        "tier": "platform",
        "functions": [
            ("Six estimators", [
                ("heritability spectral", "Laplacian-eigenspectrum h²"),
                ("heritability conductance", "Graph-conductance h²"),
                ("heritability flow", "Pathway-decomposed h²"),
                ("heritability gnn", "GNN-prediction h²"),
                ("heritability multi", "Variant/gene/pathway-scale h²"),
                ("heritability report", "Unified h² report"),
            ]),
        ],
    },
    {
        "name": "Cross-trait\n& pleiotropy\n(platform)",
        "color": "#00695C", "bg": "#B2DFDB",
        "tier": "platform",
        "functions": [
            ("Genetic correlation", [
                ("multivariate correlation", "r_G between two traits"),
                ("multivariate g-matrix", "T×T G-matrix"),
                ("multivariate multi-resolution", "r_G at variant/gene/pathway"),
            ]),
            ("Pleiotropy", [
                ("multivariate pleiotropy", "Shared genes between runs"),
                ("multivariate psi", "Gene Pleiotropy Strength Index"),
                ("multivariate coherence", "Frequency-resolved coherence"),
                ("multivariate report", "Full multivariate report"),
            ]),
        ],
    },
    {
        "name": "G × E / MET\n(platform)",
        "color": "#558B2F", "bg": "#DCEDC8",
        "tier": "platform",
        "functions": [
            ("Environment analysis", [
                ("met status", "Trial summary"),
                ("met env-similarity", "Genetic r_G across environments"),
                ("met mega-env", "Mega-environment discovery"),
                ("met gxe", "G×E variance decomposition"),
                ("met diagnostics", "Design graph diagnostics"),
            ]),
            ("Reaction norms", [
                ("met reaction-norm", "Variant-level reaction norms"),
                ("met impute", "Graph-neighbour missing-cell imputation"),
                ("met report", "Comprehensive MET report"),
            ]),
        ],
    },
    {
        "name": "PRS & MR\n(platform)",
        "color": "#AD1457", "bg": "#F8BBD0",
        "tier": "platform",
        "functions": [
            ("Polygenic risk score", [
                ("prs classical", "Weighted-sum PRS"),
                ("prs pathway", "Pathway-partitioned PRS"),
                ("prs report", "Compare all PRS methods"),
            ]),
            ("Mendelian randomisation", [
                ("mr report", "IVW + Egger + weighted median"),
            ]),
        ],
    },
    {
        "name": "GNN & agent\n(platform)",
        "color": "#4E342E", "bg": "#D7CCC8",
        "tier": "platform",
        "functions": [
            ("Graph neural net", [
                ("gnn export", "Export subgraph to PyTorch Geometric"),
            ]),
            ("LLM interface", [
                ("agent", "LangGraph natural-language interface"),
                ("interpret", "Deterministic post-hoc interpretation"),
            ]),
        ],
    },
    {
        "name": "Interfaces\n& IO",
        "color": "#424242", "bg": "#E0E0E0",
        "tier": "infrastructure",
        "functions": [
            ("Server", [
                ("serve", "FastAPI REST server (37 endpoints)"),
                ("mcp", "MCP server (16 tools, MCP-compatible clients)"),
            ]),
            ("Plot", [
                ("plot manhattan", "Manhattan plot from results"),
                ("plot qq", "QQ plot with λ_GC"),
            ]),
            ("Results", [
                ("results export", "TSV export (GNExT / PheWeb compatible)"),
                ("results delete", "Delete all results for a GWAS run"),
            ]),
        ],
    },
]


TIER_LEGEND = [
    ("benchmarked", "Rigorously benchmarked in this paper", "#C62828"),
    ("preview", "Previewed here; full benchmark in paper #2", "#EF6C00"),
    ("platform", "Implemented; awaiting dedicated benchmarks", "#6A1B9A"),
    ("infrastructure", "Shared infrastructure", "#424242"),
]


# ===================================================================
# Layout
# ===================================================================

def main() -> None:
    # Geometry (enlarged from the 6pt original)
    line_h = 0.32
    func_gap = 0.14
    domain_gap = 0.36
    left_col_w = 3.0
    cmd_col_w = 4.8
    desc_col_w = 6.0

    # Compute total height
    total_h = 0.0
    for d in DOMAINS:
        n_lines = sum(len(cmds) for _, cmds in d["functions"])
        n_funcs = len(d["functions"])
        total_h += n_lines * line_h + (n_funcs - 1) * func_gap + domain_gap
    total_h += 2.6  # title + legend (enlarged layout)

    fig_w = left_col_w + cmd_col_w + desc_col_w + 0.8
    fig = plt.figure(figsize=(fig_w, total_h))
    ax = fig.add_axes((0.01, 0.01, 0.98, 0.98))
    ax.set_xlim(0, fig_w)
    ax.set_ylim(0, total_h)
    ax.axis("off")

    # Title
    n_cmds = sum(len(c) for d in DOMAINS for _, cmds in d["functions"] for c in [cmds])
    n_funcs = sum(len(d["functions"]) for d in DOMAINS)
    ax.text(fig_w / 2, total_h - 0.25, "GraphGWAS CLI command hierarchy",
            ha="center", va="top", fontsize=17, fontweight="bold", color="#222")
    ax.text(fig_w / 2, total_h - 0.65,
            f"{len(DOMAINS)} functional domains  ·  {n_funcs} groups  ·  "
            f"{n_cmds} commands  ·  entry point: graphgwas …",
            ha="center", va="top", fontsize=11, color="#666")

    # Legend — horizontal row beneath the title, not in corner
    legend_y = total_h - 1.4
    n_legend = len(TIER_LEGEND)
    legend_spacing = fig_w / (n_legend + 1)
    for i, (tier, label, color) in enumerate(TIER_LEGEND):
        legend_x = (i + 1) * legend_spacing
        ax.plot(legend_x - 0.25, legend_y, "s", color=color, markersize=10)
        ax.text(legend_x - 0.1, legend_y, label, fontsize=9.5, va="center", color="#333")

    # Domain rendering
    y = total_h - 2.1
    for d in DOMAINS:
        n_lines = sum(len(cmds) for _, cmds in d["functions"])
        n_funcs = len(d["functions"])
        domain_h = n_lines * line_h + (n_funcs - 1) * func_gap
        block_top = y
        block_bot = y - domain_h

        # Left-column domain box
        box = FancyBboxPatch(
            (0.2, block_bot - 0.05), left_col_w - 0.3, domain_h + 0.1,
            boxstyle="round,pad=0.02",
            linewidth=0.8, edgecolor=d["color"], facecolor=d["bg"],
            alpha=0.55, zorder=1,
        )
        ax.add_patch(box)
        ax.text(0.2 + (left_col_w - 0.3) / 2, (block_top + block_bot) / 2,
                d["name"], ha="center", va="center", fontsize=10.5,
                fontweight="bold", color=d["color"])

        # Function groups
        cur_y = block_top
        for f_name, cmds in d["functions"]:
            # Function label (subtle, above commands)
            ax.text(left_col_w + 0.05, cur_y - 0.05, f_name,
                    ha="left", va="top", fontsize=9.5, color=d["color"],
                    fontweight="bold", style="italic")
            # Commands
            for i, (cmd, desc) in enumerate(cmds):
                cmd_y = cur_y - 0.05 - (i + 1) * line_h
                ax.text(left_col_w + 0.15, cmd_y, cmd,
                        ha="left", va="center",
                        fontsize=9, family="monospace", color="#222")
                ax.text(left_col_w + cmd_col_w + 0.1, cmd_y, desc,
                        ha="left", va="center", fontsize=8.5, color="#555")
            cur_y = cur_y - 0.05 - len(cmds) * line_h - func_gap

        y = block_bot - domain_gap

    # Footnote
    ax.text(fig_w / 2, 0.25,
            "All commands share the same graph-database / BGEN / Pan-UKB input "
            "interface and write outputs as typed graph nodes "
            "(AssociationResult, CredibleSet) queryable via the database's "
            "native query language.",
            ha="center", va="bottom", fontsize=9, color="#666", style="italic")

    out_dir = Path("/mnt/data/GraphGWAS/paper/finemapping_v1/figures")
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out_dir / f"figS3_cli_tree.{ext}",
                    bbox_inches="tight", pad_inches=0.12)
    # Mirror to results/benchmark_v2
    alt_dir = Path("/mnt/data/GraphGWAS/results/benchmark_v2/paper_figures")
    alt_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(alt_dir / f"figS3_cli_tree.{ext}",
                    bbox_inches="tight", pad_inches=0.12)
    print(f"Wrote {out_dir}/figS3_cli_tree.{{png,pdf}}")


if __name__ == "__main__":
    main()
