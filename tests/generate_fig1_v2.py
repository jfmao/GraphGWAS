"""Generate Fig 1: flat-vs-relational annotation priors in fine-mapping.

Two-panel figure contrasting the per-variant scalar prior used by
flat-prior fine-mappers (panel a) with the typed biological knowledge
graph used by GraphGWAS (panel b). Panel (b) is deliberately designed
to emphasise the *data model* (node types, edge types, real labels),
not the *algorithm* — Fig 2 handles the HBP message-passing schematic.

Panel (b) thus shows:
  - a heterogeneous knowledge graph: Variant, Gene, Tissue,
    Pathway, RegulatoryElement, PPI-partner nodes
  - real labels (FTO, IRX3, Adipose_Subcutaneous, energy_metabolism)
  - a highlighted traversal from a query variant out to its
    multi-omics neighbourhood
  - the credible-set output shown as a co-queryable graph object
    (with gene/tissue/pathway attached), not just a PIP bar chart

Writes to paper/manuscript_v1/figures/fig1_architecture.{png,pdf}.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle


# ===================================================================
# Styling
# ===================================================================

VARIANT_COLOR_CAUSAL = "#d62728"
VARIANT_COLOR_NEUTRAL = "#b0b0b0"
GENE_COLOR = "#2ca02c"
PATHWAY_COLOR = "#9467bd"
PPI_COLOR = "#8c564b"
EQTL_COLOR = "#ff7f0e"
TISSUE_COLOR = "#17becf"
CCRE_COLOR = "#e377c2"
ANNOT_COLOR = "#7f7f7f"
METHOD_BOX_COLOR = "#1f77b4"
OUTPUT_COLOR = "#d62728"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
})


# ===================================================================
# Panel A: flat-prior fine-mapping (PolyFun-style) — unchanged
# ===================================================================

def _draw_variant_row(ax, y, n=10, causal_idx=4, x_start=0.08, x_end=0.92,
                      label="variants"):
    xs = np.linspace(x_start, x_end, n)
    for i, x in enumerate(xs):
        c = VARIANT_COLOR_CAUSAL if i == causal_idx else VARIANT_COLOR_NEUTRAL
        ax.plot(x, y, "o", color=c, markersize=8, markeredgecolor="black",
                markeredgewidth=0.7, zorder=3)
    if label:
        ax.text(x_start - 0.03, y, label, ha="right", va="center",
                fontsize=8, color="#444")
    return xs


def _draw_method_box(ax, x, y, w, h, text, color=METHOD_BOX_COLOR):
    box = FancyBboxPatch((x - w/2, y - h/2), w, h,
                         boxstyle="round,pad=0.01",
                         linewidth=1.2, edgecolor=color,
                         facecolor=color, alpha=0.12, zorder=2)
    ax.add_patch(box)
    ax.text(x, y, text, ha="center", va="center",
            fontsize=9, color=color, fontweight="bold")


def draw_panel_a(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("a  Flat-prior fine-mapping (PolyFun / SuSiE family)",
                 loc="left", fontsize=11, fontweight="bold", pad=8)

    xs = _draw_variant_row(ax, y=0.82, label="variants")

    rng = np.random.default_rng(0)
    annot_heights = rng.uniform(0.2, 0.9, len(xs))
    annot_heights[4] = 0.7
    for x, h in zip(xs, annot_heights):
        ax.bar(x, 0.08 * h, width=0.03, bottom=0.67,
               color=ANNOT_COLOR, edgecolor="black", linewidth=0.3, zorder=3)
    ax.text(0.05, 0.71, "scalar\nprior $w_i$", ha="right", va="center",
            fontsize=8, color="#444")

    ax.annotate("", xy=(0.5, 0.52), xytext=(0.5, 0.64),
                arrowprops=dict(arrowstyle="->", lw=1.2, color="#444"))

    _draw_method_box(ax, 0.5, 0.46, 0.7, 0.10,
                     "Bayesian regression on  ($z$, $R$, $w$)",
                     color=METHOD_BOX_COLOR)
    ax.text(0.5, 0.38,
            "context collapsed to a single scalar per variant",
            ha="center", va="center", fontsize=8, color="#666",
            style="italic")

    ax.annotate("", xy=(0.5, 0.22), xytext=(0.5, 0.36),
                arrowprops=dict(arrowstyle="->", lw=1.2, color="#444"))

    xs_out = xs
    pips = np.array([0.10, 0.08, 0.18, 0.22, 0.45, 0.38,
                     0.18, 0.10, 0.06, 0.05])
    for i, (x, p) in enumerate(zip(xs_out, pips)):
        color = OUTPUT_COLOR if i == 4 else "#606060"
        ax.bar(x, p * 0.10, width=0.03, bottom=0.08, color=color,
               edgecolor="black", linewidth=0.4, zorder=3)
    ax.text(0.05, 0.13, "PIPs", ha="right", va="center", fontsize=8, color="#444")

    ax.text(0.5, 0.02,
            "credible set = top-PIP variants only\n"
            "(post-hoc enrichment required to attach biology)",
            ha="center", va="bottom", fontsize=7.5, color="#666",
            style="italic")


# ===================================================================
# Panel B: relational-prior fine-mapping — biological knowledge graph
# ===================================================================

def _node(ax, x, y, shape, color, size=14, label=None, sub=None,
          label_offset=(0, -0.035), halo=False, zorder=4):
    """Draw a typed graph node with optional label + sublabel."""
    if halo:
        ax.plot(x, y, shape, color=color, markersize=size + 8,
                markerfacecolor="none", markeredgecolor=color,
                markeredgewidth=1.2, alpha=0.35, zorder=zorder - 1)
    ax.plot(x, y, shape, color=color, markersize=size,
            markeredgecolor="black", markeredgewidth=0.8, zorder=zorder)
    if label:
        ax.text(x + label_offset[0], y + label_offset[1], label,
                ha="center", va="top", fontsize=8, fontweight="bold",
                color="#222", zorder=zorder + 1)
    if sub:
        ax.text(x + label_offset[0], y + label_offset[1] - 0.03, sub,
                ha="center", va="top", fontsize=6.5,
                style="italic", color="#666", zorder=zorder + 1)


def _edge(ax, p1, p2, color, lw, label=None, label_pos=0.5,
          style="-", rad=0.0, alpha=0.85, zorder=2):
    """Draw a typed edge with optional label at fractional position."""
    if rad != 0.0:
        arrow = FancyArrowPatch(p1, p2, arrowstyle="-",
                                connectionstyle=f"arc3,rad={rad}",
                                linewidth=lw, color=color, alpha=alpha,
                                zorder=zorder)
        ax.add_patch(arrow)
    else:
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], style, color=color,
                linewidth=lw, alpha=alpha, zorder=zorder)
    if label:
        mx = p1[0] + label_pos * (p2[0] - p1[0])
        my = p1[1] + label_pos * (p2[1] - p1[1])
        ax.text(mx, my + 0.012, label, ha="center", va="bottom",
                fontsize=6.5, style="italic", color=color,
                zorder=zorder + 2,
                bbox=dict(boxstyle="round,pad=0.1", facecolor="white",
                          edgecolor="none", alpha=0.85))


def draw_panel_b(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("b  Relational-prior fine-mapping — multi-omics knowledge graph",
                 loc="left", fontsize=11, fontweight="bold", pad=8)

    # Query variant — highlighted at centre-left
    qx, qy = 0.18, 0.55
    _node(ax, qx, qy, "o", VARIANT_COLOR_CAUSAL, size=14,
          label="query variant", sub="chr16:53,767,042:T/C",
          halo=True, zorder=5)

    # Gene node FTO (coding neighbour)
    fto_x, fto_y = 0.42, 0.72
    _node(ax, fto_x, fto_y, "s", GENE_COLOR, size=13,
          label="FTO", sub="gene")
    _edge(ax, (qx, qy), (fto_x, fto_y), "#888", 1.0,
          label="HAS_CONSEQUENCE", label_pos=0.55)

    # Gene node IRX3 (eQTL target in a specific tissue) — shows the
    # tissue as a sub-attribute of the eQTL edge
    irx3_x, irx3_y = 0.42, 0.38
    _node(ax, irx3_x, irx3_y, "s", GENE_COLOR, size=13,
          label="IRX3", sub="gene (trans-eQTL)")
    _edge(ax, (qx, qy), (irx3_x, irx3_y), EQTL_COLOR, 2.0,
          label="eQTL  [Adipose_Subcutaneous]", label_pos=0.52)

    # Regulatory element (cCRE from ENCODE)
    ccre_x, ccre_y = 0.22, 0.80
    _node(ax, ccre_x, ccre_y, "d", CCRE_COLOR, size=11,
          label="cCRE", sub="enhancer", label_offset=(0, 0.05))
    _edge(ax, (qx, qy), (ccre_x, ccre_y), "#888", 1.0,
          label="IN_REGULATORY", label_pos=0.5)

    # PPI partners of FTO
    ppi_partners = [
        (0.68, 0.82, "AKT1"),
        (0.72, 0.72, "ADIPOQ"),
        (0.66, 0.62, "RPGRIP1L"),
    ]
    for px, py, name in ppi_partners:
        _node(ax, px, py, "s", GENE_COLOR, size=10,
              label=name, sub=None)
        _edge(ax, (fto_x, fto_y), (px, py), PPI_COLOR, 1.2, alpha=0.75)
    ax.text(0.70, 0.92, "INTERACTS_WITH (STRING $\\geq 700$)",
            ha="center", va="center", fontsize=7,
            color=PPI_COLOR, style="italic")

    # Pathway nodes
    pw1_x, pw1_y = 0.60, 0.50
    _node(ax, pw1_x, pw1_y, "^", PATHWAY_COLOR, size=12,
          label="fatty acid\nmetabolism", sub="pathway",
          label_offset=(0, -0.04))
    _edge(ax, (fto_x, fto_y), (pw1_x, pw1_y), "#999", 1.0,
          label="IN_PATHWAY", label_pos=0.55)
    _edge(ax, (irx3_x, irx3_y), (pw1_x, pw1_y), "#999", 0.8, alpha=0.6)

    # Tissue expression annotation node (for eQTL context)
    tiss_x, tiss_y = 0.30, 0.22
    _node(ax, tiss_x, tiss_y, "h", TISSUE_COLOR, size=11,
          label="Adipose_Subcutaneous", sub="GTEx tissue",
          label_offset=(0, -0.045))
    _edge(ax, (irx3_x, irx3_y), (tiss_x, tiss_y), EQTL_COLOR, 1.0,
          alpha=0.55)

    # Output box: credible set IS a graph object
    out_box_x, out_box_y, out_w, out_h = 0.76, 0.28, 0.42, 0.22
    box = FancyBboxPatch(
        (out_box_x - out_w/2, out_box_y - out_h/2), out_w, out_h,
        boxstyle="round,pad=0.01",
        linewidth=1.3, edgecolor=OUTPUT_COLOR,
        facecolor=OUTPUT_COLOR, alpha=0.08, zorder=1,
    )
    ax.add_patch(box)
    ax.text(out_box_x, out_box_y + 0.07,
            "CredibleSet (graph object)",
            ha="center", va="center", fontsize=8.5, fontweight="bold",
            color=OUTPUT_COLOR)
    ax.text(out_box_x, out_box_y + 0.02,
            "chr16:53,767,042:T/C  PIP = 0.87",
            ha="center", va="center", fontsize=7.5,
            family="monospace", color="#333")
    ax.text(out_box_x, out_box_y - 0.02,
            "+ linked to: FTO, IRX3, cCRE,",
            ha="center", va="center", fontsize=7.2, color="#555")
    ax.text(out_box_x, out_box_y - 0.05,
            "  Adipose_Subcutaneous, fatty-acid metabolism",
            ha="center", va="center", fontsize=7.2, color="#555")
    ax.text(out_box_x, out_box_y - 0.09,
            "all co-queryable in one traversal",
            ha="center", va="center", fontsize=6.8, color="#888",
            style="italic")

    # Arrow from graph region to credible-set output box
    ax.annotate("", xy=(out_box_x - out_w/2 - 0.015, out_box_y),
                xytext=(fto_x + 0.05, 0.55),
                arrowprops=dict(arrowstyle="->", lw=1.4,
                                color=OUTPUT_COLOR, alpha=0.9,
                                connectionstyle="arc3,rad=-0.2"))

    # Legend / node-type key (bottom-left, compact)
    legend_y = 0.02
    legend_items = [
        ("o", VARIANT_COLOR_CAUSAL, "Variant"),
        ("s", GENE_COLOR, "Gene"),
        ("^", PATHWAY_COLOR, "Pathway"),
        ("d", CCRE_COLOR, "RegElement"),
        ("h", TISSUE_COLOR, "Tissue"),
    ]
    lx = 0.04
    for marker, color, label in legend_items:
        ax.plot(lx, legend_y, marker, color=color, markersize=8,
                markeredgecolor="black", markeredgewidth=0.5)
        ax.text(lx + 0.015, legend_y, label, va="center",
                fontsize=7, color="#333")
        lx += 0.125

    # Method-box-equivalent at top (minimal — full algorithm is Fig 2)
    ax.text(0.98, 0.98, "$\\rightarrow$ HBP / GAFM message passing "
                         "on this graph (Fig.\\,2)",
            ha="right", va="top", fontsize=7.5, style="italic",
            color=METHOD_BOX_COLOR,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                      edgecolor=METHOD_BOX_COLOR, linewidth=0.8,
                      alpha=0.9))

    # 1KG loaded counts — inline on right edge (kept from earlier)
    counts_text = ("1KG-loaded schema:\n"
                   "Variant: 70.7 M\n"
                   "Gene: 20,092 (GENCODE v47)\n"
                   "eQTL edges: 43.2 M (GTEx v8, 49 tissues)\n"
                   "PPI edges: 230,850 (STRING $\\geq$ 700)\n"
                   "RegElement: 370,000 (ENCODE cCRE v4)")
    ax.text(0.99, 0.11, counts_text,
            ha="right", va="bottom", fontsize=6.5, color="#444",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                      edgecolor="#aaa", linewidth=0.5), zorder=5)


# ===================================================================
# Compose figure
# ===================================================================

def main() -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 6.4))
    draw_panel_a(ax1)
    draw_panel_b(ax2)

    fig.suptitle(
        "Figure 1 — Flat-prior vs relational-prior fine-mapping",
        fontsize=12, fontweight="bold", y=0.995,
    )

    plt.tight_layout(rect=[0, 0, 1, 0.97])

    out_dir = Path("/mnt/data/GraphGWAS/paper/manuscript_v1/figures")
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out_dir / f"fig1_architecture.{ext}", dpi=200,
                    bbox_inches="tight")

    alt_dir = Path("/mnt/data/GraphGWAS/results/benchmark_v2/paper_figures")
    for ext in ("png", "pdf"):
        fig.savefig(alt_dir / f"fig1_architecture.{ext}", dpi=200,
                    bbox_inches="tight")

    print(f"Wrote {out_dir}/fig1_architecture.{{png,pdf}}")
    print(f"Wrote {alt_dir}/fig1_architecture.{{png,pdf}}")


if __name__ == "__main__":
    main()
