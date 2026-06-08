#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from pathlib import Path

outdir = Path("figures/conceptual_model")
outdir.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(13, 8.5))
ax.set_xlim(0, 13)
ax.set_ylim(0, 8.5)
ax.axis("off")

def box(x, y, w, h, text, fc="#F7F7F7", ec="#333333",
        fontsize=10, weight="normal", lw=1.2):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.06,rounding_size=0.08",
        linewidth=lw,
        edgecolor=ec,
        facecolor=fc
    )
    ax.add_patch(patch)
    ax.text(
        x + w / 2,
        y + h / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        fontweight=weight,
        wrap=True
    )

def arrow(x1, y1, x2, y2, color="#333333", lw=1.4, rad=0.0):
    arr = FancyArrowPatch(
        (x1, y1), (x2, y2),
        arrowstyle="-|>",
        mutation_scale=14,
        linewidth=lw,
        color=color,
        connectionstyle=f"arc3,rad={rad}"
    )
    ax.add_patch(arr)

# Title
ax.text(
    6.5, 8.05,
    "Shared osmoadaptive transcriptional response under salt stress",
    ha="center",
    va="center",
    fontsize=17,
    fontweight="bold"
)

# Top salt stress box
box(
    5.05, 7.05, 2.9, 0.6,
    "Salt stress",
    fc="#D9EAF7",
    fontsize=13,
    weight="bold"
)

# Species boxes
box(
    0.65, 5.95, 3.45, 0.75,
    "Bacillus velezensis\nControl vs salt",
    fc="#FFF2CC",
    fontsize=11.5,
    weight="bold"
)

box(
    8.9, 5.95, 3.45, 0.75,
    "Bacillus paralicheniformis\nControl vs salt",
    fc="#FFF2CC",
    fontsize=11.5,
    weight="bold"
)

# RBH integration
box(
    4.55, 5.75, 3.9, 0.85,
    "RBH orthology integration\n2,963 ortholog pairs",
    fc="#E2F0D9",
    fontsize=11.5,
    weight="bold"
)

# Arrows from salt stress
arrow(6.1, 7.05, 2.4, 6.7, rad=0.05)
arrow(6.9, 7.05, 10.6, 6.7, rad=-0.05)
arrow(6.5, 7.05, 6.5, 6.6)

# Arrows from species to RBH
arrow(4.10, 6.32, 4.55, 6.18)
arrow(8.90, 6.32, 8.45, 6.18)

# Response classes
box(
    0.45, 4.55, 3.55, 0.9,
    "Conserved induction\n249 orthologs",
    fc="#F4CCCC",
    fontsize=11.5,
    weight="bold"
)

box(
    4.75, 4.55, 3.55, 0.9,
    "Conserved repression\n256 orthologs",
    fc="#CFE2F3",
    fontsize=11.5,
    weight="bold"
)

box(
    9.05, 4.55, 3.55, 0.9,
    "Opposite responses\n275 orthologs",
    fc="#D9EAD3",
    fontsize=11.5,
    weight="bold"
)

# Arrows from RBH to response classes
arrow(5.35, 5.75, 2.25, 5.45, rad=0.02)
arrow(6.5, 5.75, 6.5, 5.45)
arrow(7.65, 5.75, 10.85, 5.45, rad=-0.02)

# Lower modules conserved induction
box(
    0.25, 3.15, 3.95, 0.9,
    "Compatible solutes\nproB, proC, betB, gbsB",
    fc="#FCE5CD",
    fontsize=10.5
)

box(
    0.25, 2.0, 3.95, 0.9,
    "Transport and ion homeostasis\nopuAA, opuAB, nhaC",
    fc="#FCE5CD",
    fontsize=10.5
)

# Lower modules conserved repression
box(
    4.55, 3.15, 3.95, 0.9,
    "Biofilm and extracellular matrix\ntasA, sipW, epsB",
    fc="#D9D2E9",
    fontsize=10.5
)

box(
    4.55, 2.0, 3.95, 0.9,
    "Cell-envelope and extracellular\nfunctions repressed",
    fc="#D9D2E9",
    fontsize=10.5
)

# Lower modules opposite
box(
    8.85, 3.15, 3.95, 0.9,
    "Regulatory divergence\nmainly accessory/unclassified genes",
    fc="#EADCF8",
    fontsize=10.5
)

box(
    8.85, 2.0, 3.95, 0.9,
    "Sodium homeostasis depleted\n0 opposite orthologs; FDR = 0.003",
    fc="#EADCF8",
    fontsize=10.5
)

# Vertical arrows within each column
arrow(2.22, 4.55, 2.22, 4.05)
arrow(2.22, 3.15, 2.22, 2.90)

arrow(6.52, 4.55, 6.52, 4.05)
arrow(6.52, 3.15, 6.52, 2.90)

arrow(10.82, 4.55, 10.82, 4.05)
arrow(10.82, 3.15, 10.82, 2.90)

# Final outcome box
box(
    3.25, 0.75, 6.5, 0.85,
    "Shared osmoadaptive core\nwith partial regulatory divergence",
    fc="#D9EAF7",
    fontsize=12,
    weight="bold"
)

# Clean arrows to final box
arrow(2.22, 2.0, 4.05, 1.60, rad=0.0)
arrow(6.52, 2.0, 6.52, 1.60)
arrow(10.82, 2.0, 8.95, 1.60, rad=0.0)

# Note outside the final box
ax.text(
    6.5, 0.22,
    "Model based on orthology-guided differential expression, functional annotation, and opposite-response depletion analysis",
    ha="center",
    va="center",
    fontsize=8.5
)

plt.tight_layout()

png = outdir / "Figure8_Conceptual_Model_Osmoadaptation.png"
pdf = outdir / "Figure8_Conceptual_Model_Osmoadaptation.pdf"

fig.savefig(png, dpi=600, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")

print(f"Saved: {png}")
print(f"Saved: {pdf}")
