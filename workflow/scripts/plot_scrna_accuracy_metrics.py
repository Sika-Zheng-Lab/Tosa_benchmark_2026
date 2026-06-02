"""
Generate grouped bar chart of per-tool scRNA-seq SJ accuracy metrics (Figure 4A).

Reads the simulation accuracy TSV produced by compare_scrna_junctions.py
and plots sensitivity, precision, F1 score, and exact match rate for each tool.

Input:
  - simulation_accuracy.tsv (columns: tool, sensitivity, precision,
    f1_score, exact_match_rate, ...)

Output:
  - Grouped bar chart PNG
"""

import os
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Snakemake I/O
# ---------------------------------------------------------------------------
in_accuracy = snakemake.input.accuracy

out_plot = snakemake.output.plot

os.makedirs(os.path.dirname(out_plot), exist_ok=True)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
df = pd.read_csv(in_accuracy, sep="\t")

METRICS = ["sensitivity", "precision", "f1_score", "exact_match_rate"]
METRIC_LABELS = ["Sensitivity", "Precision", "F1 Score", "Exact Match Rate"]
TOOL_ORDER = ["STARsolo", "regtools", "Tosa"]
TOOLS = [t for t in TOOL_ORDER if t in df["tool"].values]

# Tool colour palette (consistent with Figure 2D-E)
TOOL_COLORS = {
    "STARsolo": "#9467bd",
    "regtools": "#ff7f0e",
    "Tosa": "#1f77b4",
}

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
n_metrics = len(METRICS)
n_tools = len(TOOLS)
x = np.arange(n_metrics)
bar_width = 0.8 / n_tools

fig, ax = plt.subplots(figsize=(6, 4))

for i, tool in enumerate(TOOLS):
    row = df[df["tool"] == tool].iloc[0]
    values = [row[m] for m in METRICS]
    offset = (i - n_tools / 2 + 0.5) * bar_width
    bars = ax.bar(
        x + offset,
        values,
        bar_width,
        label=tool,
        color=TOOL_COLORS.get(tool, "#999999"),
        edgecolor="black",
        linewidth=0.5,
    )
    # Add value labels on top of bars
    for bar, val in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            bar.get_height() + 0.01,
            f"{val:.2f}",
            ha="center",
            va="bottom",
            fontsize=10,
        )

ax.set_xticks(x)
ax.set_xticklabels(METRIC_LABELS, fontsize=12)
ax.set_ylabel("Score", fontsize=14)
ax.set_title("scRNA-seq SJ Counting Accuracy (Simulation)", fontsize=15)
ax.set_ylim(0, 1.15)
ax.legend(fontsize=11)
ax.tick_params(labelsize=12)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
fig.savefig(out_plot, dpi=300, bbox_inches="tight")
plt.close(fig)
