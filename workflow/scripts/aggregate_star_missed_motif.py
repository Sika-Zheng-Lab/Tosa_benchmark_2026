"""Aggregate STAR-missed motif analysis across all samples.

Reads per-sample star_missed_motif.tsv files produced by
analyze_star_missed_motif.py and generates:

  1. star_missed_motif_summary.tsv  — mean junction counts and fractions
     per motif category × group across all samples
  2. star_missed_motif_barplot.png  — stacked bar chart comparing motif
     category proportions between tosa_only and both groups
  3. star_missed_count_distribution.png — box plot of log10(tosa_count) for
     tosa_only vs both groups, stratified by motif category
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Snakemake I/O
# ---------------------------------------------------------------------------
in_files = snakemake.input.motif_tsvs
samples = snakemake.params.samples

out_summary = snakemake.output.motif_summary
out_barplot = snakemake.output.motif_barplot
out_countplot = snakemake.output.count_distribution

for p in [out_summary, out_barplot, out_countplot]:
    os.makedirs(os.path.dirname(p), exist_ok=True)

# ---------------------------------------------------------------------------
# Load and concatenate all per-sample data
# ---------------------------------------------------------------------------
frames = []
for f, sample in zip(in_files, samples):
    df = pd.read_csv(f, sep="\t")
    df["sample"] = sample
    frames.append(df)

all_data = pd.concat(frames, ignore_index=True)

# ---------------------------------------------------------------------------
# Compute per-sample summary, then average across samples
# ---------------------------------------------------------------------------
MOTIF_ORDER = ["canonical", "semi_canonical", "non_canonical"]
GROUP_ORDER = ["both", "tosa_only"]

per_sample = (
    all_data.groupby(["sample", "group", "motif_category"])
    .size()
    .reset_index(name="count")
)

# Compute fraction within each sample × group
totals = per_sample.groupby(["sample", "group"])["count"].transform("sum")
per_sample["fraction"] = per_sample["count"] / totals

# Mean across samples
mean_summary = (
    per_sample.groupby(["group", "motif_category"])[["count", "fraction"]]
    .mean()
    .reset_index()
)
mean_summary.to_csv(out_summary, sep="\t", index=False)

# ---------------------------------------------------------------------------
# Plot 1: Stacked bar chart — motif category proportions per group
# ---------------------------------------------------------------------------
plt.rcParams["font.family"] = "Arial"

fig, ax = plt.subplots(figsize=(4, 5))

COLORS = {
    "canonical": "#4C72B0",
    "semi_canonical": "#DD8452",
    "non_canonical": "#55A868",
}

pivot = mean_summary.pivot(index="group", columns="motif_category", values="fraction")
pivot = pivot.reindex(index=GROUP_ORDER, columns=MOTIF_ORDER).fillna(0)

bottom = np.zeros(len(GROUP_ORDER))
for motif in MOTIF_ORDER:
    vals = pivot[motif].values
    ax.bar(
        GROUP_ORDER,
        vals,
        bottom=bottom,
        color=COLORS[motif],
        label=motif.replace("_", " "),
        edgecolor="white",
        linewidth=0.5,
    )
    bottom += vals

ax.set_ylabel("Fraction of junctions")
ax.set_ylim(0, 1.05)
ax.set_title("Splice-site motif composition")
ax.legend(title="Motif", bbox_to_anchor=(1.02, 1), loc="upper left")
ax.set_xticklabels(["Detected by\nboth", "Tosa only\n(STAR missed)"])

fig.tight_layout()
fig.savefig(out_barplot, dpi=300, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# Plot 2: Count distribution — tosa_only vs both, by motif category
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6, 5))

# Combine group and motif for x-axis categories
plot_data = all_data[["tosa_count", "group", "motif_category"]].copy()
plot_data["log10_count"] = np.log10(plot_data["tosa_count"].clip(lower=1))

positions = []
labels = []
box_data = []
colors_list = []

pos = 0
for group in GROUP_ORDER:
    for motif in MOTIF_ORDER:
        subset = plot_data[
            (plot_data["group"] == group) & (plot_data["motif_category"] == motif)
        ]["log10_count"]
        if len(subset) > 0:
            box_data.append(subset.values)
        else:
            box_data.append([0])
        positions.append(pos)
        labels.append(f"{motif.replace('_', ' ')}")
        colors_list.append(COLORS[motif])
        pos += 1
    pos += 0.5  # gap between groups

bp = ax.boxplot(
    box_data,
    positions=positions,
    widths=0.7,
    patch_artist=True,
    showfliers=False,
    medianprops=dict(color="black", linewidth=1.5),
)
for patch, color in zip(bp["boxes"], colors_list):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

ax.set_xticks(positions)
ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
ax.set_ylabel("log10(Tosa read count)")
ax.set_title("Read count distribution by motif and detection group")

# Add group labels
group_centers = [
    np.mean(positions[:len(MOTIF_ORDER)]),
    np.mean(positions[len(MOTIF_ORDER):]),
]
y_min = ax.get_ylim()[0]
for center, group_label in zip(group_centers, ["Detected by both", "Tosa only\n(STAR missed)"]):
    ax.text(center, y_min - 0.3, group_label, ha="center", va="top", fontsize=9, fontweight="bold")

ax.set_xlim(positions[0] - 0.5, positions[-1] + 0.5)
fig.tight_layout()
fig.subplots_adjust(bottom=0.25)
fig.savefig(out_countplot, dpi=300, bbox_inches="tight")
plt.close(fig)
