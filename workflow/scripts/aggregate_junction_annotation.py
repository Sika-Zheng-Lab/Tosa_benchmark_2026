"""Aggregate per-sample junction annotation classification results.

Generates:
- annotation_summary.tsv          — per-tool mean fractions across samples
- star_missed_summary.tsv         — category breakdown of junctions missed by STAR
- annotation_proportion.png       — stacked bar: annotated/genic_novel/intergenic
                                    proportion per tool (mean across samples)
- star_missed_annotation.png      — stacked bar: category breakdown of junctions
                                    detected by each other tool but missed by STAR
"""

import os
import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Snakemake I/O
# ---------------------------------------------------------------------------
in_summaries = snakemake.input.summaries
in_classifications = snakemake.input.classifications

out_annotation_summary = snakemake.output.annotation_summary
out_star_missed_summary = snakemake.output.star_missed_summary
out_annotation_plot = snakemake.output.annotation_plot
out_star_missed_plot = snakemake.output.star_missed_plot

samples = snakemake.params.samples

for out_path in [out_annotation_summary, out_star_missed_summary,
                 out_annotation_plot, out_star_missed_plot]:
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

TOOLS = ["STAR", "regtools", "featureCounts", "Tosa"]
CATEGORIES = ["annotated", "genic_novel", "intergenic"]
CAT_LABELS = {"annotated": "Annotated", "genic_novel": "Genic novel",
              "intergenic": "Intergenic"}
CAT_COLORS = {"annotated": "#4C72B0", "genic_novel": "#DD8452",
              "intergenic": "#55A868"}

# ---------------------------------------------------------------------------
# 1. Aggregate per-tool annotation proportions across samples
# ---------------------------------------------------------------------------
frames = []
for path in in_summaries:
    frames.append(pd.read_csv(path, sep="\t"))
df_all = pd.concat(frames, ignore_index=True)

# Mean fraction per tool × category across samples
df_mean = (
    df_all.groupby(["tool", "category"])["fraction"]
    .mean()
    .reset_index()
    .rename(columns={"fraction": "mean_fraction"})
)
# Also compute SD
df_sd = (
    df_all.groupby(["tool", "category"])["fraction"]
    .std()
    .reset_index()
    .rename(columns={"fraction": "sd_fraction"})
)
df_summary = df_mean.merge(df_sd, on=["tool", "category"])

# Mean count
df_mean_count = (
    df_all.groupby(["tool", "category"])["count"]
    .mean()
    .reset_index()
    .rename(columns={"count": "mean_count"})
)
df_summary = df_summary.merge(df_mean_count, on=["tool", "category"])
df_summary.to_csv(out_annotation_summary, sep="\t", index=False)

# ---------------------------------------------------------------------------
# 2. Analyse junctions missed by STAR
# ---------------------------------------------------------------------------
# For each sample, load classification.tsv and identify junctions where
# STAR == 0 but each other tool == 1, then tally categories.
star_missed_rows = []

for path in in_classifications:
    df_cls = pd.read_csv(path, sep="\t")
    # Extract sample name from path
    sample = os.path.basename(os.path.dirname(path))
    for tool in ["Tosa", "regtools", "featureCounts"]:
        missed = df_cls[(df_cls["STAR"] == 0) & (df_cls[tool] == 1)]
        total_missed = len(missed)
        for cat in CATEGORIES:
            n = (missed["category"] == cat).sum()
            frac = n / total_missed if total_missed > 0 else 0.0
            star_missed_rows.append({
                "sample": sample,
                "tool": tool,
                "category": cat,
                "count": n,
                "fraction": frac,
                "total_missed": total_missed,
            })

df_star_missed = pd.DataFrame(star_missed_rows)

# Mean across samples
df_sm_mean = (
    df_star_missed.groupby(["tool", "category"])
    .agg(mean_fraction=("fraction", "mean"),
         sd_fraction=("fraction", "std"),
         mean_count=("count", "mean"))
    .reset_index()
)
df_sm_mean.to_csv(out_star_missed_summary, sep="\t", index=False)

# ---------------------------------------------------------------------------
# 3. Plot: Annotation proportion per tool (stacked bar)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(5, 5))

# Pivot to tools × categories
pivot = df_summary.pivot(index="tool", columns="category",
                         values="mean_fraction").reindex(TOOLS)[CATEGORIES]

bottom = np.zeros(len(TOOLS))
for cat in CATEGORIES:
    vals = pivot[cat].values
    ax.bar(TOOLS, vals, bottom=bottom, label=CAT_LABELS[cat],
           color=CAT_COLORS[cat], edgecolor="white", linewidth=0.5)
    # Add percentage labels
    for i, (v, b) in enumerate(zip(vals, bottom)):
        if v > 0.03:  # only label if visible
            ax.text(i, b + v / 2, f"{v:.1%}", ha="center", va="center",
                    fontsize=9, color="white", fontweight="bold")
    bottom += vals

ax.set_ylabel("Fraction of detected junctions", fontsize=11)
ax.set_title("Junction annotation status per tool\n(mean across samples)",
             fontsize=12)
ax.set_ylim(0, 1.05)
ax.legend(loc="upper right", fontsize=10)
ax.tick_params(labelsize=11)

plt.tight_layout()
plt.savefig(out_annotation_plot, dpi=300)
plt.close()

# ---------------------------------------------------------------------------
# 4. Plot: STAR-missed junctions — category breakdown (stacked bar)
# ---------------------------------------------------------------------------
OTHER_TOOLS = ["regtools", "featureCounts", "Tosa"]

fig, ax = plt.subplots(figsize=(5, 5))

pivot_sm = df_sm_mean.pivot(index="tool", columns="category",
                            values="mean_fraction").reindex(OTHER_TOOLS)[CATEGORIES]

bottom = np.zeros(len(OTHER_TOOLS))
for cat in CATEGORIES:
    vals = pivot_sm[cat].values
    ax.bar(OTHER_TOOLS, vals, bottom=bottom, label=CAT_LABELS[cat],
           color=CAT_COLORS[cat], edgecolor="white", linewidth=0.5)
    for i, (v, b) in enumerate(zip(vals, bottom)):
        if v > 0.03:
            ax.text(i, b + v / 2, f"{v:.1%}", ha="center", va="center",
                    fontsize=9, color="white", fontweight="bold")
    bottom += vals

ax.set_ylabel("Fraction of STAR-missed junctions", fontsize=11)
ax.set_title("Annotation status of junctions\ndetected by tool but missed by STAR\n"
             "(mean across samples)", fontsize=12)
ax.set_ylim(0, 1.05)
ax.legend(loc="upper right", fontsize=10)
ax.tick_params(labelsize=11)

plt.tight_layout()
plt.savefig(out_star_missed_plot, dpi=300)
plt.close()
