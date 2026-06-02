#!/usr/bin/env python3
"""
Aggregate large-scale benchmark results and generate comparison plots.

Reads Snakemake benchmark TSV files and per-sample read counts, then produces:
  1. Per-sample benchmark table (time, memory, read count for each tool)
  2. Summary table with total wall-clock time per tool
  3. Bar plot: total execution time comparison
  4. Box plot: per-sample execution time distribution
  5. Scatter plot: read count vs execution time
  6. Scatter plot: read count vs memory usage
  7. Box plot: memory usage distribution
  8. Speed-up comparison (1-thread vs 8-thread)
"""

import os
import sys
import logging

import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from scipy import stats

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Parameters from Snakemake
# ---------------------------------------------------------------------------

SAMPLES = snakemake.params.samples
TOOLS = snakemake.params.tools
BENCH_DIR = snakemake.params.benchmark_dir
RC_DIR = snakemake.params.read_count_dir

# Tool display names
TOOL_LABELS = {
    "tosa_1t": "Tosa (1 thread)",
    "tosa_8t": "Tosa (8 threads)",
    "regtools_1t": "regtools (1 thread)",
    "featurecounts_1t": "featureCounts (1 thread)",
    "featurecounts_8t": "featureCounts (8 threads)",
}

# Color palette
TOOL_COLORS = {
    "Tosa (1 thread)": "#1f77b4",
    "Tosa (8 threads)": "#aec7e8",
    "regtools (1 thread)": "#ff7f0e",
    "featureCounts (1 thread)": "#2ca02c",
    "featureCounts (8 threads)": "#98df8a",
}

# ---------------------------------------------------------------------------
# 1. Read benchmark data
# ---------------------------------------------------------------------------

logger.info("Reading benchmark files ...")

records = []
for sample in SAMPLES:
    for tool in TOOLS:
        bench_file = os.path.join(BENCH_DIR, sample, f"{tool}.tsv")
        try:
            df = pd.read_csv(bench_file, sep="\t")
            # Snakemake benchmark columns: s, h:m:s, max_rss, max_vms, ...
            row = {
                "sample": sample,
                "tool": tool,
                "tool_label": TOOL_LABELS[tool],
                "wall_time_s": df["s"].mean(),
                "max_rss_mb": df["max_rss"].mean(),
                "max_vms_mb": df["max_vms"].mean(),
            }
            records.append(row)
        except Exception as e:
            logger.warning(f"Failed to read {bench_file}: {e}")

bench_df = pd.DataFrame(records)
logger.info(f"Read {len(bench_df)} benchmark records for {len(SAMPLES)} samples x {len(TOOLS)} tools.")

# ---------------------------------------------------------------------------
# 2. Read per-sample read counts
# ---------------------------------------------------------------------------

logger.info("Reading read count files ...")

read_counts = {}
for sample in SAMPLES:
    rc_file = os.path.join(RC_DIR, f"{sample}.txt")
    try:
        with open(rc_file) as f:
            read_counts[sample] = int(f.read().strip())
    except Exception as e:
        logger.warning(f"Failed to read {rc_file}: {e}")
        read_counts[sample] = np.nan

bench_df["read_count"] = bench_df["sample"].map(read_counts)

# ---------------------------------------------------------------------------
# 3. Save per-sample table
# ---------------------------------------------------------------------------

out_per_sample = snakemake.output.per_sample
bench_df.to_csv(out_per_sample, sep="\t", index=False)
logger.info(f"Saved per-sample table: {out_per_sample}")

# ---------------------------------------------------------------------------
# 4. Compute summary (total wall-clock time per tool)
# ---------------------------------------------------------------------------

summary = (
    bench_df.groupby(["tool", "tool_label"])
    .agg(
        total_wall_time_s=("wall_time_s", "sum"),
        mean_wall_time_s=("wall_time_s", "mean"),
        median_wall_time_s=("wall_time_s", "median"),
        mean_max_rss_mb=("max_rss_mb", "mean"),
        median_max_rss_mb=("max_rss_mb", "median"),
        max_max_rss_mb=("max_rss_mb", "max"),
        n_samples=("sample", "count"),
    )
    .reset_index()
    .sort_values("total_wall_time_s")
)
summary["total_wall_time_min"] = summary["total_wall_time_s"] / 60.0
summary["total_wall_time_h"] = summary["total_wall_time_s"] / 3600.0

out_summary = snakemake.output.summary
summary.to_csv(out_summary, sep="\t", index=False)
logger.info(f"Saved summary table: {out_summary}")

# ---------------------------------------------------------------------------
# Plot settings
# ---------------------------------------------------------------------------

sns.set_theme(style="whitegrid", font_scale=1.4)
plt.rcParams["figure.dpi"] = 300

# Order for consistent plotting
TOOL_ORDER = [
    "regtools (1 thread)",
    "featureCounts (1 thread)",
    "featureCounts (8 threads)",
    "Tosa (1 thread)",
    "Tosa (8 threads)",
]

# ---------------------------------------------------------------------------
# Plot 1: Bar plot — total execution time
# ---------------------------------------------------------------------------

logger.info("Generating bar plot: total execution time ...")

fig, ax = plt.subplots(figsize=(10, 6))
plot_data = summary.set_index("tool_label").reindex(TOOL_ORDER)
bars = ax.bar(
    plot_data.index,
    plot_data["total_wall_time_min"],
    color=[TOOL_COLORS[t] for t in plot_data.index],
    edgecolor="black",
    linewidth=0.5,
)
# Add value labels on bars
for bar in bars:
    height = bar.get_height()
    ax.text(
        bar.get_x() + bar.get_width() / 2.0,
        height,
        f"{height:.1f}",
        ha="center",
        va="bottom",
        fontsize=10,
    )
ax.set_ylabel("Total wall-clock time (minutes)")
ax.set_title(f"Total execution time across {len(SAMPLES)} samples")
ax.set_xticklabels(plot_data.index, rotation=30, ha="right")
plt.tight_layout()
fig.savefig(snakemake.output.barplot, bbox_inches="tight")
plt.close(fig)
logger.info(f"Saved: {snakemake.output.barplot}")

# ---------------------------------------------------------------------------
# Plot 2: Box plot — per-sample execution time
# ---------------------------------------------------------------------------

logger.info("Generating box plot: per-sample execution time ...")

fig, ax = plt.subplots(figsize=(10, 6))
sns.boxplot(
    data=bench_df,
    x="tool_label",
    y="wall_time_s",
    order=TOOL_ORDER,
    palette=TOOL_COLORS,
    ax=ax,
    showfliers=True,
    flierprops=dict(marker="o", markersize=3, alpha=0.5),
)
ax.set_ylabel("Wall-clock time per sample (seconds)")
ax.set_xlabel("")
ax.set_title("Per-sample execution time distribution")
ax.set_xticklabels(TOOL_ORDER, rotation=30, ha="right")
plt.tight_layout()
fig.savefig(snakemake.output.boxplot_time, bbox_inches="tight")
plt.close(fig)
logger.info(f"Saved: {snakemake.output.boxplot_time}")

# ---------------------------------------------------------------------------
# Plot 3: Scatter — read count vs execution time
# ---------------------------------------------------------------------------

logger.info("Generating scatter plot: read count vs execution time ...")

fig, ax = plt.subplots(figsize=(5, 4.5))
for tool_label in TOOL_ORDER:
    subset = bench_df[bench_df["tool_label"] == tool_label].dropna(subset=["read_count"])
    ax.scatter(
        subset["read_count"] / 1e6,
        subset["wall_time_s"],
        label=tool_label,
        color=TOOL_COLORS[tool_label],
        alpha=0.6,
        s=20,
        edgecolors="none",
    )
    # Add regression line
    if len(subset) > 2:
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            subset["read_count"] / 1e6, subset["wall_time_s"]
        )
        x_range = np.linspace(subset["read_count"].min() / 1e6, subset["read_count"].max() / 1e6, 100)
        ax.plot(x_range, slope * x_range + intercept, color=TOOL_COLORS[tool_label], linewidth=1.5, alpha=0.8)

ax.set_xlabel("Read count (millions)")
ax.set_ylabel("Wall-clock time (seconds)")
ax.set_title("Read count vs execution time")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_box_aspect(1)
ax.legend(fontsize=9, bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
fig.savefig(snakemake.output.scatter_time, dpi=300, bbox_inches="tight")
plt.close(fig)
logger.info(f"Saved: {snakemake.output.scatter_time}")

# ---------------------------------------------------------------------------
# Plot 4: Scatter — read count vs memory usage (max_rss)
# ---------------------------------------------------------------------------

logger.info("Generating scatter plot: read count vs memory usage ...")

fig, ax = plt.subplots(figsize=(5, 4.5))
for tool_label in TOOL_ORDER:
    subset = bench_df[bench_df["tool_label"] == tool_label].dropna(subset=["read_count"])
    ax.scatter(
        subset["read_count"] / 1e6,
        subset["max_rss_mb"],
        label=tool_label,
        color=TOOL_COLORS[tool_label],
        alpha=0.6,
        s=20,
        edgecolors="none",
    )
    # Add regression line
    if len(subset) > 2:
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            subset["read_count"] / 1e6, subset["max_rss_mb"]
        )
        x_range = np.linspace(subset["read_count"].min() / 1e6, subset["read_count"].max() / 1e6, 100)
        ax.plot(x_range, slope * x_range + intercept, color=TOOL_COLORS[tool_label], linewidth=1.5, alpha=0.8)

ax.set_xlabel("Read count (millions)")
ax.set_ylabel("Max RSS memory (MB)")
ax.set_title("Read count vs memory usage")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_box_aspect(1)
ax.legend(fontsize=9, bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
fig.savefig(snakemake.output.scatter_memory, dpi=300, bbox_inches="tight")
plt.close(fig)
logger.info(f"Saved: {snakemake.output.scatter_memory}")

# ---------------------------------------------------------------------------
# Plot 5: Box plot — memory usage distribution
# ---------------------------------------------------------------------------

logger.info("Generating box plot: memory usage ...")

fig, ax = plt.subplots(figsize=(10, 6))
sns.boxplot(
    data=bench_df,
    x="tool_label",
    y="max_rss_mb",
    order=TOOL_ORDER,
    palette=TOOL_COLORS,
    ax=ax,
    showfliers=True,
    flierprops=dict(marker="o", markersize=3, alpha=0.5),
)
ax.set_ylabel("Max RSS memory (MB)")
ax.set_xlabel("")
ax.set_title("Per-sample memory usage distribution")
ax.set_xticklabels(TOOL_ORDER, rotation=30, ha="right")
plt.tight_layout()
fig.savefig(snakemake.output.boxplot_memory, bbox_inches="tight")
plt.close(fig)
logger.info(f"Saved: {snakemake.output.boxplot_memory}")

# ---------------------------------------------------------------------------
# Plot 6: Speed-up comparison (1t → 8t)
# ---------------------------------------------------------------------------

logger.info("Generating speed-up comparison ...")

# Pivot: sample × tool → wall_time_s
pivot = bench_df.pivot_table(index="sample", columns="tool", values="wall_time_s")

speedup_records = []
for sample in SAMPLES:
    if sample not in pivot.index:
        continue
    row = pivot.loc[sample]
    # Tosa speed-up
    if "tosa_1t" in row and "tosa_8t" in row and row["tosa_8t"] > 0:
        speedup_records.append({
            "sample": sample,
            "tool": "Tosa",
            "speedup": row["tosa_1t"] / row["tosa_8t"],
            "read_count": read_counts.get(sample, np.nan),
        })
    # featureCounts speed-up
    if "featurecounts_1t" in row and "featurecounts_8t" in row and row["featurecounts_8t"] > 0:
        speedup_records.append({
            "sample": sample,
            "tool": "featureCounts",
            "speedup": row["featurecounts_1t"] / row["featurecounts_8t"],
            "read_count": read_counts.get(sample, np.nan),
        })

speedup_df = pd.DataFrame(speedup_records)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Left: box plot of speed-up
if not speedup_df.empty:
    sns.boxplot(
        data=speedup_df,
        x="tool",
        y="speedup",
        palette={"Tosa": "#1f77b4", "featureCounts": "#2ca02c"},
        ax=axes[0],
    )
    axes[0].axhline(y=1, color="red", linestyle="--", alpha=0.5, label="No speed-up")
    axes[0].axhline(y=8, color="gray", linestyle=":", alpha=0.5, label="Linear (8x)")
    axes[0].set_ylabel("Speed-up (1t / 8t)")
    axes[0].set_title("Multi-thread speed-up (1→8 threads)")
    axes[0].legend(fontsize=9)

    # Right: scatter — read count vs speed-up
    for tool_name, color in [("Tosa", "#1f77b4"), ("featureCounts", "#2ca02c")]:
        subset = speedup_df[speedup_df["tool"] == tool_name].dropna(subset=["read_count"])
        axes[1].scatter(
            subset["read_count"] / 1e6,
            subset["speedup"],
            label=tool_name,
            color=color,
            alpha=0.6,
            s=20,
            edgecolors="none",
        )
    axes[1].axhline(y=1, color="red", linestyle="--", alpha=0.5)
    axes[1].axhline(y=8, color="gray", linestyle=":", alpha=0.5)
    axes[1].set_xlabel("Read count (millions)")
    axes[1].set_ylabel("Speed-up (1t / 8t)")
    axes[1].set_title("Read count vs speed-up")
    axes[1].legend(fontsize=9)

plt.tight_layout()
fig.savefig(snakemake.output.speedup, bbox_inches="tight")
plt.close(fig)
logger.info(f"Saved: {snakemake.output.speedup}")

logger.info("All done.")
