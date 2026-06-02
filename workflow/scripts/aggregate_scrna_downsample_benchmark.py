"""
Aggregate scRNA-seq downsample benchmark results across read-depth fractions
(Figure 4A-B).

For each fraction, loads Tosa (1 thread) / Tosa (multi-threaded) / regtools
benchmark TSVs and plots wall-clock time and peak RSS memory vs read count as
scatter+line plots with mean ± SD error bars.

Color palette matches the bulk benchmark (aggregate_benchmark.py).
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
fractions = snakemake.params.fractions       # e.g. [1.0, 0.75, 0.5, 0.25]
frac_pcts = snakemake.params.frac_pcts       # e.g. ["100", "075", "050", "025"]
out_plot = snakemake.output.plot
out_tsv = snakemake.output.summary_tsv

os.makedirs(os.path.dirname(out_plot), exist_ok=True)

tosa_threads = snakemake.config.get("scrna_benchmark", {}).get("tosa_threads", 8)

# Tool order & colours — consistent with bulk aggregate_benchmark.py
TOOL_ORDER = ["regtools (1 thread)", "Tosa (1 thread)", f"Tosa ({tosa_threads} threads)"]
TOOL_COLORS = {
    "Tosa (1 thread)": "#1f77b4",
    f"Tosa ({tosa_threads} threads)": "#aec7e8",
    "regtools (1 thread)": "#ff7f0e",
}

# ---------------------------------------------------------------------------
# Load benchmark data + read counts
# ---------------------------------------------------------------------------
rows = []
for frac, pct in zip(fractions, frac_pcts):
    # Read count for this fraction
    rc_path = getattr(snakemake.input, f"read_count_{pct}")
    with open(rc_path) as fh:
        read_count = int(fh.read().strip())

    for tool_label, key_prefix in [
        ("Tosa (1 thread)", "bench_tosa_1t_"),
        (f"Tosa ({tosa_threads} threads)", "bench_tosa_mt_"),
        ("regtools (1 thread)", "bench_regtools_"),
    ]:
        path = getattr(snakemake.input, f"{key_prefix}{pct}")
        df_b = pd.read_csv(path, sep="\t")
        for _, row in df_b.iterrows():
            rows.append({
                "read_count": read_count,
                "read_count_M": read_count / 1e6,
                "tool": tool_label,
                "wall_time_s": row["s"],
                "max_rss_mb": row["max_rss"],
            })

df = pd.DataFrame(rows)

# Sort by read count for consistent x-axis
read_counts_sorted = sorted(df["read_count_M"].unique())

# ---------------------------------------------------------------------------
# Plot: 2-panel scatter+line (wall time + peak memory vs read count)
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

# --- Panel A: Wall-clock time ---
for tool in TOOL_ORDER:
    sub = df[df["tool"] == tool]
    means = sub.groupby("read_count_M")["wall_time_s"].mean().reindex(read_counts_sorted)
    sds = sub.groupby("read_count_M")["wall_time_s"].std().reindex(read_counts_sorted).fillna(0)
    axes[0].errorbar(
        means.index, means.values, yerr=sds.values,
        marker="o", markersize=5, capsize=3, linewidth=1.5,
        color=TOOL_COLORS[tool], label=tool,
    )

axes[0].set_xlabel("Read count (millions)", fontsize=14)
axes[0].set_ylabel("Wall-clock time (seconds)", fontsize=14)
axes[0].set_title("Execution Time", fontsize=15)
axes[0].spines["top"].set_visible(False)
axes[0].spines["right"].set_visible(False)
axes[0].tick_params(labelsize=12)
axes[0].legend(fontsize=9, loc="upper left")

# --- Panel B: Peak memory ---
for tool in TOOL_ORDER:
    sub = df[df["tool"] == tool]
    means = sub.groupby("read_count_M")["max_rss_mb"].mean().reindex(read_counts_sorted)
    sds = sub.groupby("read_count_M")["max_rss_mb"].std().reindex(read_counts_sorted).fillna(0)
    axes[1].errorbar(
        means.index, means.values, yerr=sds.values,
        marker="o", markersize=5, capsize=3, linewidth=1.5,
        color=TOOL_COLORS[tool], label=tool,
    )

axes[1].set_xlabel("Read count (millions)", fontsize=14)
axes[1].set_ylabel("Peak RSS memory (MB)", fontsize=14)
axes[1].set_title("Peak Memory", fontsize=15)
axes[1].spines["top"].set_visible(False)
axes[1].spines["right"].set_visible(False)
axes[1].tick_params(labelsize=12)
axes[1].legend(fontsize=9, loc="upper left")

plt.tight_layout()
fig.savefig(out_plot, dpi=300, bbox_inches="tight")
plt.close(fig)
plt.close()

# ---------------------------------------------------------------------------
# Export summary TSV
# ---------------------------------------------------------------------------
summary_rows = []
for tool in TOOL_ORDER:
    sub = df[df["tool"] == tool]
    for rc in read_counts_sorted:
        rc_sub = sub[sub["read_count_M"] == rc]
        summary_rows.append({
            "tool": tool,
            "read_count": int(rc_sub["read_count"].iloc[0]),
            "read_count_M": rc,
            "wall_time_mean_s": rc_sub["wall_time_s"].mean(),
            "wall_time_sd_s": rc_sub["wall_time_s"].std(),
            "max_rss_mean_mb": rc_sub["max_rss_mb"].mean(),
            "max_rss_sd_mb": rc_sub["max_rss_mb"].std(),
            "n_replicates": len(rc_sub),
        })

df_summary = pd.DataFrame(summary_rows)
df_summary.to_csv(out_tsv, sep="\t", index=False, float_format="%.2f")

print(f"Downsample benchmark plot saved to {out_plot}")
print(f"Downsample benchmark summary saved to {out_tsv}")
