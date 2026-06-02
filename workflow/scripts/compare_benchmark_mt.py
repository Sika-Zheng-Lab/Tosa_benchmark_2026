"""
Compare multi-threaded benchmark results: Tosa (8 threads), featureCounts (8 threads),
and regtools (1 thread, single-threaded tool).

Generates a bar plot with mean ± 95% CI error bars for wall-clock time and peak RSS memory.
"""

import os
import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# ---------------------------------------------------------------------------
# Snakemake I/O
# ---------------------------------------------------------------------------
in_bench_tosa_mt = snakemake.input.bench_tosa_mt
in_bench_featurecounts_mt = snakemake.input.bench_featurecounts_mt
in_bench_regtools = snakemake.input.bench_regtools

out_benchmark = snakemake.output.benchmark

os.makedirs(os.path.dirname(out_benchmark), exist_ok=True)

# ---------------------------------------------------------------------------
# Load benchmark data (repeat produces multiple rows per TSV)
# ---------------------------------------------------------------------------
TOOL_ORDER = ["regtools (1 thread)", "featureCounts (8 threads)", "Tosa (8 threads)"]

bench_frames = []
for tool_label, path in [
    ("regtools (1 thread)", in_bench_regtools),
    ("featureCounts (8 threads)", in_bench_featurecounts_mt),
    ("Tosa (8 threads)", in_bench_tosa_mt),
]:
    df_b = pd.read_csv(path, sep="\t")
    df_b["tool"] = tool_label
    df_b["wall_time_s"] = df_b["s"]
    df_b["max_rss_mb"] = df_b["max_rss"]
    bench_frames.append(df_b[["tool", "wall_time_s", "max_rss_mb"]])

df_bench = pd.concat(bench_frames, ignore_index=True)

# ---------------------------------------------------------------------------
# Bar plot: wall time & peak memory (mean ± 95% CI)
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Wall time
sns.barplot(
    data=df_bench, x="tool", y="wall_time_s", ax=axes[0],
    palette="Set2", errorbar=("ci", 95), capsize=0.1, order=TOOL_ORDER,
)
axes[0].set_title("Wall Clock Time (multi-threaded)")
axes[0].set_ylabel("Time (seconds)")
axes[0].set_xlabel("")
axes[0].tick_params(axis="x", rotation=15)

# Memory
sns.barplot(
    data=df_bench, x="tool", y="max_rss_mb", ax=axes[1],
    palette="Set2", errorbar=("ci", 95), capsize=0.1, order=TOOL_ORDER,
)
axes[1].set_title("Peak Memory / RSS (multi-threaded)")
axes[1].set_ylabel("Memory (MB)")
axes[1].set_xlabel("")
axes[1].tick_params(axis="x", rotation=15)

plt.tight_layout()
plt.savefig(out_benchmark, dpi=150)
plt.close()
