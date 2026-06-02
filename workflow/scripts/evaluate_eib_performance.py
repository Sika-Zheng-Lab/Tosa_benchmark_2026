"""
Evaluate Tosa performance across thread counts and anchor lengths (EIB benchmark).

Loads benchmark TSV files for all combinations of:
  anchor_length × library_type × strandedness × threads

Produces:
  - Thread scaling plot (wall time)
  - Thread scaling plot (peak memory / RSS)
  - Speedup factor plot (T1 / Tn vs ideal linear)
  - Performance summary table
"""

import os
import re
import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# ---------------------------------------------------------------------------
# Snakemake I/O
# ---------------------------------------------------------------------------
in_benchmark_files = snakemake.input.benchmarks      # list of benchmark TSV files
benchmark_labels = snakemake.params.benchmark_labels  # list of "b{AL}_{LT}_{SD}_{threads}"
thread_values = snakemake.params.thread_values        # list of thread counts as strings
anchor_lengths = snakemake.params.anchor_lengths      # list of anchor length strings

out_time_plot = snakemake.output.time_plot
out_memory_plot = snakemake.output.memory_plot
out_speedup_plot = snakemake.output.speedup_plot
out_summary = snakemake.output.summary

os.makedirs(os.path.dirname(out_time_plot), exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Load benchmark data
# ---------------------------------------------------------------------------
frames = []
for i, label in enumerate(benchmark_labels):
    # label format: b{anchor}_{LT}_{SD}_{threads}, e.g., b1_PE_XS_4
    m = re.match(r"^b(\d+)_([^_]+)_([^_]+)_(\d+)$", label)
    if not m:
        raise ValueError(f"Cannot parse benchmark label: {label}")
    anchor = m.group(1)
    library_type = m.group(2)
    strandedness = m.group(3)
    threads = m.group(4)

    df_b = pd.read_csv(in_benchmark_files[i], sep="\t")
    df_b["anchor_length"] = int(anchor)
    df_b["library_type"] = library_type
    df_b["strandedness"] = strandedness
    df_b["threads"] = int(threads)
    df_b["dataset"] = f"b{anchor}_{library_type}_{strandedness}"
    df_b["wall_time_s"] = df_b["s"]
    df_b["max_rss_mb"] = df_b["max_rss"]
    frames.append(
        df_b[["anchor_length", "library_type", "strandedness", "threads",
              "dataset", "wall_time_s", "max_rss_mb"]]
    )

df = pd.concat(frames, ignore_index=True)
datasets = sorted(df["dataset"].unique())
thread_ints = sorted(df["threads"].unique())

# ---------------------------------------------------------------------------
# 2. Thread scaling — Wall Time
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 6))

palette = sns.color_palette("Set2", len(datasets))
for idx, ds in enumerate(datasets):
    sub = df[df["dataset"] == ds]
    medians = sub.groupby("threads")["wall_time_s"].median().reset_index()
    ci = sub.groupby("threads")["wall_time_s"].agg(
        lo=lambda x: np.percentile(x, 2.5),
        hi=lambda x: np.percentile(x, 97.5),
    ).reset_index()
    merged = medians.merge(ci, on="threads")

    ax.errorbar(
        merged["threads"], merged["wall_time_s"],
        yerr=[merged["wall_time_s"] - merged["lo"], merged["hi"] - merged["wall_time_s"]],
        marker="o", capsize=4, label=ds, color=palette[idx % len(palette)],
    )

ax.set_xlabel("Threads")
ax.set_ylabel("Wall Clock Time (seconds)")
ax.set_title("Tosa Thread Scaling — Wall Time")
ax.set_xticks(thread_ints)
ax.legend(title="Dataset", bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(out_time_plot, dpi=150, bbox_inches="tight")
plt.close()

# ---------------------------------------------------------------------------
# 3. Thread scaling — Peak Memory (RSS)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 6))

for idx, ds in enumerate(datasets):
    sub = df[df["dataset"] == ds]
    medians = sub.groupby("threads")["max_rss_mb"].median().reset_index()
    ci = sub.groupby("threads")["max_rss_mb"].agg(
        lo=lambda x: np.percentile(x, 2.5),
        hi=lambda x: np.percentile(x, 97.5),
    ).reset_index()
    merged = medians.merge(ci, on="threads")

    ax.errorbar(
        merged["threads"], merged["max_rss_mb"],
        yerr=[merged["max_rss_mb"] - merged["lo"], merged["hi"] - merged["max_rss_mb"]],
        marker="s", capsize=4, label=ds, color=palette[idx % len(palette)],
    )

ax.set_xlabel("Threads")
ax.set_ylabel("Peak Memory / RSS (MB)")
ax.set_title("Tosa Thread Scaling — Peak Memory")
ax.set_xticks(thread_ints)
ax.legend(title="Dataset", bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(out_memory_plot, dpi=150, bbox_inches="tight")
plt.close()

# ---------------------------------------------------------------------------
# 4. Speedup factor plot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 6))

for idx, ds in enumerate(datasets):
    sub = df[df["dataset"] == ds]
    medians = sub.groupby("threads")["wall_time_s"].median().reset_index()
    medians = medians.sort_values("threads")

    t1 = medians.loc[medians["threads"] == min(thread_ints), "wall_time_s"].values
    if len(t1) == 0:
        continue
    t1 = t1[0]

    medians["speedup"] = t1 / medians["wall_time_s"]

    ax.plot(
        medians["threads"], medians["speedup"],
        marker="o", label=ds, color=palette[idx % len(palette)],
    )

# Ideal linear speedup
ax.plot(thread_ints, thread_ints, "--", color="gray", alpha=0.5, label="Ideal linear")
ax.set_xlabel("Threads")
ax.set_ylabel("Speedup (T₁ / Tₙ)")
ax.set_title("Tosa Thread Scaling — Speedup Factor")
ax.set_xticks(thread_ints)
ax.legend(title="Dataset", bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(out_speedup_plot, dpi=150, bbox_inches="tight")
plt.close()

# ---------------------------------------------------------------------------
# 5. Performance summary table
# ---------------------------------------------------------------------------
summary = df.groupby(["dataset", "anchor_length", "library_type", "strandedness", "threads"]).agg(
    median_time_s=("wall_time_s", "median"),
    mean_time_s=("wall_time_s", "mean"),
    std_time_s=("wall_time_s", "std"),
    median_rss_mb=("max_rss_mb", "median"),
    mean_rss_mb=("max_rss_mb", "mean"),
).reset_index()

# Add speedup column
speedup_rows = []
for ds in datasets:
    sub = summary[summary["dataset"] == ds].sort_values("threads")
    t1 = sub.loc[sub["threads"] == min(thread_ints), "median_time_s"].values
    if len(t1) == 0:
        for _, row in sub.iterrows():
            speedup_rows.append(np.nan)
        continue
    t1 = t1[0]
    for _, row in sub.iterrows():
        speedup_rows.append(round(t1 / row["median_time_s"], 2) if row["median_time_s"] > 0 else np.nan)

summary["speedup"] = speedup_rows

# Round numeric columns
for col in ["median_time_s", "mean_time_s", "std_time_s", "median_rss_mb", "mean_rss_mb"]:
    summary[col] = summary[col].round(4)

summary.to_csv(out_summary, sep="\t", index=False)

print("=== EIB Performance Summary ===")
print(summary.to_string(index=False))
