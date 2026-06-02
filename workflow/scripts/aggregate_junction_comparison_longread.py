"""Aggregate per-sample junction comparison results across all samples
for the long-read benchmark.

Generates:
- Per-sample Jaccard table (concatenated)
- Set-based Jaccard heatmap (mean across samples)
- Count-based Jaccard heatmap (mean across samples)
- Boxplot of per-sample Jaccard distributions
- UpSet plot of junction detection overlap (union across all samples)
- Junction summary table (per-tool stats)

Note: No STAR — long-read data is aligned with minimap2.
"""

import os
import itertools
import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
import matplotlib.pyplot as plt
import seaborn as sns
from upsetplot import from_memberships, UpSet

# ---------------------------------------------------------------------------
# Snakemake I/O
# ---------------------------------------------------------------------------
in_jaccards = snakemake.input.jaccards
in_memberships = snakemake.input.memberships
in_norm_tosa = snakemake.input.norm_tosa
in_norm_regtools = snakemake.input.norm_regtools
in_norm_featurecounts = snakemake.input.norm_featurecounts

out_per_sample = snakemake.output.per_sample_jaccard
out_set_heatmap = snakemake.output.jaccard_set_heatmap
out_cnt_heatmap = snakemake.output.jaccard_count_heatmap
out_boxplot = snakemake.output.jaccard_boxplot
out_upset = snakemake.output.upset_plot
out_summary = snakemake.output.junction_summary

samples = snakemake.params.samples

for out_path in [out_per_sample, out_set_heatmap, out_cnt_heatmap,
                 out_boxplot, out_upset, out_summary]:
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

TOOLS = ["regtools", "featureCounts", "Tosa"]

# ---------------------------------------------------------------------------
# 1. Concatenate per-sample Jaccard tables
# ---------------------------------------------------------------------------
frames = []
for path in in_jaccards:
    frames.append(pd.read_csv(path, sep="\t"))
df_all = pd.concat(frames, ignore_index=True)
df_all.to_csv(out_per_sample, sep="\t", index=False)

# ---------------------------------------------------------------------------
# 2. Set-based Jaccard heatmap (mean across samples)
# ---------------------------------------------------------------------------
n_tools = len(TOOLS)

jaccard_set_matrix = pd.DataFrame(np.eye(n_tools), index=TOOLS, columns=TOOLS)
jaccard_cnt_matrix = pd.DataFrame(np.eye(n_tools), index=TOOLS, columns=TOOLS)

for (t1, t2), grp in df_all.groupby(["tool1", "tool2"]):
    mean_set = grp["set_jaccard"].mean()
    mean_cnt = grp["count_jaccard"].mean()
    jaccard_set_matrix.loc[t1, t2] = mean_set
    jaccard_set_matrix.loc[t2, t1] = mean_set
    jaccard_cnt_matrix.loc[t1, t2] = mean_cnt
    jaccard_cnt_matrix.loc[t2, t1] = mean_cnt

# Data-driven colour limits
_off_diag_set = jaccard_set_matrix.values[~np.eye(n_tools, dtype=bool)]
_off_diag_cnt = jaccard_cnt_matrix.values[~np.eye(n_tools, dtype=bool)]
_vmin_set = max(0, np.floor(_off_diag_set.min() * 20) / 20 - 0.05)
_vmin_cnt = max(0, np.floor(_off_diag_cnt.min() * 20) / 20 - 0.05)

fig, ax = plt.subplots(figsize=(5.5, 5))
sns.heatmap(
    jaccard_set_matrix, annot=True, fmt=".4f", cmap="YlOrRd",
    vmin=_vmin_set, vmax=1, square=True, linewidths=0.5, ax=ax,
    cbar_kws={"label": "Jaccard Index"},
    annot_kws={"size": 12},
)
ax.set_title("Set-based Jaccard Index\n(mean across samples, junction presence/absence)", fontsize=11)
ax.tick_params(labelsize=12)
plt.tight_layout()
plt.savefig(out_set_heatmap, dpi=300)
plt.close()

fig, ax = plt.subplots(figsize=(5.5, 5))
sns.heatmap(
    jaccard_cnt_matrix, annot=True, fmt=".4f", cmap="YlOrRd",
    vmin=_vmin_cnt, vmax=1, square=True, linewidths=0.5, ax=ax,
    cbar_kws={"label": "Generalized Jaccard Index"},
    annot_kws={"size": 12},
)
ax.set_title("Count-based Jaccard Index\n(mean across samples, Σmin / Σmax)", fontsize=11)
ax.tick_params(labelsize=12)
plt.tight_layout()
plt.savefig(out_cnt_heatmap, dpi=300)
plt.close()

# ---------------------------------------------------------------------------
# 3. Boxplot of per-sample Jaccard distributions
# ---------------------------------------------------------------------------
df_all["pair"] = df_all["tool1"] + " vs " + df_all["tool2"]

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

sns.boxplot(
    data=df_all, x="pair", y="set_jaccard", ax=axes[0],
    palette="Set2", order=sorted(df_all["pair"].unique()),
)
axes[0].set_title("Set-based Jaccard Index per Sample (long-read)")
axes[0].set_ylabel("Jaccard Index")
axes[0].set_xlabel("")
axes[0].tick_params(axis="x", rotation=45)

sns.boxplot(
    data=df_all, x="pair", y="count_jaccard", ax=axes[1],
    palette="Set2", order=sorted(df_all["pair"].unique()),
)
axes[1].set_title("Count-based Jaccard Index per Sample (long-read)")
axes[1].set_ylabel("Generalized Jaccard Index")
axes[1].set_xlabel("")
axes[1].tick_params(axis="x", rotation=45)

plt.tight_layout()
plt.savefig(out_boxplot, dpi=150)
plt.close()

# ---------------------------------------------------------------------------
# 4. UpSet plot (union of junctions across all samples per tool)
# ---------------------------------------------------------------------------
tool_junctions = {tool: set() for tool in TOOLS}
tool_inputs = {
    "Tosa": in_norm_tosa,
    "regtools": in_norm_regtools,
    "featureCounts": in_norm_featurecounts,
}

for tool, paths in tool_inputs.items():
    for path in paths:
        df = pd.read_csv(path, sep="\t")
        tool_junctions[tool].update(df["junction"].tolist())

all_junctions = set()
for s in tool_junctions.values():
    all_junctions |= s

memberships = []
for junc in all_junctions:
    member = [t for t in TOOLS if junc in tool_junctions[t]]
    memberships.append(member)

upset_data = from_memberships(memberships)
fig = plt.figure(figsize=(10, 6))
upset = UpSet(
    upset_data,
    subset_size="count",
    show_counts=True,
    sort_by="cardinality",
    sort_categories_by="cardinality",
)
upset.plot(fig=fig)
fig.suptitle(
    f"Junction Detection Overlap — Long-read\n(union across {len(samples)} samples)",
    fontsize=14, y=1.02,
)
fig.savefig(out_upset, dpi=150, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# 5. Junction summary table
# ---------------------------------------------------------------------------
summary_rows = []
for tool, paths in tool_inputs.items():
    counts_per_sample = []
    all_juncs = set()
    for path in paths:
        df = pd.read_csv(path, sep="\t")
        counts_per_sample.append(len(df))
        all_juncs.update(df["junction"].tolist())
    summary_rows.append({
        "tool": tool,
        "mean_junctions_per_sample": np.mean(counts_per_sample),
        "median_junctions_per_sample": np.median(counts_per_sample),
        "min_junctions_per_sample": np.min(counts_per_sample),
        "max_junctions_per_sample": np.max(counts_per_sample),
        "total_unique_junctions": len(all_juncs),
        "n_samples": len(paths),
    })

df_summary = pd.DataFrame(summary_rows)
df_summary.to_csv(out_summary, sep="\t", index=False)
