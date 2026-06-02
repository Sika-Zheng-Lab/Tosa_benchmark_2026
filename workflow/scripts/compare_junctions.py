"""
Compare junction counts across STAR, Tosa, regtools, and featureCounts.

Generates:
- Merged junction count table
- Pairwise scatter plots with correlation
- UpSet plot of junction detection overlap
- Correlation table
- Jaccard index bar plot
- Benchmark bar plot (time & memory)
"""

import os
import itertools
import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from upsetplot import from_memberships, UpSet
from scipy import stats

# ---------------------------------------------------------------------------
# Snakemake I/O
# ---------------------------------------------------------------------------
in_star = snakemake.input.star
in_tosa = snakemake.input.tosa
in_regtools = snakemake.input.regtools
in_featurecounts = snakemake.input.featurecounts
in_bench_tosa = snakemake.input.bench_tosa
in_bench_regtools = snakemake.input.bench_regtools
in_bench_featurecounts = snakemake.input.bench_featurecounts

out_counts = snakemake.output.counts
out_scatter = snakemake.output.scatter
out_venn = snakemake.output.venn
out_correlation = snakemake.output.correlation
out_jaccard = snakemake.output.jaccard
out_benchmark = snakemake.output.benchmark

os.makedirs(os.path.dirname(out_counts), exist_ok=True)

# ---------------------------------------------------------------------------
# Tool names (display labels)
# ---------------------------------------------------------------------------
TOOLS = ["STAR", "regtools", "featureCounts", "Tosa"]

# ---------------------------------------------------------------------------
# 1. Load normalized data
# ---------------------------------------------------------------------------
df_star = pd.read_csv(in_star, sep="\t")
df_tosa = pd.read_csv(in_tosa, sep="\t")
df_regtools = pd.read_csv(in_regtools, sep="\t")
df_featurecounts = pd.read_csv(in_featurecounts, sep="\t")

# ---------------------------------------------------------------------------
# 2. Merge into a single table (outer join on junction)
# ---------------------------------------------------------------------------
merged = (
    df_star[["junction", "count"]]
    .rename(columns={"count": "STAR"})
    .merge(
        df_tosa[["junction", "count"]].rename(columns={"count": "Tosa"}),
        on="junction",
        how="outer",
    )
    .merge(
        df_regtools[["junction", "count"]].rename(columns={"count": "regtools"}),
        on="junction",
        how="outer",
    )
    .merge(
        df_featurecounts[["junction", "count"]].rename(
            columns={"count": "featureCounts"}
        ),
        on="junction",
        how="outer",
    )
)
merged = merged.fillna(0).astype(
    {"STAR": int, "Tosa": int, "regtools": int, "featureCounts": int}
)
merged.to_csv(out_counts, sep="\t", index=False)

# ---------------------------------------------------------------------------
# 3. Pairwise scatter plots — 4×4 matrix layout
# ---------------------------------------------------------------------------
pairs = list(itertools.combinations(TOOLS, 2))
n_tools = len(TOOLS)
fig, axes = plt.subplots(n_tools, n_tools, figsize=(20, 20))

for i, row_tool in enumerate(TOOLS):
    for j, col_tool in enumerate(TOOLS):
        ax = axes[i][j]
        if i == j:
            # Diagonal: histogram of junction counts for this tool
            vals = merged[row_tool]
            vals_pos = vals[vals > 0]
            if len(vals_pos) > 0:
                ax.hist(
                    np.log10(vals_pos + 1), bins=50, color="steelblue", alpha=0.7
                )
                ax.set_xlabel("log10(count + 1)", fontsize=13)
                ax.set_ylabel("Frequency", fontsize=13)
            ax.set_title(row_tool, fontsize=16, fontweight="bold")
            ax.tick_params(axis="both", labelsize=11)
        else:
            # Off-diagonal: scatter plot (col_tool on x-axis, row_tool on y-axis)
            x = merged[col_tool]
            y = merged[row_tool]
            ax.scatter(x, y, alpha=0.3, s=8, edgecolors="none")

            # Diagonal line
            max_val = max(x.max(), y.max())
            if max_val > 0:
                ax.plot(
                    [0, max_val], [0, max_val], "r--", linewidth=0.8
                )

            # Correlation (only where both > 0)
            mask = (x > 0) & (y > 0)
            if mask.sum() > 2:
                r_pearson, _ = stats.pearsonr(x[mask], y[mask])
                r_spearman, _ = stats.spearmanr(x[mask], y[mask])
                ax.set_title(
                    f"r={r_pearson:.3f}  ρ={r_spearman:.3f}", fontsize=13
                )
            ax.tick_params(axis="both", labelsize=11)

        # Axis labels on edges only
        if i == n_tools - 1:
            ax.set_xlabel(col_tool, fontsize=14)
        else:
            ax.set_xlabel("")
        if j == 0:
            if i != j:
                ax.set_ylabel(row_tool, fontsize=14)
        else:
            ax.set_ylabel("")

plt.suptitle("Pairwise Junction Count Comparison", fontsize=20, y=1.01)
plt.tight_layout()
plt.savefig(out_scatter, dpi=150, bbox_inches="tight")
plt.close()

# ---------------------------------------------------------------------------
# 4. UpSet plot (replaces Venn diagram for 4-set overlap)
# ---------------------------------------------------------------------------
set_star = set(df_star["junction"])
set_tosa = set(df_tosa["junction"])
set_regtools = set(df_regtools["junction"])
set_featurecounts = set(df_featurecounts["junction"])

all_junctions = set_star | set_tosa | set_regtools | set_featurecounts

memberships = []
for junc in all_junctions:
    member = []
    if junc in set_star:
        member.append("STAR")
    if junc in set_tosa:
        member.append("Tosa")
    if junc in set_regtools:
        member.append("regtools")
    if junc in set_featurecounts:
        member.append("featureCounts")
    memberships.append(member)

upset_data = from_memberships(memberships)
fig = plt.figure(figsize=(12, 7))
upset = UpSet(
    upset_data,
    subset_size="count",
    show_counts=True,
    sort_by="cardinality",
    sort_categories_by="cardinality",
)
upset.plot(fig=fig)
fig.suptitle("Junction Detection Overlap", fontsize=14, y=1.02)
fig.savefig(out_venn, dpi=150, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# 5. Correlation table — all 6 pairs
# ---------------------------------------------------------------------------
corr_rows = []
for x_tool, y_tool in pairs:
    x = merged[x_tool]
    y = merged[y_tool]
    mask = (x > 0) & (y > 0)
    n_common = mask.sum()
    if n_common > 2:
        r_pearson, p_pearson = stats.pearsonr(x[mask], y[mask])
        r_spearman, p_spearman = stats.spearmanr(x[mask], y[mask])
    else:
        r_pearson = p_pearson = r_spearman = p_spearman = float("nan")

    # Jaccard (set-based: presence/absence of junctions)
    s1 = set(merged.loc[merged[x_tool] > 0, "junction"])
    s2 = set(merged.loc[merged[y_tool] > 0, "junction"])
    jaccard = len(s1 & s2) / len(s1 | s2) if len(s1 | s2) > 0 else float("nan")

    # Generalized Jaccard (count-based: sum(min) / sum(max))
    min_sum = np.minimum(x, y).sum()
    max_sum = np.maximum(x, y).sum()
    jaccard_count = min_sum / max_sum if max_sum > 0 else float("nan")

    corr_rows.append(
        {
            "pair": f"{x_tool}_vs_{y_tool}",
            "n_common_junctions": n_common,
            "pearson_r": r_pearson,
            "pearson_p": p_pearson,
            "spearman_r": r_spearman,
            "spearman_p": p_spearman,
            "jaccard_index": jaccard,
            "jaccard_count": jaccard_count,
        }
    )

df_corr = pd.DataFrame(corr_rows)
df_corr.to_csv(out_correlation, sep="\t", index=False)

# ---------------------------------------------------------------------------
# 6. Jaccard index heatmaps
# ---------------------------------------------------------------------------

# Build symmetric matrices from pairwise values
jaccard_set_matrix = pd.DataFrame(
    np.eye(n_tools), index=TOOLS, columns=TOOLS
)
jaccard_cnt_matrix = pd.DataFrame(
    np.eye(n_tools), index=TOOLS, columns=TOOLS
)

for _, row in df_corr.iterrows():
    pair = row["pair"]
    t1, t2 = pair.split("_vs_")
    jaccard_set_matrix.loc[t1, t2] = row["jaccard_index"]
    jaccard_set_matrix.loc[t2, t1] = row["jaccard_index"]
    jaccard_cnt_matrix.loc[t1, t2] = row["jaccard_count"]
    jaccard_cnt_matrix.loc[t2, t1] = row["jaccard_count"]

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Compute data-driven colour limits for better contrast
_off_diag_set = jaccard_set_matrix.values[~np.eye(n_tools, dtype=bool)]
_off_diag_cnt = jaccard_cnt_matrix.values[~np.eye(n_tools, dtype=bool)]
_vmin_set = max(0, np.floor(_off_diag_set.min() * 20) / 20 - 0.05)
_vmin_cnt = max(0, np.floor(_off_diag_cnt.min() * 20) / 20 - 0.05)

# (a) Set-based Jaccard heatmap
sns.heatmap(
    jaccard_set_matrix, annot=True, fmt=".4f", cmap="YlOrRd",
    vmin=_vmin_set, vmax=1, square=True, linewidths=0.5, ax=axes[0],
    cbar_kws={"label": "Jaccard Index"},
)
axes[0].set_title("Set-based Jaccard\n(junction presence/absence)")

# (b) Count-based Generalized Jaccard heatmap
sns.heatmap(
    jaccard_cnt_matrix, annot=True, fmt=".4f", cmap="YlOrRd",
    vmin=_vmin_cnt, vmax=1, square=True, linewidths=0.5, ax=axes[1],
    cbar_kws={"label": "Generalized Jaccard Index"},
)
axes[1].set_title("Count-based Jaccard\n(Σmin / Σmax of counts)")

plt.tight_layout()
plt.savefig(out_jaccard, dpi=150)
plt.close()

# ---------------------------------------------------------------------------
# 7. Benchmark bar plot (time & memory) — repeat-aware with error bars
# ---------------------------------------------------------------------------
BENCH_ORDER = ["regtools", "featureCounts", "Tosa"]

bench_frames = []
for tool, path in [
    ("regtools", in_bench_regtools),
    ("featureCounts", in_bench_featurecounts),
    ("Tosa", in_bench_tosa),
]:
    df_b = pd.read_csv(path, sep="\t")
    df_b["tool"] = tool
    df_b["wall_time_s"] = df_b["s"]
    df_b["max_rss_mb"] = df_b["max_rss"]
    bench_frames.append(df_b[["tool", "wall_time_s", "max_rss_mb"]])

df_bench = pd.concat(bench_frames, ignore_index=True)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Wall time
sns.barplot(
    data=df_bench, x="tool", y="wall_time_s", ax=axes[0],
    palette="Set2", errorbar=("ci", 95), capsize=0.1, order=BENCH_ORDER,
)
axes[0].set_title("Wall Clock Time (1 thread)")
axes[0].set_ylabel("Time (seconds)")
axes[0].set_xlabel("")

# Memory
sns.barplot(
    data=df_bench, x="tool", y="max_rss_mb", ax=axes[1],
    palette="Set2", errorbar=("ci", 95), capsize=0.1, order=BENCH_ORDER,
)
axes[1].set_title("Peak Memory / RSS (1 thread)")
axes[1].set_ylabel("Memory (MB)")
axes[1].set_xlabel("")

plt.tight_layout()
plt.savefig(out_benchmark, dpi=150)
plt.close()
