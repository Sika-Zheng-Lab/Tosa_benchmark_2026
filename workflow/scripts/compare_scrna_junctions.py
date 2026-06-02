"""
Compare scRNA-seq junction counts across STARsolo, Tosa single, and regtools.

Supports two modes:
  1. Simulation mode: Compare each tool against ground truth
  2. Cross-tool mode: Pairwise comparison between tools

Input: Normalized sparse TSV files (barcode\tjunction\tcount) from each tool.
Ground truth (simulation only): expected_junction_counts.tsv

Generates:
- Per-cell sensitivity/specificity (simulation mode)
- Cross-tool correlation at cell × junction level
- Junction detection overlap (UpSet plot)
- Per-cell concordance distribution
- Benchmark performance bar plots
"""

import os
import sys
import itertools

import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# ---------------------------------------------------------------------------
# Snakemake I/O
# ---------------------------------------------------------------------------
in_starsolo = snakemake.input.starsolo
in_tosa = snakemake.input.tosa
in_regtools = snakemake.input.regtools

# Optional ground truth (only in simulation mode)
in_expected = snakemake.input.get("expected", None)

# Benchmarks (optional)
in_bench_starsolo = snakemake.input.get("bench_starsolo", None)
in_bench_tosa = snakemake.input.get("bench_tosa", None)
in_bench_regtools = snakemake.input.get("bench_regtools", None)

out_counts = snakemake.output.counts
out_scatter = snakemake.output.scatter
out_venn = snakemake.output.venn
out_correlation = snakemake.output.correlation
out_per_cell = snakemake.output.per_cell
out_benchmark = snakemake.output.get("benchmark", None)

# Simulation-specific outputs
out_accuracy = snakemake.output.get("accuracy", None)

os.makedirs(os.path.dirname(out_counts), exist_ok=True)

TOOLS = ["STARsolo", "Tosa", "regtools"]

# ---------------------------------------------------------------------------
# 1. Load normalized sparse data
# ---------------------------------------------------------------------------

def load_sparse_tsv(path):
    """Load sparse TSV (barcode, junction, count) into a DataFrame."""
    df = pd.read_csv(path, sep="\t")
    df.columns = ["barcode", "junction", "count"]
    return df


df_starsolo = load_sparse_tsv(in_starsolo)
df_tosa = load_sparse_tsv(in_tosa)
df_regtools = load_sparse_tsv(in_regtools)

tool_dfs = {
    "STARsolo": df_starsolo,
    "Tosa": df_tosa,
    "regtools": df_regtools,
}

# ---------------------------------------------------------------------------
# 2. Build unified cell × junction matrix (pseudo-bulk + per-cell)
# ---------------------------------------------------------------------------

# Pseudo-bulk: sum across all cells per junction
def to_pseudobulk(df):
    return df.groupby("junction")["count"].sum().reset_index()


pseudobulk = {}
for tool, df in tool_dfs.items():
    pb = to_pseudobulk(df)
    pseudobulk[tool] = pb.set_index("junction")["count"]

# Merge pseudo-bulk counts
merged = pd.DataFrame(pseudobulk).fillna(0).astype(int)
merged.index.name = "junction"
merged.to_csv(out_counts, sep="\t")

print(f"Total junctions detected: {len(merged)}", file=sys.stderr)
for tool in TOOLS:
    n_detected = (merged[tool] > 0).sum()
    print(f"  {tool}: {n_detected} junctions", file=sys.stderr)


# ---------------------------------------------------------------------------
# 3. Pairwise scatter plots (pseudo-bulk level)
# ---------------------------------------------------------------------------
n_tools = len(TOOLS)
fig, axes = plt.subplots(n_tools, n_tools, figsize=(15, 15))

for i, row_tool in enumerate(TOOLS):
    for j, col_tool in enumerate(TOOLS):
        ax = axes[i][j]
        if i == j:
            vals = merged[row_tool]
            vals_pos = vals[vals > 0]
            if len(vals_pos) > 0:
                ax.hist(np.log10(vals_pos + 1), bins=30, color="steelblue", alpha=0.7)
            ax.set_title(row_tool, fontsize=10)
            ax.set_xlabel("log10(count + 1)")
        elif i > j:
            x = merged[col_tool]
            y = merged[row_tool]
            mask = (x > 0) | (y > 0)
            ax.scatter(
                np.log10(x[mask] + 1),
                np.log10(y[mask] + 1),
                s=3,
                alpha=0.3,
            )
            if mask.sum() > 2:
                r, _ = stats.pearsonr(x[mask], y[mask])
                rho, _ = stats.spearmanr(x[mask], y[mask])
                ax.set_title(f"r={r:.3f} ρ={rho:.3f}", fontsize=8)
            ax.set_xlabel(col_tool, fontsize=8)
            ax.set_ylabel(row_tool, fontsize=8)
        else:
            ax.axis("off")

fig.suptitle("Pseudo-bulk Junction Count Comparison (scRNA-seq)", fontsize=14)
fig.tight_layout()
fig.savefig(out_scatter, dpi=150, bbox_inches="tight")
plt.close(fig)


# ---------------------------------------------------------------------------
# 4. UpSet plot — junction detection overlap
# ---------------------------------------------------------------------------
try:
    from upsetplot import UpSet

    # Build a Boolean indicator DataFrame for UpSet: one row per junction,
    # one column per tool indicating whether the junction was detected.
    indicators = pd.DataFrame(
        {t: merged[t] > 0 for t in TOOLS},
        index=merged.index,
    )
    # UpSet expects a MultiIndex of booleans
    upset_series = indicators.set_index(TOOLS).assign(
        count=merged[TOOLS].max(axis=1)
    )["count"]
    fig_upset = plt.figure(figsize=(10, 6))
    upset = UpSet(upset_series, show_counts=True, sort_by="cardinality")
    upset.plot(fig=fig_upset)
    fig_upset.savefig(out_venn, dpi=150, bbox_inches="tight")
    plt.close(fig_upset)
except (ImportError, Exception) as e:
    print(f"UpSet plot fallback: {e}", file=sys.stderr)
    # Fallback: simple Venn-style bar chart
    fig, ax = plt.subplots(figsize=(8, 5))
    for t in TOOLS:
        n = (merged[t] > 0).sum()
        ax.bar(t, n)
    ax.set_ylabel("Number of junctions detected")
    ax.set_title("Junction Detection per Tool (scRNA-seq)")
    fig.savefig(out_venn, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# 5. Correlation table
# ---------------------------------------------------------------------------
pairs = list(itertools.combinations(TOOLS, 2))
corr_rows = []
for t1, t2 in pairs:
    x = merged[t1]
    y = merged[t2]
    mask = (x > 0) | (y > 0)
    if mask.sum() > 2:
        r, _ = stats.pearsonr(x[mask], y[mask])
        rho, _ = stats.spearmanr(x[mask], y[mask])
    else:
        r, rho = np.nan, np.nan
    n_both = ((x > 0) & (y > 0)).sum()
    n_union = mask.sum()
    jaccard = n_both / n_union if n_union > 0 else 0
    corr_rows.append({
        "tool1": t1,
        "tool2": t2,
        "pearson_r": round(r, 4),
        "spearman_rho": round(rho, 4),
        "jaccard": round(jaccard, 4),
        "n_both": n_both,
        "n_union": n_union,
    })

df_corr = pd.DataFrame(corr_rows)
df_corr.to_csv(out_correlation, sep="\t", index=False)


# ---------------------------------------------------------------------------
# 6. Per-cell concordance
# ---------------------------------------------------------------------------
# For each cell, compute the pairwise Pearson correlation of junction counts

# Build per-cell junction vectors for each tool
all_barcodes = set()
all_junctions = set()
for df in tool_dfs.values():
    all_barcodes.update(df["barcode"])
    all_junctions.update(df["junction"])

all_barcodes = sorted(all_barcodes)
all_junctions = sorted(all_junctions)

# Create junction index
junc_idx = {j: i for i, j in enumerate(all_junctions)}
n_junc = len(all_junctions)

per_cell_corrs = []

for bc in all_barcodes:
    vectors = {}
    for tool, df in tool_dfs.items():
        vec = np.zeros(n_junc)
        sub = df[df["barcode"] == bc]
        for _, row in sub.iterrows():
            idx = junc_idx.get(row["junction"])
            if idx is not None:
                vec[idx] = row["count"]
        vectors[tool] = vec

    for t1, t2 in itertools.combinations(TOOLS, 2):
        v1 = vectors[t1]
        v2 = vectors[t2]
        mask = (v1 > 0) | (v2 > 0)
        if mask.sum() >= 2:
            r, _ = stats.pearsonr(v1[mask], v2[mask])
        else:
            r = np.nan
        per_cell_corrs.append({
            "barcode": bc,
            "tool1": t1,
            "tool2": t2,
            "pearson_r": r,
        })

df_per_cell = pd.DataFrame(per_cell_corrs)

# Pair colours — consistent across figures
PAIR_COLORS = {
    ("STARsolo", "Tosa"): "#9467bd",
    ("STARsolo", "regtools"): "#ff7f0e",
    ("Tosa", "regtools"): "#1f77b4",
}

matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial"]

fig, ax = plt.subplots(figsize=(5, 4.5))
for (t1, t2), grp in df_per_cell.groupby(["tool1", "tool2"]):
    vals = grp["pearson_r"].dropna()
    color = PAIR_COLORS.get((t1, t2), "#999999")
    if len(vals) > 0:
        val_range = vals.max() - vals.min()
        if val_range == 0 or np.isclose(val_range, 0):
            center = vals.iloc[0]
            bin_edges = [center - 0.005, center + 0.005]
            ax.hist(vals, bins=bin_edges, alpha=0.6, label=f"{t1} vs {t2}",
                    color=color, edgecolor="black", linewidth=0.5)
        else:
            n_bins = min(20, max(1, len(vals.unique())))
            ax.hist(vals, bins=n_bins, alpha=0.6, label=f"{t1} vs {t2}",
                    color=color, edgecolor="black", linewidth=0.5)
ax.set_xlabel("Per-cell Pearson r", fontsize=14)
ax.set_ylabel("Number of cells", fontsize=14)
ax.set_title("Per-cell Junction Count Concordance", fontsize=15)
ax.legend(fontsize=11)
ax.tick_params(labelsize=12)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig(out_per_cell, dpi=300, bbox_inches="tight")
plt.close(fig)


# ---------------------------------------------------------------------------
# 7. Simulation accuracy (if ground truth provided)
# ---------------------------------------------------------------------------
if in_expected and out_accuracy:
    df_expected = pd.read_csv(in_expected, sep="\t")
    df_expected.columns = ["barcode", "junction", "count"]

    accuracy_rows = []
    for tool, df_tool in tool_dfs.items():
        # Merge with expected on (barcode, junction)
        merged_gt = df_expected.merge(
            df_tool,
            on=["barcode", "junction"],
            how="outer",
            suffixes=("_expected", "_observed"),
        )
        merged_gt = merged_gt.fillna(0)

        expected_counts = merged_gt["count_expected"].values
        observed_counts = merged_gt["count_observed"].values

        # Sensitivity: fraction of expected (count > 0) that are detected
        expected_positive = expected_counts > 0
        observed_positive = observed_counts > 0

        tp = (expected_positive & observed_positive).sum()
        fn = (expected_positive & ~observed_positive).sum()
        fp = (~expected_positive & observed_positive).sum()

        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        f1 = 2 * precision * sensitivity / (precision + sensitivity) if (precision + sensitivity) > 0 else 0

        # Exact match rate (among expected positives)
        exact_match = (expected_counts == observed_counts)[expected_positive].mean()

        # RMSE
        rmse = np.sqrt(((expected_counts - observed_counts) ** 2).mean())

        # Correlation
        mask = expected_positive | observed_positive
        if mask.sum() > 2:
            pearson_r, _ = stats.pearsonr(expected_counts[mask], observed_counts[mask])
            spearman_rho, _ = stats.spearmanr(expected_counts[mask], observed_counts[mask])
        else:
            pearson_r, spearman_rho = np.nan, np.nan

        accuracy_rows.append({
            "tool": tool,
            "true_positives": tp,
            "false_negatives": fn,
            "false_positives": fp,
            "sensitivity": round(sensitivity, 4),
            "precision": round(precision, 4),
            "f1_score": round(f1, 4),
            "exact_match_rate": round(exact_match, 4),
            "rmse": round(rmse, 4),
            "pearson_r": round(pearson_r, 4),
            "spearman_rho": round(spearman_rho, 4),
        })

    df_accuracy = pd.DataFrame(accuracy_rows)
    df_accuracy.to_csv(out_accuracy, sep="\t", index=False)
    print("Simulation accuracy:", file=sys.stderr)
    print(df_accuracy.to_string(index=False), file=sys.stderr)


# ---------------------------------------------------------------------------
# 8. Benchmark performance plots (if benchmark files provided)
# ---------------------------------------------------------------------------
if out_benchmark:
    bench_data = {}
    for tool, path in [("STARsolo", in_bench_starsolo),
                        ("Tosa", in_bench_tosa),
                        ("regtools", in_bench_regtools)]:
        if path:
            df_b = pd.read_csv(path, sep="\t")
            bench_data[tool] = {
                "time_mean": df_b["s"].mean(),
                "time_std": df_b["s"].std(),
                "mem_mean": df_b["max_rss"].mean(),
                "mem_std": df_b["max_rss"].std(),
            }

    if bench_data:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        tools = list(bench_data.keys())
        colors = {"STARsolo": "#4C72B0", "Tosa": "#DD8452", "regtools": "#55A868"}
        bar_colors = [colors.get(t, "#999999") for t in tools]
        times = [bench_data[t]["time_mean"] for t in tools]
        time_errs = [bench_data[t]["time_std"] for t in tools]
        mems = [bench_data[t]["mem_mean"] for t in tools]
        mem_errs = [bench_data[t]["mem_std"] for t in tools]

        ax1.bar(tools, times, yerr=time_errs, capsize=5, color=bar_colors)
        ax1.set_ylabel("Wall-clock time (s)")
        ax1.set_title("Runtime Comparison (scRNA-seq)")

        ax2.bar(tools, mems, yerr=mem_errs, capsize=5, color=bar_colors)
        ax2.set_ylabel("Peak RSS (MB)")
        ax2.set_title("Memory Comparison (scRNA-seq)")

        fig.tight_layout()
        fig.savefig(out_benchmark, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        # Create placeholder to satisfy Snakemake output
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No benchmark data available", ha="center", va="center")
        fig.savefig(out_benchmark, dpi=150, bbox_inches="tight")
        plt.close(fig)
