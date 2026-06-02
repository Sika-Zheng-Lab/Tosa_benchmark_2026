"""
Evaluate Tosa single-cell boundary counting accuracy against ground truth
(scRNA-seq EIB benchmark).

Compares Tosa boundary output (from the 1-thread run) with expected boundary
counts at the cell × boundary level for each dataset
(anchor_length x strandedness).

Produces:
  - accuracy_summary.tsv: overall metrics per dataset
  - per_cell_accuracy.tsv: per-cell metrics
  - accuracy_by_boundary_type.png
  - accuracy_by_anchor_length.png
  - per_cell_accuracy.png
  - accuracy_heatmap.png
"""

import os
import re
import sys

import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# ---------------------------------------------------------------------------
# Snakemake I/O
# ---------------------------------------------------------------------------
in_tosa_files = snakemake.input.tosa_boundaries
in_expected_files = snakemake.input.expected_boundaries  # list, one per anchor length
in_gene_metadata = snakemake.input.gene_metadata

dataset_labels = snakemake.params.dataset_labels
anchor_lengths = snakemake.params.anchor_lengths
strandedness_modes = snakemake.params.strandedness_modes

out_accuracy_summary = snakemake.output.accuracy_summary
out_per_cell_accuracy = snakemake.output.per_cell_accuracy
out_by_boundary_type = snakemake.output.by_boundary_type
out_by_anchor = snakemake.output.by_anchor_length
out_per_cell_plot = snakemake.output.per_cell_plot
out_heatmap = snakemake.output.heatmap

os.makedirs(os.path.dirname(out_accuracy_summary), exist_ok=True)

# Build lookup: anchor_length -> expected counts file
_expected_lookup = {}
for idx_al, al in enumerate(anchor_lengths):
    _expected_lookup[al] = in_expected_files[idx_al]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_label(label):
    """Parse label like 'b1_XS' -> (anchor_length, strandedness)."""
    m = re.match(r"^b(\d+)_(.+)$", label)
    if not m:
        raise ValueError(f"Cannot parse dataset label: {label}")
    return m.group(1), m.group(2)


def boundary_center(boundary_str):
    """Extract (chrom, center_pos) from boundary string like 'chr:start-end'."""
    chrom, coords = boundary_str.rsplit(":", 1)
    start, end = coords.split("-")
    center = (int(start) + int(end)) // 2
    return chrom, center


def add_boundary_center(df, boundary_col="boundary"):
    """Add 'boundary_center' column (chrom:center) to dataframe."""
    centers = []
    for b in df[boundary_col]:
        chrom, center = boundary_center(str(b))
        centers.append(f"{chrom}:{center}")
    df = df.copy()
    df["boundary_center"] = centers
    return df


# ---------------------------------------------------------------------------
# 1. Load gene metadata
# ---------------------------------------------------------------------------
df_gene_meta = pd.read_csv(in_gene_metadata, sep="\t")

# ---------------------------------------------------------------------------
# 2. Process each dataset
# ---------------------------------------------------------------------------

all_accuracy_rows = []
all_per_cell_rows = []

for i, label in enumerate(dataset_labels):
    anchor_len, strandedness = parse_label(label)

    # Load anchor-length-specific expected boundary counts
    df_expected = pd.read_csv(_expected_lookup[anchor_len], sep="\t")

    # Load Tosa boundary output
    # tosa single -g outputs {prefix}_boundary_barcodes.tsv.gz with columns:
    # Feature, Strand, Barcode, Count  (Feature = boundary coordinate)
    tosa_path = in_tosa_files[i]
    df_tosa = pd.read_csv(tosa_path, sep="\t", compression="gzip")

    # Normalize column names
    # Tosa boundary_barcodes_detail.tsv.gz columns: Boundary, Type, Strand, Barcode, Count
    col_map = {}
    for c in df_tosa.columns:
        cl = c.lower().strip()
        if cl == "barcode":
            col_map[c] = "barcode"
        elif cl == "boundary":
            col_map[c] = "boundary"
        elif cl == "type":
            col_map[c] = "boundary_type"
        elif cl == "strand":
            col_map[c] = "strand"
        elif cl == "count":
            col_map[c] = "count"
    df_tosa = df_tosa.rename(columns=col_map)

    # Merge with expected using boundary center position (not exact span
    # coordinates), because Tosa -b <anchor_length> changes the output span
    # width while the boundary center remains identical.
    df_exp_c = add_boundary_center(df_expected, "boundary")
    df_tosa_c = add_boundary_center(df_tosa, "boundary")

    common_keys = ["barcode", "boundary_center"]
    if "boundary_type" in df_exp_c.columns and "boundary_type" in df_tosa_c.columns:
        common_keys = ["barcode", "boundary_center", "boundary_type"]

    merged = df_exp_c.merge(
        df_tosa_c[common_keys + ["count"]],
        on=common_keys,
        how="outer",
        suffixes=("_expected", "_observed"),
    )
    merged = merged.fillna(0)

    exp = merged["count_expected"].values
    obs = merged["count_observed"].values

    # Overall metrics
    exp_pos = exp > 0
    obs_pos = obs > 0
    tp = (exp_pos & obs_pos).sum()
    fn = (exp_pos & ~obs_pos).sum()
    fp = (~exp_pos & obs_pos).sum()

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    f1 = 2 * precision * sensitivity / (precision + sensitivity) if (precision + sensitivity) > 0 else 0

    exact_match = (exp == obs)[exp_pos].mean() if exp_pos.sum() > 0 else 0
    rmse = np.sqrt(((exp - obs) ** 2).mean())

    mask = exp_pos | obs_pos
    if mask.sum() > 2:
        pearson_r, _ = stats.pearsonr(exp[mask], obs[mask])
        spearman_rho, _ = stats.spearmanr(exp[mask], obs[mask])
    else:
        pearson_r, spearman_rho = np.nan, np.nan

    all_accuracy_rows.append({
        "dataset": label,
        "anchor_length": anchor_len,
        "strandedness": strandedness,
        "true_positives": tp,
        "false_negatives": fn,
        "false_positives": fp,
        "sensitivity": round(sensitivity, 4),
        "precision": round(precision, 4),
        "f1_score": round(f1, 4),
        "exact_match_rate": round(exact_match, 4),
        "rmse": round(rmse, 4),
        "pearson_r": round(pearson_r, 4) if not np.isnan(pearson_r) else np.nan,
        "spearman_rho": round(spearman_rho, 4) if not np.isnan(spearman_rho) else np.nan,
    })

    # Per-cell accuracy
    for bc in merged["barcode"].unique():
        cell_data = merged[merged["barcode"] == bc]
        cell_exp = cell_data["count_expected"].values
        cell_obs = cell_data["count_observed"].values
        cell_exp_pos = cell_exp > 0
        cell_obs_pos = cell_obs > 0

        cell_tp = (cell_exp_pos & cell_obs_pos).sum()
        cell_fn = (cell_exp_pos & ~cell_obs_pos).sum()
        cell_sensitivity = cell_tp / (cell_tp + cell_fn) if (cell_tp + cell_fn) > 0 else 0
        cell_exact = (cell_exp == cell_obs)[cell_exp_pos].mean() if cell_exp_pos.sum() > 0 else 0

        all_per_cell_rows.append({
            "dataset": label,
            "barcode": bc,
            "sensitivity": round(cell_sensitivity, 4),
            "exact_match_rate": round(cell_exact, 4),
            "n_expected": cell_exp_pos.sum(),
            "n_detected": cell_obs_pos.sum(),
        })

# ---------------------------------------------------------------------------
# 4. Write summary tables
# ---------------------------------------------------------------------------
df_accuracy = pd.DataFrame(all_accuracy_rows)
df_accuracy.to_csv(out_accuracy_summary, sep="\t", index=False)

df_per_cell = pd.DataFrame(all_per_cell_rows)
df_per_cell.to_csv(out_per_cell_accuracy, sep="\t", index=False)

print("Accuracy summary:", file=sys.stderr)
print(df_accuracy.to_string(index=False), file=sys.stderr)

# ---------------------------------------------------------------------------
# 5. Plots
# ---------------------------------------------------------------------------

# 5a. Accuracy by anchor length — simple bar chart (XS only, exact match rate)
sub = df_accuracy[df_accuracy["strandedness"] == "XS"].sort_values(
    "anchor_length", key=lambda x: x.astype(int)
).copy()

fig, ax = plt.subplots(figsize=(4, 3.5))
bars = ax.bar(
    sub["anchor_length"].astype(str),
    sub["exact_match_rate"],
    color="#1f77b4",
    width=0.5,
)
for bar, val in zip(bars, sub["exact_match_rate"]):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + 0.005,
        f"{val:.3f}",
        ha="center",
        va="bottom",
        fontsize=12,
    )
ax.set_xlabel("Boundary Anchor Length (bp)", fontsize=14)
ax.set_ylabel("Exact Match Rate", fontsize=14)
ax.set_ylim(0, 1.12)
ax.set_title("scRNA EIB Counting Accuracy (XS)", fontsize=15)
ax.tick_params(labelsize=12)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()
fig.savefig(out_by_anchor, dpi=300, bbox_inches="tight")
plt.close(fig)

# 5b. Accuracy by boundary type (if available in expected)
# Use the first anchor length's expected counts for type distribution
df_expected_first = pd.read_csv(_expected_lookup[anchor_lengths[0]], sep="\t")
if "boundary_type" in df_expected_first.columns:
    fig, ax = plt.subplots(figsize=(8, 5))
    for btype in df_expected_first["boundary_type"].unique():
        sub = df_expected_first[df_expected_first["boundary_type"] == btype]
        ax.bar(btype, len(sub), alpha=0.7)
    ax.set_ylabel("Count")
    ax.set_title("Boundary Type Distribution")
    fig.savefig(out_by_boundary_type, dpi=150, bbox_inches="tight")
    plt.close(fig)
else:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.text(0.5, 0.5, "No boundary type information", ha="center", va="center", transform=ax.transAxes)
    fig.savefig(out_by_boundary_type, dpi=150, bbox_inches="tight")
    plt.close(fig)

# 5c. Per-cell accuracy distribution (exact match rate only, anchor=1 bp)
fig, ax = plt.subplots(figsize=(4, 3.5))
# Use first dataset (anchor_length=1) for per-cell distribution
first_label = dataset_labels[0] if dataset_labels else None
if first_label:
    sub = df_per_cell[df_per_cell["dataset"] == first_label]
    if len(sub) > 0:
        ax.hist(sub["exact_match_rate"], bins=20, color="#1f77b4", edgecolor="black", linewidth=0.5)
ax.set_xlabel("Exact Match Rate", fontsize=14)
ax.set_ylabel("Number of Cells", fontsize=14)
ax.set_title("Per-cell Exact Match Rate\n(Anchor Length 1 bp)", fontsize=15)
ax.tick_params(labelsize=12)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()
fig.savefig(out_per_cell_plot, dpi=300, bbox_inches="tight")
plt.close(fig)

# 5d. Heatmap (anchor_length × strandedness)
if len(df_accuracy) > 1:
    pivot = df_accuracy.pivot_table(
        index="anchor_length",
        columns="strandedness",
        values="sensitivity",
        aggfunc="mean",
    )
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(pivot, annot=True, fmt=".3f", cmap="YlOrRd", ax=ax, vmin=0, vmax=1)
    ax.set_title("Sensitivity: Anchor Length × Strandedness (scRNA-seq EIB)")
    fig.tight_layout()
    fig.savefig(out_heatmap, dpi=150, bbox_inches="tight")
    plt.close(fig)
else:
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.text(0.5, 0.5, "Insufficient data for heatmap", ha="center", va="center", transform=ax.transAxes)
    fig.savefig(out_heatmap, dpi=150, bbox_inches="tight")
    plt.close(fig)
