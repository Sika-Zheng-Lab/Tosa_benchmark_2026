"""
Evaluate Tosa boundary counting accuracy against ground truth (EIB benchmark).

Compares Tosa boundary output (from the 1-thread run) with expected boundary
counts for each dataset (anchor_length x library_type x strandedness).
Produces stratified analysis by:
  - Exon count (2, 3, 5, 10)
  - Intron size (binned: small/medium/large)
  - Coverage level (5, 20, 100)
  - Boundary type (5p vs 3p)
  - Boundary anchor length
  - Library type x Strandedness (per anchor length)

Outputs:
  - accuracy_summary.tsv: overall and per-dataset metrics
  - strand_accuracy.tsv: strand assignment accuracy
  - accuracy_by_exon_count.png
  - accuracy_by_intron_size.png
  - accuracy_by_coverage.png
  - accuracy_by_boundary_type.png
  - accuracy_by_anchor_length.png
  - accuracy_heatmap.png (library_type x strandedness, faceted by anchor)
"""

import os
import re
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
in_tosa_files = snakemake.input.tosa_boundaries       # list of Tosa boundary output files
in_expected_files = snakemake.input.expected_counts     # list (anchor_length x library_type)
in_read_meta_files = snakemake.input.read_metadata      # list per library_type (PE, SE)
in_gene_metadata = snakemake.input.gene_metadata        # single gene metadata file

# Dataset labels like "b1_PE_XS", "b4_PE_RF", ..., "b8_SE_FR"
dataset_labels = snakemake.params.dataset_labels
anchor_lengths = snakemake.params.anchor_lengths        # list of anchor length strings
library_types = snakemake.params.library_types           # ["PE", "SE"]
strandedness_modes = snakemake.params.strandedness_modes # ["XS", "RF", "FR"]

out_accuracy_summary = snakemake.output.accuracy_summary
out_strand_accuracy = snakemake.output.strand_accuracy
out_by_exon = snakemake.output.by_exon_count
out_by_intron = snakemake.output.by_intron_size
out_by_coverage = snakemake.output.by_coverage
out_by_boundary_type = snakemake.output.by_boundary_type
out_by_anchor = snakemake.output.by_anchor_length
out_heatmap = snakemake.output.heatmap

os.makedirs(os.path.dirname(out_accuracy_summary), exist_ok=True)


# ---------------------------------------------------------------------------
# Helpers — parse dataset label
# ---------------------------------------------------------------------------

def parse_label(label):
    """Parse dataset label like 'b1_PE_XS' -> (anchor_length, library_type, strandedness)."""
    m = re.match(r"^b(\d+)_([^_]+)_(.+)$", label)
    if not m:
        raise ValueError(f"Cannot parse dataset label: {label}")
    return m.group(1), m.group(2), m.group(3)


# ---------------------------------------------------------------------------
# 1. Load gene metadata
# ---------------------------------------------------------------------------
df_gene_meta = pd.read_csv(in_gene_metadata, sep="\t")


def mean_intron_size(row):
    if row["intron_sizes"] == "NA" or pd.isna(row["intron_sizes"]):
        return 0
    sizes = [int(x) for x in str(row["intron_sizes"]).split(",")]
    return np.mean(sizes)


df_gene_meta["mean_intron_size"] = df_gene_meta.apply(mean_intron_size, axis=1)


def bin_intron_size(val):
    if val <= 1000:
        return "small (<=1kb)"
    elif val <= 10000:
        return "medium (1-10kb)"
    else:
        return "large (>10kb)"


df_gene_meta["intron_size_bin"] = df_gene_meta["mean_intron_size"].apply(bin_intron_size)

# ---------------------------------------------------------------------------
# 2. Build lookups for expected counts and read metadata
# ---------------------------------------------------------------------------
# expected_counts files are ordered by expand(anchor_lengths x library_types)
# => (a1, PE), (a1, SE), (a4, PE), (a4, SE), ...
_expected_lookup = {}  # (anchor_length_str, library_type) -> filepath
idx = 0
for al in anchor_lengths:
    for lt in library_types:
        _expected_lookup[(al, lt)] = in_expected_files[idx]
        idx += 1

# read_metadata files are per library_type only (anchor-independent)
_readmeta_lookup = {}
for idx_lt, lt in enumerate(library_types):
    _readmeta_lookup[lt] = in_read_meta_files[idx_lt]

# ---------------------------------------------------------------------------
# 3. Load and compare per dataset
# ---------------------------------------------------------------------------

all_results = []
strand_results = []

for i, label in enumerate(dataset_labels):
    al, library_type, strandedness = parse_label(label)

    # Load Tosa boundary output (columns: Boundary, Type, Strand, Count)
    df_tosa = pd.read_csv(in_tosa_files[i], sep="\t", compression="gzip")
    df_tosa.columns = ["boundary", "type_observed", "strand_observed", "count"]

    # Load expected boundary counts for this anchor_length + library_type
    df_exp = pd.read_csv(_expected_lookup[(al, library_type)], sep="\t")

    # Merge on boundary coordinate
    df_merged = df_exp.merge(
        df_tosa[["boundary", "type_observed", "strand_observed", "count"]],
        on="boundary",
        how="left",
    )
    df_merged["count"] = df_merged["count"].fillna(0).astype(int)
    df_merged["type_observed"] = df_merged["type_observed"].fillna(".")
    df_merged["strand_observed"] = df_merged["strand_observed"].fillna(".")

    # Add gene metadata
    df_merged = df_merged.merge(
        df_gene_meta[["gene_id", "n_exons", "mean_intron_size", "intron_size_bin", "coverage"]],
        on="gene_id",
        how="left",
    )

    # Accuracy metrics
    df_merged["exact_match"] = df_merged["expected_count"] == df_merged["count"]
    df_merged["diff"] = df_merged["count"] - df_merged["expected_count"]
    df_merged["dataset"] = label
    df_merged["anchor_length"] = int(al)
    df_merged["library_type"] = library_type
    df_merged["strandedness"] = strandedness

    all_results.append(df_merged)

    # Strand accuracy (compare gene strand from GTF vs Tosa's assigned strand)
    df_merged["strand_correct"] = df_merged["strand"] == df_merged["strand_observed"]
    strand_results.append(
        df_merged[
            ["boundary", "gene_id", "dataset", "anchor_length",
             "strand", "strand_observed", "strand_correct"]
        ]
    )

    # Check for false-positive boundaries (in Tosa but not in expected)
    expected_boundaries = set(df_exp["boundary"])
    tosa_boundaries = set(df_tosa["boundary"])
    false_positives = tosa_boundaries - expected_boundaries
    if false_positives:
        for fp in false_positives:
            fp_row = df_tosa[df_tosa["boundary"] == fp].iloc[0]
            all_results[-1] = pd.concat(
                [
                    all_results[-1],
                    pd.DataFrame(
                        [
                            {
                                "boundary": fp,
                                "type": fp_row["type_observed"],
                                "strand": ".",
                                "expected_count": 0,
                                "gene_id": "unknown",
                                "count": int(fp_row["count"]),
                                "type_observed": fp_row["type_observed"],
                                "strand_observed": fp_row["strand_observed"],
                                "n_exons": np.nan,
                                "mean_intron_size": np.nan,
                                "intron_size_bin": "unknown",
                                "coverage": np.nan,
                                "exact_match": False,
                                "diff": int(fp_row["count"]),
                                "dataset": label,
                                "anchor_length": int(al),
                                "library_type": library_type,
                                "strandedness": strandedness,
                                "strand_correct": False,
                            }
                        ]
                    ),
                ],
                ignore_index=True,
            )

df_all = pd.concat(all_results, ignore_index=True)
df_strand_all = pd.concat(strand_results, ignore_index=True)

# ---------------------------------------------------------------------------
# 4. Summary table
# ---------------------------------------------------------------------------
summary_rows = []
for label in dataset_labels:
    al, lt, sd = parse_label(label)
    sub = df_all[df_all["dataset"] == label]
    n_total = len(sub)
    n_exact = sub["exact_match"].sum()
    n_detected = (sub["count"] > 0).sum()
    n_expected = (sub["expected_count"] > 0).sum()
    sensitivity = n_detected / n_expected if n_expected > 0 else 0

    valid = sub[sub["expected_count"] > 0]
    if len(valid) > 1:
        pearson_r, _ = stats.pearsonr(valid["expected_count"], valid["count"])
        spearman_r, _ = stats.spearmanr(valid["expected_count"], valid["count"])
        rmse = np.sqrt(np.mean(valid["diff"] ** 2))
    else:
        pearson_r = spearman_r = rmse = np.nan

    summary_rows.append(
        {
            "dataset": label,
            "anchor_length": int(al),
            "library_type": lt,
            "strandedness": sd,
            "n_expected_boundaries": n_expected,
            "n_detected": n_detected,
            "sensitivity": round(sensitivity, 4),
            "n_exact_match": int(n_exact),
            "exact_match_rate": round(n_exact / n_total if n_total > 0 else 0, 4),
            "rmse": round(rmse, 4) if not np.isnan(rmse) else "NA",
            "pearson_r": round(pearson_r, 4) if not np.isnan(pearson_r) else "NA",
            "spearman_r": round(spearman_r, 4) if not np.isnan(spearman_r) else "NA",
            "n_false_positives": int((sub["expected_count"] == 0).sum()),
        }
    )

df_summary = pd.DataFrame(summary_rows)
df_summary.to_csv(out_accuracy_summary, sep="\t", index=False)

# ---------------------------------------------------------------------------
# 5. Strand accuracy table
# ---------------------------------------------------------------------------
strand_summary_rows = []
for label in dataset_labels:
    al, lt, sd = parse_label(label)
    sub = df_strand_all[df_strand_all["dataset"] == label]
    n_total = len(sub)
    n_correct = sub["strand_correct"].sum()
    strand_summary_rows.append(
        {
            "dataset": label,
            "anchor_length": int(al),
            "library_type": lt,
            "strandedness": sd,
            "n_boundaries": n_total,
            "n_strand_correct": int(n_correct),
            "strand_accuracy": round(n_correct / n_total if n_total > 0 else 0, 4),
        }
    )

df_strand_summary = pd.DataFrame(strand_summary_rows)
df_strand_summary.to_csv(out_strand_accuracy, sep="\t", index=False)

# ---------------------------------------------------------------------------
# 6. Stratified plots
# ---------------------------------------------------------------------------
plt.rcParams.update({"font.size": 11})

# 6a. By exon count
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
valid = df_all[df_all["n_exons"].notna()].copy()
valid["n_exons"] = valid["n_exons"].astype(int)

grouped = valid.groupby("n_exons").agg(
    exact_match_rate=("exact_match", "mean"),
    rmse=("diff", lambda x: np.sqrt(np.mean(x**2))),
).reset_index()

axes[0].bar(grouped["n_exons"].astype(str), grouped["exact_match_rate"], color="#4C72B0")
axes[0].set_xlabel("Number of Exons")
axes[0].set_ylabel("Exact Match Rate")
axes[0].set_title("Accuracy by Exon Count")
axes[0].set_ylim(0, 1.05)

axes[1].bar(grouped["n_exons"].astype(str), grouped["rmse"], color="#DD8452")
axes[1].set_xlabel("Number of Exons")
axes[1].set_ylabel("RMSE")
axes[1].set_title("RMSE by Exon Count")

plt.tight_layout()
plt.savefig(out_by_exon, dpi=150)
plt.close()

# 6b. By intron size bin
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
valid = df_all[df_all["intron_size_bin"] != "unknown"].copy()
bin_order = ["small (<=1kb)", "medium (1-10kb)", "large (>10kb)"]

grouped = valid.groupby("intron_size_bin").agg(
    exact_match_rate=("exact_match", "mean"),
    rmse=("diff", lambda x: np.sqrt(np.mean(x**2))),
).reindex(bin_order).reset_index()

axes[0].bar(grouped["intron_size_bin"], grouped["exact_match_rate"], color="#4C72B0")
axes[0].set_xlabel("Intron Size")
axes[0].set_ylabel("Exact Match Rate")
axes[0].set_title("Accuracy by Intron Size")
axes[0].set_ylim(0, 1.05)

axes[1].bar(grouped["intron_size_bin"], grouped["rmse"], color="#DD8452")
axes[1].set_xlabel("Intron Size")
axes[1].set_ylabel("RMSE")
axes[1].set_title("RMSE by Intron Size")

plt.tight_layout()
plt.savefig(out_by_intron, dpi=150)
plt.close()

# 6c. By coverage level
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
valid = df_all[df_all["coverage"].notna()].copy()
valid["coverage"] = valid["coverage"].astype(int)

grouped = valid.groupby("coverage").agg(
    exact_match_rate=("exact_match", "mean"),
    rmse=("diff", lambda x: np.sqrt(np.mean(x**2))),
).reset_index()

axes[0].bar(grouped["coverage"].astype(str), grouped["exact_match_rate"], color="#4C72B0")
axes[0].set_xlabel("Coverage Level")
axes[0].set_ylabel("Exact Match Rate")
axes[0].set_title("Accuracy by Coverage Level")
axes[0].set_ylim(0, 1.05)

axes[1].bar(grouped["coverage"].astype(str), grouped["rmse"], color="#DD8452")
axes[1].set_xlabel("Coverage Level")
axes[1].set_ylabel("RMSE")
axes[1].set_title("RMSE by Coverage Level")

plt.tight_layout()
plt.savefig(out_by_coverage, dpi=150)
plt.close()

# 6d. By boundary type (5p vs 3p)
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
valid = df_all[df_all["type"].isin(["5p", "3p"])].copy()

grouped = valid.groupby("type").agg(
    exact_match_rate=("exact_match", "mean"),
    rmse=("diff", lambda x: np.sqrt(np.mean(x**2))),
).reindex(["5p", "3p"]).reset_index()

axes[0].bar(grouped["type"], grouped["exact_match_rate"], color="#4C72B0")
axes[0].set_xlabel("Boundary Type")
axes[0].set_ylabel("Exact Match Rate")
axes[0].set_title("Accuracy by Boundary Type")
axes[0].set_ylim(0, 1.05)

axes[1].bar(grouped["type"], grouped["rmse"], color="#DD8452")
axes[1].set_xlabel("Boundary Type")
axes[1].set_ylabel("RMSE")
axes[1].set_title("RMSE by Boundary Type")

plt.tight_layout()
plt.savefig(out_by_boundary_type, dpi=150)
plt.close()

# 6e. By anchor length
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

grouped = df_all.groupby("anchor_length").agg(
    exact_match_rate=("exact_match", "mean"),
    rmse=("diff", lambda x: np.sqrt(np.mean(x**2))),
).reset_index()

axes[0].bar(grouped["anchor_length"].astype(str), grouped["exact_match_rate"], color="#4C72B0")
axes[0].set_xlabel("Boundary Anchor Length")
axes[0].set_ylabel("Exact Match Rate")
axes[0].set_title("Accuracy by Anchor Length")
axes[0].set_ylim(0, 1.05)

axes[1].bar(grouped["anchor_length"].astype(str), grouped["rmse"], color="#DD8452")
axes[1].set_xlabel("Boundary Anchor Length")
axes[1].set_ylabel("RMSE")
axes[1].set_title("RMSE by Anchor Length")

plt.tight_layout()
plt.savefig(out_by_anchor, dpi=150)
plt.close()

# 6f. Bar chart: exact match rate by anchor length (PE + XS only)
sub = df_summary[
    (df_summary["library_type"] == "PE") & (df_summary["strandedness"] == "XS")
].sort_values("anchor_length").copy()

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
        bar.get_height() + 0.002,
        f"{val:.3f}",
        ha="center",
        va="bottom",
        fontsize=12,
    )
ax.set_xlabel("Boundary Anchor Length (bp)", fontsize=14)
ax.set_ylabel("Exact Match Rate", fontsize=14)
ax.set_ylim(0.98, 1.005)
ax.set_title("EIB Counting Accuracy (PE, XS)", fontsize=15)
ax.tick_params(labelsize=12)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig(out_heatmap, dpi=300)
plt.close()

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
print("=== EIB Boundary Accuracy Evaluation ===")
print(df_summary.to_string(index=False))
print("\n=== Strand Accuracy ===")
print(df_strand_summary.to_string(index=False))
