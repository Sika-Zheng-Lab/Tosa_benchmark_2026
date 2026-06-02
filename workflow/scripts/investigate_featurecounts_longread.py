"""Investigate why featureCounts shows extremely low Jaccard index
compared to Tosa and regtools in the long-read RNA-seq benchmark.

Analyses performed
------------------
1. Basic statistics: junction counts per tool per sample
2. Overlap analysis: shared / unique / missing junctions
3. Coordinate offset check: do featureCounts-unique junctions become
   matches when shifted by ±1 or ±2?
4. Annotation status: annotated vs novel fraction per tool
5. Intron size distributions by overlap category
6. Raw .jcounts sanity checks (Chr_SP1 ≠ Chr_SP2, low-count junctions)
7. Diagnostic plots and summary report

Inputs (Snakemake)
------------------
- norm_tosa       : list of normalized tosa.tsv (one per sample)
- norm_regtools   : list of normalized regtools.tsv
- norm_featurecounts : list of normalized featurecounts.tsv
- raw_jcounts     : one raw .jcounts file (first sample) for sanity checks
- gtf             : GTF annotation file

Outputs
-------
- report          : TSV diagnostic report
- plot_counts     : junction count bar chart
- plot_offset     : coordinate offset histogram
- plot_annotation : annotated vs novel stacked bar
- plot_intron_size: intron size violin plot
- plot_overlap    : overlap bar chart per sample
"""

import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
import matplotlib.pyplot as plt
import seaborn as sns

# ---------------------------------------------------------------------------
# Snakemake I/O
# ---------------------------------------------------------------------------
in_norm_tosa = snakemake.input.norm_tosa
in_norm_regtools = snakemake.input.norm_regtools
in_norm_featurecounts = snakemake.input.norm_featurecounts
in_raw_jcounts = snakemake.input.raw_jcounts
in_gtf = snakemake.input.gtf

out_report = snakemake.output.report
out_plot_counts = snakemake.output.plot_counts
out_plot_offset = snakemake.output.plot_offset
out_plot_annotation = snakemake.output.plot_annotation
out_plot_intron_size = snakemake.output.plot_intron_size
out_plot_overlap = snakemake.output.plot_overlap

samples = snakemake.params.samples

for out_path in [out_report, out_plot_counts, out_plot_offset,
                 out_plot_annotation, out_plot_intron_size, out_plot_overlap]:
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

TOOLS = ["Tosa", "regtools", "featureCounts"]

# ===================================================================
# Phase 1: Load data & basic statistics
# ===================================================================

print("Phase 1: Loading normalized junctions ...")

tool_inputs = {
    "Tosa": in_norm_tosa,
    "regtools": in_norm_regtools,
    "featureCounts": in_norm_featurecounts,
}

# Per-sample junction sets and counts
per_sample_data = {}  # sample → {tool: DataFrame}
per_sample_counts = []  # rows for summary

for i, sample in enumerate(samples):
    per_sample_data[sample] = {}
    for tool, paths in tool_inputs.items():
        df = pd.read_csv(paths[i], sep="\t")
        per_sample_data[sample][tool] = df
        per_sample_counts.append({
            "sample": sample,
            "tool": tool,
            "n_junctions": len(df),
            "total_count": int(df["count"].sum()),
        })

df_counts = pd.DataFrame(per_sample_counts)

# ===================================================================
# Phase 2: Overlap analysis (per sample)
# ===================================================================

print("Phase 2: Overlap analysis ...")

overlap_rows = []
all_fc_unique = set()
all_fc_missing = set()
all_tosa_regtools_union = set()

for sample in samples:
    s_tosa = set(per_sample_data[sample]["Tosa"]["junction"])
    s_regtools = set(per_sample_data[sample]["regtools"]["junction"])
    s_fc = set(per_sample_data[sample]["featureCounts"]["junction"])

    shared_all = s_tosa & s_regtools & s_fc
    tosa_regtools_only = (s_tosa & s_regtools) - s_fc
    fc_unique = s_fc - s_tosa - s_regtools
    fc_missing = (s_tosa & s_regtools) - s_fc

    all_fc_unique.update(fc_unique)
    all_fc_missing.update(fc_missing)
    all_tosa_regtools_union.update(s_tosa | s_regtools)

    overlap_rows.append({
        "sample": sample,
        "shared_all_3": len(shared_all),
        "tosa_regtools_only": len(tosa_regtools_only),
        "fc_unique": len(fc_unique),
        "fc_missing": len(fc_missing),
        "tosa_total": len(s_tosa),
        "regtools_total": len(s_regtools),
        "fc_total": len(s_fc),
    })

df_overlap = pd.DataFrame(overlap_rows)

# ===================================================================
# Phase 3: Coordinate offset check
# ===================================================================

print("Phase 3: Coordinate offset check ...")

_junc_re = re.compile(r"^(.+):(\d+)-(\d+)$")


def parse_junction(junc_str):
    m = _junc_re.match(junc_str)
    if not m:
        return None, None, None
    return m.group(1), int(m.group(2)), int(m.group(3))


def shift_junction(junc_str, d_start, d_end):
    chrom, start, end = parse_junction(junc_str)
    if chrom is None:
        return None
    return f"{chrom}:{start + d_start}-{end + d_end}"


# Build reference set from Tosa + regtools (union across all samples)
reference_junctions = set()
for sample in samples:
    reference_junctions.update(
        per_sample_data[sample]["Tosa"]["junction"]
    )
    reference_junctions.update(
        per_sample_data[sample]["regtools"]["junction"]
    )

# Check each offset combination for featureCounts-unique junctions
offsets = [
    (0, 0), (-1, 0), (1, 0), (0, -1), (0, 1),
    (-1, -1), (1, 1), (-1, 1), (1, -1),
    (-2, 0), (2, 0), (0, -2), (0, 2),
    (-2, -2), (2, 2),
]

offset_match_counts = {o: 0 for o in offsets}
for junc in all_fc_unique:
    for d_s, d_e in offsets:
        shifted = shift_junction(junc, d_s, d_e)
        if shifted and shifted in reference_junctions:
            offset_match_counts[(d_s, d_e)] += 1

offset_rows = []
for (d_s, d_e), count in sorted(offset_match_counts.items(),
                                  key=lambda x: -x[1]):
    offset_rows.append({
        "start_offset": d_s,
        "end_offset": d_e,
        "n_fc_unique_rescued": count,
        "pct_fc_unique_rescued": (
            100.0 * count / len(all_fc_unique) if all_fc_unique else 0
        ),
    })

df_offsets = pd.DataFrame(offset_rows)

# ===================================================================
# Phase 4: Annotation status analysis
# ===================================================================

print("Phase 4: Parsing GTF for known introns ...")

transcript_exons = defaultdict(list)
_re_tid = re.compile(r'transcript_id "([^"]+)"')

with open(in_gtf) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        if parts[2] != "exon":
            continue
        m = _re_tid.search(parts[8])
        if m:
            tid = m.group(1)
            transcript_exons[tid].append(
                (parts[0], int(parts[3]), int(parts[4]))
            )

known_introns = set()
for tid, exons in transcript_exons.items():
    exons_sorted = sorted(exons, key=lambda x: x[1])
    for i in range(len(exons_sorted) - 1):
        chrom = exons_sorted[i][0]
        intron_start = exons_sorted[i][2] + 1
        intron_end = exons_sorted[i + 1][1] - 1
        if intron_start <= intron_end:
            known_introns.add(f"{chrom}:{intron_start}-{intron_end}")

del transcript_exons

print(f"  Known introns from GTF: {len(known_introns):,}")

# Classify per tool (union across samples)
tool_all_junctions = {tool: set() for tool in TOOLS}
for sample in samples:
    for tool in TOOLS:
        tool_all_junctions[tool].update(
            per_sample_data[sample][tool]["junction"]
        )

annotation_rows = []
for tool in TOOLS:
    juncs = tool_all_junctions[tool]
    n_annotated = sum(1 for j in juncs if j in known_introns)
    n_novel = len(juncs) - n_annotated
    annotation_rows.append({
        "tool": tool,
        "n_annotated": n_annotated,
        "n_novel": n_novel,
        "pct_annotated": 100.0 * n_annotated / len(juncs) if juncs else 0,
        "total": len(juncs),
    })

# Also classify featureCounts-unique and featureCounts-missing
n_fc_unique_ann = sum(1 for j in all_fc_unique if j in known_introns)
n_fc_missing_ann = sum(1 for j in all_fc_missing if j in known_introns)
annotation_rows.append({
    "tool": "fc_unique",
    "n_annotated": n_fc_unique_ann,
    "n_novel": len(all_fc_unique) - n_fc_unique_ann,
    "pct_annotated": (
        100.0 * n_fc_unique_ann / len(all_fc_unique) if all_fc_unique else 0
    ),
    "total": len(all_fc_unique),
})
annotation_rows.append({
    "tool": "fc_missing",
    "n_annotated": n_fc_missing_ann,
    "n_novel": len(all_fc_missing) - n_fc_missing_ann,
    "pct_annotated": (
        100.0 * n_fc_missing_ann / len(all_fc_missing) if all_fc_missing else 0
    ),
    "total": len(all_fc_missing),
})

df_annotation = pd.DataFrame(annotation_rows)

# ===================================================================
# Phase 5: Intron size analysis
# ===================================================================

print("Phase 5: Intron size analysis ...")


def get_intron_size(junc_str):
    _, start, end = parse_junction(junc_str)
    if start is None:
        return None
    return end - start + 1


# Categorise junctions
# Union of first sample for convenience
sample0 = samples[0]
s_tosa_0 = set(per_sample_data[sample0]["Tosa"]["junction"])
s_regtools_0 = set(per_sample_data[sample0]["regtools"]["junction"])
s_fc_0 = set(per_sample_data[sample0]["featureCounts"]["junction"])

categories = {
    "shared_all_3": s_tosa_0 & s_regtools_0 & s_fc_0,
    "fc_unique": s_fc_0 - s_tosa_0 - s_regtools_0,
    "fc_missing": (s_tosa_0 & s_regtools_0) - s_fc_0,
}

intron_size_rows = []
for cat, juncs in categories.items():
    for j in juncs:
        sz = get_intron_size(j)
        if sz is not None:
            intron_size_rows.append({"category": cat, "intron_size": sz})

df_intron_size = pd.DataFrame(intron_size_rows)

# ===================================================================
# Phase 6: Raw .jcounts sanity checks
# ===================================================================

print("Phase 6: Raw .jcounts sanity checks ...")

df_raw = pd.read_csv(in_raw_jcounts, sep="\t", comment="#", low_memory=False)
count_col_raw = df_raw.columns[-1]

# Check Chr_SP1 ≠ Chr_SP2
n_diff_chr = int((df_raw["Chr_SP1"] != df_raw["Chr_SP2"]).sum())

# Check Strand_SP1 ≠ Strand_SP2
n_diff_strand = int((
    df_raw["Strand_SP1"].astype(str) != df_raw["Strand_SP2"].astype(str)
).sum())

# Count distribution of featureCounts-unique junctions in raw data
raw_junctions = set()
for _, row in df_raw.iterrows():
    chrom = row["Chr_SP1"]
    site1 = int(row["Location_SP1"])
    site2 = int(row["Location_SP2"])
    count = int(row[count_col_raw])
    if count == 0:
        continue
    intron_start = site1 + 1
    intron_end = site2 - 1
    raw_junctions.add(f"{chrom}:{intron_start}-{intron_end}")

# How many fc_unique junctions from first sample exist in raw data?
fc_unique_in_raw = all_fc_unique & raw_junctions
fc_unique_not_in_raw = all_fc_unique - raw_junctions

# Low-count junctions in raw file
raw_count_vals = df_raw[count_col_raw].astype(int)
n_count_1 = int((raw_count_vals == 1).sum())
n_count_le5 = int((raw_count_vals <= 5).sum())
n_total_raw = len(df_raw)

# Example discrepant junctions (first sample)
# Top 20 featureCounts-unique junctions sorted by count
df_fc_0 = per_sample_data[sample0]["featureCounts"]
fc_unique_list = df_fc_0[df_fc_0["junction"].isin(s_fc_0 - s_tosa_0 - s_regtools_0)]
fc_unique_top = fc_unique_list.nlargest(20, "count") if len(fc_unique_list) > 0 else fc_unique_list

# Top 20 featureCounts-missing junctions (shared by Tosa+regtools)
df_tosa_0 = per_sample_data[sample0]["Tosa"]
fc_missing_set = (s_tosa_0 & s_regtools_0) - s_fc_0
fc_missing_df = df_tosa_0[df_tosa_0["junction"].isin(fc_missing_set)]
fc_missing_top = fc_missing_df.nlargest(20, "count") if len(fc_missing_df) > 0 else fc_missing_df

# ===================================================================
# Phase 7: Write report
# ===================================================================

print("Phase 7: Writing diagnostic report ...")

with open(out_report, "w") as fh:
    fh.write("=" * 80 + "\n")
    fh.write("featureCounts Long-read Junction Diagnostic Report\n")
    fh.write("=" * 80 + "\n\n")

    # --- Section 1: Basic counts ---
    fh.write("## 1. Junction counts per tool per sample\n\n")
    fh.write(df_counts.to_string(index=False))
    fh.write("\n\n")

    # --- Section 2: Overlap ---
    fh.write("## 2. Junction overlap analysis\n\n")
    fh.write(df_overlap.to_string(index=False))
    fh.write("\n\n")

    # --- Section 3: Coordinate offset ---
    fh.write("## 3. Coordinate offset check\n")
    fh.write(f"Total featureCounts-unique junctions (union): {len(all_fc_unique)}\n")
    fh.write(f"Total featureCounts-missing junctions (union): {len(all_fc_missing)}\n\n")
    fh.write("Offset → rescued matches (featureCounts-unique → Tosa/regtools):\n")
    fh.write(df_offsets.to_string(index=False))
    fh.write("\n\n")

    # --- Section 4: Annotation status ---
    fh.write("## 4. Annotation status (annotated vs novel)\n\n")
    fh.write(f"Known introns from GTF: {len(known_introns):,}\n\n")
    fh.write(df_annotation.to_string(index=False))
    fh.write("\n\n")

    # --- Section 5: Intron size ---
    fh.write("## 5. Intron size statistics (first sample)\n\n")
    if not df_intron_size.empty:
        for cat in ["shared_all_3", "fc_unique", "fc_missing"]:
            subset = df_intron_size[df_intron_size["category"] == cat]["intron_size"]
            if len(subset) > 0:
                fh.write(f"  {cat}: n={len(subset):,}, "
                         f"median={subset.median():.0f}, "
                         f"mean={subset.mean():.0f}, "
                         f"min={subset.min():.0f}, max={subset.max():.0f}\n")
    fh.write("\n")

    # --- Section 6: Raw .jcounts checks ---
    fh.write("## 6. Raw .jcounts sanity checks\n\n")
    fh.write(f"Total rows in .jcounts: {n_total_raw:,}\n")
    fh.write(f"Rows with Chr_SP1 != Chr_SP2: {n_diff_chr:,}\n")
    fh.write(f"Rows with Strand_SP1 != Strand_SP2: {n_diff_strand:,}\n")
    fh.write(f"Rows with count == 1: {n_count_1:,}\n")
    fh.write(f"Rows with count <= 5: {n_count_le5:,}\n\n")

    # --- Section 7: Example junctions ---
    fh.write("## 7. Top featureCounts-unique junctions (first sample, by count)\n\n")
    if len(fc_unique_top) > 0:
        fh.write(fc_unique_top.to_string(index=False))
    else:
        fh.write("(none)")
    fh.write("\n\n")

    fh.write("## 8. Top featureCounts-missing junctions "
             "(Tosa+regtools shared, first sample, by Tosa count)\n\n")
    if len(fc_missing_top) > 0:
        fh.write(fc_missing_top.to_string(index=False))
    else:
        fh.write("(none)")
    fh.write("\n\n")

    fh.write("=" * 80 + "\n")
    fh.write("End of report\n")

# ===================================================================
# Phase 8: Diagnostic plots
# ===================================================================

print("Phase 8: Generating diagnostic plots ...")

# --- Plot 1: Junction count bar chart ---
fig, ax = plt.subplots(figsize=(10, 5))
sns.barplot(data=df_counts, x="sample", y="n_junctions", hue="tool",
            palette="Set2", ax=ax)
ax.set_title("Junction count per tool per sample (long-read)")
ax.set_ylabel("Number of junctions")
ax.set_xlabel("")
ax.tick_params(axis="x", rotation=45)
ax.legend(title="Tool")
plt.tight_layout()
plt.savefig(out_plot_counts, dpi=150)
plt.close()

# --- Plot 2: Coordinate offset histogram ---
fig, ax = plt.subplots(figsize=(10, 5))
labels = [f"({int(r['start_offset']):+d},{int(r['end_offset']):+d})"
          for _, r in df_offsets.iterrows()]
values = df_offsets["n_fc_unique_rescued"].tolist()
colors = ["#e74c3c" if l == "(+0,+0)" else "#3498db" for l in labels]
bars = ax.bar(labels, values, color=colors, edgecolor="black", linewidth=0.5)
ax.set_title(
    f"Coordinate offset rescue: featureCounts-unique junctions\n"
    f"(n={len(all_fc_unique):,} unique junctions across all samples)"
)
ax.set_ylabel("Number rescued")
ax.set_xlabel("(start_offset, end_offset)")
ax.tick_params(axis="x", rotation=45)
# Add value labels
for bar, val in zip(bars, values):
    if val > 0:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                str(val), ha="center", va="bottom", fontsize=8)
plt.tight_layout()
plt.savefig(out_plot_offset, dpi=150)
plt.close()

# --- Plot 3: Annotated vs novel stacked bar ---
df_ann_plot = df_annotation[df_annotation["tool"].isin(
    TOOLS + ["fc_unique", "fc_missing"]
)].copy()
fig, ax = plt.subplots(figsize=(8, 5))
x = np.arange(len(df_ann_plot))
w = 0.5
ax.bar(x, df_ann_plot["n_annotated"], w, label="Annotated", color="#2ecc71")
ax.bar(x, df_ann_plot["n_novel"], w, bottom=df_ann_plot["n_annotated"],
       label="Novel", color="#e74c3c")
ax.set_xticks(x)
ax.set_xticklabels(df_ann_plot["tool"], rotation=30, ha="right")
ax.set_ylabel("Number of junctions")
ax.set_title("Annotation status: annotated vs novel junctions")
ax.legend()
# Add percentage annotation
for i, (_, row) in enumerate(df_ann_plot.iterrows()):
    total = row["n_annotated"] + row["n_novel"]
    if total > 0:
        ax.text(i, total + total * 0.01,
                f'{row["pct_annotated"]:.1f}%\nann.',
                ha="center", va="bottom", fontsize=8)
plt.tight_layout()
plt.savefig(out_plot_annotation, dpi=150)
plt.close()

# --- Plot 4: Intron size violin plot ---
if not df_intron_size.empty:
    # Cap at 500kb for better visualization
    df_is_plot = df_intron_size.copy()
    df_is_plot["intron_size_log10"] = np.log10(
        df_is_plot["intron_size"].clip(lower=1)
    )
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.violinplot(
        data=df_is_plot, x="category", y="intron_size_log10",
        palette="Set2", ax=ax, inner="quartile",
        order=["shared_all_3", "fc_unique", "fc_missing"],
    )
    ax.set_title("Intron size distribution by overlap category (first sample)")
    ax.set_ylabel("log10(intron size)")
    ax.set_xlabel("")
    plt.tight_layout()
    plt.savefig(out_plot_intron_size, dpi=150)
    plt.close()
else:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(0.5, 0.5, "No data", ha="center", va="center",
            fontsize=14, transform=ax.transAxes)
    ax.set_axis_off()
    plt.savefig(out_plot_intron_size, dpi=150)
    plt.close()

# --- Plot 5: Overlap bar chart per sample ---
fig, ax = plt.subplots(figsize=(12, 5))
x = np.arange(len(samples))
width = 0.2
ax.bar(x - width, df_overlap["shared_all_3"], width,
       label="Shared all 3", color="#2ecc71")
ax.bar(x, df_overlap["tosa_regtools_only"], width,
       label="Tosa+regtools only", color="#3498db")
ax.bar(x + width, df_overlap["fc_unique"], width,
       label="featureCounts unique", color="#e74c3c")
ax.set_xticks(x)
ax.set_xticklabels(samples, rotation=45, ha="right")
ax.set_ylabel("Number of junctions")
ax.set_title("Junction overlap categories per sample")
ax.legend()
plt.tight_layout()
plt.savefig(out_plot_overlap, dpi=150)
plt.close()

print("Done.")
