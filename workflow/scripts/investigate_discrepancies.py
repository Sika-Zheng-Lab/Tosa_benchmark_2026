"""
Investigate junction count discrepancies between Tosa, STAR, regtools,
and featureCounts.

For every junction where any two tools disagree on count, this script
extracts the supporting reads from the BAM and classifies the likely
cause of the difference.

Cause categories
----------------
- multi_mapping    : reads with NH > 1 support the junction
- secondary_supp   : secondary (0x100) or supplementary (0x800) alignments
                     present at the junction
- anchor_length    : some reads have anchor < 8 bp (regtools -a 8 threshold)
- pe_mate_dedup    : both PE mates span the same junction (PE only)
- intron_size      : intron size < 70 bp or > 500 000 bp (regtools -m/-M)
- unknown          : none of the above explains the difference

Outputs
-------
- discrepancy_table.tsv  : per-junction diagnostic table
- discrepancy_summary.png: overview plots
"""

import os
import re
import pysam
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# ---------------------------------------------------------------------------
# Snakemake I/O
# ---------------------------------------------------------------------------
in_counts = snakemake.input.counts
in_bam = snakemake.input.bam
in_bai = snakemake.input.bai
in_sj = snakemake.input.sj

out_table = snakemake.output.table
out_plot = snakemake.output.plot

is_paired = snakemake.params.is_paired

os.makedirs(os.path.dirname(out_table), exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Load merged counts & identify discrepant junctions
# ---------------------------------------------------------------------------
df = pd.read_csv(in_counts, sep="\t")
# Columns: junction, STAR, Tosa, regtools, featureCounts

tools = ["Tosa", "STAR", "regtools", "featureCounts"]
# A junction is discrepant if any pair of tools differs
import itertools as _itertools_max
df["max_diff"] = df[tools].apply(
    lambda row: max(
        abs(row[t1] - row[t2])
        for t1, t2 in _itertools_max.combinations(tools, 2)
    ),
    axis=1,
)
disc = df[df["max_diff"] >= 1].copy()

if disc.empty:
    # No discrepancies — write empty outputs
    disc_out = pd.DataFrame(
        columns=[
            "junction", "Tosa", "STAR", "regtools", "featureCounts",
            "intron_size", "n_total_reads", "n_unique", "n_multi",
            "n_secondary_supp", "min_anchor", "n_anchor_lt8",
            "n_pe_both_mates", "star_multi_reads", "suspected_causes",
        ]
    )
    disc_out.to_csv(out_table, sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(0.5, 0.5, "No discrepancies found", ha="center", va="center",
            fontsize=14, transform=ax.transAxes)
    ax.set_axis_off()
    plt.savefig(out_plot, dpi=150)
    plt.close()
    import sys
    sys.exit(0)

# ---------------------------------------------------------------------------
# 2. Load STAR SJ.out.tab for multi_reads column
# ---------------------------------------------------------------------------
star_multi = {}
with open(in_sj) as fh:
    for line in fh:
        parts = line.strip().split("\t")
        chrom, start, end = parts[0], parts[1], parts[2]
        junction = f"{chrom}:{start}-{end}"
        multi_reads = int(parts[7])
        star_multi[junction] = star_multi.get(junction, 0) + multi_reads

# ---------------------------------------------------------------------------
# 3. Parse junction string → genomic coords
# ---------------------------------------------------------------------------
_junc_re = re.compile(r"^(.+):(\d+)-(\d+)$")


def parse_junction(junc_str):
    """Return (chrom, start_1based, end_1based) from 'chr:start-end'."""
    m = _junc_re.match(junc_str)
    return m.group(1), int(m.group(2)), int(m.group(3))


# ---------------------------------------------------------------------------
# 4. Extract reads from BAM for each discrepant junction
# ---------------------------------------------------------------------------
def get_junction_reads(bamfile, chrom, intron_start_1, intron_end_1):
    """
    Fetch reads whose CIGAR has an N-operation matching the given intron.

    Parameters
    ----------
    intron_start_1 : int   1-based first position of the intron
    intron_end_1   : int   1-based last position of the intron

    Returns a list of dicts with per-read information.
    """
    # pysam uses 0-based half-open coordinates
    intron_start_0 = intron_start_1 - 1  # 0-based start
    intron_end_0 = intron_end_1          # 0-based exclusive end

    results = []
    # Fetch reads that overlap the intron region
    try:
        fetch_start = max(0, intron_start_0 - 500)
        fetch_end = intron_end_0 + 500
        reads = bamfile.fetch(chrom, fetch_start, fetch_end)
    except ValueError:
        return results

    for read in reads:
        if read.is_unmapped:
            continue
        # Walk the CIGAR to find N operations
        ref_pos = read.reference_start  # 0-based
        query_consumed = 0
        for op, length in read.cigartuples:
            if op == 3:  # N = intron skip
                n_start = ref_pos          # 0-based intron start
                n_end = ref_pos + length   # 0-based exclusive intron end
                if n_start == intron_start_0 and n_end == intron_end_0:
                    # This read spans our junction
                    # Compute left and right anchor lengths
                    left_anchor = n_start - read.reference_start
                    right_anchor = read.reference_end - n_end
                    # Recompute left anchor = total M before this N
                    # We use the simple approach: position diff
                    nh = read.get_tag("NH") if read.has_tag("NH") else 1
                    is_secondary = read.is_secondary
                    is_supplementary = read.is_supplementary
                    results.append({
                        "read_name": read.query_name,
                        "nh": nh,
                        "is_secondary": is_secondary,
                        "is_supplementary": is_supplementary,
                        "left_anchor": left_anchor,
                        "right_anchor": right_anchor,
                        "is_read1": read.is_read1 if read.is_paired else None,
                    })
                    break  # found the matching N, no need to check further
            # Advance reference position for reference-consuming ops
            if op in (0, 2, 3, 7, 8):  # M, D, N, =, X
                ref_pos += length

    return results


# ---------------------------------------------------------------------------
# 5. Analyse each discrepant junction
# ---------------------------------------------------------------------------
bamfile = pysam.AlignmentFile(in_bam, "rb")

rows = []
for _, row in disc.iterrows():
    junc = row["junction"]
    chrom, start_1, end_1 = parse_junction(junc)
    intron_size = end_1 - start_1 + 1

    read_infos = get_junction_reads(bamfile, chrom, start_1, end_1)

    n_total = len(read_infos)
    n_unique = sum(1 for r in read_infos if r["nh"] == 1)
    n_multi = sum(1 for r in read_infos if r["nh"] > 1)
    n_sec_supp = sum(
        1 for r in read_infos if r["is_secondary"] or r["is_supplementary"]
    )

    anchors = [min(r["left_anchor"], r["right_anchor"]) for r in read_infos]
    min_anchor = min(anchors) if anchors else -1
    n_anchor_lt8 = sum(1 for a in anchors if a < 8)

    # PE mate dedup: count read names where both mates span the junction
    n_pe_both = 0
    if is_paired and read_infos:
        from collections import Counter
        name_counts = Counter(r["read_name"] for r in read_infos)
        n_pe_both = sum(1 for c in name_counts.values() if c >= 2)

    s_multi = star_multi.get(junc, 0)

    # Classify causes
    causes = []
    if n_multi > 0 and (row["regtools"] > row["Tosa"] or row["regtools"] > row["STAR"]):
        causes.append("multi_mapping")
    if n_sec_supp > 0:
        causes.append("secondary_supp")
    if n_anchor_lt8 > 0:
        causes.append("anchor_length")
    if n_pe_both > 0 and is_paired:
        causes.append("pe_mate_dedup")
    if intron_size < 70 or intron_size > 500000:
        causes.append("intron_size")
    if not causes:
        causes.append("unknown")

    rows.append({
        "junction": junc,
        "Tosa": int(row["Tosa"]),
        "STAR": int(row["STAR"]),
        "regtools": int(row["regtools"]),
        "featureCounts": int(row["featureCounts"]),
        "intron_size": intron_size,
        "n_total_reads": n_total,
        "n_unique": n_unique,
        "n_multi": n_multi,
        "n_secondary_supp": n_sec_supp,
        "min_anchor": min_anchor,
        "n_anchor_lt8": n_anchor_lt8,
        "n_pe_both_mates": n_pe_both,
        "star_multi_reads": s_multi,
        "suspected_causes": ",".join(causes),
    })

bamfile.close()

df_out = pd.DataFrame(rows)
df_out.to_csv(out_table, sep="\t", index=False)

print(f"Discrepant junctions: {len(df_out)}")
print(f"Cause distribution:")
all_causes = []
for c in df_out["suspected_causes"]:
    all_causes.extend(c.split(","))
cause_counts = pd.Series(all_causes).value_counts()
print(cause_counts.to_string())

# ---------------------------------------------------------------------------
# 6. Summary plots
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(14, 11))

# (a) Cause category counts
cause_list = []
for causes_str in df_out["suspected_causes"]:
    cause_list.extend(causes_str.split(","))
cause_series = pd.Series(cause_list).value_counts()
sns.barplot(x=cause_series.values, y=cause_series.index, ax=axes[0, 0],
            palette="Set2", orient="h")
axes[0, 0].set_xlabel("Number of junctions")
axes[0, 0].set_title("Discrepancy causes\n(non-exclusive)")
for i, v in enumerate(cause_series.values):
    axes[0, 0].text(v + 0.3, i, str(v), va="center")

# (b) Count difference histograms (pairwise)
import itertools as _itertools
pairs = list(_itertools.combinations(tools, 2))
colors = sns.color_palette("Set2", len(pairs))
for idx, (t1, t2) in enumerate(pairs):
    diff = df_out[t1] - df_out[t2]
    axes[0, 1].hist(diff, bins=50, alpha=0.5, label=f"{t1} − {t2}",
                     color=colors[idx])
axes[0, 1].set_xlabel("Count difference")
axes[0, 1].set_ylabel("Number of junctions")
axes[0, 1].set_title("Distribution of count differences")
axes[0, 1].legend(fontsize=6)
axes[0, 1].axvline(0, color="black", linewidth=0.8, linestyle="--")

# (c) Multi-mapping reads vs count diff (Tosa vs regtools)
diff_tosa_reg = df_out["regtools"] - df_out["Tosa"]
axes[1, 0].scatter(df_out["n_multi"], diff_tosa_reg, alpha=0.4, s=12,
                    edgecolors="none")
axes[1, 0].set_xlabel("Number of multi-mapping reads (NH > 1)")
axes[1, 0].set_ylabel("regtools − Tosa")
axes[1, 0].set_title("Multi-mapping vs count difference")
axes[1, 0].axhline(0, color="black", linewidth=0.8, linestyle="--")

# (d) Anchor length vs count diff (Tosa vs regtools)
axes[1, 1].scatter(df_out["n_anchor_lt8"], diff_tosa_reg, alpha=0.4, s=12,
                    edgecolors="none", color="orange")
axes[1, 1].set_xlabel("Number of reads with anchor < 8 bp")
axes[1, 1].set_ylabel("regtools − Tosa")
axes[1, 1].set_title("Short anchor reads vs count difference")
axes[1, 1].axhline(0, color="black", linewidth=0.8, linestyle="--")

plt.tight_layout()
plt.savefig(out_plot, dpi=150)
plt.close()
