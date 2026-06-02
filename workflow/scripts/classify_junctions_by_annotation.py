"""Classify junctions by GTF annotation status.

For each junction detected by any of the four tools (STAR, Tosa, regtools,
featureCounts), classify it into one of three categories:

  - annotated:    matches a known intron in the GTF transcript models
  - genic_novel:  at least one end falls within a gene region but the
                  junction is not a known intron
  - intergenic:   both ends lie outside any annotated gene region

Strand information is ignored for classification — only coordinates are used.

Outputs:
  classification.tsv — per-junction: junction, category, per-tool detection flags
  summary.tsv        — per-tool counts and fractions for each category
"""

import os
import re
from bisect import bisect_right
from collections import defaultdict

import pandas as pd

# ---------------------------------------------------------------------------
# Snakemake I/O
# ---------------------------------------------------------------------------
in_star = snakemake.input.star
in_tosa = snakemake.input.tosa
in_regtools = snakemake.input.regtools
in_featurecounts = snakemake.input.featurecounts
in_gtf = snakemake.input.gtf

out_classification = snakemake.output.classification
out_summary = snakemake.output.summary

for out_path in [out_classification, out_summary]:
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

TOOLS = ["STAR", "Tosa", "regtools", "featureCounts"]

# ---------------------------------------------------------------------------
# 1. Parse GTF — extract known introns and gene regions
# ---------------------------------------------------------------------------

# transcript_id → list of (chrom, exon_start, exon_end)  (1-based inclusive)
transcript_exons = defaultdict(list)
# gene regions per chromosome: list of (start, end)  (1-based inclusive)
gene_regions_raw = defaultdict(list)

_re_transcript_id = re.compile(r'transcript_id "([^"]+)"')

with open(in_gtf) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        feature = parts[2]
        chrom = parts[0]
        start = int(parts[3])  # 1-based inclusive
        end = int(parts[4])    # 1-based inclusive

        if feature == "gene":
            gene_regions_raw[chrom].append((start, end))
        elif feature == "exon":
            m = _re_transcript_id.search(parts[8])
            if m:
                tid = m.group(1)
                transcript_exons[tid].append((chrom, start, end))

# --- Build known intron set ---
# For each transcript, sort exons by start and derive introns from
# consecutive exon pairs: intron_start = exon_i_end + 1,
#                          intron_end   = exon_(i+1)_start - 1
known_introns = set()
for tid, exons in transcript_exons.items():
    exons_sorted = sorted(exons, key=lambda x: x[1])
    for i in range(len(exons_sorted) - 1):
        chrom = exons_sorted[i][0]
        intron_start = exons_sorted[i][2] + 1
        intron_end = exons_sorted[i + 1][1] - 1
        if intron_start <= intron_end:
            known_introns.add(f"{chrom}:{intron_start}-{intron_end}")

del transcript_exons  # free memory

# --- Build sorted gene-region index for bisect-based overlap queries ---
# For each chromosome, merge overlapping gene regions, then store sorted
# start and end arrays for binary search.
gene_starts = {}  # chrom → sorted list of merged-region starts
gene_ends = {}    # chrom → sorted list of merged-region ends

for chrom, regions in gene_regions_raw.items():
    regions.sort()
    merged = []
    for s, e in regions:
        if merged and s <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    gene_starts[chrom] = [m[0] for m in merged]
    gene_ends[chrom] = [m[1] for m in merged]

del gene_regions_raw


def is_in_gene_region(chrom, pos):
    """Check whether a 1-based position falls within any merged gene region."""
    if chrom not in gene_starts:
        return False
    starts = gene_starts[chrom]
    ends = gene_ends[chrom]
    # Find the rightmost region whose start <= pos
    idx = bisect_right(starts, pos) - 1
    if idx < 0:
        return False
    return pos <= ends[idx]


# ---------------------------------------------------------------------------
# 2. Load normalised junction data from the four tools
# ---------------------------------------------------------------------------
dfs = {}
for tool, path in zip(TOOLS, [in_star, in_tosa, in_regtools, in_featurecounts]):
    dfs[tool] = pd.read_csv(path, sep="\t")

# Outer merge on junction column
merged = dfs["STAR"][["junction"]].assign(STAR=1)
for tool in TOOLS[1:]:
    merged = merged.merge(
        dfs[tool][["junction"]].assign(**{tool: 1}),
        on="junction",
        how="outer",
    )
merged = merged.fillna(0).astype({t: int for t in TOOLS})

# ---------------------------------------------------------------------------
# 3. Classify each junction
# ---------------------------------------------------------------------------

def classify_junction(junc_str):
    """Return 'annotated', 'genic_novel', or 'intergenic'."""
    if junc_str in known_introns:
        return "annotated"
    # Parse coordinates
    chrom, coords = junc_str.rsplit(":", 1)
    start_str, end_str = coords.split("-")
    start = int(start_str)
    end = int(end_str)
    if is_in_gene_region(chrom, start) or is_in_gene_region(chrom, end):
        return "genic_novel"
    return "intergenic"


merged["category"] = merged["junction"].apply(classify_junction)

# ---------------------------------------------------------------------------
# 4. Output: classification table
# ---------------------------------------------------------------------------
merged.to_csv(out_classification, sep="\t", index=False)

# ---------------------------------------------------------------------------
# 5. Output: summary table (per-tool × category counts and fractions)
# ---------------------------------------------------------------------------
summary_rows = []
for tool in TOOLS:
    detected = merged[merged[tool] == 1]
    total = len(detected)
    for cat in ["annotated", "genic_novel", "intergenic"]:
        n = (detected["category"] == cat).sum()
        frac = n / total if total > 0 else 0.0
        summary_rows.append({
            "tool": tool,
            "category": cat,
            "count": n,
            "fraction": frac,
            "total": total,
        })

df_summary = pd.DataFrame(summary_rows)
df_summary.to_csv(out_summary, sep="\t", index=False)
