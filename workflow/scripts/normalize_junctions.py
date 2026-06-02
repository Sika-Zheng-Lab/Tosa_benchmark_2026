"""Normalize junction outputs from STAR, Tosa, regtools, and featureCounts
into a common format.

Common output format (TSV, no header in individual files):
  junction    strand    count

Junction coordinates: 1-based inclusive (chr:start-end)
"""

import gzip
import os
import pandas as pd

# ---------------------------------------------------------------------------
# Snakemake I/O
# ---------------------------------------------------------------------------
in_star = snakemake.input.star
in_tosa = snakemake.input.tosa
in_regtools = snakemake.input.regtools
in_featurecounts = snakemake.input.featurecounts

out_star = snakemake.output.star
out_tosa = snakemake.output.tosa
out_regtools = snakemake.output.regtools
out_featurecounts = snakemake.output.featurecounts


# ---------------------------------------------------------------------------
# 1. STAR SJ.out.tab
# ---------------------------------------------------------------------------
# Columns: chr, start(1-based), end(1-based), strand(0/1/2), motif, annotated,
#          unique_reads, multi_reads, max_overhang
# Already 1-based inclusive — no conversion needed.
# strand: 0=undefined, 1=+, 2=-

strand_map_star = {0: ".", 1: "+", 2: "-"}

rows_star = []
with open(in_star) as fh:
    for line in fh:
        parts = line.strip().split("\t")
        chrom = parts[0]
        start = parts[1]
        end = parts[2]
        strand_code = int(parts[3])
        unique_count = int(parts[6])
        if unique_count == 0:
            continue
        junction = f"{chrom}:{start}-{end}"
        strand = strand_map_star.get(strand_code, ".")
        rows_star.append((junction, strand, unique_count))

df_star = pd.DataFrame(rows_star, columns=["junction", "strand", "count"])
os.makedirs(os.path.dirname(out_star), exist_ok=True)
df_star.to_csv(out_star, sep="\t", index=False)


# ---------------------------------------------------------------------------
# 2. Tosa junction.tsv.gz
# ---------------------------------------------------------------------------
# Columns: Junction, Strand, Count (develop branch)
# Already 1-based inclusive — no conversion needed.

df_tosa = pd.read_csv(in_tosa, sep="\t", compression="gzip")
df_tosa.columns = ["junction", "strand", "count"]
df_tosa.to_csv(out_tosa, sep="\t", index=False)


# ---------------------------------------------------------------------------
# 3. regtools junctions BED12
# ---------------------------------------------------------------------------
# BED12 columns: chrom, chromStart, chromEnd, name, score, strand,
#                thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
# Junction coords (0-based):
#   junction_start = chromStart + blockSizes[0]
#   junction_end   = chromEnd - blockSizes[1]
# Convert to 1-based inclusive:
#   start_1based = junction_start + 1
#   end_1based   = junction_end       (since half-open end == inclusive end in 1-based)

rows_regtools = []
with open(in_regtools) as fh:
    for line in fh:
        if line.startswith("#") or line.startswith("track"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 12:
            continue
        chrom = parts[0]
        chrom_start = int(parts[1])
        chrom_end = int(parts[2])
        score = int(parts[4])
        strand = parts[5]
        block_sizes = parts[10].rstrip(",").split(",")
        block_size_left = int(block_sizes[0])
        block_size_right = int(block_sizes[1])

        junction_start_0 = chrom_start + block_size_left
        junction_end_0 = chrom_end - block_size_right

        # Convert to 1-based inclusive
        start_1 = junction_start_0 + 1
        end_1 = junction_end_0

        junction = f"{chrom}:{start_1}-{end_1}"
        rows_regtools.append((junction, strand, score))

df_regtools = pd.DataFrame(rows_regtools, columns=["junction", "strand", "count"])
# regtools may list the same junction coords from different records —
# group and sum just in case
df_regtools = df_regtools.groupby(["junction", "strand"], as_index=False)["count"].sum()
df_regtools.to_csv(out_regtools, sep="\t", index=False)


# ---------------------------------------------------------------------------
# 4. featureCounts .jcounts
# ---------------------------------------------------------------------------
# The .jcounts file produced by featureCounts -J has columns:
#   Gene_SP1, Gene_SP2, Transcript, Status, Donor, Acceptor,
#   Chr_SP1, Location_SP1, Strand_SP1,
#   Chr_SP2, Location_SP2, Strand_SP2,
#   NearestExonBoundary_SP1, NearestExonBoundary_SP2, <BAM_file_name>
# Location_SP1 = last base of upstream exon (1-based)
# Location_SP2 = first base of downstream exon (1-based)
# Intron in 1-based inclusive coords:
#   start = Location_SP1 + 1
#   end   = Location_SP2 - 1
# The count column is the last column (named after the BAM path).

df_fc = pd.read_csv(in_featurecounts, sep="\t", comment="#")

# The count column is the last column (named after the BAM file path)
count_col = df_fc.columns[-1]

rows_fc = []
for _, row in df_fc.iterrows():
    chrom = row["Chr_SP1"]
    site1 = int(row["Location_SP1"])
    site2 = int(row["Location_SP2"])
    strand = row["Strand_SP1"]
    count = int(row[count_col])
    if count == 0:
        continue
    # Convert exon boundary positions to intron coordinates
    intron_start = site1 + 1
    intron_end = site2 - 1
    junction = f"{chrom}:{intron_start}-{intron_end}"
    # Strand may be "NA" for unstranded data
    if str(strand) == "NA" or str(strand) == "nan":
        strand = "."
    rows_fc.append((junction, strand, count))

df_featurecounts = pd.DataFrame(rows_fc, columns=["junction", "strand", "count"])
# Group and sum in case of duplicate junction entries
df_featurecounts = df_featurecounts.groupby(
    ["junction", "strand"], as_index=False
)["count"].sum()
os.makedirs(os.path.dirname(out_featurecounts), exist_ok=True)
df_featurecounts.to_csv(out_featurecounts, sep="\t", index=False)
