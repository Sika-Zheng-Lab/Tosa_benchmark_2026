"""
Normalize scRNA-seq junction outputs from STARsolo, Tosa single, and regtools
into a common sparse cell × junction matrix format.

Each tool's output is converted to:
  barcode\tjunction\tcount

Junction coordinates: 1-based inclusive (chr:start-end), matching the bulk
normalization convention.

## Output formats of each tool:

### STARsolo SJ
  Solo.out/SJ/raw/ contains:
    - features.tsv: junction features (chr, start_0based, end_0based_exclusive, strand)
    - barcodes.tsv: cell barcode list
    - matrix.mtx: MatrixMarket sparse matrix (feature × barcode)
  Junction coords in features.tsv are 0-based half-open → convert to 1-based inclusive.

### Tosa single
  {prefix}_junction_barcodes.tsv.gz: Feature, Strand, Barcode, Count
  Feature column is the junction coordinate string (chr:start-end), already
  1-based inclusive (same as bulk Tosa).

### regtools -b
  BED12 file: aggregate junction counts (same format as bulk)
  Barcode file (same line order as BED12): <num_barcodes>\t<bc1:count1,bc2:count2,...>
  Junction coords from BED12: 0-based → convert to 1-based inclusive.
  Barcode file lines are 1:1 matched with BED12 lines.
"""

import gzip
import os
import sys

import pandas as pd
from scipy import io as spio
from scipy.sparse import csr_matrix
import numpy as np

# ---------------------------------------------------------------------------
# Snakemake I/O
# ---------------------------------------------------------------------------
in_starsolo_sj_raw = snakemake.input.starsolo_sj_raw
in_tosa = snakemake.input.tosa
in_regtools_bed = snakemake.input.regtools_bed
in_regtools_barcodes = snakemake.input.regtools_barcodes
in_whitelist = snakemake.input.whitelist

out_starsolo = snakemake.output.starsolo
out_tosa = snakemake.output.tosa
out_regtools = snakemake.output.regtools

# ---------------------------------------------------------------------------
# Load valid cell barcodes (whitelist or STARsolo filtered barcodes)
# ---------------------------------------------------------------------------
with open(in_whitelist) as fh:
    valid_barcodes = set(line.strip().split("\t")[0] for line in fh if line.strip())

print(f"Loaded {len(valid_barcodes)} valid cell barcodes", file=sys.stderr)


def write_sparse_tsv(records, outpath):
    """Write list of (barcode, junction, count) tuples to TSV."""
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, "w") as fh:
        fh.write("barcode\tjunction\tcount\n")
        for bc, junc, count in sorted(records):
            fh.write(f"{bc}\t{junc}\t{count}\n")


# ---------------------------------------------------------------------------
# 1. STARsolo SJ (Solo.out/SJ/raw/)
# ---------------------------------------------------------------------------

sj_dir = in_starsolo_sj_raw
features_path = os.path.join(sj_dir, "features.tsv.gz")
barcodes_path = os.path.join(sj_dir, "barcodes.tsv.gz")
matrix_path = os.path.join(sj_dir, "matrix.mtx.gz")

# Also try uncompressed if gzipped versions don't exist
if not os.path.exists(features_path):
    features_path = os.path.join(sj_dir, "features.tsv")
if not os.path.exists(barcodes_path):
    barcodes_path = os.path.join(sj_dir, "barcodes.tsv")
if not os.path.exists(matrix_path):
    matrix_path = os.path.join(sj_dir, "matrix.mtx")

# Load features (junction definitions)
# STARsolo SJ features.tsv has the same format as SJ.out.tab:
#   chr\tstart\tend\tstrand\tmotif\tannotated\tunique_count\tmulti_count\toverhang
# start and end are 1-based intron start/end (same as SJ.out.tab convention).
# We build junction strings as "chr:start-end" (1-based inclusive).
features = []
open_fn = gzip.open if features_path.endswith(".gz") else open
with open_fn(features_path, "rt") as fh:
    for line in fh:
        parts = line.strip().split("\t")
        chrom = parts[0]
        start = parts[1]
        end = parts[2]
        junction_str = f"{chrom}:{start}-{end}"
        features.append(junction_str)

# Load barcodes
barcodes_star = []
open_fn = gzip.open if barcodes_path.endswith(".gz") else open
with open_fn(barcodes_path, "rt") as fh:
    for line in fh:
        barcodes_star.append(line.strip())

# Load sparse matrix (MatrixMarket format: feature × barcode)
open_fn = gzip.open if matrix_path.endswith(".gz") else open
with open_fn(matrix_path, "rt") as fh:
    mtx = spio.mmread(fh)

# Convert to CSR for efficient row slicing
mtx_csr = csr_matrix(mtx)

records_starsolo = []
for bc_idx in range(mtx_csr.shape[1]):
    bc = barcodes_star[bc_idx]
    if bc not in valid_barcodes:
        continue
    col = mtx_csr.getcol(bc_idx).toarray().ravel()
    for feat_idx in range(len(col)):
        if col[feat_idx] > 0:
            records_starsolo.append((bc, features[feat_idx], int(col[feat_idx])))

write_sparse_tsv(records_starsolo, out_starsolo)
print(f"STARsolo: {len(records_starsolo)} non-zero entries", file=sys.stderr)


# ---------------------------------------------------------------------------
# 2. Tosa single — {prefix}_junction_barcodes.tsv.gz
# ---------------------------------------------------------------------------
# Columns: Feature, Strand, Barcode, Count
# Feature is the junction coordinate string (chr:start-end), already 1-based inclusive

df_tosa = pd.read_csv(in_tosa, sep="\t", compression="gzip")
# Normalize column names (Tosa develop branch)
col_map = {}
for c in df_tosa.columns:
    cl = c.lower()
    if "barcode" in cl or "cell" in cl:
        col_map[c] = "barcode"
    elif "feature" in cl or "junction" in cl:
        col_map[c] = "junction"
    elif "strand" in cl:
        col_map[c] = "strand"
    elif "count" in cl:
        col_map[c] = "count"
df_tosa = df_tosa.rename(columns=col_map)

records_tosa = []
for _, row in df_tosa.iterrows():
    bc = str(row["barcode"])
    if bc not in valid_barcodes:
        continue
    junction = str(row["junction"])
    count = int(row["count"])
    if count > 0:
        records_tosa.append((bc, junction, count))

write_sparse_tsv(records_tosa, out_tosa)
print(f"Tosa single: {len(records_tosa)} non-zero entries", file=sys.stderr)


# ---------------------------------------------------------------------------
# 3. regtools BED12 + barcode file
# ---------------------------------------------------------------------------
# BED12 and barcode file have matching line order.
# BED12 → junction coords (0-based → 1-based inclusive)
# Barcode file line format: <num_barcodes>\t<bc1:count1,bc2:count2,...>

records_regtools = []

bed_lines = []
with open(in_regtools_bed) as fh:
    for line in fh:
        if line.startswith("#") or line.startswith("track"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 12:
            continue
        bed_lines.append(parts)

barcode_lines = []
with open(in_regtools_barcodes) as fh:
    for line in fh:
        barcode_lines.append(line.strip())

if len(bed_lines) != len(barcode_lines):
    print(
        f"WARNING: BED12 has {len(bed_lines)} junctions but barcode file "
        f"has {len(barcode_lines)} lines. Using minimum.",
        file=sys.stderr,
    )

for i in range(min(len(bed_lines), len(barcode_lines))):
    parts = bed_lines[i]
    chrom = parts[0]
    chrom_start = int(parts[1])
    chrom_end = int(parts[2])
    strand = parts[5]
    block_sizes = parts[10].rstrip(",").split(",")
    block_size_left = int(block_sizes[0])
    block_size_right = int(block_sizes[1])

    # Convert BED12 to 1-based inclusive junction coords
    junction_start_0 = chrom_start + block_size_left
    junction_end_0 = chrom_end - block_size_right
    start_1 = junction_start_0 + 1
    end_1 = junction_end_0
    junction_str = f"{chrom}:{start_1}-{end_1}"

    # Parse barcode line
    bc_line = barcode_lines[i]
    if "\t" not in bc_line:
        continue
    bc_parts = bc_line.split("\t")
    if len(bc_parts) < 2:
        continue
    # num_barcodes = int(bc_parts[0])  # not strictly needed
    bc_entries = bc_parts[1].split(",")

    for entry in bc_entries:
        if ":" not in entry:
            continue
        bc, count_str = entry.rsplit(":", 1)
        bc = bc.strip()
        try:
            count = int(count_str.strip())
        except ValueError:
            continue
        if bc not in valid_barcodes:
            continue
        if count > 0:
            records_regtools.append((bc, junction_str, count))

write_sparse_tsv(records_regtools, out_regtools)
print(f"regtools: {len(records_regtools)} non-zero entries", file=sys.stderr)
