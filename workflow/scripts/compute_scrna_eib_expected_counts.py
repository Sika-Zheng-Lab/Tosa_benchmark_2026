"""
Compute expected scRNA boundary counts per cell for a given anchor length.

Analogous to compute_eib_expected_counts.py but tracks cell barcodes.
A read is counted for a boundary only if its aligned segment fully covers
the boundary interval [splice_site - anchor_len, splice_site + anchor_len).

Fragment-level deduplication: since both R1 (cDNA) and R2 can hit a boundary,
we take the union of boundary hits per read-pair (fragment), consistent with
Tosa's fragment-aware counting.

Inputs:
  - gene_metadata.tsv
  - read_positions.tsv (with barcode column)
  - genome FASTA (for chromosome name)

Parameters:
  - boundary_anchor_length

Outputs:
  - expected_boundary_counts.tsv (barcode, boundary, boundary_type, count)
"""

import bisect
import os

in_genome = snakemake.input.genome
in_metadata = snakemake.input.metadata
in_read_positions = snakemake.input.read_positions

anchor_len = int(snakemake.params.boundary_anchor_length)

out_expected = snakemake.output.expected


def read_chrom_name(fasta_path):
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                return line[1:].strip().split()[0]
    raise ValueError(f"No FASTA header found in {fasta_path}")


CHROM = read_chrom_name(in_genome)

# ---------------------------------------------------------------------------
# 1. Load gene metadata and reconstruct intron coordinates
# ---------------------------------------------------------------------------

genes = []
with open(in_metadata) as fh:
    header = fh.readline()
    for line in fh:
        parts = line.strip().split("\t")
        gene_id = parts[0]
        strand = parts[2]
        n_exons = int(parts[3])
        exon_sizes = [int(x) for x in parts[4].split(",")]
        intron_sizes = [int(x) for x in parts[5].split(",")] if parts[5] != "NA" else []
        gene_start = int(parts[6])

        exons = []
        pos = gene_start
        for j in range(n_exons):
            ex_start = pos
            ex_end = pos + exon_sizes[j] - 1
            exons.append((ex_start, ex_end))
            pos = ex_end + 1
            if j < n_exons - 1:
                pos += intron_sizes[j]

        genes.append({
            "gene_id": gene_id,
            "strand": strand,
            "n_exons": n_exons,
            "exons": exons,
        })

# ---------------------------------------------------------------------------
# 2. Compute boundary intervals at the specified anchor length
# ---------------------------------------------------------------------------

all_boundaries = []

for gene in genes:
    for j in range(gene["n_exons"] - 1):
        intron_start_0 = gene["exons"][j][1]
        intron_end_0 = gene["exons"][j + 1][0] - 1

        b5_start = intron_start_0 - anchor_len
        b5_end = intron_start_0 + anchor_len
        b5_id = f"{CHROM}:{b5_start}-{b5_end}"
        all_boundaries.append({
            "boundary_id": b5_id,
            "start": b5_start,
            "end": b5_end,
            "gene_id": gene["gene_id"],
            "strand": gene["strand"],
            "type": "5p",
        })

        b3_start = intron_end_0 - anchor_len
        b3_end = intron_end_0 + anchor_len
        b3_id = f"{CHROM}:{b3_start}-{b3_end}"
        all_boundaries.append({
            "boundary_id": b3_id,
            "start": b3_start,
            "end": b3_end,
            "gene_id": gene["gene_id"],
            "strand": gene["strand"],
            "type": "3p",
        })

boundary_starts = sorted(set(b["start"] for b in all_boundaries))
boundaries_by_start = {}
for b in all_boundaries:
    boundaries_by_start.setdefault(b["start"], []).append(b)

# ---------------------------------------------------------------------------
# 3. Containment helpers
# ---------------------------------------------------------------------------


def find_contained_boundaries(seg_start, seg_len):
    seg_end = seg_start + seg_len
    result = set()
    lo = bisect.bisect_left(boundary_starts, seg_start)
    hi = bisect.bisect_right(boundary_starts, seg_end - 1)
    for idx in range(lo, hi):
        bs = boundary_starts[idx]
        for b in boundaries_by_start[bs]:
            if b["end"] <= seg_end:
                result.add(b["boundary_id"])
    return result


# ---------------------------------------------------------------------------
# 4. Load read positions and compute per-cell boundary counts
# ---------------------------------------------------------------------------

# (barcode, boundary_id) -> count
cell_boundary_counts = {}

with open(in_read_positions) as fh:
    header = fh.readline()
    for line in fh:
        parts = line.strip().split("\t")
        rid = parts[0]
        barcode = parts[1]
        r1_start = int(parts[2])
        r1_len = int(parts[3])
        r2_start = int(parts[4])
        r2_len = int(parts[5])

        hits = find_contained_boundaries(r1_start, r1_len)
        if r2_start >= 0 and r2_len > 0:
            hits |= find_contained_boundaries(r2_start, r2_len)

        for bid in hits:
            key = (barcode, bid)
            cell_boundary_counts[key] = cell_boundary_counts.get(key, 0) + 1

# ---------------------------------------------------------------------------
# 5. Build boundary type lookup
# ---------------------------------------------------------------------------

boundary_type_map = {b["boundary_id"]: b["type"] for b in all_boundaries}

# ---------------------------------------------------------------------------
# 6. Write expected boundary counts
# ---------------------------------------------------------------------------

os.makedirs(os.path.dirname(out_expected), exist_ok=True)

with open(out_expected, "w") as fh:
    fh.write("barcode\tboundary\tboundary_type\tcount\n")
    for (barcode, bid), count in sorted(cell_boundary_counts.items()):
        btype = boundary_type_map.get(bid, ".")
        fh.write(f"{barcode}\t{bid}\t{btype}\t{count}\n")

n_entries = len(cell_boundary_counts)
print(f"Anchor length: {anchor_len}")
print(f"Total cell x boundary entries with reads: {n_entries}")
