"""
Compute expected boundary counts for the EIB benchmark at a given anchor length.

This script is decoupled from FASTQ read generation so that multiple boundary
anchor lengths can be evaluated without re-generating reads or re-aligning.

Boundary intervals follow the new Tosa develop-branch convention:
  - 5' boundary: [intron_start - anchor_len, intron_start + anchor_len)
  - 3' boundary: [intron_end   - anchor_len, intron_end   + anchor_len)

A read (or read-pair) is counted for a boundary when its aligned segment
**fully contains** the boundary interval.

Inputs:
  - gene_metadata.tsv: per-gene information (exon/intron coordinates)
  - read_positions.tsv: per-read alignment positions (from generate_eib_reads.py)

Parameters:
  - boundary_anchor_length: bases on each side of the splice site

Outputs:
  - expected_boundary_counts.tsv: ground-truth boundary counts
"""

import bisect
import os

# ---------------------------------------------------------------------------
# Snakemake variables
# ---------------------------------------------------------------------------
in_genome = snakemake.input.genome
in_metadata = snakemake.input.metadata
in_read_positions = snakemake.input.read_positions

anchor_len = int(snakemake.params.boundary_anchor_length)

out_expected = snakemake.output.expected


def read_chrom_name(fasta_path):
    """Read just the chromosome name from a FASTA file header."""
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
        # Columns: gene_id, transcript_id, strand, n_exons, exon_sizes,
        #          intron_sizes, gene_start, gene_end, coverage
        gene_id = parts[0]
        strand = parts[2]
        n_exons = int(parts[3])
        exon_sizes = [int(x) for x in parts[4].split(",")]
        intron_sizes = [int(x) for x in parts[5].split(",")] if parts[5] != "NA" else []
        gene_start = int(parts[6])  # 1-based

        # Reconstruct exon coordinates (1-based inclusive)
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
        # 0-based intron coordinates (matching Tosa's internal representation)
        intron_start_0 = gene["exons"][j][1]
        intron_end_0 = gene["exons"][j + 1][0] - 1

        # 5' boundary: straddles the exon->intron splice site
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

        # 3' boundary: straddles the intron->exon splice site
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

# Build a sorted index of boundary intervals for efficient containment checks
boundary_starts = sorted(set(b["start"] for b in all_boundaries))
boundaries_by_start = {}
for b in all_boundaries:
    boundaries_by_start.setdefault(b["start"], []).append(b)

# ---------------------------------------------------------------------------
# 3. Load read positions and compute containment
# ---------------------------------------------------------------------------

def segment_contains_boundary(seg_start, seg_len, b_start, b_end):
    """Check if aligned segment [seg_start, seg_start+seg_len) fully contains [b_start, b_end)."""
    return seg_start <= b_start and (seg_start + seg_len) >= b_end


def find_contained_boundaries(seg_start, seg_len):
    """Find all boundary IDs whose interval is fully contained within the segment."""
    seg_end = seg_start + seg_len
    result = set()
    # Only need to check boundaries whose start >= seg_start and end <= seg_end
    lo = bisect.bisect_left(boundary_starts, seg_start)
    # Upper bound: boundary start must be < seg_end (since b_end > b_start)
    hi = bisect.bisect_right(boundary_starts, seg_end - 1)
    for idx in range(lo, hi):
        bs = boundary_starts[idx]
        for b in boundaries_by_start[bs]:
            if b["end"] <= seg_end:
                result.add(b["boundary_id"])
    return result


# Per-read boundary containment (union of R1 and R2 for PE)
read_boundary_hits = {}  # read_id -> set of boundary_ids

with open(in_read_positions) as fh:
    header = fh.readline()
    for line in fh:
        parts = line.strip().split("\t")
        rid = parts[0]
        r1_start = int(parts[1])
        r1_len = int(parts[2])
        r2_start = int(parts[3])
        r2_len = int(parts[4])

        hits = find_contained_boundaries(r1_start, r1_len)
        if r2_start >= 0 and r2_len > 0:
            hits |= find_contained_boundaries(r2_start, r2_len)

        if hits:
            read_boundary_hits[rid] = hits

# ---------------------------------------------------------------------------
# 4. Aggregate expected counts per boundary
# ---------------------------------------------------------------------------

expected_counts = {}
for rid, boundaries in read_boundary_hits.items():
    for bid in boundaries:
        expected_counts[bid] = expected_counts.get(bid, 0) + 1

# Build output table with metadata for all known boundaries
expected_rows = []
for b in all_boundaries:
    bid = b["boundary_id"]
    count = expected_counts.get(bid, 0)
    expected_rows.append({
        "boundary": bid,
        "type": b["type"],
        "strand": b["strand"],
        "expected_count": count,
        "gene_id": b["gene_id"],
    })

# ---------------------------------------------------------------------------
# 5. Write expected boundary counts
# ---------------------------------------------------------------------------
os.makedirs(os.path.dirname(out_expected), exist_ok=True)

with open(out_expected, "w") as fh:
    fh.write("boundary\ttype\tstrand\texpected_count\tgene_id\n")
    for row in expected_rows:
        fh.write(
            f"{row['boundary']}\t{row['type']}\t{row['strand']}\t"
            f"{row['expected_count']}\t{row['gene_id']}\n"
        )

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
n_total = len(expected_rows)
n_with_reads = sum(1 for r in expected_rows if r["expected_count"] > 0)
total_reads_hitting = len(read_boundary_hits)

print(f"Anchor length: {anchor_len}")
print(f"Total boundaries: {n_total} ({n_with_reads} with reads)")
print(f"Reads hitting at least one boundary: {total_reads_hitting}")
