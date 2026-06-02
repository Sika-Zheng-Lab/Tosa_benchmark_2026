"""
Generate synthetic FASTQ reads for the EIB benchmark.

For each gene with multiple exons, generates reads that span exon-intron
splice sites (simulating intron retention events). These reads align as
simple matches (e.g., 100M) because their sequence is taken directly from
the reference genome at positions crossing splice sites.

Read positions are recorded in ``read_positions.tsv`` so that expected
boundary counts can be computed separately for each boundary anchor length
(see ``compute_eib_expected_counts.py``).

Also generates negative-control reads (entirely within exons) that should
not be counted as boundary reads.

Outputs:
  - FASTQ files (PE: R1+R2, SE: reads)
  - read_positions.tsv: per-read alignment positions for downstream expected
    count computation
  - read_metadata.tsv: per-read metadata for stratified analysis
"""

import gzip
import os
import random

# ---------------------------------------------------------------------------
# Snakemake variables
# ---------------------------------------------------------------------------
in_genome = snakemake.input.genome
in_gtf = snakemake.input.gtf
in_metadata = snakemake.input.metadata

library_type = snakemake.params.library_type   # "PE" or "SE"
read_length = snakemake.params.read_length
fragment_length = snakemake.params.fragment_length
fragment_length_sd = snakemake.params.fragment_length_sd
seed = snakemake.params.seed

out_read_positions = snakemake.output.read_positions
out_read_meta = snakemake.output.read_metadata

if library_type == "PE":
    out_r1 = snakemake.output.r1
    out_r2 = snakemake.output.r2
else:
    out_reads = snakemake.output.reads

random.seed(seed)

# ---------------------------------------------------------------------------
# 1. Load genome sequence
# ---------------------------------------------------------------------------

def load_fasta(path):
    """Load single-chromosome FASTA, return (chrom_name, sequence)."""
    chrom = None
    seqs = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                chrom = line[1:].split()[0]
            else:
                seqs.append(line)
    return chrom, "".join(seqs)


CHROM, genome_seq = load_fasta(in_genome)
GENOME_LEN = len(genome_seq)
L = read_length

# ---------------------------------------------------------------------------
# 2. Load gene metadata
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
        gene_start = int(parts[6])  # 1-based
        coverage = int(parts[8])

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
            "exon_sizes": exon_sizes,
            "intron_sizes": intron_sizes,
            "exons": exons,
            "coverage": coverage,
        })

# ---------------------------------------------------------------------------
# 3. Compute splice-site positions (anchor-length agnostic)
# ---------------------------------------------------------------------------
# For each intron we record two splice-site points (0-based):
#   intron_start_0 = 1-based exon end  (first intron base)
#   intron_end_0   = 1-based next exon start - 1  (0-based half-open intron end)
# Reads must straddle these points.

splice_sites = []  # list of (gene dict, intron_start_0, intron_end_0, btype)

for gene in genes:
    for j in range(gene["n_exons"] - 1):
        intron_start_0 = gene["exons"][j][1]
        intron_end_0 = gene["exons"][j + 1][0] - 1
        splice_sites.append((gene, intron_start_0, intron_end_0, "5p"))
        splice_sites.append((gene, intron_start_0, intron_end_0, "3p"))

# ---------------------------------------------------------------------------
# 4. Helpers
# ---------------------------------------------------------------------------

def rev_comp(seq):
    """Reverse complement of a DNA sequence."""
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def make_qual(length):
    """High-quality PHRED string."""
    return "I" * length


# ---------------------------------------------------------------------------
# 5. Generate splice-site-spanning reads
# ---------------------------------------------------------------------------

reads_r1 = []
reads_r2 = []
reads_se = []
read_meta_rows = []
read_position_rows = []

read_id_counter = 0

for gene, intron_start_0, intron_end_0, btype in splice_sites:
    coverage = gene["coverage"]

    # The splice-site point that reads must straddle
    if btype == "5p":
        SP = intron_start_0   # exon->intron junction point
    else:
        SP = intron_end_0     # intron->exon junction point

    # Read at position P with length L covers [P, P+L).
    # To straddle the splice-site SP, the read must cover at least
    # 1 base on each side: P <= SP-1  and  P+L >= SP+1
    #   => P <= SP - 1  and  P >= SP + 1 - L  =>  SP + 1 - L <= P <= SP - 1
    p_min = max(0, SP + 1 - L)
    p_max = min(SP - 1, GENOME_LEN - L)

    if p_min > p_max:
        continue

    for _ in range(coverage):
        read_id_counter += 1
        rid = f"eib_br_{read_id_counter:08d}"

        P = random.randint(p_min, p_max)
        seq = genome_seq[P : P + L]
        assert len(seq) == L, f"Read at {P} has length {len(seq)} != {L}"

        if library_type == "PE":
            frag_len = max(L + 10, int(random.gauss(fragment_length, fragment_length_sd)))
            P2 = P + frag_len - L
            P2 = max(0, min(P2, GENOME_LEN - L))

            seq2 = genome_seq[P2 : P2 + L]
            assert len(seq2) == L

            reads_r1.append((f"@{rid}", seq, make_qual(L)))
            reads_r2.append((f"@{rid}", rev_comp(seq2), make_qual(L)))

            read_position_rows.append({
                "read_id": rid,
                "R1_start": P,
                "R1_length": L,
                "R2_start": P2,
                "R2_length": L,
            })
        else:
            reads_se.append((f"@{rid}", seq, make_qual(L)))

            read_position_rows.append({
                "read_id": rid,
                "R1_start": P,
                "R1_length": L,
                "R2_start": -1,
                "R2_length": 0,
            })

        read_meta_rows.append({
            "read_id": rid,
            "gene_id": gene["gene_id"],
            "target_splice_site": SP,
            "boundary_type": btype,
            "read_pos_0based": P,
            "coverage_category": coverage,
        })

# -----------------------------------------------------------------------
# Negative control reads (entirely within exons -- should NOT be counted)
# -----------------------------------------------------------------------
for gene in genes:
    coverage = gene["coverage"]
    best_exon_idx = max(range(gene["n_exons"]), key=lambda k: gene["exon_sizes"][k])
    best_exon_size = gene["exon_sizes"][best_exon_idx]

    if best_exon_size >= L:
        ex_start_0 = gene["exons"][best_exon_idx][0] - 1  # 0-based
        ex_end_0 = gene["exons"][best_exon_idx][1]         # 0-based half-open

        n_neg = max(1, coverage // 5)
        for _ in range(n_neg):
            read_id_counter += 1
            rid = f"eib_neg_{read_id_counter:08d}"

            p_max_neg = ex_end_0 - L
            if p_max_neg < ex_start_0:
                break

            P = random.randint(ex_start_0, p_max_neg)
            seq = genome_seq[P : P + L]
            assert len(seq) == L

            if library_type == "PE":
                frag_len_neg = min(
                    max(L + 10, int(random.gauss(fragment_length, fragment_length_sd))),
                    ex_end_0 - P,
                )
                P2 = P + frag_len_neg - L
                P2 = max(ex_start_0, min(P2, ex_end_0 - L))

                seq2 = genome_seq[P2 : P2 + L]
                assert len(seq2) == L

                reads_r1.append((f"@{rid}", seq, make_qual(L)))
                reads_r2.append((f"@{rid}", rev_comp(seq2), make_qual(L)))

                read_position_rows.append({
                    "read_id": rid,
                    "R1_start": P,
                    "R1_length": L,
                    "R2_start": P2,
                    "R2_length": L,
                })
            else:
                reads_se.append((f"@{rid}", seq, make_qual(L)))

                read_position_rows.append({
                    "read_id": rid,
                    "R1_start": P,
                    "R1_length": L,
                    "R2_start": -1,
                    "R2_length": 0,
                })


# ---------------------------------------------------------------------------
# 6. Write FASTQ files
# ---------------------------------------------------------------------------

if library_type == "PE":
    os.makedirs(os.path.dirname(out_r1), exist_ok=True)

    with gzip.open(out_r1, "wt") as fh:
        for name, seq, qual in reads_r1:
            fh.write(f"{name}/1\n{seq}\n+\n{qual}\n")

    with gzip.open(out_r2, "wt") as fh:
        for name, seq, qual in reads_r2:
            fh.write(f"{name}/2\n{seq}\n+\n{qual}\n")

    print(f"Generated {len(reads_r1)} PE pairs ({library_type})")
else:
    os.makedirs(os.path.dirname(out_reads), exist_ok=True)

    with gzip.open(out_reads, "wt") as fh:
        for name, seq, qual in reads_se:
            fh.write(f"{name}\n{seq}\n+\n{qual}\n")

    print(f"Generated {len(reads_se)} SE reads ({library_type})")


# ---------------------------------------------------------------------------
# 7. Write read positions
# ---------------------------------------------------------------------------
os.makedirs(os.path.dirname(out_read_positions), exist_ok=True)

with open(out_read_positions, "w") as fh:
    fh.write("read_id\tR1_start\tR1_length\tR2_start\tR2_length\n")
    for row in read_position_rows:
        fh.write(
            f"{row['read_id']}\t{row['R1_start']}\t{row['R1_length']}\t"
            f"{row['R2_start']}\t{row['R2_length']}\n"
        )


# ---------------------------------------------------------------------------
# 8. Write read metadata
# ---------------------------------------------------------------------------
os.makedirs(os.path.dirname(out_read_meta), exist_ok=True)

with open(out_read_meta, "w") as fh:
    fh.write("read_id\tgene_id\ttarget_splice_site\tboundary_type\tread_pos_0based\tcoverage_category\n")
    for row in read_meta_rows:
        fh.write(
            f"{row['read_id']}\t{row['gene_id']}\t{row['target_splice_site']}\t"
            f"{row['boundary_type']}\t{row['read_pos_0based']}\t{row['coverage_category']}\n"
        )


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
total_boundary_reads = sum(1 for r in read_position_rows if r["read_id"].startswith("eib_br_"))
total_neg_reads = sum(1 for r in read_position_rows if r["read_id"].startswith("eib_neg_"))

print(f"Boundary-spanning reads: {total_boundary_reads}")
print(f"Negative control reads: {total_neg_reads}")
print(f"Total read position entries: {len(read_position_rows)}")
print(f"Read metadata entries: {len(read_meta_rows)}")
