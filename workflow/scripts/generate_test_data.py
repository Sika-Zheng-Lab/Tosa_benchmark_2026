"""
Generate synthetic test data for Tosa benchmark.

Creates:
- A minimal genome with 3 exons and 2 introns
- Paired-end FASTQ reads with controlled junction-spanning patterns
- Expected junction count table per tool

Genome layout (chr_test, 1600 bp):
  Exon1:   1 -  200  (200 bp)
  Intron1: 201 - 700  (500 bp)
  Exon2:   701 - 900  (200 bp)
  Intron2: 901 - 1400 (500 bp)
  Exon3:   1401 - 1600 (200 bp)

mRNA transcript (600 bp):
  Exon1(1-200) | Exon2(201-400) | Exon3(401-600)
  Junction1 at mRNA 200/201  => genomic intron 201-700  (1-based inclusive)
  Junction2 at mRNA 400/401  => genomic intron 901-1400 (1-based inclusive)

Read cases (all target Junction1 = chr_test:201-700):
  Case A: Both mates span junction1
  Case B: Only R1 spans junction1
  Case C: Neither mate spans a junction (negative control)
"""

import gzip
import os
import random

# ---------------------------------------------------------------------------
# Access Snakemake variables
# ---------------------------------------------------------------------------
read_length = snakemake.params.read_length
num_both = snakemake.params.num_both
num_one = snakemake.params.num_one
num_none = snakemake.params.num_none

out_genome = snakemake.output.genome
out_fai = snakemake.output.fai
out_gtf = snakemake.output.gtf
out_r1 = snakemake.output.r1
out_r2 = snakemake.output.r2
out_expected = snakemake.output.expected

# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------
random.seed(42)

# ---------------------------------------------------------------------------
# 1. Build synthetic genome
# ---------------------------------------------------------------------------
EXON_LEN = 200
INTRON_LEN = 500
CHROM = "chr_test"

def random_seq(length):
    """Generate a random DNA sequence avoiding long poly-N or low complexity."""
    return "".join(random.choice("ACGT") for _ in range(length))

exon1 = random_seq(EXON_LEN)
# Introns MUST start with GT (donor) and end with AG (acceptor)
# so that STAR recognises the canonical splice signals.
intron1 = "GT" + random_seq(INTRON_LEN - 4) + "AG"
exon2 = random_seq(EXON_LEN)
intron2 = "GT" + random_seq(INTRON_LEN - 4) + "AG"
exon3 = random_seq(EXON_LEN)

genome_seq = exon1 + intron1 + exon2 + intron2 + exon3
genome_len = len(genome_seq)  # 1600

# mRNA (spliced transcript)
mrna = exon1 + exon2 + exon3  # 600 bp

os.makedirs(os.path.dirname(out_genome), exist_ok=True)

with open(out_genome, "w") as fh:
    fh.write(f">{CHROM}\n")
    for i in range(0, genome_len, 80):
        fh.write(genome_seq[i : i + 80] + "\n")

# FAI: name\tlength\toffset\tlinebases\tlinewidth
header_len = len(f">{CHROM}\n")
with open(out_fai, "w") as fh:
    fh.write(f"{CHROM}\t{genome_len}\t{header_len}\t80\t81\n")

# GTF annotation — needed by featureCounts -J
# gene_id "test_gene"; transcript_id "test_tx";
os.makedirs(os.path.dirname(out_gtf), exist_ok=True)
with open(out_gtf, "w") as fh:
    attr_gene = 'gene_id "test_gene";'
    attr_tx = 'gene_id "test_gene"; transcript_id "test_tx";'
    # Gene
    fh.write(f"{CHROM}\tsynthetic\tgene\t1\t{genome_len}\t.\t+\t.\t{attr_gene}\n")
    # Transcript
    fh.write(f"{CHROM}\tsynthetic\ttranscript\t1\t{genome_len}\t.\t+\t.\t{attr_tx}\n")
    # Exon1: 1-200
    fh.write(f"{CHROM}\tsynthetic\texon\t1\t{EXON_LEN}\t.\t+\t.\t{attr_tx}\n")
    # Exon2: 701-900
    exon2_start = EXON_LEN + INTRON_LEN + 1  # 701
    exon2_end = exon2_start + EXON_LEN - 1   # 900
    fh.write(f"{CHROM}\tsynthetic\texon\t{exon2_start}\t{exon2_end}\t.\t+\t.\t{attr_tx}\n")
    # Exon3: 1401-1600
    exon3_start = 2 * EXON_LEN + 2 * INTRON_LEN + 1  # 1401
    exon3_end = exon3_start + EXON_LEN - 1           # 1600
    fh.write(f"{CHROM}\tsynthetic\texon\t{exon3_start}\t{exon3_end}\t.\t+\t.\t{attr_tx}\n")

# ---------------------------------------------------------------------------
# 2. Helper: convert mRNA coords to genomic CIGAR-aware reads
# ---------------------------------------------------------------------------

def mrna_to_genomic_cigar(mrna_start, mrna_end):
    """
    Given 0-based mRNA start/end, return (genomic_start_0based, cigar_string, seq).

    Exon1: mRNA 0-199   -> genome 0-199
    Exon2: mRNA 200-399 -> genome 700-899
    Exon3: mRNA 400-599 -> genome 1400-1599
    Intron1: genome 200-699 (500 bp)
    Intron2: genome 900-1399 (500 bp)
    """
    exon_bounds_mrna = [(0, 200), (200, 400), (400, 600)]
    exon_starts_genome = [0, 700, 1400]
    intron_sizes = [500, 500]  # between exon pairs

    segments = []  # list of (genome_start, length)
    pos = mrna_start
    while pos < mrna_end:
        for idx, (es, ee) in enumerate(exon_bounds_mrna):
            if es <= pos < ee:
                seg_end = min(mrna_end, ee)
                offset = pos - es
                gstart = exon_starts_genome[idx] + offset
                length = seg_end - pos
                segments.append((gstart, length, idx))
                pos = seg_end
                break

    # Build CIGAR
    cigar = ""
    genomic_start = segments[0][0]
    for i, (gs, ln, eidx) in enumerate(segments):
        if i > 0:
            prev_eidx = segments[i - 1][2]
            if eidx != prev_eidx:
                # Insert N (intron skip) for each skipped intron
                for j in range(prev_eidx, eidx):
                    cigar += f"{intron_sizes[j]}N"
        cigar += f"{ln}M"

    seq = mrna[mrna_start:mrna_end]
    return genomic_start, cigar, seq


def rev_comp(seq):
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def make_qual(length):
    return "I" * length


# ---------------------------------------------------------------------------
# 3. Generate FASTQ reads
# ---------------------------------------------------------------------------

reads_r1 = []
reads_r2 = []
pair_id = 0

# Case A: Both R1 and R2 span junction1 (mRNA pos 200)
# Fragment on mRNA: 151-270 (120 bp insert)
# R1: mRNA 151-250 (100 bp), spans junction at 200
# R2: mRNA 171-270 (100 bp), spans junction at 200 (reversed)
for i in range(num_both):
    pair_id += 1
    name = f"@read_both_{pair_id}"

    r1_start, r1_end = 151, 151 + read_length  # 151-250
    r2_start, r2_end = 171, 171 + read_length  # 171-270

    g1_start, cigar1, seq1 = mrna_to_genomic_cigar(r1_start, r1_end)
    g2_start, cigar2, seq2 = mrna_to_genomic_cigar(r2_start, r2_end)

    # R1 forward, R2 reverse complement
    reads_r1.append((name, seq1, make_qual(len(seq1))))
    reads_r2.append((name, rev_comp(seq2), make_qual(len(seq2))))

# Case B: Only R1 spans junction1
# R1: mRNA 151-250 (spans junction at 200)
# R2: mRNA 301-400 (entirely in exon2, no junction)
for i in range(num_one):
    pair_id += 1
    name = f"@read_one_{pair_id}"

    r1_start, r1_end = 151, 151 + read_length  # 151-250
    r2_start, r2_end = 301, 301 + read_length  # 301-400

    g1_start, cigar1, seq1 = mrna_to_genomic_cigar(r1_start, r1_end)
    g2_start, cigar2, seq2 = mrna_to_genomic_cigar(r2_start, r2_end)

    reads_r1.append((name, seq1, make_qual(len(seq1))))
    reads_r2.append((name, rev_comp(seq2), make_qual(len(seq2))))

# Case C: No junction — both mates in exon1
# R1: mRNA 10-109
# R2: mRNA 50-149
for i in range(num_none):
    pair_id += 1
    name = f"@read_none_{pair_id}"

    r1_start, r1_end = 10, 10 + read_length   # 10-109
    r2_start, r2_end = 50, 50 + read_length   # 50-149

    g1_start, cigar1, seq1 = mrna_to_genomic_cigar(r1_start, r1_end)
    g2_start, cigar2, seq2 = mrna_to_genomic_cigar(r2_start, r2_end)

    reads_r1.append((name, seq1, make_qual(len(seq1))))
    reads_r2.append((name, rev_comp(seq2), make_qual(len(seq2))))

# Write FASTQ files
os.makedirs(os.path.dirname(out_r1), exist_ok=True)

with gzip.open(out_r1, "wt") as fh:
    for name, seq, qual in reads_r1:
        fh.write(f"{name}/1\n{seq}\n+\n{qual}\n")

with gzip.open(out_r2, "wt") as fh:
    for name, seq, qual in reads_r2:
        fh.write(f"{name}/2\n{seq}\n+\n{qual}\n")

# ---------------------------------------------------------------------------
# 4. Write expected counts
# ---------------------------------------------------------------------------
# Junction1 in 1-based inclusive coords: chr_test:201-700
junction1 = f"{CHROM}:201-700"

# Expected counts (deduplicated — count each PE pair once regardless of tool):
# Tool           | Case A (both mates) | Case B (one mate) | Case C (none) | Total
# ---------------|---------------------|-------------------|---------------|------
# Tosa           | num_both            | num_one           | 0             | num_both + num_one
# STAR           | num_both            | num_one           | 0             | num_both + num_one
# regtools       | num_both            | num_one           | 0             | num_both + num_one
# featurecounts  | num_both            | num_one           | 0             | num_both + num_one

os.makedirs(os.path.dirname(out_expected), exist_ok=True)

with open(out_expected, "w") as fh:
    fh.write("junction\ttool\texpected_count\tcaseA_count\tcaseB_count\tcaseC_count\n")
    for tool in ["tosa", "star", "regtools", "featurecounts"]:
        fh.write(f"{junction1}\t{tool}\t{num_both + num_one}\t{num_both}\t{num_one}\t0\n")

print(f"Generated synthetic genome: {genome_len} bp")
print(f"Generated {pair_id} read pairs ({num_both} both, {num_one} one, {num_none} none)")
print(f"Expected junction: {junction1}")
