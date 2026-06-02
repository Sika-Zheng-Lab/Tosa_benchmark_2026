"""
Generate synthetic scRNA-seq FASTQ reads for the EIB benchmark.

Extends the bulk EIB read generation to include cell barcodes and UMIs.
Each cell gets a different coverage level for boundary-spanning reads.
The cDNA reads (R2) are placed to continuously span exon-intron boundaries
(like the bulk version), while R1 carries CB(16bp) + UMI(12bp).

Uses the SAME EIB genome/GTF from the bulk EIB benchmark.

Outputs:
  - FASTQ files (R1 = CB+UMI, R2 = cDNA)
  - Whitelist of valid cell barcodes
  - Ground truth: cell × boundary count matrix (sparse TSV)
  - Read positions TSV for expected count computation
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

n_cells = snakemake.params.n_cells
coverage_per_cell = snakemake.params.coverage_per_cell
read_length = snakemake.params.read_length
fragment_length = snakemake.params.fragment_length
fragment_length_sd = snakemake.params.fragment_length_sd
cb_length = snakemake.params.cb_length
umi_length = snakemake.params.umi_length
seed = snakemake.params.seed

out_r1 = snakemake.output.r1
out_r2 = snakemake.output.r2
out_whitelist = snakemake.output.whitelist
out_read_positions = snakemake.output.read_positions
out_expected_boundaries = snakemake.output.expected_boundaries

random.seed(seed)

BASES = "ACGT"
L = read_length

# ---------------------------------------------------------------------------
# 1. Generate cell barcodes
# ---------------------------------------------------------------------------

def generate_unique_barcodes(n, length):
    barcodes = set()
    while len(barcodes) < n:
        bc = "".join(random.choice(BASES) for _ in range(length))
        barcodes.add(bc)
    return sorted(barcodes)


cell_barcodes = generate_unique_barcodes(n_cells, cb_length)

os.makedirs(os.path.dirname(out_whitelist), exist_ok=True)
with open(out_whitelist, "w") as fh:
    for bc in cell_barcodes:
        fh.write(bc + "\n")


# ---------------------------------------------------------------------------
# 2. Load genome sequence
# ---------------------------------------------------------------------------

def load_fasta(path):
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


# ---------------------------------------------------------------------------
# 3. Load gene metadata
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
        coverage = int(parts[8])

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
# 4. Compute splice-site positions (same as bulk EIB)
# ---------------------------------------------------------------------------

splice_sites = []
for gene in genes:
    for j in range(gene["n_exons"] - 1):
        intron_start_0 = gene["exons"][j][1]
        intron_end_0 = gene["exons"][j + 1][0] - 1
        splice_sites.append((gene, intron_start_0, intron_end_0, "5p"))
        splice_sites.append((gene, intron_start_0, intron_end_0, "3p"))


# ---------------------------------------------------------------------------
# 5. Helpers
# ---------------------------------------------------------------------------

def rev_comp(seq):
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def make_qual(length):
    return "I" * length


def generate_umi(length):
    return "".join(random.choice(BASES) for _ in range(length))


# ---------------------------------------------------------------------------
# 6. Generate reads per cell (streaming writes to avoid high memory usage)
# ---------------------------------------------------------------------------

boundary_counts = {}  # (barcode, boundary_str, boundary_type) -> count
read_id_counter = 0

os.makedirs(os.path.dirname(out_r1), exist_ok=True)
os.makedirs(os.path.dirname(out_read_positions), exist_ok=True)

with gzip.open(out_r1, "wt") as fh_r1, \
     gzip.open(out_r2, "wt") as fh_r2, \
     open(out_read_positions, "w") as fh_pos:
    fh_pos.write("read_id\tbarcode\tR1_start\tR1_length\tR2_start\tR2_length\n")

    for cell_idx, barcode in enumerate(cell_barcodes):
        # Each cell gets a random coverage level
        cell_cov = random.choice(coverage_per_cell)
        if cell_cov == 0:
            continue

        for gene, intron_start_0, intron_end_0, btype in splice_sites:
            # Boundary coordinate string — must match Tosa's convention:
            # boundary = [splice_site - anchor_len, splice_site + anchor_len]
            # With default anchor_len=1, this is a 3bp span.
            if btype == "5p":
                SP = intron_start_0
                boundary_str = f"{CHROM}:{intron_start_0 - 1}-{intron_start_0 + 1}"
            else:
                SP = intron_end_0
                boundary_str = f"{CHROM}:{intron_end_0 - 1}-{intron_end_0 + 1}"

            p_min = max(0, SP + 1 - L)
            p_max = min(SP - 1, GENOME_LEN - L)

            if p_min > p_max:
                continue

            for _ in range(cell_cov):
                read_id_counter += 1
                rid = f"scrna_eib_{read_id_counter:08d}"

                P = random.randint(p_min, p_max)
                seq = genome_seq[P : P + L]
                assert len(seq) == L

                umi = generate_umi(umi_length)

                # R1: CB + UMI
                r1_seq = barcode + umi

                # R2: cDNA (reverse complement of genomic sequence)
                # For PE, generate second read
                frag_len = max(L + 10, int(random.gauss(fragment_length, fragment_length_sd)))
                P2 = P + frag_len - L
                P2 = max(0, min(P2, GENOME_LEN - L))
                seq2 = genome_seq[P2 : P2 + L]
                assert len(seq2) == L

                r1_qual = make_qual(len(r1_seq))
                r2_qual = make_qual(L)
                fh_r1.write(f"@{rid}/1\n{r1_seq}\n+\n{r1_qual}\n")
                fh_r2.write(f"@{rid}/2\n{rev_comp(seq)}\n+\n{r2_qual}\n")

                fh_pos.write(
                    f"{rid}\t{barcode}\t{P}\t{L}\t{P2}\t{L}\n"
                )

            # Track boundary counts for this cell
            key = (barcode, boundary_str, btype)
            boundary_counts[key] = boundary_counts.get(key, 0) + cell_cov

print(f"Generated {read_id_counter} PE pairs for {n_cells} cells")

# ---------------------------------------------------------------------------
# 7. Write expected boundary counts (sparse: barcode × boundary)
# ---------------------------------------------------------------------------
os.makedirs(os.path.dirname(out_expected_boundaries), exist_ok=True)
with open(out_expected_boundaries, "w") as fh:
    fh.write("barcode\tboundary\tboundary_type\tcount\n")
    for (bc, boundary_str, btype), count in sorted(boundary_counts.items()):
        fh.write(f"{bc}\t{boundary_str}\t{btype}\t{count}\n")

print(f"Ground truth: {len(boundary_counts)} non-zero cell×boundary entries")
