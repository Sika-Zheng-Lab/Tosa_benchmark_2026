"""
Generate synthetic scRNA-seq test data for benchmarking.

Creates:
- A synthetic genome with multiple genes (each with introns)
- A set of cell barcodes (10x Chromium v3 style: 16 bp)
- Paired-end FASTQ reads: R1 = CB(16bp) + UMI(12bp), R2 = cDNA (junction-spanning)
- Ground truth: cell × junction count matrix
- Whitelist file of valid cell barcodes

Genome layout:
  Multiple genes on a single chromosome, each with variable exon/intron
  counts. Genes are placed sequentially with 1 kb spacers.

Read generation:
  For each cell, for each gene/junction, generate a random number of
  reads (drawn from coverage_per_cell options). Each cDNA read (R2) is
  derived from spliced mRNA so STAR maps it with N-cigar (splice junction).
  R1 carries the cell barcode + UMI.
"""

import gzip
import itertools
import os
import random
import string

# ---------------------------------------------------------------------------
# Snakemake variables
# ---------------------------------------------------------------------------
out_genome = snakemake.output.genome
out_fai = snakemake.output.fai
out_gtf = snakemake.output.gtf
out_r1 = snakemake.output.r1
out_r2 = snakemake.output.r2
out_whitelist = snakemake.output.whitelist
out_expected_junctions = snakemake.output.expected_junctions
out_gene_metadata = snakemake.output.gene_metadata

n_cells = snakemake.params.n_cells
n_genes = snakemake.params.n_genes
exon_counts = snakemake.params.exon_counts
exon_size_range = snakemake.params.exon_size_range
intron_size_range = snakemake.params.intron_size_range
coverage_per_cell = snakemake.params.coverage_per_cell
read_length = snakemake.params.read_length
cb_length = snakemake.params.cb_length
umi_length = snakemake.params.umi_length
seed = snakemake.params.seed

random.seed(seed)

# ---------------------------------------------------------------------------
# 1. Generate cell barcodes
# ---------------------------------------------------------------------------
BASES = "ACGT"


def generate_unique_barcodes(n, length):
    """Generate n unique random DNA barcodes of given length."""
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
# 2. Build synthetic genome
# ---------------------------------------------------------------------------
CHROM = "chr_scrna"
SPACER_LEN = 1000


def random_seq(length):
    return "".join(random.choice(BASES) for _ in range(length))


genes = []
genome_parts = []
current_pos = 0  # 0-based position on chromosome

for g in range(n_genes):
    gene_id = f"scrna_gene_{g + 1:03d}"
    n_exons = random.choice(exon_counts)
    strand = random.choice(["+", "-"])

    exon_sizes = [random.randint(*exon_size_range) for _ in range(n_exons)]
    intron_sizes = [random.randint(*intron_size_range) for _ in range(n_exons - 1)]

    # Gene start in 1-based
    gene_start_1 = current_pos + 1

    # Build gene sequence
    exon_coords = []  # list of (start_1based, end_1based)
    gene_seq_parts = []
    pos = current_pos

    for j in range(n_exons):
        ex_start = pos
        ex_end = pos + exon_sizes[j] - 1
        exon_coords.append((ex_start + 1, ex_end + 1))  # 1-based inclusive
        gene_seq_parts.append(random_seq(exon_sizes[j]))
        pos = ex_end + 1

        if j < n_exons - 1:
            # Intron: GT...AG canonical splice sites
            isize = intron_sizes[j]
            intron_seq = "GT" + random_seq(max(0, isize - 4)) + "AG"
            gene_seq_parts.append(intron_seq)
            pos += isize

    gene_end_1 = pos  # 1-based inclusive
    gene_seq = "".join(gene_seq_parts)

    # Compute junction coordinates (1-based inclusive intron boundaries)
    junctions = []
    for j in range(n_exons - 1):
        # intron_start = exon_end + 1 (1-based)
        # intron_end = next_exon_start - 1 (1-based)
        intron_start = exon_coords[j][1] + 1
        intron_end = exon_coords[j + 1][0] - 1
        junctions.append((intron_start, intron_end))

    # Build mRNA (spliced transcript)
    mrna_parts = []
    for j in range(n_exons):
        es0 = exon_coords[j][0] - 1 - current_pos  # offset in gene_seq would be complex
        mrna_parts.append(gene_seq_parts[j * 2] if j * 2 < len(gene_seq_parts) else "")
    # Actually: gene_seq_parts has exon, intron, exon, intron, ..., exon
    mrna = "".join(gene_seq_parts[i] for i in range(0, len(gene_seq_parts), 2))

    genes.append({
        "gene_id": gene_id,
        "strand": strand,
        "n_exons": n_exons,
        "exon_sizes": exon_sizes,
        "intron_sizes": intron_sizes,
        "exon_coords": exon_coords,
        "junctions": junctions,
        "gene_start_1": gene_start_1,
        "gene_end_1": gene_end_1,
        "gene_offset_0": current_pos,
        "mrna": mrna,
        "exon_seqs": [gene_seq_parts[i] for i in range(0, len(gene_seq_parts), 2)],
    })

    genome_parts.append(gene_seq)
    current_pos = pos

    # Add spacer between genes
    if g < n_genes - 1:
        genome_parts.append(random_seq(SPACER_LEN))
        current_pos += SPACER_LEN

genome_seq = "".join(genome_parts)
genome_len = len(genome_seq)

# Write genome FASTA
os.makedirs(os.path.dirname(out_genome), exist_ok=True)
with open(out_genome, "w") as fh:
    fh.write(f">{CHROM}\n")
    for i in range(0, genome_len, 80):
        fh.write(genome_seq[i : i + 80] + "\n")

# Write FAI
header_len = len(f">{CHROM}\n")
with open(out_fai, "w") as fh:
    fh.write(f"{CHROM}\t{genome_len}\t{header_len}\t80\t81\n")

# Write GTF
os.makedirs(os.path.dirname(out_gtf), exist_ok=True)
with open(out_gtf, "w") as fh:
    for gene in genes:
        gid = gene["gene_id"]
        strand = gene["strand"]
        attr_gene = f'gene_id "{gid}";'
        attr_tx = f'gene_id "{gid}"; transcript_id "{gid}_tx";'
        # Gene
        fh.write(
            f"{CHROM}\tsynthetic\tgene\t{gene['gene_start_1']}\t{gene['gene_end_1']}"
            f"\t.\t{strand}\t.\t{attr_gene}\n"
        )
        # Transcript
        fh.write(
            f"{CHROM}\tsynthetic\ttranscript\t{gene['gene_start_1']}\t{gene['gene_end_1']}"
            f"\t.\t{strand}\t.\t{attr_tx}\n"
        )
        # Exons
        for ex_start, ex_end in gene["exon_coords"]:
            fh.write(
                f"{CHROM}\tsynthetic\texon\t{ex_start}\t{ex_end}"
                f"\t.\t{strand}\t.\t{attr_tx}\n"
            )

# Write gene metadata
os.makedirs(os.path.dirname(out_gene_metadata), exist_ok=True)
with open(out_gene_metadata, "w") as fh:
    fh.write(
        "gene_id\tstrand\tn_exons\texon_sizes\tintron_sizes\t"
        "gene_start\tgene_end\tjunctions\n"
    )
    for gene in genes:
        fh.write(
            f"{gene['gene_id']}\t{gene['strand']}\t{gene['n_exons']}\t"
            f"{','.join(map(str, gene['exon_sizes']))}\t"
            f"{','.join(map(str, gene['intron_sizes'])) if gene['intron_sizes'] else 'NA'}\t"
            f"{gene['gene_start_1']}\t{gene['gene_end_1']}\t"
            f"{';'.join(f'{CHROM}:{s}-{e}' for s, e in gene['junctions'])}\n"
        )

# ---------------------------------------------------------------------------
# 3. Generate FASTQ reads with cell barcodes + UMIs
# ---------------------------------------------------------------------------

def rev_comp(seq):
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def make_qual(length):
    return "I" * length


def mrna_to_genomic_read(gene, mrna_start, mrna_end):
    """
    Given 0-based mRNA coords, extract the cDNA read sequence from the mRNA.
    Returns the read sequence (already spliced/mRNA-derived).
    """
    return gene["mrna"][mrna_start:mrna_end]


def generate_umi(length):
    return "".join(random.choice(BASES) for _ in range(length))


reads_r1_count = 0
read_id_counter = 0

# Ground truth: dict of (barcode, junction_str) -> count
junction_counts = {}
# Also track UMI for deduplication ground truth
# (barcode, junction_str, umi) -> number of reads with this exact combo
umi_tracker = {}

os.makedirs(os.path.dirname(out_r1), exist_ok=True)

with gzip.open(out_r1, "wt") as fh_r1, gzip.open(out_r2, "wt") as fh_r2:
    for cell_idx, barcode in enumerate(cell_barcodes):
        for gene in genes:
            if not gene["junctions"]:
                continue

            mrna_len = len(gene["mrna"])

            for junc_idx, (junc_start, junc_end) in enumerate(gene["junctions"]):
                # Random coverage for this cell × junction
                cov = random.choice(coverage_per_cell)
                if cov == 0:
                    continue

                # Junction position in mRNA coordinates
                # Exon boundaries in mRNA space: cumulative exon sizes
                mrna_junc_pos = sum(gene["exon_sizes"][: junc_idx + 1])

                junction_str = f"{CHROM}:{junc_start}-{junc_end}"

                # Generate unique UMIs for this cell × junction
                umis = [generate_umi(umi_length) for _ in range(cov)]

                for umi in umis:
                    read_id_counter += 1
                    rid = f"scrna_{read_id_counter:08d}"

                    # Place cDNA read to span the junction
                    # Read should overlap both sides of the mRNA junction point
                    # mrna_junc_pos is the first base of the next exon in mRNA space
                    min_anchor = 8
                    r2_start = max(0, mrna_junc_pos - read_length + min_anchor)
                    r2_start = min(r2_start, max(0, mrna_junc_pos - min_anchor))
                    # Randomize position within valid range
                    valid_start_min = max(0, mrna_junc_pos - read_length + min_anchor)
                    valid_start_max = min(mrna_len - read_length, mrna_junc_pos - min_anchor)

                    if valid_start_min > valid_start_max:
                        continue

                    r2_start = random.randint(valid_start_min, valid_start_max)
                    r2_end = r2_start + read_length

                    cDNA_seq = mrna_to_genomic_read(gene, r2_start, r2_end)
                    if len(cDNA_seq) < read_length:
                        continue

                    # R1: CB + UMI (padded to 28bp for v3, rest is polyT/adapter)
                    r1_seq = barcode + umi
                    r1_len = len(r1_seq)

                    # R2: cDNA read (reverse complement as per 10x convention)
                    r2_seq = rev_comp(cDNA_seq)

                    r1_qual = make_qual(r1_len)
                    r2_qual = make_qual(read_length)
                    fh_r1.write(f"@{rid}/1\n{r1_seq}\n+\n{r1_qual}\n")
                    fh_r2.write(f"@{rid}/2\n{r2_seq}\n+\n{r2_qual}\n")
                    reads_r1_count += 1

                    # Track ground truth (after UMI dedup, each unique UMI = 1 count)
                    key = (barcode, junction_str)
                    umi_key = (barcode, junction_str, umi)
                    if umi_key not in umi_tracker:
                        umi_tracker[umi_key] = 0
                        # First time seeing this UMI for this cell×junction
                        junction_counts[key] = junction_counts.get(key, 0) + 1
                    umi_tracker[umi_key] += 1

                # Also generate duplicate UMI reads (PCR duplicates) for some junctions
                # to test UMI deduplication
                if cov >= 3:
                    dup_umi = umis[0]  # reuse first UMI
                    for _ in range(2):  # 2 PCR duplicates
                        read_id_counter += 1
                        rid = f"scrna_{read_id_counter:08d}"

                        valid_start_min = max(0, mrna_junc_pos - read_length + min_anchor)
                        valid_start_max = min(mrna_len - read_length, mrna_junc_pos - min_anchor)
                        if valid_start_min > valid_start_max:
                            continue

                        r2_start = random.randint(valid_start_min, valid_start_max)
                        r2_end = r2_start + read_length

                        cDNA_seq = mrna_to_genomic_read(gene, r2_start, r2_end)
                        if len(cDNA_seq) < read_length:
                            continue

                        r1_seq = barcode + dup_umi
                        r2_seq = rev_comp(cDNA_seq)

                        r1_qual = make_qual(len(r1_seq))
                        r2_qual = make_qual(read_length)
                        fh_r1.write(f"@{rid}/1\n{r1_seq}\n+\n{r1_qual}\n")
                        fh_r2.write(f"@{rid}/2\n{r2_seq}\n+\n{r2_qual}\n")
                        reads_r1_count += 1

                        # These are PCR duplicates — same UMI, should NOT increase count
                        umi_key = (barcode, junction_str, dup_umi)
                        umi_tracker[umi_key] += 1

print(f"Generated {reads_r1_count} read pairs for {n_cells} cells × {n_genes} genes")
print(f"Unique cell×junction entries with count > 0: {len(junction_counts)}")

# ---------------------------------------------------------------------------
# 4. Write ground truth junction count matrix
# ---------------------------------------------------------------------------
# Format: barcode\tjunction\tcount
# (sparse representation — only non-zero entries)

os.makedirs(os.path.dirname(out_expected_junctions), exist_ok=True)

# Collect all junctions
all_junctions = set()
for gene in genes:
    for junc_start, junc_end in gene["junctions"]:
        all_junctions.add(f"{CHROM}:{junc_start}-{junc_end}")
all_junctions = sorted(all_junctions)

with open(out_expected_junctions, "w") as fh:
    fh.write("barcode\tjunction\tcount\n")
    for barcode in cell_barcodes:
        for junc in all_junctions:
            count = junction_counts.get((barcode, junc), 0)
            if count > 0:
                fh.write(f"{barcode}\t{junc}\t{count}\n")

print(f"Ground truth written: {sum(junction_counts.values())} total UMI-deduplicated counts")
