"""
Generate a synthetic genome and GTF annotation for the EIB benchmark.

Creates a multi-gene synthetic genome with diverse exon structures:
- Variable exon counts (2, 3, 5, 10)
- Variable exon sizes (100–2000 bp)
- Variable intron sizes (100–100,000 bp) with canonical GT-AG splice sites
- Genes on both + and - strands
- Coverage levels assigned per gene for downstream read generation

Outputs:
  genome.fa      — FASTA genome
  genome.fa.fai  — FASTA index
  annotation.gtf — GTF with gene/transcript/exon features
  gene_metadata.tsv — per-gene metadata for stratified analysis
"""

import math
import os
import random

# ---------------------------------------------------------------------------
# Snakemake variables
# ---------------------------------------------------------------------------
n_genes = snakemake.params.n_genes
exon_counts = snakemake.params.exon_counts
exon_size_range = snakemake.params.exon_size_range
intron_size_range = snakemake.params.intron_size_range
coverage_levels = snakemake.params.coverage_levels
seed = snakemake.params.seed

out_genome = snakemake.output.genome
out_fai = snakemake.output.fai
out_gtf = snakemake.output.gtf
out_metadata = snakemake.output.metadata

random.seed(seed)

# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def random_seq(length):
    """Generate a random DNA sequence."""
    return "".join(random.choice("ACGT") for _ in range(length))


# ---------------------------------------------------------------------------
# 1. Design gene structures
# ---------------------------------------------------------------------------

CHROM = "chr_eib"
INTER_GENE_GAP = 1000  # bp between genes

genes = []
current_pos = 1  # 1-based genomic coordinate

for i in range(n_genes):
    gene_id = f"eib_gene_{i + 1:04d}"
    tx_id = f"eib_tx_{i + 1:04d}"
    n_exons = random.choice(exon_counts)
    strand = random.choice(["+", "-"])
    coverage = random.choice(coverage_levels)

    # Generate exon and intron sizes
    exon_sizes = [random.randint(exon_size_range[0], exon_size_range[1]) for _ in range(n_exons)]
    intron_sizes = [random.randint(intron_size_range[0], intron_size_range[1]) for _ in range(n_exons - 1)]

    # Compute exon coordinates (1-based inclusive)
    exons = []
    pos = current_pos
    for j in range(n_exons):
        exon_start = pos
        exon_end = pos + exon_sizes[j] - 1
        exons.append((exon_start, exon_end))
        pos = exon_end + 1
        if j < n_exons - 1:
            pos += intron_sizes[j]  # skip intron

    gene_start = exons[0][0]
    gene_end = exons[-1][1]

    genes.append({
        "gene_id": gene_id,
        "tx_id": tx_id,
        "strand": strand,
        "n_exons": n_exons,
        "exon_sizes": exon_sizes,
        "intron_sizes": intron_sizes,
        "exons": exons,
        "gene_start": gene_start,
        "gene_end": gene_end,
        "coverage": coverage,
    })

    current_pos = gene_end + INTER_GENE_GAP + 1

genome_length = current_pos - INTER_GENE_GAP  # remove trailing gap

# ---------------------------------------------------------------------------
# 2. Build genome FASTA
# ---------------------------------------------------------------------------
# Build full sequence: exons = random, introns = GT..AG + random middle

genome_parts = []
prev_end = 0  # 0-based end of previous feature

for gene in genes:
    # Fill gap before gene with random sequence
    gap_before = gene["gene_start"] - 1 - prev_end
    if gap_before > 0:
        genome_parts.append(random_seq(gap_before))

    for j, (ex_start, ex_end) in enumerate(gene["exons"]):
        exon_len = ex_end - ex_start + 1
        genome_parts.append(random_seq(exon_len))

        if j < len(gene["intron_sizes"]):
            intron_len = gene["intron_sizes"][j]
            if intron_len >= 4:
                intron_seq = "GT" + random_seq(intron_len - 4) + "AG"
            else:
                intron_seq = random_seq(intron_len)
            genome_parts.append(intron_seq)

    prev_end = gene["gene_end"]

genome_seq = "".join(genome_parts)
actual_genome_length = len(genome_seq)

# Write FASTA
os.makedirs(os.path.dirname(out_genome), exist_ok=True)

with open(out_genome, "w") as fh:
    fh.write(f">{CHROM}\n")
    for i in range(0, actual_genome_length, 80):
        fh.write(genome_seq[i : i + 80] + "\n")

# Write FAI
header_len = len(f">{CHROM}\n")
with open(out_fai, "w") as fh:
    fh.write(f"{CHROM}\t{actual_genome_length}\t{header_len}\t80\t81\n")

# ---------------------------------------------------------------------------
# 3. Write GTF annotation
# ---------------------------------------------------------------------------
os.makedirs(os.path.dirname(out_gtf), exist_ok=True)

with open(out_gtf, "w") as fh:
    for gene in genes:
        gid = gene["gene_id"]
        tid = gene["tx_id"]
        strand = gene["strand"]
        attr_gene = f'gene_id "{gid}";'
        attr_tx = f'gene_id "{gid}"; transcript_id "{tid}";'

        # Gene line
        fh.write(f"{CHROM}\tsynthetic\tgene\t{gene['gene_start']}\t{gene['gene_end']}\t.\t{strand}\t.\t{attr_gene}\n")
        # Transcript line
        fh.write(f"{CHROM}\tsynthetic\ttranscript\t{gene['gene_start']}\t{gene['gene_end']}\t.\t{strand}\t.\t{attr_tx}\n")
        # Exon lines
        for j, (ex_start, ex_end) in enumerate(gene["exons"]):
            exon_attr = f'gene_id "{gid}"; transcript_id "{tid}"; exon_number "{j + 1}";'
            fh.write(f"{CHROM}\tsynthetic\texon\t{ex_start}\t{ex_end}\t.\t{strand}\t.\t{exon_attr}\n")

# ---------------------------------------------------------------------------
# 4. Write gene metadata
# ---------------------------------------------------------------------------
os.makedirs(os.path.dirname(out_metadata), exist_ok=True)

with open(out_metadata, "w") as fh:
    fh.write("gene_id\ttranscript_id\tstrand\tn_exons\texon_sizes\tintron_sizes\tgene_start\tgene_end\tcoverage\n")
    for gene in genes:
        exon_sizes_str = ",".join(str(s) for s in gene["exon_sizes"])
        intron_sizes_str = ",".join(str(s) for s in gene["intron_sizes"]) if gene["intron_sizes"] else "NA"
        fh.write(
            f"{gene['gene_id']}\t{gene['tx_id']}\t{gene['strand']}\t{gene['n_exons']}\t"
            f"{exon_sizes_str}\t{intron_sizes_str}\t{gene['gene_start']}\t{gene['gene_end']}\t"
            f"{gene['coverage']}\n"
        )

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
print(f"Generated synthetic genome: {actual_genome_length:,} bp on {CHROM}")
print(f"Number of genes: {n_genes}")
print(f"Genome FASTA index Nbases suggestion: {min(14, int(math.log2(actual_genome_length) / 2 - 1))}")
