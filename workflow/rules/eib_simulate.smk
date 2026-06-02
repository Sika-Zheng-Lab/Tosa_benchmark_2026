## =============================================================================
## EIB Benchmark — synthetic genome & read generation rules
## =============================================================================

import math

EIB_CFG = config.get("eib_benchmark", {})
EIB_LIBRARY_TYPES = EIB_CFG.get("library_types", ["PE"])
EIB_STRANDEDNESS = EIB_CFG.get("strandedness", ["XS"])
EIB_THREADS = [str(t) for t in EIB_CFG.get("threads", [1])]
EIB_ANCHOR_LENGTHS = [str(a) for a in EIB_CFG.get("boundary_anchor_lengths", [1])]


rule eib_generate_genome:
    """Generate synthetic genome, GTF, and gene metadata for EIB benchmark."""
    output:
        genome="eib_benchmark/genome/genome.fa",
        fai="eib_benchmark/genome/genome.fa.fai",
        gtf="eib_benchmark/genome/annotation.gtf",
        metadata="eib_benchmark/genome/gene_metadata.tsv",
    params:
        n_genes=EIB_CFG.get("n_genes", 200),
        exon_counts=EIB_CFG.get("exon_counts", [2, 3, 5, 10]),
        exon_size_range=EIB_CFG.get("exon_size_range", [100, 2000]),
        intron_size_range=EIB_CFG.get("intron_size_range", [100, 100000]),
        coverage_levels=EIB_CFG.get("coverage_levels", [5, 20, 100]),
        seed=EIB_CFG.get("seed", 42),
    log:
        "log/eib/generate_genome.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/generate_eib_genome.py"


rule eib_star_index:
    """Build STAR genome index for the EIB synthetic genome."""
    input:
        genome="eib_benchmark/genome/genome.fa",
        gtf="eib_benchmark/genome/annotation.gtf",
    output:
        directory("eib_benchmark/genome/star_index"),
    params:
        # Compute appropriate --genomeSAindexNbases for synthetic genome
        # Actual size computed at run time; use a safe small value
        sa_bases=6,
    threads: 1
    log:
        "log/eib/star_index.log",
    container:
        config["containers"]["star"]
    shell:
        """
        GENOME_LEN=$(awk '{{sum += $2}} END {{print sum}}' {input.genome}.fai)
        SA_BASES=$(awk -v gl="$GENOME_LEN" 'BEGIN {{v=log(gl)/log(2)/2-1; v=int(v); if(v>14) v=14; if(v<1) v=1; print v}}')
        mkdir -p {output}
        STAR \
            --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.gtf} \
            --genomeSAindexNbases $SA_BASES \
            --runThreadN {threads} \
            > {log} 2>&1
        """


rule eib_generate_reads_pe:
    """Generate paired-end FASTQ reads for EIB benchmark (shared across strandedness and anchor lengths)."""
    input:
        genome="eib_benchmark/genome/genome.fa",
        gtf="eib_benchmark/genome/annotation.gtf",
        metadata="eib_benchmark/genome/gene_metadata.tsv",
    output:
        r1="eib_benchmark/PE/fastq/eib_R1.fastq.gz",
        r2="eib_benchmark/PE/fastq/eib_R2.fastq.gz",
        read_positions="eib_benchmark/PE/read_positions.tsv",
        read_metadata="eib_benchmark/PE/read_metadata.tsv",
    params:
        library_type="PE",
        read_length=EIB_CFG.get("read_length", 100),
        fragment_length=EIB_CFG.get("fragment_length", 300),
        fragment_length_sd=EIB_CFG.get("fragment_length_sd", 50),
        seed=EIB_CFG.get("seed", 42),
    log:
        "log/eib/generate_reads_PE.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/generate_eib_reads.py"


rule eib_generate_reads_se:
    """Generate single-end FASTQ reads for EIB benchmark (shared across strandedness and anchor lengths)."""
    input:
        genome="eib_benchmark/genome/genome.fa",
        gtf="eib_benchmark/genome/annotation.gtf",
        metadata="eib_benchmark/genome/gene_metadata.tsv",
    output:
        reads="eib_benchmark/SE/fastq/eib_reads.fastq.gz",
        read_positions="eib_benchmark/SE/read_positions.tsv",
        read_metadata="eib_benchmark/SE/read_metadata.tsv",
    params:
        library_type="SE",
        read_length=EIB_CFG.get("read_length", 100),
        fragment_length=EIB_CFG.get("fragment_length", 300),
        fragment_length_sd=EIB_CFG.get("fragment_length_sd", 50),
        seed=EIB_CFG.get("seed", 42),
    log:
        "log/eib/generate_reads_SE.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/generate_eib_reads.py"


rule eib_compute_expected_counts:
    """Compute expected boundary counts for a given anchor length (lightweight, no FASTQ regen)."""
    input:
        genome="eib_benchmark/genome/genome.fa",
        metadata="eib_benchmark/genome/gene_metadata.tsv",
        read_positions="eib_benchmark/{library_type}/read_positions.tsv",
    output:
        expected="eib_benchmark/anchor_{anchor_length}/{library_type}/expected_boundary_counts.tsv",
    params:
        boundary_anchor_length=lambda wildcards: int(wildcards.anchor_length),
    log:
        "log/eib/compute_expected_{library_type}_anchor{anchor_length}.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/compute_eib_expected_counts.py"
