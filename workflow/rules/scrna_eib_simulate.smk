## =============================================================================
## scRNA-seq EIB Benchmark — simulation rules
## =============================================================================
## Reuses the bulk EIB genome (eib_benchmark/genome/) and adds scRNA-seq
## specific read generation with cell barcodes + UMIs.
## =============================================================================

SCRNA_EIB_CFG = config.get("scrna_eib_benchmark", {})
SCRNA_CFG = config.get("scrna_benchmark", {})
SCRNA_EIB_ANCHOR_LENGTHS = [str(a) for a in SCRNA_EIB_CFG.get("boundary_anchor_lengths", [1, 4, 8])]
SCRNA_EIB_THREADS = [str(t) for t in SCRNA_EIB_CFG.get("threads", [1])]
SCRNA_EIB_STRANDEDNESS = SCRNA_EIB_CFG.get("strandedness", ["XS"])


rule scrna_eib_generate_reads:
    """Generate scRNA-seq FASTQ reads with cell barcodes for EIB benchmark."""
    input:
        genome="eib_benchmark/genome/genome.fa",
        gtf="eib_benchmark/genome/annotation.gtf",
        metadata="eib_benchmark/genome/gene_metadata.tsv",
    output:
        r1="scrna_eib_benchmark/fastq/scrna_eib_R1.fastq.gz",
        r2="scrna_eib_benchmark/fastq/scrna_eib_R2.fastq.gz",
        whitelist="scrna_eib_benchmark/reference/whitelist.txt",
        read_positions="scrna_eib_benchmark/read_positions.tsv",
        expected_boundaries="scrna_eib_benchmark/expected_boundary_counts.tsv",
    params:
        n_cells=SCRNA_EIB_CFG.get("n_cells", 30),
        coverage_per_cell=SCRNA_EIB_CFG.get("coverage_per_cell", [1, 5, 20]),
        read_length=SCRNA_EIB_CFG.get("read_length", 100),
        fragment_length=SCRNA_EIB_CFG.get("fragment_length", 300),
        fragment_length_sd=SCRNA_EIB_CFG.get("fragment_length_sd", 50),
        cb_length=SCRNA_CFG.get("cb_length", 16),
        umi_length=SCRNA_CFG.get("umi_length", 12),
        seed=SCRNA_EIB_CFG.get("seed", 42),
    log:
        "log/scrna_eib/generate_reads.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/generate_scrna_eib_reads.py"


rule scrna_eib_starsolo_align:
    """Align scRNA-seq EIB reads with STARsolo."""
    input:
        r1="scrna_eib_benchmark/fastq/scrna_eib_R1.fastq.gz",
        r2="scrna_eib_benchmark/fastq/scrna_eib_R2.fastq.gz",
        index="eib_benchmark/genome/star_index",
        whitelist="scrna_eib_benchmark/reference/whitelist.txt",
    output:
        bam="scrna_eib_benchmark/starsolo/scrna_eib.Aligned.sortedByCoord.out.bam",
        sj="scrna_eib_benchmark/starsolo/scrna_eib.SJ.out.tab",
    params:
        prefix="scrna_eib_benchmark/starsolo/scrna_eib.",
        cb_len=SCRNA_CFG.get("cb_length", 16),
        umi_start=SCRNA_CFG.get("cb_length", 16) + 1,
        umi_len=SCRNA_CFG.get("umi_length", 12),
    threads:
        SCRNA_EIB_CFG.get("star_threads", EIB_CFG.get("star_threads", 32))
    log:
        "log/scrna_eib/starsolo_align.log",
    benchmark:
        "benchmarks/scrna_eib/starsolo_align.tsv"
    container:
        config["containers"]["star"]
    shell:
        """
        ulimit -n 65535 2>/dev/null || true
        STAR \
            --runMode alignReads \
            --genomeDir {input.index} \
            --readFilesIn {input.r2} {input.r1} \
            --readFilesCommand zcat \
            --soloType CB_UMI_Simple \
            --soloCBwhitelist {input.whitelist} \
            --soloCBstart 1 \
            --soloCBlen {params.cb_len} \
            --soloUMIstart {params.umi_start} \
            --soloUMIlen {params.umi_len} \
            --soloFeatures Gene SJ \
            --soloCellFilter None \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS NM XS CB UB \
            --outFilterMultimapNmax 1 \
            --outBAMsortingBinsN 50 \
            --alignIntronMin 20 \
            --outFileNamePrefix {params.prefix} \
            --runThreadN {threads} \
            > {log} 2>&1
        """


rule scrna_eib_samtools_index:
    """Index scRNA-seq EIB BAM file."""
    input:
        bam="scrna_eib_benchmark/starsolo/scrna_eib.Aligned.sortedByCoord.out.bam",
    output:
        bai="scrna_eib_benchmark/starsolo/scrna_eib.Aligned.sortedByCoord.out.bam.bai",
    log:
        "log/scrna_eib/samtools_index.log",
    container:
        config["containers"]["samtools"]
    shell:
        """
        samtools index {input.bam} > {log} 2>&1
        """
