## =============================================================================
## EIB Benchmark — STAR alignment rules
## =============================================================================


def _eib_star_input_fastqs(wildcards):
    """Return FASTQ file list for EIB STAR alignment."""
    lt = wildcards.library_type
    base = f"eib_benchmark/{lt}/fastq"
    if lt == "PE":
        return [f"{base}/eib_R1.fastq.gz", f"{base}/eib_R2.fastq.gz"]
    else:
        return [f"{base}/eib_reads.fastq.gz"]


def _eib_star_readfiles_param(wildcards):
    """Build --readFilesIn parameter string for EIB STAR alignment."""
    lt = wildcards.library_type
    base = f"eib_benchmark/{lt}/fastq"
    if lt == "PE":
        return f"{base}/eib_R1.fastq.gz {base}/eib_R2.fastq.gz"
    else:
        return f"{base}/eib_reads.fastq.gz"


rule eib_star_align:
    """Align EIB synthetic FASTQ reads with STAR (one BAM per library_type)."""
    input:
        fastqs=_eib_star_input_fastqs,
        index="eib_benchmark/genome/star_index",
    output:
        bam="eib_benchmark/{library_type}/star/eib.Aligned.sortedByCoord.out.bam",
        sj="eib_benchmark/{library_type}/star/eib.SJ.out.tab",
    params:
        readfiles=_eib_star_readfiles_param,
        prefix="eib_benchmark/{library_type}/star/eib.",
    threads:
        EIB_CFG.get("star_threads", 32)
    log:
        "log/eib/star_align_{library_type}.log",
    benchmark:
        "benchmarks/eib/{library_type}/star_align.tsv"
    container:
        config["containers"]["star"]
    shell:
        """
        ulimit -n 65535 2>/dev/null || true
        STAR \
            --runMode alignReads \
            --genomeDir {input.index} \
            --readFilesIn {params.readfiles} \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS NM XS \
            --outFilterMultimapNmax 1 \
            --outBAMsortingBinsN 50 \
            --alignIntronMin 20 \
            --outFileNamePrefix {params.prefix} \
            --runThreadN {threads} \
            > {log} 2>&1
        """


rule eib_samtools_index:
    """Index EIB BAM file."""
    input:
        bam="eib_benchmark/{library_type}/star/eib.Aligned.sortedByCoord.out.bam",
    output:
        bai="eib_benchmark/{library_type}/star/eib.Aligned.sortedByCoord.out.bam.bai",
    log:
        "log/eib/samtools_index_{library_type}.log",
    container:
        config["containers"]["samtools"]
    shell:
        """
        samtools index {input.bam} > {log} 2>&1
        """
