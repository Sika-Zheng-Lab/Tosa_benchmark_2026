## =============================================================================
## STAR alignment rules
## =============================================================================


def get_star_input_fastqs(wildcards):
    """Return list of FASTQ files for a sample."""
    s = wildcards.sample
    if s == "simulation":
        fqs = ["simulation/fastq/sim_R1.fastq.gz",
               "simulation/fastq/sim_R2.fastq.gz"]
    else:
        fqs = [config["samples"][s]["R1"]]
        if is_paired(s):
            fqs.append(config["samples"][s]["R2"])
    return fqs


def get_star_readfiles_param(wildcards):
    """Build --readFilesIn parameter string."""
    s = wildcards.sample
    if s == "simulation":
        return "simulation/fastq/sim_R1.fastq.gz simulation/fastq/sim_R2.fastq.gz"
    r1 = config["samples"][s]["R1"]
    if is_paired(s):
        return f"{r1} {config['samples'][s]['R2']}"
    return r1


def get_star_index_input(wildcards):
    """Return STAR index dir (pre-built or simulation-built)."""
    return get_star_index(wildcards.sample)


rule star_align:
    """Align FASTQ reads with STAR."""
    input:
        fastqs=get_star_input_fastqs,
        index=get_star_index_input,
    output:
        bam="{sample}/star/{sample}.Aligned.sortedByCoord.out.bam",
        sj="{sample}/star/{sample}.SJ.out.tab",
    params:
        readfiles=get_star_readfiles_param,
        prefix="{sample}/star/{sample}.",
    threads:
        config["star_threads"]
    log:
        "log/star_align/{sample}.log",
    benchmark:
        "benchmarks/{sample}/star_align.tsv"
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


rule samtools_index:
    """Index a BAM file."""
    input:
        bam="{sample}/star/{sample}.Aligned.sortedByCoord.out.bam",
    output:
        bai="{sample}/star/{sample}.Aligned.sortedByCoord.out.bam.bai",
    log:
        "log/samtools_index/{sample}.log",
    container:
        config["containers"]["samtools"]
    shell:
        """
        samtools index {input.bam} > {log} 2>&1
        """
