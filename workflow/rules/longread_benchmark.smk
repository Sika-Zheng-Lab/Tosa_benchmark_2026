## =============================================================================
## Long-read benchmark rules
## Tosa (1t, 8t), regtools (1t), featureCounts (1t, 8t)
## =============================================================================


# ---------------------------------------------------------------------------
# Read count per BAM (samtools flagstat)
# ---------------------------------------------------------------------------

rule count_reads:
    """Count total reads in BAM file using samtools flagstat."""
    input:
        bam=get_bam,
    output:
        read_count="read_counts/{sample}.txt",
    log:
        "log/count_reads/{sample}.log",
    container:
        config["containers"]["samtools"]
    shell:
        """
        samtools flagstat {input.bam} 2> {log} | head -1 | awk '{{print $1}}' > {output.read_count}
        """


# ---------------------------------------------------------------------------
# Tosa — 1 thread
# ---------------------------------------------------------------------------

rule tosa_1t:
    """Count junction reads with Tosa (1 thread)."""
    input:
        bam=get_bam,
        bai=get_bai,
    output:
        junction="tosa_1t/{sample}/{sample}_junction.tsv.gz",
    params:
        strand=TOSA_STRAND,
        prefix="tosa_1t/{sample}/{sample}",
    log:
        "log/tosa_1t/{sample}.log",
    benchmark:
        repeat("benchmarks/{sample}/tosa_1t.tsv", config["benchmark_repeat"])
    threads: 1
    container:
        config["containers"]["tosa"]
    shell:
        """
        tosa bulk \
            -m 20 \
            -p 1 \
            {params.strand} \
            {input.bam} \
            {params.prefix} \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Tosa — 8 threads
# ---------------------------------------------------------------------------

rule tosa_8t:
    """Count junction reads with Tosa (8 threads)."""
    input:
        bam=get_bam,
        bai=get_bai,
    output:
        junction="tosa_8t/{sample}/{sample}_junction.tsv.gz",
    params:
        strand=TOSA_STRAND,
        prefix="tosa_8t/{sample}/{sample}",
    log:
        "log/tosa_8t/{sample}.log",
    benchmark:
        repeat("benchmarks/{sample}/tosa_8t.tsv", config["benchmark_repeat"])
    threads: config["tosa_threads"]
    container:
        config["containers"]["tosa"]
    shell:
        """
        tosa bulk \
            -m 20 \
            -p {threads} \
            {params.strand} \
            {input.bam} \
            {params.prefix} \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# regtools — 1 thread (regtools is single-threaded only)
# Uses intron-motif method with genome FASTA for strand detection,
# which works with minimap2 long-read BAMs (no XS tag required).
# ---------------------------------------------------------------------------

rule regtools_1t:
    """Extract junctions with regtools (1 thread).
    minimap2 uses 'ts' tag (not 'XS') for strand — pass -t ts so regtools reads it.
    """
    input:
        bam=get_bam,
        bai=get_bai,
    output:
        bed="regtools_1t/{sample}/{sample}_junctions.bed",
    log:
        "log/regtools_1t/{sample}.log",
    benchmark:
        repeat("benchmarks/{sample}/regtools_1t.tsv", config["benchmark_repeat"])
    threads: 1
    container:
        config["containers"]["regtools"]
    shell:
        """
        regtools junctions extract \
            -s XS \
            -t ts \
            -a 8 \
            -m 20 \
            -M 500000 \
            -o {output.bed} \
            {input.bam} \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# featureCounts — 1 thread
# Uses -L for long-read mode.
# ---------------------------------------------------------------------------

rule featurecounts_1t:
    """Count junction reads with featureCounts (1 thread)."""
    input:
        bam=get_bam,
        bai=get_bai,
    output:
        main="featurecounts_1t/{sample}/{sample}_featurecounts",
        summary="featurecounts_1t/{sample}/{sample}_featurecounts.summary",
        jcounts="featurecounts_1t/{sample}/{sample}_featurecounts.jcounts",
    params:
        strand=FC_STRAND,
        gtf=config["gtf"],
        tmpdir="/tmp/featurecounts_1t_{sample}",
    log:
        "log/featurecounts_1t/{sample}.log",
    benchmark:
        repeat("benchmarks/{sample}/featurecounts_1t.tsv", config["benchmark_repeat"])
    threads: 1
    resources:
        featurecounts_io=1,
    container:
        config["containers"]["featurecounts"]
    shell:
        """
        mkdir -p {params.tmpdir}
        featureCounts \
            -a {params.gtf} \
            -J \
            -L \
            -s {params.strand} \
            -T 1 \
            --tmpDir {params.tmpdir} \
            -o {output.main} \
            {input.bam} \
            > {log} 2>&1
        rm -rf {params.tmpdir}
        """


# ---------------------------------------------------------------------------
# featureCounts — 8 threads
# ---------------------------------------------------------------------------

rule featurecounts_8t:
    """Count junction reads with featureCounts (8 threads)."""
    input:
        bam=get_bam,
        bai=get_bai,
    output:
        main="featurecounts_8t/{sample}/{sample}_featurecounts",
        summary="featurecounts_8t/{sample}/{sample}_featurecounts.summary",
        jcounts="featurecounts_8t/{sample}/{sample}_featurecounts.jcounts",
    params:
        strand=FC_STRAND,
        gtf=config["gtf"],
        tmpdir="/tmp/featurecounts_8t_{sample}",
    log:
        "log/featurecounts_8t/{sample}.log",
    benchmark:
        repeat("benchmarks/{sample}/featurecounts_8t.tsv", config["benchmark_repeat"])
    threads: config["featurecounts_threads"]
    resources:
        featurecounts_io=1,
    container:
        config["containers"]["featurecounts"]
    shell:
        """
        mkdir -p {params.tmpdir}
        featureCounts \
            -a {params.gtf} \
            -J \
            -L \
            -s {params.strand} \
            -T {threads} \
            --tmpDir {params.tmpdir} \
            -o {output.main} \
            {input.bam} \
            > {log} 2>&1
        rm -rf {params.tmpdir}
        """


# ---------------------------------------------------------------------------
# Aggregate benchmark results and generate plots
# ---------------------------------------------------------------------------

rule aggregate_results:
    """Aggregate all benchmark results and generate comparison plots."""
    input:
        benchmarks=expand(
            "benchmarks/{sample}/{tool}.tsv",
            sample=SAMPLES,
            tool=TOOLS,
        ),
        read_counts=expand("read_counts/{sample}.txt", sample=SAMPLES),
    output:
        per_sample="longread_results/per_sample_benchmark.tsv",
        summary="longread_results/summary_total_time.tsv",
        barplot="longread_results/barplot_total_time.png",
        boxplot_time="longread_results/boxplot_time_per_sample.png",
        scatter_time="longread_results/scatter_reads_vs_time.png",
        scatter_memory="longread_results/scatter_reads_vs_memory.png",
        boxplot_memory="longread_results/boxplot_memory.png",
        speedup="longread_results/speedup_comparison.png",
    params:
        samples=SAMPLES,
        tools=TOOLS,
        benchmark_dir="benchmarks",
        read_count_dir="read_counts",
    log:
        "log/aggregate_results.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/aggregate_benchmark_longread.py"
