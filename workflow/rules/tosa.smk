## =============================================================================
## Tosa junction counting rules
## =============================================================================


rule tosa_count:
    """Count junction reads with Tosa (bulk mode)."""
    input:
        bam="{sample}/star/{sample}.Aligned.sortedByCoord.out.bam",
        bai="{sample}/star/{sample}.Aligned.sortedByCoord.out.bam.bai",
    output:
        junction="{sample}/tosa/{sample}_junction.tsv.gz",
    params:
        strand=get_strandedness_tosa,
        prefix="{sample}/tosa/{sample}",
    log:
        "log/tosa_count/{sample}.log",
    benchmark:
        repeat("benchmarks/{sample}/tosa.tsv", config["benchmark_repeat"])
    container:
        config["containers"]["tosa"]
    shell:
        """
        tosa bulk \
        -m 20 \
        {params.strand} \
        {input.bam} \
        {params.prefix} \
        > {log} 2>&1
        """


rule tosa_count_mt:
    """Count junction reads with Tosa (bulk mode, multi-threaded)."""
    input:
        bam="{sample}/star/{sample}.Aligned.sortedByCoord.out.bam",
        bai="{sample}/star/{sample}.Aligned.sortedByCoord.out.bam.bai",
    output:
        junction="{sample}/tosa_mt/{sample}_junction.tsv.gz",
    params:
        strand=get_strandedness_tosa,
        prefix="{sample}/tosa_mt/{sample}",
    log:
        "log/tosa_count_mt/{sample}.log",
    benchmark:
        repeat("benchmarks/{sample}/tosa_mt.tsv", config["benchmark_repeat"])
    threads:
        config["tosa_threads"]
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
