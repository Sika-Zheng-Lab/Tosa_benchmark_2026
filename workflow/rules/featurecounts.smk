## =============================================================================
## featureCounts junction counting rules
## =============================================================================


rule featurecounts_count:
    """Count junction reads with featureCounts (-J option)."""
    input:
        bam="{sample}/star/{sample}.Aligned.sortedByCoord.out.bam",
        bai="{sample}/star/{sample}.Aligned.sortedByCoord.out.bam.bai",
        gtf=lambda wc: get_gtf(wc.sample),
    output:
        main="{sample}/featurecounts/{sample}_featurecounts",
        summary="{sample}/featurecounts/{sample}_featurecounts.summary",
        jcounts="{sample}/featurecounts/{sample}_featurecounts.jcounts",
    params:
        strand=get_strandedness_featurecounts,
        pe_flag=lambda wc: "-p --countReadPairs" if is_paired(wc.sample) else "",
    log:
        "log/featurecounts/{sample}.log",
    benchmark:
        repeat("benchmarks/{sample}/featurecounts.tsv", config["benchmark_repeat"])
    container:
        config["containers"]["featurecounts"]
    shell:
        """
        featureCounts \
            -a {input.gtf} \
            -J \
            -s {params.strand} \
            {params.pe_flag} \
            -T 1 \
            -o {output.main} \
            {input.bam} \
            > {log} 2>&1
        """


rule featurecounts_count_mt:
    """Count junction reads with featureCounts (-J option, multi-threaded)."""
    input:
        bam="{sample}/star/{sample}.Aligned.sortedByCoord.out.bam",
        bai="{sample}/star/{sample}.Aligned.sortedByCoord.out.bam.bai",
        gtf=lambda wc: get_gtf(wc.sample),
    output:
        main="{sample}/featurecounts_mt/{sample}_featurecounts",
        summary="{sample}/featurecounts_mt/{sample}_featurecounts.summary",
        jcounts="{sample}/featurecounts_mt/{sample}_featurecounts.jcounts",
    params:
        strand=get_strandedness_featurecounts,
        pe_flag=lambda wc: "-p --countReadPairs" if is_paired(wc.sample) else "",
    log:
        "log/featurecounts_count_mt/{sample}.log",
    benchmark:
        repeat("benchmarks/{sample}/featurecounts_mt.tsv", config["benchmark_repeat"])
    threads:
        config["featurecounts_threads"]
    container:
        config["containers"]["featurecounts"]
    shell:
        """
        featureCounts \
            -a {input.gtf} \
            -J \
            -s {params.strand} \
            {params.pe_flag} \
            -T {threads} \
            -o {output.main} \
            {input.bam} \
            > {log} 2>&1
        """
