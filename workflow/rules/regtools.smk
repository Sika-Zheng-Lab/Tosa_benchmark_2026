## =============================================================================
## regtools junction extraction rules
## =============================================================================


rule regtools_extract:
    """Extract junctions with regtools."""
    input:
        bam="{sample}/star/{sample}.Aligned.sortedByCoord.out.bam",
        bai="{sample}/star/{sample}.Aligned.sortedByCoord.out.bam.bai",
    output:
        bed="{sample}/regtools/{sample}_junctions.bed",
    params:
        strand=get_strandedness_regtools,
    log:
        "log/regtools_extract/{sample}.log",
    benchmark:
        repeat("benchmarks/{sample}/regtools.tsv", config["benchmark_repeat"])
    container:
        config["containers"]["regtools"]
    shell:
        """
        regtools junctions extract \
            -s {params.strand} \
            -a 8 \
            -m 20 \
            -M 500000 \
            -o {output.bed} \
            {input.bam} \
            > {log} 2>&1
        """
