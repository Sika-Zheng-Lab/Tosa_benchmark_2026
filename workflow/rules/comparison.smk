## =============================================================================
## Comparison / analysis rules
## =============================================================================


rule normalize_junctions:
    """Normalise outputs from STAR, Tosa, regtools, featureCounts into a common format."""
    input:
        star="{sample}/star/{sample}.SJ.out.tab",
        tosa="{sample}/tosa/{sample}_junction.tsv.gz",
        regtools="{sample}/regtools/{sample}_junctions.bed",
        featurecounts="{sample}/featurecounts/{sample}_featurecounts.jcounts",
    output:
        star="{sample}/normalized/star.tsv",
        tosa="{sample}/normalized/tosa.tsv",
        regtools="{sample}/normalized/regtools.tsv",
        featurecounts="{sample}/normalized/featurecounts.tsv",
    log:
        "log/normalize_junctions/{sample}.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/normalize_junctions.py"


rule compare_junctions:
    """Compare junction counts across tools and generate plots."""
    input:
        star="{sample}/normalized/star.tsv",
        tosa="{sample}/normalized/tosa.tsv",
        regtools="{sample}/normalized/regtools.tsv",
        featurecounts="{sample}/normalized/featurecounts.tsv",
        bench_tosa="benchmarks/{sample}/tosa.tsv",
        bench_regtools="benchmarks/{sample}/regtools.tsv",
        bench_featurecounts="benchmarks/{sample}/featurecounts.tsv",
    output:
        counts="comparison/{sample}_junction_counts.tsv",
        scatter="comparison/{sample}_scatter.png",
        venn="comparison/{sample}_venn.png",
        correlation="comparison/{sample}_correlation.tsv",
        jaccard="comparison/{sample}_jaccard.png",
        benchmark="comparison/{sample}_benchmark.png",
    log:
        "log/compare_junctions/{sample}.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/compare_junctions.py"


rule validate_simulation:
    """Validate simulation counts against expected values."""
    input:
        star="simulation/normalized/star.tsv",
        tosa="simulation/normalized/tosa.tsv",
        regtools="simulation/normalized/regtools.tsv",
        featurecounts="simulation/normalized/featurecounts.tsv",
        expected="simulation/expected_counts.tsv",
    output:
        table="comparison/simulation_validation.tsv",
        plot="comparison/simulation_validation.png",
    log:
        "log/validate_simulation/simulation.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/validate_simulation.py"


rule investigate_discrepancies:
    """Investigate per-junction count discrepancies between tools."""
    input:
        counts="comparison/{sample}_junction_counts.tsv",
        bam="{sample}/star/{sample}.Aligned.sortedByCoord.out.bam",
        bai="{sample}/star/{sample}.Aligned.sortedByCoord.out.bam.bai",
        sj="{sample}/star/{sample}.SJ.out.tab",
    output:
        table="comparison/{sample}_discrepancy_table.tsv",
        plot="comparison/{sample}_discrepancy_summary.png",
    params:
        is_paired=lambda wc: is_paired(wc.sample),
    log:
        "log/investigate_discrepancies/{sample}.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/investigate_discrepancies.py"


rule compare_benchmark_mt:
    """Compare multi-threaded benchmark results (Tosa 8t, featureCounts 8t, regtools 1t)."""
    input:
        bench_tosa_mt="benchmarks/{sample}/tosa_mt.tsv",
        bench_featurecounts_mt="benchmarks/{sample}/featurecounts_mt.tsv",
        bench_regtools="benchmarks/{sample}/regtools.tsv",
    output:
        benchmark="comparison/{sample}_benchmark_mt.png",
    log:
        "log/compare_benchmark_mt/{sample}.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/compare_benchmark_mt.py"
