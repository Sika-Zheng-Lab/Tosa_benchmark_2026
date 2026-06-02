## =============================================================================
## EIB Benchmark — Tosa junction counting rules (multi-thread sweep)
## =============================================================================


def _get_eib_strandedness_tosa(wildcards):
    """Return Tosa strand flag for the given strandedness wildcard."""
    s = wildcards.strandedness
    if s:
        return f"--strand {s}"
    return ""


rule eib_tosa_count:
    """Count boundary reads with Tosa at variable thread counts and anchor lengths."""
    input:
        bam="eib_benchmark/{library_type}/star/eib.Aligned.sortedByCoord.out.bam",
        bai="eib_benchmark/{library_type}/star/eib.Aligned.sortedByCoord.out.bam.bai",
        gtf="eib_benchmark/genome/annotation.gtf",
    output:
        boundary="eib_benchmark/anchor_{anchor_length}/{library_type}_{strandedness}/tosa/{threads}/eib_boundary.tsv.gz",
        junction="eib_benchmark/anchor_{anchor_length}/{library_type}_{strandedness}/tosa/{threads}/eib_junction.tsv.gz",
    params:
        strand=_get_eib_strandedness_tosa,
        prefix="eib_benchmark/anchor_{anchor_length}/{library_type}_{strandedness}/tosa/{threads}/eib",
    threads:
        lambda wildcards: int(wildcards.threads)
    log:
        "log/eib/tosa_{library_type}_{strandedness}_anchor{anchor_length}_{threads}.log",
    benchmark:
        repeat(
            "benchmarks/eib/anchor_{anchor_length}/{library_type}_{strandedness}/tosa_{threads}.tsv",
            EIB_CFG.get("benchmark_repeat", 10),
        )
    container:
        config["containers"]["tosa"]
    shell:
        """
        tosa bulk \
        -m 20 \
        -p {threads} \
        -b {wildcards.anchor_length} \
        {params.strand} \
        -g {input.gtf} \
        {input.bam} \
        {params.prefix} \
        > {log} 2>&1
        """
