## =============================================================================
## scRNA-seq EIB Benchmark — Tosa single-cell mode with boundary counting
## =============================================================================

SCRNA_EIB_CFG = config.get("scrna_eib_benchmark", {})
SCRNA_CFG = config.get("scrna_benchmark", {})


def _scrna_eib_tosa_strandedness(wildcards):
    """Return Tosa strand flag for the given strandedness wildcard."""
    s = wildcards.strandedness
    if s:
        return f"--strand {s}"
    return ""


rule scrna_eib_tosa_count:
    """Count boundary reads with Tosa single-cell mode at variable thread counts and anchor lengths."""
    input:
        bam="scrna_eib_benchmark/starsolo/scrna_eib.Aligned.sortedByCoord.out.bam",
        bai="scrna_eib_benchmark/starsolo/scrna_eib.Aligned.sortedByCoord.out.bam.bai",
        gtf="eib_benchmark/genome/annotation.gtf",
        whitelist="scrna_eib_benchmark/reference/whitelist.txt",
    output:
        boundary_detail="scrna_eib_benchmark/anchor_{anchor_length}/{strandedness}/tosa/{threads}/scrna_eib_boundary_barcodes_detail.tsv.gz",
        junction_barcodes="scrna_eib_benchmark/anchor_{anchor_length}/{strandedness}/tosa/{threads}/scrna_eib_junction_barcodes.tsv.gz",
    params:
        strand=_scrna_eib_tosa_strandedness,
        prefix="scrna_eib_benchmark/anchor_{anchor_length}/{strandedness}/tosa/{threads}/scrna_eib",
    threads:
        lambda wildcards: int(wildcards.threads)
    log:
        "log/scrna_eib/tosa_{strandedness}_anchor{anchor_length}_{threads}.log",
    benchmark:
        repeat(
            "benchmarks/scrna_eib/anchor_{anchor_length}/{strandedness}/tosa_{threads}.tsv",
            SCRNA_EIB_CFG.get("benchmark_repeat", 10),
        )
    container:
        config["containers"]["tosa"]
    shell:
        """
        tosa single \
        -m 20 \
        -p {threads} \
        -b {wildcards.anchor_length} \
        -c {input.whitelist} \
        {params.strand} \
        -g {input.gtf} \
        {input.bam} \
        {params.prefix} \
        > {log} 2>&1
        """
