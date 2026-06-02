## =============================================================================
## scRNA-seq Tosa rules — Tosa single-cell mode
## =============================================================================

import re as _re

SCRNA_CFG = config.get("scrna_benchmark", {})


def _parse_ds_sample(sample_name):
    """Parse a ds sample name. Returns (base_sample, frac_pct) or (name, None)."""
    m = _re.match(r"^(.+)_ds(\d{3})$", sample_name)
    if m:
        return m.group(1), m.group(2)
    return sample_name, None


def _scrna_tosa_strandedness(wildcards):
    """Return strand flag for Tosa in scRNA-seq mode."""
    s = SCRNA_CFG.get("strandedness", "XS")
    if s:
        return f"--strand {s}"
    return ""


def _scrna_tosa_whitelist(wildcards):
    """Return cell barcode whitelist for Tosa single mode."""
    s = wildcards.scrna_sample
    base, _ = _parse_ds_sample(s)
    if base == "scrna_simulation":
        return "scrna_simulation/reference/whitelist.txt"
    # For real data (including ds), use STARsolo filtered barcodes from the sample's own dir
    return f"{s}/starsolo/{s}.Solo.out/Gene/filtered/barcodes.tsv"


rule scrna_tosa_count:
    """Count junction reads with Tosa single-cell mode (1 thread)."""
    input:
        bam="{scrna_sample}/starsolo/{scrna_sample}.Aligned.sortedByCoord.out.bam",
        bai="{scrna_sample}/starsolo/{scrna_sample}.Aligned.sortedByCoord.out.bam.bai",
        whitelist=_scrna_tosa_whitelist,
    output:
        junction_barcodes="{scrna_sample}/tosa_sc/{scrna_sample}_junction_barcodes.tsv.gz",
    params:
        strand=_scrna_tosa_strandedness,
        prefix="{scrna_sample}/tosa_sc/{scrna_sample}",
    log:
        "log/scrna_tosa_count/{scrna_sample}.log",
    benchmark:
        repeat(
            "benchmarks/{scrna_sample}/tosa_sc.tsv",
            SCRNA_CFG.get("benchmark_repeat", 10),
        )
    container:
        config["containers"]["tosa"]
    shell:
        """
        tosa single \
        -m 20 \
        -c {input.whitelist} \
        {params.strand} \
        {input.bam} \
        {params.prefix} \
        > {log} 2>&1
        """


rule scrna_tosa_count_mt:
    """Count junction reads with Tosa single-cell mode (multi-threaded)."""
    input:
        bam="{scrna_sample}/starsolo/{scrna_sample}.Aligned.sortedByCoord.out.bam",
        bai="{scrna_sample}/starsolo/{scrna_sample}.Aligned.sortedByCoord.out.bam.bai",
        whitelist=_scrna_tosa_whitelist,
    output:
        junction_barcodes="{scrna_sample}/tosa_sc_mt/{scrna_sample}_junction_barcodes.tsv.gz",
    params:
        strand=_scrna_tosa_strandedness,
        prefix="{scrna_sample}/tosa_sc_mt/{scrna_sample}",
    threads:
        SCRNA_CFG.get("tosa_threads", 8)
    log:
        "log/scrna_tosa_count_mt/{scrna_sample}.log",
    benchmark:
        repeat(
            "benchmarks/{scrna_sample}/tosa_sc_mt.tsv",
            SCRNA_CFG.get("benchmark_repeat", 10),
        )
    container:
        config["containers"]["tosa"]
    shell:
        """
        tosa single \
        -m 20 \
        -p {threads} \
        -c {input.whitelist} \
        {params.strand} \
        {input.bam} \
        {params.prefix} \
        > {log} 2>&1
        """
