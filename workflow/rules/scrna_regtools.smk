## =============================================================================
## scRNA-seq regtools rules — junction extraction with cell barcode output
## =============================================================================
##
## regtools -b option:
##   -b <file>  specifies an OUTPUT file path for barcode information.
##   regtools reads the CB tag from each BAM record and, for each junction,
##   outputs: <num_barcodes>\t<bc1:count1,bc2:count2,...>
##   The barcode file has the same line ordering as the BED12 output.
##
## NOTE: -b does NOT filter by a whitelist — it writes per-junction barcode
## counts for ALL barcodes found in the BAM. Filtering to valid barcodes
## is done downstream in the normalization script.
## =============================================================================

SCRNA_CFG = config.get("scrna_benchmark", {})


def _scrna_regtools_strandedness(wildcards):
    """Return strand flag for regtools."""
    s = SCRNA_CFG.get("strandedness", "XS")
    if s:
        return s
    return "XS"


rule scrna_regtools_extract:
    """Extract junctions with regtools, outputting per-junction barcode counts."""
    input:
        bam="{scrna_sample}/starsolo/{scrna_sample}.Aligned.sortedByCoord.out.bam",
        bai="{scrna_sample}/starsolo/{scrna_sample}.Aligned.sortedByCoord.out.bam.bai",
    output:
        bed="{scrna_sample}/regtools_sc/{scrna_sample}_junctions.bed",
        barcodes="{scrna_sample}/regtools_sc/{scrna_sample}_barcodes.tsv",
    params:
        strand=_scrna_regtools_strandedness,
    log:
        "log/scrna_regtools_extract/{scrna_sample}.log",
    benchmark:
        repeat(
            "benchmarks/{scrna_sample}/regtools_sc.tsv",
            SCRNA_CFG.get("benchmark_repeat", 10),
        )
    container:
        config["containers"]["regtools"]
    shell:
        """
        regtools junctions extract \
            -s {params.strand} \
            -a 8 \
            -m 20 \
            -M 500000 \
            -b {output.barcodes} \
            -o {output.bed} \
            {input.bam} \
            > {log} 2>&1
        """
