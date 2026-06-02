## =============================================================================
## STARsolo alignment rules — scRNA-seq alignment
## =============================================================================

import re as _star_re

SCRNA_CFG = config.get("scrna_benchmark", {})


def _parse_ds_sample_star(sample_name):
    """Parse a ds sample name. Returns (base_sample, frac_pct) or (name, None)."""
    m = _star_re.match(r"^(.+)_ds(\d{3})$", sample_name)
    if m:
        return m.group(1), m.group(2)
    return sample_name, None


def _scrna_get_star_index(sample):
    """Return STAR index directory for scRNA-seq sample."""
    if sample == "scrna_simulation":
        return "scrna_simulation/star_index"
    return config["star_index"]


def _scrna_get_whitelist(sample):
    """Return cell barcode whitelist path."""
    if sample == "scrna_simulation":
        return "scrna_simulation/reference/whitelist.txt"
    wl = SCRNA_CFG.get("whitelist")
    if wl:
        return wl
    raise ValueError(f"No whitelist specified for scRNA-seq sample {sample}")


def _starsolo_input_fastqs(wildcards):
    """Return FASTQ files for STARsolo alignment."""
    s = wildcards.scrna_sample
    if s == "scrna_simulation":
        return {
            "r1": "scrna_simulation/fastq/scrna_R1.fastq.gz",
            "r2": "scrna_simulation/fastq/scrna_R2.fastq.gz",
        }
    base, pct = _parse_ds_sample_star(s)
    if pct is not None:
        return {
            "r1": f"{base}/downsampled/frac{pct}/R1.fastq.gz",
            "r2": f"{base}/downsampled/frac{pct}/R2.fastq.gz",
        }
    samples = SCRNA_CFG.get("samples", {})
    return {
        "r1": samples[s]["R1"],
        "r2": samples[s]["R2"],
    }


def _starsolo_readfiles_param(wildcards):
    """Build --readFilesIn parameter: R2 (cDNA) first, then R1 (CB+UMI) for STARsolo."""
    fqs = _starsolo_input_fastqs(wildcards)
    r1 = fqs["r1"]
    r2 = fqs["r2"]
    if isinstance(r2, list):
        return f"{','.join(r2)} {','.join(r1)}"
    return f"{r2} {r1}"


def _starsolo_index_input(wildcards):
    return _scrna_get_star_index(wildcards.scrna_sample)


def _starsolo_whitelist_input(wildcards):
    return _scrna_get_whitelist(wildcards.scrna_sample)


rule starsolo_align:
    """Align scRNA-seq reads with STARsolo (10x Chromium v3)."""
    input:
        r1=lambda wc: _starsolo_input_fastqs(wc)["r1"],
        r2=lambda wc: _starsolo_input_fastqs(wc)["r2"],
        index=_starsolo_index_input,
        whitelist=_starsolo_whitelist_input,
    output:
        bam="{scrna_sample}/starsolo/{scrna_sample}.Aligned.sortedByCoord.out.bam",
        sj="{scrna_sample}/starsolo/{scrna_sample}.SJ.out.tab",
        solo_sj_raw=directory("{scrna_sample}/starsolo/{scrna_sample}.Solo.out/SJ/raw"),
        filtered_barcodes="{scrna_sample}/starsolo/{scrna_sample}.Solo.out/Gene/filtered/barcodes.tsv",
    params:
        readfiles=_starsolo_readfiles_param,
        prefix="{scrna_sample}/starsolo/{scrna_sample}.",
        cb_len=SCRNA_CFG.get("cb_length", 16),
        umi_start=SCRNA_CFG.get("cb_length", 16) + 1,
        umi_len=SCRNA_CFG.get("umi_length", 12),
        features=SCRNA_CFG.get("starsolo_features", "Gene SJ"),
        cell_filter=SCRNA_CFG.get("starsolo_cell_filter", "EmptyDrops_CR"),
    threads:
        SCRNA_CFG.get("starsolo_threads", 32)
    log:
        "log/starsolo_align/{scrna_sample}.log",
    benchmark:
        "benchmarks/{scrna_sample}/starsolo_align.tsv"
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
            --soloType CB_UMI_Simple \
            --soloCBwhitelist {input.whitelist} \
            --soloCBstart 1 \
            --soloCBlen {params.cb_len} \
            --soloUMIstart {params.umi_start} \
            --soloUMIlen {params.umi_len} \
            --soloBarcodeReadLength 0 \
            --soloFeatures {params.features} \
            --soloCellFilter {params.cell_filter} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS NM XS CB UB \
            --outFilterMultimapNmax 1 \
            --outBAMsortingBinsN 50 \
            --alignIntronMin 20 \
            --outFileNamePrefix {params.prefix} \
            --runThreadN {threads} \
            > {log} 2>&1
        """


rule starsolo_samtools_index:
    """Index STARsolo BAM file."""
    input:
        bam="{scrna_sample}/starsolo/{scrna_sample}.Aligned.sortedByCoord.out.bam",
    output:
        bai="{scrna_sample}/starsolo/{scrna_sample}.Aligned.sortedByCoord.out.bam.bai",
    log:
        "log/starsolo_samtools_index/{scrna_sample}.log",
    container:
        config["containers"]["samtools"]
    shell:
        """
        samtools index {input.bam} > {log} 2>&1
        """
