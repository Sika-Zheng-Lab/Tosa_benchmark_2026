## =============================================================================
## scRNA-seq Downsampling Benchmark rules
## =============================================================================
## For each real scRNA-seq sample, merge multi-lane FASTQs, downsample at
## several fractions with seqtk, then run the full STARsolo → Tosa/regtools
## pipeline for each fraction to measure runtime/memory vs read depth.
## =============================================================================

import re as _re

SCRNA_CFG = config.get("scrna_benchmark", {})
SCRNA_SAMPLES = list(SCRNA_CFG.get("samples", {}).keys()) if SCRNA_CFG.get("samples") else []

# Downsample fractions (e.g. [1.0, 0.75, 0.5, 0.25])
_DS_FRACTIONS = SCRNA_CFG.get("downsample_fractions", [1.0, 0.75, 0.5, 0.25])
_DS_SEED = SCRNA_CFG.get("downsample_seed", 42)
_DS_BENCH_REPEAT = SCRNA_CFG.get("downsample_benchmark_repeat", 3)

# Build fraction percent tags: 1.0 → "100", 0.75 → "075", 0.5 → "050", 0.25 → "025"
def _frac_to_pct(f):
    return f"{int(round(f * 100)):03d}"

DS_FRAC_PCTS = [_frac_to_pct(f) for f in _DS_FRACTIONS]
DS_PCT_TO_FRAC = {_frac_to_pct(f): f for f in _DS_FRACTIONS}

# Downsampled sample names: e.g. scrna_sample1_ds100, scrna_sample1_ds075, ...
SCRNA_DS_SAMPLES = [
    f"{s}_ds{pct}" for s in SCRNA_SAMPLES for pct in DS_FRAC_PCTS
]

# ---------------------------------------------------------------------------
# Helper: get config FASTQ paths for a base sample (as lists)
# ---------------------------------------------------------------------------

def _get_sample_fastqs(sample):
    """Return (r1_list, r2_list) for a base scRNA-seq sample."""
    samples = SCRNA_CFG.get("samples", {})
    r1 = samples[sample]["R1"]
    r2 = samples[sample]["R2"]
    if isinstance(r1, str):
        r1 = [r1]
    if isinstance(r2, str):
        r2 = [r2]
    return r1, r2


# ---------------------------------------------------------------------------
# Rule: Merge multi-lane FASTQs into a single file per read per sample
# ---------------------------------------------------------------------------

def _merge_fastq_input(wildcards):
    """Return the list of FASTQ files for merging."""
    r1, r2 = _get_sample_fastqs(wildcards.scrna_sample)
    return {"r1": r1, "r2": r2}


rule scrna_merge_fastq:
    """Concatenate multi-lane FASTQs into a single R1 and R2 per sample."""
    input:
        unpack(_merge_fastq_input),
    output:
        r1="{scrna_sample}/merged_fastq/R1.fastq.gz",
        r2="{scrna_sample}/merged_fastq/R2.fastq.gz",
    log:
        "log/scrna_merge_fastq/{scrna_sample}.log",
    shell:
        """
        cat {input.r1} > {output.r1} 2> {log}
        cat {input.r2} >> {output.r2} 2>> {log}
        """


# ---------------------------------------------------------------------------
# Rule: Downsample FASTQs with seqtk
# ---------------------------------------------------------------------------

rule scrna_downsample_fastq:
    """Downsample merged FASTQs using seqtk sample.

    For fraction 1.0 (frac_pct=100), symlink the merged file to avoid
    redundant I/O.  For other fractions, run seqtk sample with matched
    seeds so R1 and R2 stay in sync.
    """
    input:
        r1="{scrna_sample}/merged_fastq/R1.fastq.gz",
        r2="{scrna_sample}/merged_fastq/R2.fastq.gz",
    output:
        r1="{scrna_sample}/downsampled/frac{frac_pct}/R1.fastq.gz",
        r2="{scrna_sample}/downsampled/frac{frac_pct}/R2.fastq.gz",
    params:
        fraction=lambda wc: DS_PCT_TO_FRAC.get(wc.frac_pct, 1.0),
        seed=_DS_SEED,
    log:
        "log/scrna_downsample/{scrna_sample}_frac{frac_pct}.log",
    container:
        config["containers"]["seqtk"]
    shell:
        """
        FRAC={params.fraction}
        if [ "$FRAC" = "1.0" ] || [ "$FRAC" = "1" ]; then
            ln -sf $(realpath {input.r1}) {output.r1}
            ln -sf $(realpath {input.r2}) {output.r2}
            echo "Fraction=1.0: symlinked merged FASTQs" > {log} 2>&1
        else
            seqtk sample -s{params.seed} {input.r1} $FRAC \
                | gzip -c > {output.r1} 2> {log}
            seqtk sample -s{params.seed} {input.r2} $FRAC \
                | gzip -c > {output.r2} 2>> {log}
        fi
        """


# ---------------------------------------------------------------------------
# Rule: Count reads in a downsampled BAM
# ---------------------------------------------------------------------------

rule scrna_ds_count_reads:
    """Count reads in a downsampled BAM using samtools flagstat."""
    input:
        bam="{scrna_sample}/starsolo/{scrna_sample}.Aligned.sortedByCoord.out.bam",
    output:
        read_count="read_counts/{scrna_sample}.txt",
    log:
        "log/scrna_ds_count_reads/{scrna_sample}.log",
    container:
        config["containers"]["samtools"]
    shell:
        """
        samtools flagstat {input.bam} 2> {log} | head -1 | awk '{{print $1}}' > {output.read_count}
        """


# ---------------------------------------------------------------------------
# Rule: Aggregate downsample benchmarks into a summary plot
# ---------------------------------------------------------------------------

def _ds_benchmark_inputs(wildcards):
    """Collect benchmark TSVs and read counts across all fractions for a base sample."""
    sample = wildcards.scrna_sample
    inputs = {}
    for pct in DS_FRAC_PCTS:
        ds_name = f"{sample}_ds{pct}"
        inputs[f"bench_tosa_1t_{pct}"] = f"benchmarks/{ds_name}/tosa_sc.tsv"
        inputs[f"bench_tosa_mt_{pct}"] = f"benchmarks/{ds_name}/tosa_sc_mt.tsv"
        inputs[f"bench_regtools_{pct}"] = f"benchmarks/{ds_name}/regtools_sc.tsv"
        inputs[f"read_count_{pct}"] = f"read_counts/{ds_name}.txt"
    return inputs


rule scrna_downsample_benchmark_summary:
    """Aggregate runtime/memory across read-depth fractions into line plot (Figure 4A-B)."""
    input:
        unpack(_ds_benchmark_inputs),
    output:
        plot="scrna_comparison/{scrna_sample}_downsample_benchmark.png",
        summary_tsv="scrna_comparison/{scrna_sample}_downsample_benchmark_summary.tsv",
    params:
        fractions=_DS_FRACTIONS,
        frac_pcts=DS_FRAC_PCTS,
    log:
        "log/scrna_downsample_benchmark/{scrna_sample}.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/aggregate_scrna_downsample_benchmark.py"
