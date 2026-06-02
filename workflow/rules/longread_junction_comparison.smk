## =============================================================================
## Long-read junction comparison rules
## Tosa vs regtools vs featureCounts (no STAR — minimap2-aligned data)
## =============================================================================


# ---------------------------------------------------------------------------
# Normalize junction outputs to common format (per sample)
# ---------------------------------------------------------------------------

rule normalize_junctions:
    """Normalize Tosa/regtools/featureCounts junctions to common format."""
    input:
        tosa="tosa_1t/{sample}/{sample}_junction.tsv.gz",
        regtools="regtools_1t/{sample}/{sample}_junctions.bed",
        featurecounts="featurecounts_1t/{sample}/{sample}_featurecounts.jcounts",
    output:
        tosa="normalized/{sample}/tosa.tsv",
        regtools="normalized/{sample}/regtools.tsv",
        featurecounts="normalized/{sample}/featurecounts.tsv",
    log:
        "log/normalize_junctions/{sample}.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/normalize_junctions_longread.py"


# ---------------------------------------------------------------------------
# Per-sample junction comparison (Jaccard + membership)
# ---------------------------------------------------------------------------

rule compare_junctions_per_sample:
    """Compute pairwise Jaccard indices and membership counts per sample."""
    input:
        tosa="normalized/{sample}/tosa.tsv",
        regtools="normalized/{sample}/regtools.tsv",
        featurecounts="normalized/{sample}/featurecounts.tsv",
    output:
        jaccard="junction_comparison/{sample}/jaccard.tsv",
        membership="junction_comparison/{sample}/membership.tsv",
    log:
        "log/compare_junctions/{sample}.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/compare_junctions_per_sample_longread.py"


# ---------------------------------------------------------------------------
# Aggregate junction comparison across all samples
# ---------------------------------------------------------------------------

rule aggregate_junction_comparison:
    """Aggregate Jaccard indices and generate summary plots."""
    input:
        jaccards=expand(
            "junction_comparison/{sample}/jaccard.tsv",
            sample=SAMPLES,
        ),
        memberships=expand(
            "junction_comparison/{sample}/membership.tsv",
            sample=SAMPLES,
        ),
        norm_tosa=expand("normalized/{sample}/tosa.tsv", sample=SAMPLES),
        norm_regtools=expand("normalized/{sample}/regtools.tsv", sample=SAMPLES),
        norm_featurecounts=expand("normalized/{sample}/featurecounts.tsv", sample=SAMPLES),
    output:
        per_sample_jaccard="longread_results/junction_comparison/per_sample_jaccard.tsv",
        jaccard_set_heatmap="longread_results/junction_comparison/jaccard_set_heatmap.png",
        jaccard_count_heatmap="longread_results/junction_comparison/jaccard_count_heatmap.png",
        jaccard_boxplot="longread_results/junction_comparison/jaccard_boxplot.png",
        upset_plot="longread_results/junction_comparison/upset_plot.png",
        junction_summary="longread_results/junction_comparison/junction_summary.tsv",
    params:
        samples=SAMPLES,
    log:
        "log/aggregate_junction_comparison.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/aggregate_junction_comparison_longread.py"


# ---------------------------------------------------------------------------
# Tosa vs regtools per-SJ scatter (single sample)
# ---------------------------------------------------------------------------

rule plot_tosa_vs_regtools_scatter:
    """Scatter plot of per-SJ counts: Tosa vs regtools for a single sample."""
    input:
        tosa="normalized/{sample}/tosa.tsv",
        regtools="normalized/{sample}/regtools.tsv",
    output:
        plot="longread_results/junction_comparison/{sample}_tosa_vs_regtools_scatter.png",
    log:
        "log/plot_tosa_vs_regtools_scatter/{sample}.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/plot_tosa_vs_regtools_scatter.py"


# ---------------------------------------------------------------------------
# Investigate featureCounts low Jaccard index
# ---------------------------------------------------------------------------

rule investigate_featurecounts:
    """Diagnose why featureCounts shows low Jaccard vs Tosa/regtools."""
    input:
        norm_tosa=expand("normalized/{sample}/tosa.tsv", sample=SAMPLES),
        norm_regtools=expand("normalized/{sample}/regtools.tsv", sample=SAMPLES),
        norm_featurecounts=expand("normalized/{sample}/featurecounts.tsv", sample=SAMPLES),
        raw_jcounts="featurecounts_1t/" + SAMPLES[0] + "/" + SAMPLES[0] + "_featurecounts.jcounts",
        gtf=config["gtf"],
    output:
        report="longread_results/fc_investigation/diagnostic_report.txt",
        plot_counts="longread_results/fc_investigation/junction_counts.png",
        plot_offset="longread_results/fc_investigation/offset_analysis.png",
        plot_annotation="longread_results/fc_investigation/annotation_status.png",
        plot_intron_size="longread_results/fc_investigation/intron_size.png",
        plot_overlap="longread_results/fc_investigation/overlap_categories.png",
    params:
        samples=SAMPLES,
    log:
        "log/investigate_featurecounts.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/investigate_featurecounts_longread.py"
