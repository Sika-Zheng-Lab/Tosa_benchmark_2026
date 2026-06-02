## =============================================================================
## Large-scale junction comparison rules
## STAR vs Tosa vs regtools vs featureCounts
## =============================================================================


# ---------------------------------------------------------------------------
# Normalize junction outputs to common format (per sample)
# ---------------------------------------------------------------------------

rule normalize_junctions:
    """Normalize STAR/Tosa/regtools/featureCounts junctions to common format."""
    input:
        star=get_star_sj,
        tosa="tosa_1t/{sample}/{sample}_junction.tsv.gz",
        regtools="regtools_1t/{sample}/{sample}_junctions.bed",
        featurecounts="featurecounts_1t/{sample}/{sample}_featurecounts.jcounts",
    output:
        star="normalized/{sample}/star.tsv",
        tosa="normalized/{sample}/tosa.tsv",
        regtools="normalized/{sample}/regtools.tsv",
        featurecounts="normalized/{sample}/featurecounts.tsv",
    log:
        "log/normalize_junctions/{sample}.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/normalize_junctions_largescale.py"


# ---------------------------------------------------------------------------
# Per-sample junction comparison (Jaccard + membership)
# ---------------------------------------------------------------------------

rule compare_junctions_per_sample:
    """Compute pairwise Jaccard indices and membership counts per sample."""
    input:
        star="normalized/{sample}/star.tsv",
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
        "../scripts/compare_junctions_per_sample.py"


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
        norm_star=expand("normalized/{sample}/star.tsv", sample=SAMPLES),
        norm_tosa=expand("normalized/{sample}/tosa.tsv", sample=SAMPLES),
        norm_regtools=expand("normalized/{sample}/regtools.tsv", sample=SAMPLES),
        norm_featurecounts=expand("normalized/{sample}/featurecounts.tsv", sample=SAMPLES),
    output:
        per_sample_jaccard="largescale_results/junction_comparison/per_sample_jaccard.tsv",
        jaccard_set_heatmap="largescale_results/junction_comparison/jaccard_set_heatmap.png",
        jaccard_count_heatmap="largescale_results/junction_comparison/jaccard_count_heatmap.png",
        jaccard_boxplot="largescale_results/junction_comparison/jaccard_boxplot.png",
        upset_plot="largescale_results/junction_comparison/upset_plot.png",
        junction_summary="largescale_results/junction_comparison/junction_summary.tsv",
    params:
        samples=SAMPLES,
    log:
        "log/aggregate_junction_comparison.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/aggregate_junction_comparison.py"


# ---------------------------------------------------------------------------
# Tosa vs regtools per-SJ scatter (single sample, for Figure 2C)
# ---------------------------------------------------------------------------

rule plot_tosa_vs_regtools_scatter:
    """Scatter plot of per-SJ counts: Tosa vs regtools for a single sample."""
    input:
        tosa="normalized/{sample}/tosa.tsv",
        regtools="normalized/{sample}/regtools.tsv",
    output:
        plot="largescale_results/junction_comparison/{sample}_tosa_vs_regtools_scatter.png",
    log:
        "log/plot_tosa_vs_regtools_scatter/{sample}.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/plot_tosa_vs_regtools_scatter.py"


# ---------------------------------------------------------------------------
# Per-sample junction classification by GTF annotation
# ---------------------------------------------------------------------------

rule classify_junctions_by_annotation:
    """Classify each junction as annotated / genic_novel / intergenic."""
    input:
        star="normalized/{sample}/star.tsv",
        tosa="normalized/{sample}/tosa.tsv",
        regtools="normalized/{sample}/regtools.tsv",
        featurecounts="normalized/{sample}/featurecounts.tsv",
        gtf=config["gtf"],
    output:
        classification="junction_annotation/{sample}/classification.tsv",
        summary="junction_annotation/{sample}/summary.tsv",
    log:
        "log/classify_junctions/{sample}.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/classify_junctions_by_annotation.py"


# ---------------------------------------------------------------------------
# Aggregate junction annotation classification across all samples
# ---------------------------------------------------------------------------

rule aggregate_junction_annotation:
    """Aggregate annotation classification and generate summary plots."""
    input:
        summaries=expand(
            "junction_annotation/{sample}/summary.tsv",
            sample=SAMPLES,
        ),
        classifications=expand(
            "junction_annotation/{sample}/classification.tsv",
            sample=SAMPLES,
        ),
    output:
        annotation_summary="largescale_results/junction_annotation/annotation_summary.tsv",
        star_missed_summary="largescale_results/junction_annotation/star_missed_summary.tsv",
        annotation_plot="largescale_results/junction_annotation/annotation_proportion.png",
        star_missed_plot="largescale_results/junction_annotation/star_missed_annotation.png",
    params:
        samples=SAMPLES,
    log:
        "log/aggregate_junction_annotation.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/aggregate_junction_annotation.py"


# ---------------------------------------------------------------------------
# Per-sample motif and count analysis for STAR-missed junctions
# ---------------------------------------------------------------------------

rule analyze_star_missed_motif:
    """Extract splice-site motifs and counts for Tosa junctions missed by STAR."""
    input:
        tosa="normalized/{sample}/tosa.tsv",
        star="normalized/{sample}/star.tsv",
        fasta=config["fasta"],
    output:
        motif_tsv="junction_annotation/{sample}/star_missed_motif.tsv",
    log:
        "log/analyze_star_missed_motif/{sample}.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/analyze_star_missed_motif.py"


# ---------------------------------------------------------------------------
# Aggregate motif analysis across all samples
# ---------------------------------------------------------------------------

rule aggregate_star_missed_motif:
    """Aggregate motif analysis and generate summary plots."""
    input:
        motif_tsvs=expand(
            "junction_annotation/{sample}/star_missed_motif.tsv",
            sample=SAMPLES,
        ),
    output:
        motif_summary="largescale_results/junction_annotation/star_missed_motif_summary.tsv",
        motif_barplot="largescale_results/junction_annotation/star_missed_motif_barplot.png",
        count_distribution="largescale_results/junction_annotation/star_missed_count_distribution.png",
    params:
        samples=SAMPLES,
    log:
        "log/aggregate_star_missed_motif.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/aggregate_star_missed_motif.py"
