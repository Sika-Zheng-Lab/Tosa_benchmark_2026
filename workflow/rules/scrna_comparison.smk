## =============================================================================
## scRNA-seq Comparison / analysis rules
## =============================================================================

import re as _re

SCRNA_CFG = config.get("scrna_benchmark", {})
SCRNA_SAMPLES = list(SCRNA_CFG.get("samples", {}).keys()) if SCRNA_CFG.get("samples") else []
ALL_SCRNA_SAMPLES = ["scrna_simulation"] + SCRNA_SAMPLES


def _parse_ds_sample(sample_name):
    """Parse a ds sample name. Returns (base_sample, frac_pct) or (name, None)."""
    m = _re.match(r"^(.+)_ds(\d{3})$", sample_name)
    if m:
        return m.group(1), m.group(2)
    return sample_name, None


def _scrna_whitelist(wildcards):
    """Return whitelist path for the given scRNA-seq sample."""
    s = wildcards.scrna_sample
    base, _ = _parse_ds_sample(s)
    if base == "scrna_simulation":
        return "scrna_simulation/reference/whitelist.txt"
    return f"{s}/starsolo/{s}.Solo.out/Gene/filtered/barcodes.tsv"


rule normalize_scrna_junctions:
    """Normalize scRNA-seq junction outputs into common sparse cell × junction format."""
    input:
        starsolo_sj_raw="{scrna_sample}/starsolo/{scrna_sample}.Solo.out/SJ/raw",
        tosa="{scrna_sample}/tosa_sc/{scrna_sample}_junction_barcodes.tsv.gz",
        regtools_bed="{scrna_sample}/regtools_sc/{scrna_sample}_junctions.bed",
        regtools_barcodes="{scrna_sample}/regtools_sc/{scrna_sample}_barcodes.tsv",
        whitelist=_scrna_whitelist,
    output:
        starsolo="{scrna_sample}/normalized_sc/starsolo.tsv",
        tosa="{scrna_sample}/normalized_sc/tosa.tsv",
        regtools="{scrna_sample}/normalized_sc/regtools.tsv",
    log:
        "log/normalize_scrna_junctions/{scrna_sample}.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/normalize_scrna_junctions.py"


rule compare_scrna_junctions_simulation:
    """Compare scRNA-seq junction counts against ground truth (simulation only)."""
    input:
        starsolo="scrna_simulation/normalized_sc/starsolo.tsv",
        tosa="scrna_simulation/normalized_sc/tosa.tsv",
        regtools="scrna_simulation/normalized_sc/regtools.tsv",
        expected="scrna_simulation/expected_junction_counts.tsv",
        bench_starsolo="benchmarks/scrna_simulation/starsolo_align.tsv",
        bench_tosa="benchmarks/scrna_simulation/tosa_sc.tsv",
        bench_regtools="benchmarks/scrna_simulation/regtools_sc.tsv",
    output:
        counts="scrna_comparison/simulation_junction_counts.tsv",
        scatter="scrna_comparison/simulation_scatter.png",
        venn="scrna_comparison/simulation_venn.png",
        correlation="scrna_comparison/simulation_correlation.tsv",
        per_cell="scrna_comparison/simulation_per_cell_concordance.png",
        accuracy="scrna_comparison/simulation_accuracy.tsv",
        benchmark="scrna_comparison/simulation_benchmark.png",
    log:
        "log/compare_scrna_junctions/simulation.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/compare_scrna_junctions.py"


rule compare_scrna_junctions_real:
    """Compare scRNA-seq junction counts across tools (real data)."""
    input:
        starsolo="{scrna_sample}/normalized_sc/starsolo.tsv",
        tosa="{scrna_sample}/normalized_sc/tosa.tsv",
        regtools="{scrna_sample}/normalized_sc/regtools.tsv",
        bench_tosa="benchmarks/{scrna_sample}/tosa_sc.tsv",
        bench_regtools="benchmarks/{scrna_sample}/regtools_sc.tsv",
    output:
        counts="scrna_comparison/{scrna_sample}_junction_counts.tsv",
        scatter="scrna_comparison/{scrna_sample}_scatter.png",
        venn="scrna_comparison/{scrna_sample}_venn.png",
        correlation="scrna_comparison/{scrna_sample}_correlation.tsv",
        per_cell="scrna_comparison/{scrna_sample}_per_cell_concordance.png",
        benchmark="scrna_comparison/{scrna_sample}_benchmark.png",
    log:
        "log/compare_scrna_junctions/{scrna_sample}.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/compare_scrna_junctions.py"


rule scrna_compare_benchmark_mt:
    """Compare regtools (1 thread) vs Tosa (multi-threaded) benchmark for scRNA-seq."""
    input:
        bench_tosa_mt="benchmarks/{scrna_sample}/tosa_sc_mt.tsv",
        bench_regtools="benchmarks/{scrna_sample}/regtools_sc.tsv",
    output:
        benchmark="scrna_comparison/{scrna_sample}_benchmark_mt.png",
    log:
        "log/scrna_compare_benchmark_mt/{scrna_sample}.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/compare_scrna_benchmark_mt.py"


rule plot_scrna_accuracy_metrics:
    """Grouped bar chart of per-tool scRNA-seq SJ accuracy metrics (Figure 4A)."""
    input:
        accuracy="scrna_comparison/simulation_accuracy.tsv",
    output:
        plot="scrna_comparison/simulation_accuracy_plot.png",
    log:
        "log/plot_scrna_accuracy_metrics.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/plot_scrna_accuracy_metrics.py"


rule plot_tosa_vs_starsolo_scatter:
    """Pseudo-bulk scatter of Tosa vs STARsolo SJ counts (Figure 4B)."""
    input:
        starsolo="{scrna_sample}/normalized_sc/starsolo.tsv",
        tosa="{scrna_sample}/normalized_sc/tosa.tsv",
    output:
        plot="scrna_comparison/{scrna_sample}_tosa_vs_starsolo_scatter.png",
    log:
        "log/plot_tosa_vs_starsolo_scatter/{scrna_sample}.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/plot_tosa_vs_starsolo_scatter.py"
