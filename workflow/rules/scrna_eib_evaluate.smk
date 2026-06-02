## =============================================================================
## scRNA-seq EIB Benchmark — evaluation rules (accuracy & performance)
## =============================================================================

SCRNA_EIB_CFG = config.get("scrna_eib_benchmark", {})
SCRNA_EIB_ANCHOR_LENGTHS = [str(a) for a in SCRNA_EIB_CFG.get("boundary_anchor_lengths", [1, 4, 8])]
SCRNA_EIB_THREADS = [str(t) for t in SCRNA_EIB_CFG.get("threads", [1])]
SCRNA_EIB_STRANDEDNESS = SCRNA_EIB_CFG.get("strandedness", ["XS"])


def _scrna_eib_dataset_labels():
    """Generate dataset labels like b1_XS, b4_XS, b8_XS."""
    labels = []
    for al in SCRNA_EIB_ANCHOR_LENGTHS:
        for sd in SCRNA_EIB_STRANDEDNESS:
            labels.append(f"b{al}_{sd}")
    return labels


def _scrna_eib_benchmark_labels():
    """Generate benchmark labels like b1_XS_1, b4_XS_4, b8_XS_8."""
    labels = []
    for al in SCRNA_EIB_ANCHOR_LENGTHS:
        for sd in SCRNA_EIB_STRANDEDNESS:
            for t in SCRNA_EIB_THREADS:
                labels.append(f"b{al}_{sd}_{t}")
    return labels


rule scrna_eib_compute_expected_counts:
    """Compute per-cell expected boundary counts at a specific anchor length."""
    input:
        genome="eib_benchmark/genome/genome.fa",
        metadata="eib_benchmark/genome/gene_metadata.tsv",
        read_positions="scrna_eib_benchmark/read_positions.tsv",
    output:
        expected="scrna_eib_benchmark/expected_boundary_counts_b{al}.tsv",
    params:
        boundary_anchor_length="{al}",
    log:
        "log/scrna_eib/compute_expected_b{al}.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/compute_scrna_eib_expected_counts.py"


rule scrna_eib_evaluate_accuracy:
    """Evaluate Tosa single-cell boundary counting accuracy against ground truth."""
    input:
        tosa_boundaries=expand(
            "scrna_eib_benchmark/anchor_{al}/{sd}/tosa/1/scrna_eib_boundary_barcodes_detail.tsv.gz",
            al=SCRNA_EIB_ANCHOR_LENGTHS,
            sd=SCRNA_EIB_STRANDEDNESS,
        ),
        expected_boundaries=expand(
            "scrna_eib_benchmark/expected_boundary_counts_b{al}.tsv",
            al=SCRNA_EIB_ANCHOR_LENGTHS,
        ),
        gene_metadata="eib_benchmark/genome/gene_metadata.tsv",
    output:
        accuracy_summary="scrna_eib_benchmark/evaluation/accuracy_summary.tsv",
        per_cell_accuracy="scrna_eib_benchmark/evaluation/per_cell_accuracy.tsv",
        by_boundary_type="scrna_eib_benchmark/evaluation/accuracy_by_boundary_type.png",
        by_anchor_length="scrna_eib_benchmark/evaluation/accuracy_by_anchor_length.png",
        per_cell_plot="scrna_eib_benchmark/evaluation/per_cell_accuracy.png",
        heatmap="scrna_eib_benchmark/evaluation/accuracy_heatmap.png",
    params:
        dataset_labels=_scrna_eib_dataset_labels(),
        anchor_lengths=SCRNA_EIB_ANCHOR_LENGTHS,
        strandedness_modes=SCRNA_EIB_STRANDEDNESS,
    log:
        "log/scrna_eib/evaluate_accuracy.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/evaluate_scrna_eib_accuracy.py"


rule scrna_eib_evaluate_performance:
    """Evaluate Tosa single-cell performance across thread counts."""
    input:
        benchmarks=expand(
            "benchmarks/scrna_eib/anchor_{al}/{sd}/tosa_{t}.tsv",
            al=SCRNA_EIB_ANCHOR_LENGTHS,
            sd=SCRNA_EIB_STRANDEDNESS,
            t=SCRNA_EIB_THREADS,
        ),
    output:
        time_plot="scrna_eib_benchmark/evaluation/thread_scaling_time.png",
        memory_plot="scrna_eib_benchmark/evaluation/thread_scaling_memory.png",
        speedup_plot="scrna_eib_benchmark/evaluation/speedup.png",
        summary="scrna_eib_benchmark/evaluation/performance_summary.tsv",
    params:
        benchmark_labels=_scrna_eib_benchmark_labels(),
        thread_values=SCRNA_EIB_THREADS,
        anchor_lengths=SCRNA_EIB_ANCHOR_LENGTHS,
    log:
        "log/scrna_eib/evaluate_performance.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/evaluate_scrna_eib_performance.py"
