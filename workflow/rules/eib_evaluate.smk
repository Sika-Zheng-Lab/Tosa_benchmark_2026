## =============================================================================
## EIB Benchmark — evaluation rules (accuracy & performance)
## =============================================================================


def _eib_dataset_labels():
    """Generate dataset labels like b1_PE_XS, b4_PE_RF, ..., b8_SE_FR."""
    labels = []
    for al in EIB_ANCHOR_LENGTHS:
        for lt in EIB_LIBRARY_TYPES:
            for sd in EIB_STRANDEDNESS:
                labels.append(f"b{al}_{lt}_{sd}")
    return labels


def _eib_benchmark_labels():
    """Generate benchmark labels like b1_PE_XS_1, b4_PE_XS_4, ..., b8_SE_FR_8."""
    labels = []
    for al in EIB_ANCHOR_LENGTHS:
        for lt in EIB_LIBRARY_TYPES:
            for sd in EIB_STRANDEDNESS:
                for t in EIB_THREADS:
                    labels.append(f"b{al}_{lt}_{sd}_{t}")
    return labels


rule eib_evaluate_accuracy:
    """Evaluate Tosa boundary counting accuracy against ground truth."""
    input:
        # Use 1-thread Tosa outputs for accuracy evaluation
        tosa_boundaries=expand(
            "eib_benchmark/anchor_{al}/{lt}_{sd}/tosa/1/eib_boundary.tsv.gz",
            al=EIB_ANCHOR_LENGTHS,
            lt=EIB_LIBRARY_TYPES,
            sd=EIB_STRANDEDNESS,
        ),
        expected_counts=expand(
            "eib_benchmark/anchor_{al}/{lt}/expected_boundary_counts.tsv",
            al=EIB_ANCHOR_LENGTHS,
            lt=EIB_LIBRARY_TYPES,
        ),
        read_metadata=expand(
            "eib_benchmark/{lt}/read_metadata.tsv",
            lt=EIB_LIBRARY_TYPES,
        ),
        gene_metadata="eib_benchmark/genome/gene_metadata.tsv",
    output:
        accuracy_summary="eib_benchmark/evaluation/accuracy_summary.tsv",
        strand_accuracy="eib_benchmark/evaluation/strand_accuracy.tsv",
        by_exon_count="eib_benchmark/evaluation/accuracy_by_exon_count.png",
        by_intron_size="eib_benchmark/evaluation/accuracy_by_intron_size.png",
        by_coverage="eib_benchmark/evaluation/accuracy_by_coverage.png",
        by_boundary_type="eib_benchmark/evaluation/accuracy_by_boundary_type.png",
        by_anchor_length="eib_benchmark/evaluation/accuracy_by_anchor_length.png",
        heatmap="eib_benchmark/evaluation/accuracy_heatmap.png",
    params:
        dataset_labels=_eib_dataset_labels(),
        anchor_lengths=EIB_ANCHOR_LENGTHS,
        library_types=EIB_LIBRARY_TYPES,
        strandedness_modes=EIB_STRANDEDNESS,
    log:
        "log/eib/evaluate_accuracy.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/evaluate_eib_accuracy.py"


rule eib_evaluate_performance:
    """Evaluate Tosa performance across thread counts and anchor lengths."""
    input:
        benchmarks=expand(
            "benchmarks/eib/anchor_{al}/{lt}_{sd}/tosa_{t}.tsv",
            al=EIB_ANCHOR_LENGTHS,
            lt=EIB_LIBRARY_TYPES,
            sd=EIB_STRANDEDNESS,
            t=EIB_THREADS,
        ),
    output:
        time_plot="eib_benchmark/evaluation/thread_scaling_time.png",
        memory_plot="eib_benchmark/evaluation/thread_scaling_memory.png",
        speedup_plot="eib_benchmark/evaluation/speedup.png",
        summary="eib_benchmark/evaluation/performance_summary.tsv",
    params:
        benchmark_labels=_eib_benchmark_labels(),
        thread_values=EIB_THREADS,
        anchor_lengths=EIB_ANCHOR_LENGTHS,
    log:
        "log/eib/evaluate_performance.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/evaluate_eib_performance.py"
