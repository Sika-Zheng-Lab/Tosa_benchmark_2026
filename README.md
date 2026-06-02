# Tosa Benchmark

Benchmark pipeline comparing junction read counting tools:

- **[Tosa](https://github.com/NaotoKubota/Tosa)** — Rust-based junction counter from BAM files
- **[regtools](https://github.com/griffithlab/regtools)** — Junction extraction from BAM files
- **[STAR](https://github.com/alexdobin/STAR)** — Splice junction counts from alignment (SJ.out.tab)
- **[featureCounts](https://subread.sourceforge.net/featureCounts.html)** — Junction counting via Subread (-J option)

## Key features

| | STAR (2.7.11a) | regtools (v1.0.0) | featureCounts (v2.1.1) | **Tosa** |
|---|---|---|---|---|
| Pipeline stage | Alignment coupled (FASTQ + genome index) | Post-alignment (BAM/CRAM) | Post-alignment (BAM) | **Post-alignment (BAM/CRAM)** |
| Multithreading | Yes | No | Yes | **Yes** |
| PE mate deduplication | Yes | No | Yes (`-p --countReadPairs`) | **Yes** |
| Exon-intron boundary read counting | No | No | No | **Yes** |
| Multimapping read handling | Yes | No | Yes | **Yes** |
| Long-read compatibility | No | Yes (`-s XS -t ts`) | Yes (`-L`) | **Yes** |
| Single-cell compatibility | Yes (STARsolo) | Yes | No | **Yes** |
| UMI deduplication | Yes (STARsolo) | No | - | **Yes** |

---

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/) ≥ 7.0
- [Singularity](https://sylabs.io/singularity/) ≥ 3.0

All bioinformatics tools run inside Singularity containers — no local installation needed.

## Containers

| Tool | Image |
|---|---|
| Tosa | `docker://naotokubota/tosa:develop` |
| regtools | `docker://naotokubota/shiba:v0.8.2` |
| featureCounts | `docker://quay.io/biocontainers/subread:2.1.1--h577a1d6_0` |
| STAR | `docker://quay.io/biocontainers/star:2.7.11a--h0033a41_0` |
| samtools | `docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1` |
| Python analysis | `docker://naotokubota/my_notebook:2.0.0` |

## Quick start

### 1. Configure

Edit `config/config.yaml`:

```yaml
# Output directory
outdir: "results"

# Pre-built STAR genome index
star_index: "/path/to/star/index"

# GTF gene annotation (required for featureCounts -J)
gtf: "/path/to/annotation.gtf"

# Strand specificity: null (unstranded), RF, FR, or XS
strandedness: null

# STAR alignment threads
star_threads: 8

# Add your samples (omit R2 for single-end)
samples:
  sample1:
    R1: "/path/to/sample1_R1.fastq.gz"
    R2: "/path/to/sample1_R2.fastq.gz"

# Directories that Singularity containers need to access
singularity_bind_paths:
  - "/path/to/data"
```

> **Note:** All paths in `config/*.yaml` (e.g., `star_index`, `gtf`, `fasta`,
> `samples[*].R1/R2`, `whitelist`, `singularity_bind_paths`) are placeholders
> (`/path/to/...`). Replace them with paths valid on your system.
>
> The large-scale and long-read benchmarks consume sample tables at
> `data/experiment_table_largescale.tsv` (Rogalska et al., ENA `ERP133481`)
> and `data/experiment_table_longread.tsv` (GEO `GSE246705`). The `bam`
> column uses relative paths of the form `bam/<accession>/.../*.bam`;
> place the corresponding BAM files (re-aligned from the public FASTQs of
> those accessions) under a `bam/` directory in your working directory, or
> rewrite the column to absolute paths on your system.

### 2. Dry run

```bash
snakemake -n --use-singularity --snakefile workflow/Snakefile
```

### 3. Run

Singularity containers need explicit bind-mounts for directories outside
the working directory. List them in `singularity_bind_paths` in the config,
then pass them via `--singularity-args`:

```bash
snakemake --cores 8 --use-singularity --snakefile workflow/Snakefile \
    --singularity-args "--bind /path/to/data"
```

> **Tip:** The Snakefile auto-collects bind paths from `outdir`, `star_index`,
> FASTQ paths, and `singularity_bind_paths` into a `SINGULARITY_ARGS` variable.
> You can print it with:
> ```bash
> snakemake --snakefile workflow/Snakefile -q --printshellcmds 2>&1 | head -1
> ```

---

## Bulk RNA-seq Benchmark

Compares junction counting across four tools (Tosa, STAR, regtools, featureCounts)
using real samples and synthetic data with known ground truth.

### What is compared

| Metric | Method |
|---|---|
| Junction count concordance | Pearson / Spearman correlation, scatter plots |
| Junction detection overlap | UpSet plot, Jaccard index |
| PE mate deduplication | Synthetic data with known expected counts |
| Performance | Wall-clock time and peak RSS memory (Snakemake `benchmark`) |

### Design

A **synthetic test dataset** is automatically generated and run through the
pipeline to validate paired-end read deduplication behaviour:

| Scenario | Tosa | STAR | regtools | featureCounts |
|---|---|---|---|---|
| Both mates span same junction | 1 | 1 | 2 | 1 |
| Only one mate spans junction | 1 | 1 | 1 | 1 |
| No junction | 0 | 0 | 0 | 0 |

### Pipeline

```
FASTQ ──► STAR align ──► BAM + SJ.out.tab
                │
                ├──► Tosa bulk ──► junction.tsv.gz
                ├──► regtools junctions extract ──► junctions.bed
                ├──► featureCounts -J ──► .jcounts
                │
                └──► Normalize ──► Compare ──► Plots & Tables
```

### Output

```
results/
├── {sample}/
│   ├── star/           # STAR BAM + SJ.out.tab
│   ├── tosa/           # Tosa junction.tsv.gz
│   ├── regtools/       # regtools junctions.bed
│   ├── featurecounts/  # featureCounts .jcounts
│   └── normalized/     # Common-format junction TSVs
├── comparison/
│   ├── {sample}_junction_counts.tsv
│   ├── {sample}_scatter.png
│   ├── {sample}_venn.png
│   ├── {sample}_correlation.tsv
│   ├── {sample}_benchmark.png
│   ├── simulation_validation.tsv
│   └── simulation_validation.png
└── benchmarks/         # Snakemake benchmark TSVs
```

---

## EIB (Exon-Intron Boundary) Benchmark — Bulk

A dedicated sub-pipeline that evaluates **Tosa's exon-intron boundary counting**
accuracy and thread-scaling performance using fully synthetic data.

Tosa counts reads that are **continuously aligned across exon-intron boundaries**
(i.e., intron retention reads), not spliced reads. When run with the `-g` GTF
flag, Tosa derives 1-base boundary intervals at each end of every intron and
reports the number of reads whose aligned segments contain each boundary:

- **5′ boundary** (exon → intron): `chr:intron_start-(intron_start+1)`
- **3′ boundary** (intron → exon): `chr:(intron_end-1)-intron_end`

### Design

1. **Synthetic genome** — 200 genes on a single chromosome with variable exon
   counts (2, 3, 5, 10), exon sizes (100–2000 bp), intron sizes (100–100 000 bp),
   and canonical GT-AG splice sites.
2. **Boundary-spanning reads** — For each exon-intron boundary, `coverage`
   reads (5, 20, or 100 per gene) are placed at random offsets so that their
   100 bp alignment continuously covers the 1-base boundary interval. The read
   sequence is taken directly from the reference, so STAR maps them as simple
   `100M` alignments (no splicing).
3. **Negative-control reads** — Reads placed entirely within exons, which should
   not be counted as boundary reads.
4. **Ground truth** — Expected boundary counts are computed analytically from
   every read's (and PE mate's) genomic position, ensuring exact correspondence
   with Tosa's coordinate system.
5. **Strandedness sweep** — The same BAM is analysed under XS, RF, and FR strand
   modes to verify consistent counting.
6. **Thread scaling** — Tosa is benchmarked at 1, 4, and 8 threads (10 repeats
   each) to measure wall-clock time, peak RSS, and speedup.

### Configuration

Configure parameters in `config/config.yaml` under `eib_benchmark`:

```yaml
eib_benchmark:
  n_genes: 200
  exon_counts: [2, 3, 5, 10]
  exon_size_range: [100, 2000]
  intron_size_range: [100, 100000]
  coverage_levels: [5, 20, 100]
  read_length: 100
  fragment_length: 300
  fragment_length_sd: 50
  threads: [1, 4, 8]
  strandedness: [XS, RF, FR]
  library_types: [PE, SE]
  benchmark_repeat: 10
  seed: 42
```

### Pipeline

```
generate_eib_genome.py
  └─► genome.fa + annotation.gtf + gene_metadata.tsv
        │
        ├─► generate_eib_reads.py (PE / SE)
        │     └─► FASTQ + expected_boundary_counts.tsv + read_metadata.tsv
        │
        ├─► STAR index ─► STAR align ─► BAM
        │
        ├─► Tosa bulk -g annotation.gtf (×3 strandedness × 3 threads)
        │     └─► *_boundary.tsv.gz
        │
        ├─► evaluate_eib_accuracy.py
        │     └─► accuracy_summary.tsv, strand_accuracy.tsv, stratified plots
        │
        └─► evaluate_eib_performance.py
              └─► performance_summary.tsv, thread scaling / speedup plots
```

### Output

```
eib_benchmark/
├── genome/                          # Synthetic genome + GTF + metadata
├── {PE,SE}/
│   ├── fastq/                       # Synthetic FASTQ
│   ├── star/                        # STAR BAM
│   ├── expected_boundary_counts.tsv # Ground-truth boundary counts
│   └── read_metadata.tsv           # Per-read metadata
├── {PE,SE}_{XS,RF,FR}/tosa/{1,4,8}/
│   ├── eib_boundary.tsv.gz         # Tosa boundary output
│   └── eib_junction.tsv.gz         # Tosa junction output
├── evaluation/
│   ├── accuracy_summary.tsv
│   ├── strand_accuracy.tsv
│   ├── accuracy_by_exon_count.png
│   ├── accuracy_by_intron_size.png
│   ├── accuracy_by_coverage.png
│   ├── accuracy_by_boundary_type.png
│   ├── accuracy_heatmap.png
│   ├── performance_summary.tsv
│   ├── thread_scaling_time.png
│   ├── thread_scaling_memory.png
│   └── speedup.png
└── benchmarks/eib/                  # Snakemake benchmark TSVs
```

### Key metrics

| Metric | Description |
|---|---|
| Sensitivity | Fraction of expected boundaries detected (count > 0) |
| Exact Match Rate | Fraction of boundaries where Tosa count = expected count |
| RMSE | Root mean squared error of count differences |
| Pearson / Spearman r | Correlation between expected and observed counts |
| False Positives | Boundaries reported by Tosa but not in ground truth |

---

## Single-cell (scRNA-seq) Benchmark

A dedicated sub-pipeline that benchmarks **cell-level junction counting** and
**exon-intron boundary (EIB) counting** for scRNA-seq data using three tools:

- **STARsolo** — STAR's built-in single-cell mode with `--soloFeatures Gene SJ`
- **Tosa single** — `tosa single -c <whitelist>` with UMI deduplication
- **regtools** — `regtools junctions extract -b` with per-barcode output

### Tools comparison

| | STARsolo | regtools `-b` | **Tosa single** |
|---|---|---|---|
| Junction detection | ✓ | ✓ | ✓ |
| Cell barcode demux | ✓ (built-in) | ✓ (`-b` option) | ✓ (`-c` option) |
| UMI deduplication | ✓ | ✗ | ✓ |
| EIB counting | ✗ | ✗ | ✓ (`-g` GTF) |
| Output format | MatrixMarket sparse | BED12 + barcode TSV | MatrixMarket sparse + TSV |

### Design

#### Junction simulation

Generates synthetic 10x v3 scRNA-seq reads with known ground truth:

- 50 cells × 20 genes with 2–5 exons each → 51 splice junctions
- R1: 28 bp (16 bp cell barcode + 12 bp UMI from 10x v3 whitelist)
- R2: 90 bp cDNA read spanning splice junctions
- Coverage per cell per junction: randomly sampled from [0, 1, 3, 10]
- Ground truth: exact (barcode, junction, UMI-deduplicated count) triplets

#### EIB simulation

Reuses the bulk EIB genome (200 genes) and generates scRNA-seq reads that
span exon-intron boundaries:

- 30 cells × 869 boundaries × coverage [1, 5, 20]
- Reads placed at random offsets across boundary intervals
- Ground truth: (barcode, boundary, boundary_type, count) per cell

#### regtools `-b` output format

The `-b` flag writes a barcode file with lines corresponding 1:1 to the
BED12 junction file. Each line has the format:

```
<n_barcodes>\t<barcode1:count1,barcode2:count2,...>
```

Where `count` is the **raw read count** (no UMI deduplication). This means
regtools junction counts are systematically higher than STARsolo/Tosa UMI
counts for cells with PCR duplicates.

### Configuration

Parameters are set in `config/config.yaml`:

```yaml
scrna_benchmark:
  chemistry: "v3"
  cb_length: 16
  umi_length: 12
  starsolo_threads: 32
  starsolo_features: "Gene SJ"
  tosa_threads: 8
  strandedness: XS

scrna_simulation:
  n_cells: 50
  n_genes: 20
  exon_counts: [2, 3, 5]
  coverage_per_cell: [0, 1, 3, 10]
  read_length: 90

scrna_eib_benchmark:
  n_cells: 30
  boundary_anchor_lengths: [1, 4, 8]
  threads: [1, 4, 8]
  strandedness: [XS]
  benchmark_repeat: 10
```

### Pipeline

```
generate_scrna_test_data.py
  └─► genome.fa + annotation.gtf + whitelist.txt + FASTQ (R1+R2)
        + expected_junction_counts.tsv
        │
        ├─► STARsolo (--soloFeatures Gene SJ)
        │     ├─► BAM (for Tosa & regtools)
        │     └─► Solo.out/SJ/raw/ (MatrixMarket)
        │
        ├─► Tosa single -c whitelist
        │     └─► *_junction_barcodes.tsv.gz
        │
        ├─► regtools junctions extract -b
        │     └─► junctions.bed + barcodes.tsv
        │
        ├─► normalize_scrna_junctions.py
        │     └─► normalized_sc/{starsolo,tosa,regtools}.tsv
        │
        └─► compare_scrna_junctions.py
              └─► accuracy, correlation, scatter, UpSet, per-cell plots

generate_scrna_eib_reads.py
  └─► FASTQ + expected_boundary_counts.tsv
        │
        ├─► STARsolo align → BAM
        │
        ├─► Tosa single -g GTF -b {1,4,8} -c whitelist
        │     └─► *_boundary_barcodes_detail.tsv.gz
        │
        ├─► evaluate_scrna_eib_accuracy.py
        │     └─► accuracy_summary.tsv, per_cell_accuracy.tsv, plots
        │
        └─► evaluate_scrna_eib_performance.py
              └─► performance_summary.tsv, thread scaling plots
```

### Output

```
scrna_simulation/
├── reference/              # Synthetic genome, GTF, whitelist
├── fastq/                  # 10x v3 FASTQ (R1: CB+UMI, R2: cDNA)
├── starsolo/               # STARsolo BAM + Solo.out/
├── tosa_sc/                # Tosa junction_barcodes.tsv.gz
├── regtools_sc/            # regtools junctions.bed + barcodes.tsv
├── normalized_sc/          # Common-format sparse TSVs
└── expected_junction_counts.tsv

scrna_comparison/
├── simulation_accuracy.tsv
├── simulation_correlation.tsv
├── simulation_junction_counts.tsv
├── simulation_scatter.png
├── simulation_venn.png
├── simulation_per_cell_concordance.png
└── simulation_benchmark.png

scrna_eib_benchmark/
├── reference/              # Whitelist
├── starsolo/               # STARsolo BAM
├── anchor_{1,4,8}/{XS}/tosa/{1,4,8}/
│   ├── scrna_eib_boundary_barcodes_detail.tsv.gz
│   └── scrna_eib_junction_barcodes.tsv.gz
├── evaluation/
│   ├── accuracy_summary.tsv
│   ├── per_cell_accuracy.tsv
│   ├── accuracy_by_boundary_type.png
│   ├── accuracy_by_anchor_length.png
│   ├── per_cell_accuracy.png
│   ├── accuracy_heatmap.png
│   ├── performance_summary.tsv
│   ├── thread_scaling_time.png
│   ├── thread_scaling_memory.png
│   └── speedup.png
├── expected_boundary_counts.tsv
└── benchmarks/scrna_eib/   # Snakemake benchmark TSVs
```

### Key metrics

#### Junction benchmark

| Metric | Description |
|---|---|
| Sensitivity / Precision / F1 | Junction detection at cell × junction level |
| Exact Match Rate | Fraction of entries where tool count = expected count |
| RMSE | Root mean squared error of count differences |
| Pearson / Spearman r | Correlation with ground truth counts |
| Cross-tool correlation | Pairwise correlation at pseudo-bulk and per-cell level |
| Junction overlap | UpSet plot of detected junction sets |

#### EIB benchmark

| Metric | Description |
|---|---|
| Sensitivity | Fraction of expected boundaries detected per cell |
| Precision | Fraction of detected boundaries that are real |
| Exact Match Rate | Count-level accuracy |
| Per-cell accuracy | Sensitivity and exact match per cell barcode |

---

## License

MIT