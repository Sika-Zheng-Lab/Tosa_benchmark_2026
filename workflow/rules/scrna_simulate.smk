## =============================================================================
## scRNA-seq Simulation rules — generate synthetic scRNA-seq data & STAR index
## =============================================================================

SCRNA_SIM_CFG = config.get("scrna_simulation", {})
SCRNA_CFG = config.get("scrna_benchmark", {})


rule generate_scrna_test_data:
    """Generate synthetic scRNA-seq genome, FASTQ reads with barcodes, and ground truth."""
    output:
        genome="scrna_simulation/reference/genome.fa",
        fai="scrna_simulation/reference/genome.fa.fai",
        gtf="scrna_simulation/reference/annotation.gtf",
        r1="scrna_simulation/fastq/scrna_R1.fastq.gz",
        r2="scrna_simulation/fastq/scrna_R2.fastq.gz",
        whitelist="scrna_simulation/reference/whitelist.txt",
        expected_junctions="scrna_simulation/expected_junction_counts.tsv",
        gene_metadata="scrna_simulation/reference/gene_metadata.tsv",
    params:
        n_cells=SCRNA_SIM_CFG.get("n_cells", 50),
        n_genes=SCRNA_SIM_CFG.get("n_genes", 20),
        exon_counts=SCRNA_SIM_CFG.get("exon_counts", [2, 3, 5]),
        exon_size_range=SCRNA_SIM_CFG.get("exon_size_range", [200, 1000]),
        intron_size_range=SCRNA_SIM_CFG.get("intron_size_range", [200, 50000]),
        coverage_per_cell=SCRNA_SIM_CFG.get("coverage_per_cell", [0, 1, 3, 10]),
        read_length=SCRNA_SIM_CFG.get("read_length", 90),
        cb_length=SCRNA_CFG.get("cb_length", 16),
        umi_length=SCRNA_CFG.get("umi_length", 12),
        seed=SCRNA_SIM_CFG.get("seed", 42),
    log:
        "log/scrna_simulation/generate_test_data.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/generate_scrna_test_data.py"


rule starsolo_index_simulation:
    """Build a STAR genome index for the scRNA-seq synthetic genome."""
    input:
        genome="scrna_simulation/reference/genome.fa",
        gtf="scrna_simulation/reference/annotation.gtf",
    output:
        directory("scrna_simulation/star_index"),
    threads: 1
    log:
        "log/scrna_simulation/star_index.log",
    container:
        config["containers"]["star"]
    shell:
        """
        GENOME_LEN=$(awk '{{sum += $2}} END {{print sum}}' {input.genome}.fai)
        SA_BASES=$(awk -v gl="$GENOME_LEN" 'BEGIN {{v=log(gl)/log(2)/2-1; v=int(v); if(v>14) v=14; if(v<1) v=1; print v}}')
        mkdir -p {output}
        STAR \
            --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.gtf} \
            --genomeSAindexNbases $SA_BASES \
            --runThreadN {threads} \
            > {log} 2>&1
        """
