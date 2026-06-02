## =============================================================================
## Simulation rules — generate synthetic test data & build STAR index
## =============================================================================


rule generate_test_data:
    """Generate synthetic genome, FASTQ reads, and expected count table."""
    output:
        genome="simulation/reference/genome.fa",
        fai="simulation/reference/genome.fa.fai",
        gtf="simulation/reference/annotation.gtf",
        r1="simulation/fastq/sim_R1.fastq.gz",
        r2="simulation/fastq/sim_R2.fastq.gz",
        expected="simulation/expected_counts.tsv",
    params:
        read_length=config["simulation"]["read_length"],
        num_both=config["simulation"]["num_pairs_both_mates"],
        num_one=config["simulation"]["num_pairs_one_mate"],
        num_none=config["simulation"]["num_pairs_no_junction"],
    log:
        "log/generate_test_data/simulation.log",
    container:
        config["containers"]["notebook"]
    script:
        "../scripts/generate_test_data.py"


rule star_index_simulation:
    """Build a STAR genome index for the synthetic genome."""
    input:
        genome="simulation/reference/genome.fa",
    output:
        directory("simulation/star_index"),
    params:
        # Small genome requires reduced SA index size
        sa_bases=4,
    threads: 1
    log:
        "log/star_index_simulation/simulation.log",
    container:
        config["containers"]["star"]
    shell:
        """
        mkdir -p {output}
        STAR \
            --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.genome} \
            --genomeSAindexNbases {params.sa_bases} \
            --runThreadN {threads} \
            > {log} 2>&1
        """
