"""Analyse splice-site motifs and read counts for STAR-missed junctions.

For each junction detected by Tosa, determine whether it was also reported
in STAR's SJ.out.tab.  Extract the donor (5') and acceptor (3') dinucleotide
motifs from the reference FASTA and classify each junction as:

  - canonical:      GT/AG or CT/AC  (reverse-complement equivalent)
  - semi_canonical:  GC/AG, CT/GC, AT/AC, GT/AT
  - non_canonical:   everything else

Outputs:
  star_missed_motif.tsv — per-junction: junction, tosa_count, donor, acceptor,
                          motif_category, group (tosa_only | both)
"""

import os

import pandas as pd
import pysam

# ---------------------------------------------------------------------------
# Snakemake I/O
# ---------------------------------------------------------------------------
in_tosa = snakemake.input.tosa
in_star = snakemake.input.star
in_fasta = snakemake.input.fasta

out_tsv = snakemake.output.motif_tsv

os.makedirs(os.path.dirname(out_tsv), exist_ok=True)

# ---------------------------------------------------------------------------
# Motif classification
# ---------------------------------------------------------------------------

CANONICAL = {("GT", "AG"), ("CT", "AC")}
SEMI_CANONICAL = {("GC", "AG"), ("CT", "GC"), ("AT", "AC"), ("GT", "AT")}


def classify_motif(donor, acceptor):
    pair = (donor, acceptor)
    if pair in CANONICAL:
        return "canonical"
    if pair in SEMI_CANONICAL:
        return "semi_canonical"
    return "non_canonical"


# ---------------------------------------------------------------------------
# Load normalised junction data
# ---------------------------------------------------------------------------
df_tosa = pd.read_csv(in_tosa, sep="\t")
df_star = pd.read_csv(in_star, sep="\t")

star_set = set(df_star["junction"])

# Tag each Tosa junction as tosa_only or both
df_tosa = df_tosa.copy()
df_tosa["group"] = df_tosa["junction"].apply(
    lambda j: "both" if j in star_set else "tosa_only"
)

# ---------------------------------------------------------------------------
# Extract motifs from reference FASTA
# ---------------------------------------------------------------------------
fa = pysam.FastaFile(in_fasta)

donors = []
acceptors = []
motif_cats = []

for junc in df_tosa["junction"]:
    # Format: chr:start-end  (1-based inclusive intron coordinates)
    chrom, coords = junc.split(":")
    start, end = coords.split("-")
    start = int(start)  # 1-based inclusive intron start
    end = int(end)      # 1-based inclusive intron end

    # pysam.fetch uses 0-based half-open coordinates
    # Donor: first 2 bases of intron  (positions start, start+1 in 1-based)
    donor = fa.fetch(chrom, start - 1, start + 1).upper()
    # Acceptor: last 2 bases of intron (positions end-1, end in 1-based)
    acceptor = fa.fetch(chrom, end - 2, end).upper()

    donors.append(donor)
    acceptors.append(acceptor)
    motif_cats.append(classify_motif(donor, acceptor))

fa.close()

df_tosa["donor"] = donors
df_tosa["acceptor"] = acceptors
df_tosa["motif_category"] = motif_cats

# Rename count column for clarity
df_out = df_tosa[["junction", "count", "donor", "acceptor", "motif_category", "group"]].copy()
df_out = df_out.rename(columns={"count": "tosa_count"})

df_out.to_csv(out_tsv, sep="\t", index=False)
