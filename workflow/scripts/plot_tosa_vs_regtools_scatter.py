"""
Generate scatter plot comparing per-SJ counts between Tosa and regtools
for a single sample (Figure 2G for bulk short-read; Figure 3E for long-read).

For bulk short-read PE data (Figure 2G), this illustrates systematic count
inflation by regtools due to PE mate double-counting: regtools consistently
reports higher counts for SJs where both mates span the same locus. For
single-end long-read data (Figure 3E), it confirms concordant counting
between Tosa and regtools in the absence of PE overcounting.

Input:
  - Normalized Tosa junction TSV (junction, strand, count)
  - Normalized regtools junction TSV (junction, strand, count)

Output:
  - Scatter plot PNG
"""

import os
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Snakemake I/O
# ---------------------------------------------------------------------------
in_tosa = snakemake.input.tosa
in_regtools = snakemake.input.regtools

out_plot = snakemake.output.plot

os.makedirs(os.path.dirname(out_plot), exist_ok=True)

# ---------------------------------------------------------------------------
# Load and merge
# ---------------------------------------------------------------------------
df_tosa = pd.read_csv(in_tosa, sep="\t")
df_regtools = pd.read_csv(in_regtools, sep="\t")

merged = df_tosa.merge(
    df_regtools,
    on="junction",
    how="inner",
    suffixes=("_tosa", "_regtools"),
)

x = merged["count_tosa"].values.astype(float)
y = merged["count_regtools"].values.astype(float)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(4.5, 4.5))

log_x = np.log10(x + 1)
log_y = np.log10(y + 1)

ax.scatter(log_x, log_y, s=5, color="#1f77b4", alpha=0.3, edgecolors="none")

# Diagonal reference (y = x)
lims = [0, max(log_x.max(), log_y.max()) * 1.05]
ax.plot(lims, lims, "r--", linewidth=1, alpha=0.7, label="y = x")

ax.set_xlabel("Tosa  log\u2081\u2080(count + 1)", fontsize=14)
ax.set_ylabel("regtools  log\u2081\u2080(count + 1)", fontsize=14)
ax.set_title("Per-SJ Count: Tosa vs regtools", fontsize=15)
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.set_aspect("equal")
ax.legend(fontsize=11)
ax.tick_params(labelsize=12)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
fig.savefig(out_plot, dpi=300, bbox_inches="tight")
plt.close(fig)
