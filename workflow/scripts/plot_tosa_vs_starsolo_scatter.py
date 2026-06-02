"""
Generate pseudo-bulk scatter plot of Tosa vs STARsolo SJ counts (Figure 4B).

Reads normalized scRNA-seq junction sparse TSV files (barcode, junction, count)
for STARsolo and Tosa, aggregates to pseudo-bulk (sum per junction), and
produces a scatter plot with Pearson r annotation.

Input:
  - Normalized STARsolo sparse TSV (barcode, junction, count)
  - Normalized Tosa sparse TSV (barcode, junction, count)

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
from scipy import stats

# ---------------------------------------------------------------------------
# Snakemake I/O
# ---------------------------------------------------------------------------
in_starsolo = snakemake.input.starsolo
in_tosa = snakemake.input.tosa

out_plot = snakemake.output.plot

os.makedirs(os.path.dirname(out_plot), exist_ok=True)

# ---------------------------------------------------------------------------
# Load and aggregate to pseudo-bulk
# ---------------------------------------------------------------------------

def load_pseudobulk(path):
    df = pd.read_csv(path, sep="\t")
    df.columns = ["barcode", "junction", "count"]
    return df.groupby("junction")["count"].sum()


pb_starsolo = load_pseudobulk(in_starsolo)
pb_tosa = load_pseudobulk(in_tosa)

# Merge on junction
merged = pd.DataFrame({
    "STARsolo": pb_starsolo,
    "Tosa": pb_tosa,
}).fillna(0).astype(int)

x = merged["Tosa"].values.astype(float)
y = merged["STARsolo"].values.astype(float)

# Correlation (non-zero union)
mask = (x > 0) | (y > 0)
if mask.sum() > 2:
    r, _ = stats.pearsonr(x[mask], y[mask])
else:
    r = float("nan")

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(4.5, 4.5))

ax.scatter(
    np.log10(x + 1),
    np.log10(y + 1),
    s=8,
    alpha=0.4,
    color="#1f77b4",
    edgecolors="none",
)

# Diagonal reference
lims = [0, max(np.log10(x.max() + 1), np.log10(y.max() + 1)) * 1.05]
ax.plot(lims, lims, "r--", linewidth=1, alpha=0.7)

ax.set_xlabel("Tosa  log₁₀(count + 1)", fontsize=14)
ax.set_ylabel("STARsolo  log₁₀(count + 1)", fontsize=14)
ax.set_title(f"Pseudo-bulk SJ Counts: Tosa vs STARsolo\n(Pearson r = {r:.4f})", fontsize=15)
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.set_aspect("equal")
ax.tick_params(labelsize=12)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
fig.savefig(out_plot, dpi=300, bbox_inches="tight")
plt.close(fig)
