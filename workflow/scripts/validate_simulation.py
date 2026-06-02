"""
Validate simulation junction counts against expected values.

Compares observed Tosa / STAR / regtools / featureCounts counts for each
junction with the expected counts from the synthetic data generator.
"""

import os
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Snakemake I/O
# ---------------------------------------------------------------------------
in_star = snakemake.input.star
in_tosa = snakemake.input.tosa
in_regtools = snakemake.input.regtools
in_featurecounts = snakemake.input.featurecounts
in_expected = snakemake.input.expected

out_table = snakemake.output.table
out_plot = snakemake.output.plot

os.makedirs(os.path.dirname(out_table), exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Load normalized observed counts
# ---------------------------------------------------------------------------
df_star = pd.read_csv(in_star, sep="\t")
df_tosa = pd.read_csv(in_tosa, sep="\t")
df_regtools = pd.read_csv(in_regtools, sep="\t")
df_featurecounts = pd.read_csv(in_featurecounts, sep="\t")

# Build lookup: junction -> count
obs = {}
for tool, df in [
    ("star", df_star),
    ("tosa", df_tosa),
    ("regtools", df_regtools),
    ("featurecounts", df_featurecounts),
]:
    obs[tool] = dict(zip(df["junction"], df["count"]))

# ---------------------------------------------------------------------------
# 2. Load expected counts
# ---------------------------------------------------------------------------
df_expected = pd.read_csv(in_expected, sep="\t")

# ---------------------------------------------------------------------------
# 3. Compare expected vs observed
# ---------------------------------------------------------------------------
rows = []
for _, row in df_expected.iterrows():
    junction = row["junction"]
    tool = row["tool"]
    expected = int(row["expected_count"])
    observed = int(obs.get(tool, {}).get(junction, 0))
    match = "PASS" if observed == expected else "FAIL"
    rows.append(
        {
            "junction": junction,
            "tool": tool,
            "expected": expected,
            "observed": observed,
            "match": match,
        }
    )

df_result = pd.DataFrame(rows)
df_result.to_csv(out_table, sep="\t", index=False)

# Print summary
n_pass = (df_result["match"] == "PASS").sum()
n_total = len(df_result)
print(f"Simulation validation: {n_pass}/{n_total} checks passed")

# ---------------------------------------------------------------------------
# 4. Grouped bar plot: expected vs observed per tool
# ---------------------------------------------------------------------------
tools = df_result["tool"].unique()
x = np.arange(len(tools))
width = 0.35

fig, ax = plt.subplots(figsize=(8, 5))

expected_vals = [df_result.loc[df_result["tool"] == t, "expected"].values[0] for t in tools]
observed_vals = [df_result.loc[df_result["tool"] == t, "observed"].values[0] for t in tools]

bars1 = ax.bar(x - width / 2, expected_vals, width, label="Expected", color="#4C72B0")
bars2 = ax.bar(x + width / 2, observed_vals, width, label="Observed", color="#DD8452")

ax.set_xlabel("Tool")
ax.set_ylabel("Junction Read Count")
ax.set_title("Simulation Validation: Expected vs Observed Junction Counts")
ax.set_xticks(x)
ax.set_xticklabels([t.capitalize() for t in tools])
ax.legend()

# Add value labels
for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        ax.annotate(
            f"{int(height)}",
            xy=(bar.get_x() + bar.get_width() / 2, height),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=10,
        )

# Add PASS/FAIL annotations
for i, t in enumerate(tools):
    match = df_result.loc[df_result["tool"] == t, "match"].values[0]
    color = "green" if match == "PASS" else "red"
    max_h = max(expected_vals[i], observed_vals[i])
    ax.text(
        i,
        max_h + max_h * 0.1,
        match,
        ha="center",
        va="bottom",
        fontsize=12,
        fontweight="bold",
        color=color,
    )

ax.set_ylim(0, max(max(expected_vals), max(observed_vals)) * 1.3)
plt.tight_layout()
plt.savefig(out_plot, dpi=150)
plt.close()
