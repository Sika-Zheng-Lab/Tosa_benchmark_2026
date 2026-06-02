"""Per-sample junction comparison across STAR, Tosa, regtools, featureCounts.

Computes pairwise Jaccard indices (set-based and count-based) and
membership counts for UpSet aggregation.
"""

import itertools
import pandas as pd
import numpy as np

# ---------------------------------------------------------------------------
# Snakemake I/O
# ---------------------------------------------------------------------------
in_star = snakemake.input.star
in_tosa = snakemake.input.tosa
in_regtools = snakemake.input.regtools
in_featurecounts = snakemake.input.featurecounts

out_jaccard = snakemake.output.jaccard
out_membership = snakemake.output.membership

sample = snakemake.wildcards.sample

# ---------------------------------------------------------------------------
# Tool names and load
# ---------------------------------------------------------------------------
TOOLS = ["STAR", "Tosa", "regtools", "featureCounts"]

dfs = {}
for tool, path in zip(TOOLS, [in_star, in_tosa, in_regtools, in_featurecounts]):
    dfs[tool] = pd.read_csv(path, sep="\t")

# ---------------------------------------------------------------------------
# Merge into a single table (outer join on junction)
# ---------------------------------------------------------------------------
merged = dfs["STAR"][["junction", "count"]].rename(columns={"count": "STAR"})
for tool in TOOLS[1:]:
    merged = merged.merge(
        dfs[tool][["junction", "count"]].rename(columns={"count": tool}),
        on="junction",
        how="outer",
    )
merged = merged.fillna(0)
for tool in TOOLS:
    merged[tool] = merged[tool].astype(int)

# ---------------------------------------------------------------------------
# Pairwise Jaccard indices
# ---------------------------------------------------------------------------
pairs = list(itertools.combinations(TOOLS, 2))
jaccard_rows = []

for t1, t2 in pairs:
    s1 = set(merged.loc[merged[t1] > 0, "junction"])
    s2 = set(merged.loc[merged[t2] > 0, "junction"])

    n_intersection = len(s1 & s2)
    n_union = len(s1 | s2)
    set_jaccard = n_intersection / n_union if n_union > 0 else float("nan")

    x = merged[t1].values.astype(float)
    y = merged[t2].values.astype(float)
    min_sum = np.minimum(x, y).sum()
    max_sum = np.maximum(x, y).sum()
    count_jaccard = min_sum / max_sum if max_sum > 0 else float("nan")

    jaccard_rows.append({
        "sample": sample,
        "tool1": t1,
        "tool2": t2,
        "set_jaccard": set_jaccard,
        "count_jaccard": count_jaccard,
        "n_intersection": n_intersection,
        "n_union": n_union,
        "n_tool1": len(s1),
        "n_tool2": len(s2),
    })

df_jaccard = pd.DataFrame(jaccard_rows)
df_jaccard.to_csv(out_jaccard, sep="\t", index=False)

# ---------------------------------------------------------------------------
# Membership counts (for UpSet aggregation)
# ---------------------------------------------------------------------------
membership_counts = {}
for _, row in merged.iterrows():
    members = tuple(sorted(t for t in TOOLS if row[t] > 0))
    if members:
        membership_counts[members] = membership_counts.get(members, 0) + 1

membership_rows = []
for members, count in membership_counts.items():
    membership_rows.append({
        "sample": sample,
        "membership": ",".join(members),
        "count": count,
    })

df_membership = pd.DataFrame(membership_rows)
df_membership.to_csv(out_membership, sep="\t", index=False)
