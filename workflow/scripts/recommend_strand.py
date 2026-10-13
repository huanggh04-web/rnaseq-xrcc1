import pandas as pd

paths = {
    0: snakemake.input.s0,
    1: snakemake.input.s1,
    2: snakemake.input.s2,
}

rows = []
for strand, path in paths.items():
    df = pd.read_csv(path, sep="\t", index_col=0)
    values = pd.to_numeric(df.iloc[:, 0], errors="coerce").fillna(0)
    total = values.sum()
    assigned = float(values.get("Assigned", 0))
    fraction = assigned / total if total > 0 else 0.0
    rows.append((strand, assigned, total, fraction))

best_strand, best_assigned, best_total, best_fraction = max(rows, key=lambda x: x[3])
sorted_rows = sorted(rows, key=lambda x: x[3], reverse=True)
second_fraction = sorted_rows[1][3] if len(sorted_rows) > 1 else 0.0

names = {
    0: "unstranded (-s 0)",
    1: "stranded forward (-s 1)",
    2: "stranded reverse (-s 2)",
}

lines = []
lines.append("featureCounts strandedness test")
lines.append("")
lines.append("strand\tassigned\ttotal\tassigned_fraction")
for strand, assigned, total, fraction in rows:
    lines.append(f"{names[strand]}\t{int(assigned)}\t{int(total)}\t{fraction:.4f}")

lines.append("")
if best_fraction > second_fraction * 1.2:
    lines.append(f"Recommendation: {names[best_strand]}")
else:
    lines.append(
        "Recommendation: ambiguous. Set featurecounts.strand in config/config.yaml "
        "after checking the library preparation protocol."
    )

with open(snakemake.output.recommendation, "w") as fh:
    fh.write("\n".join(lines) + "\n")
