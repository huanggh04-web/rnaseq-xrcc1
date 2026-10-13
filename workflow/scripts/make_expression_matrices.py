import os
import numpy as np
import pandas as pd

infile = snakemake.input.counts
df = pd.read_csv(infile, sep="\t", comment="#")

required = ["Geneid", "Length"]
for col in required:
    if col not in df.columns:
        raise ValueError(f"Missing column {col!r} in {infile}")

sample_cols = list(df.columns[6:])
if not sample_cols:
    raise ValueError("No BAM/count columns found in featureCounts output")

counts = df.set_index("Geneid")[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
lengths = pd.to_numeric(df.set_index("Geneid")["Length"], errors="coerce")
counts.index = counts.index.astype(str)
lengths.index = counts.index

clean_names = [
    os.path.basename(c).replace(".Aligned.sortedByCoord.out.bam", "")
    for c in counts.columns
]
counts.columns = clean_names

# Raw counts
counts.to_csv(snakemake.output.raw, sep="\t")

# CPM
lib_sizes = counts.sum(axis=0)
lib_sizes = lib_sizes.replace(0, np.nan)
cpm = counts.div(lib_sizes, axis=1) * 1e6
cpm = cpm.fillna(0)
cpm.to_csv(snakemake.output.cpm, sep="\t")

# FPKM and TPM
length_kb = lengths / 1000.0
rate = counts.div(length_kb, axis=0)
rate = rate.replace([np.inf, -np.inf], np.nan).fillna(0)

fpkm = rate.div(rate.sum(axis=0).replace(0, np.nan), axis=1) * 1e6
fpkm = fpkm.fillna(0)
fpkm.to_csv(snakemake.output.fpkm, sep="\t")

tpm = rate.div(rate.sum(axis=0).replace(0, np.nan), axis=1) * 1e6
tpm = tpm.fillna(0)
tpm.to_csv(snakemake.output.tpm, sep="\t")

with open(snakemake.output.summary, "w") as fh:
    fh.write(f"Genes: {counts.shape[0]}\n")
    fh.write(f"Samples: {counts.shape[1]}\n")
    fh.write("\nLibrary sizes (raw counts):\n")
    fh.write(counts.sum(axis=0).to_string())
    fh.write("\n")
