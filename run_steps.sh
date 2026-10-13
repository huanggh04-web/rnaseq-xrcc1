#!/usr/bin/env bash
set -euo pipefail

# Activate the workflow environment.
conda activate rnaseq

# 1. Check the workflow.
snakemake -n -p

# 2. Build the STAR index and run the strandedness test.
snakemake --cores 32 results/qc/strand_test/recommendation.txt

echo
echo "Strandedness recommendation:"
cat results/qc/strand_test/recommendation.txt
echo
echo "Edit config/config.yaml if needed, then run:"
echo "snakemake --cores 32 -p"
