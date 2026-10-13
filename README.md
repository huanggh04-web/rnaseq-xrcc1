# RNA-seq Snakemake workflow (mouse GRCm39)

This workflow performs:

1. FastQC on raw paired-end FASTQ
2. fastp adapter/quality trimming
3. STAR alignment to mouse GRCm39
4. featureCounts gene-level quantification
5. TPM, FPKM, and CPM matrices
6. MultiQC report
7. PCA, sample correlation, and hierarchical clustering
8. Strandedness test using featureCounts

Differential expression is intentionally not included yet because the
sample-to-group mapping for the 18 RNA-seq samples has not been found.

## 1. Install the environment

```bash
conda env create -f envs/rnaseq.yaml
conda activate rnaseq
```

If the environment already exists:

```bash
conda env update -f envs/rnaseq.yaml --prune
conda activate rnaseq
```

## 2. Configure the raw-data path

Edit `config/config.yaml`:

```yaml
rawdata_dir: /rawdata
```

Set this to the absolute path of the directory that directly contains:

```text
<sample>.R1.raw.fastq.gz
<sample>.R2.raw.fastq.gz
```

Check it:

```bash
ls -lh /rawdata/*.R1.raw.fastq.gz
ls -lh /rawdata/*.R2.raw.fastq.gz
```

If there are 18 files of each type, the configuration is correct.

## 3. Dry run

```bash
snakemake -n -p
```

## 4. Run the strandedness test first

```bash
snakemake --cores 32 results/qc/strand_test/recommendation.txt
cat results/qc/strand_test/recommendation.txt
```

Set the recommended value in `config/config.yaml`:

```yaml
featurecounts:
  strand: 0   # or 1 or 2
```

If you change the strand value after counts have already been produced,
force the quantification and downstream steps to rerun:

```bash
snakemake --cores 32 --forcerun featurecounts
```

## 5. Run the full workflow

```bash
snakemake --cores 32 -p
```

The first run downloads the reference and builds a STAR index. This
can take 1-2 hours and about 30 GB of disk space for the index alone.

## 6. Main results

```text
results/qc/fastqc/
results/fastp/
results/multiqc/multiqc_report.html
results/alignment/
results/counts/gene_counts.txt
results/expression/gene_tpm.tsv
results/expression/gene_fpkm.tsv
results/expression/gene_cpm.tsv
results/exploratory/pca.pdf
results/exploratory/sample_correlation.pdf
results/exploratory/sample_clustering.pdf
```

## 7. Important notes

- Default reference: mouse GRCm39 + GENCODE M36.
- STAR index is built locally.
- The workflow does not require group labels for the steps above.
- If the sequencing library is stranded, use the strand recommendation
  before producing the final count matrix.
- Keep the raw FASTQ files unchanged.
