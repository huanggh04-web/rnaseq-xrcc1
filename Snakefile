import os
from snakemake.io import glob_wildcards

configfile: "config/config.yaml"

RAWDIR = config["rawdata_dir"]

# Discover samples from files named:
#   <sample>.R1.raw.fastq.gz
SAMPLES, = glob_wildcards(os.path.join(RAWDIR, "{sample}.R1.raw.fastq.gz"))
SAMPLES = sorted(SAMPLES)

if not SAMPLES:
    raise ValueError(
        f"No samples found in {RAWDIR!r} using pattern "
        "'<sample>.R1.raw.fastq.gz'. Check rawdata_dir in config/config.yaml."
    )

STRAND_TEST_SAMPLE = SAMPLES[0]

def raw_r1(wildcards):
    return os.path.join(RAWDIR, f"{wildcards.sample}.R1.raw.fastq.gz")

def raw_r2(wildcards):
    return os.path.join(RAWDIR, f"{wildcards.sample}.R2.raw.fastq.gz")

rule all:
    input:
        expand("results/qc/fastqc/{sample}.R1_fastqc.html", sample=SAMPLES),
        expand("results/qc/fastqc/{sample}.R2_fastqc.html", sample=SAMPLES),
        expand("results/fastp/{sample}.fastp.html", sample=SAMPLES),
        expand("results/alignment/{sample}.Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES),
        "results/qc/strand_test/recommendation.txt",
        "results/counts/gene_counts.txt",
        "results/expression/gene_counts.tsv",
        "results/expression/gene_tpm.tsv",
        "results/expression/gene_fpkm.tsv",
        "results/multiqc/multiqc_report.html",
        "results/exploratory/pca.pdf",
        "results/exploratory/exploratory_summary.txt"

# ------------------------------------------------------------
# 1. Raw read QC
# ------------------------------------------------------------
rule fastqc:
    input:
        r1=raw_r1,
        r2=raw_r2
    output:
        r1_html="results/qc/fastqc/{sample}.R1_fastqc.html",
        r1_zip="results/qc/fastqc/{sample}.R1_fastqc.zip",
        r2_html="results/qc/fastqc/{sample}.R2_fastqc.html",
        r2_zip="results/qc/fastqc/{sample}.R2_fastqc.zip"
    threads: config["threads"]["fastqc"]
    shell:
        r"""
        mkdir -p results/qc/fastqc
        fastqc -t {threads} -o results/qc/fastqc {input.r1} {input.r2}
        r1base=$(basename {input.r1} .fastq.gz)
        r2base=$(basename {input.r2} .fastq.gz)
        mv results/qc/fastqc/${{r1base}}_fastqc.html {output.r1_html}
        mv results/qc/fastqc/${{r1base}}_fastqc.zip  {output.r1_zip}
        mv results/qc/fastqc/${{r2base}}_fastqc.html {output.r2_html}
        mv results/qc/fastqc/${{r2base}}_fastqc.zip  {output.r2_zip}
        """

rule fastp:
    input:
        r1=raw_r1,
        r2=raw_r2
    output:
        r1="results/fastp/{sample}_R1.fastq.gz",
        r2="results/fastp/{sample}_R2.fastq.gz",
        json="results/fastp/{sample}.fastp.json",
        html="results/fastp/{sample}.fastp.html"
    threads: config["threads"]["fastp"]
    shell:
        r"""
        mkdir -p results/fastp
        fastp \
          -i {input.r1} -I {input.r2} \
          -o {output.r1} -O {output.r2} \
          --detect_adapter_for_pe \
          --trim_poly_g \
          --qualified_quality_phred 20 \
          --length_required 30 \
          --thread {threads} \
          --json {output.json} \
          --html {output.html}
        """

# ------------------------------------------------------------
# 2. Reference and STAR index
# ------------------------------------------------------------
rule download_genome:
    output:
        "resources/reference/downloads/GRCm39.primary_assembly.genome.fa.gz"
    params:
        url=config["reference"]["genome_url"]
    shell:
        r"""
        mkdir -p resources/reference/downloads
        wget -c -O {output} {params.url}
        """

rule decompress_genome:
    input:
        "resources/reference/downloads/GRCm39.primary_assembly.genome.fa.gz"
    output:
        config["reference"]["genome"]
    shell:
        "mkdir -p resources/reference && zcat {input} > {output}"

rule download_gtf:
    output:
        "resources/reference/downloads/gencode.vM36.annotation.gtf.gz"
    params:
        url=config["reference"]["gtf_url"]
    shell:
        r"""
        mkdir -p resources/reference/downloads
        wget -c -O {output} {params.url}
        """

rule decompress_gtf:
    input:
        "resources/reference/downloads/gencode.vM36.annotation.gtf.gz"
    output:
        config["reference"]["gtf"]
    shell:
        "mkdir -p resources/reference && zcat {input} > {output}"

rule star_index:
    input:
        genome=config["reference"]["genome"],
        gtf=config["reference"]["gtf"]
    output:
        done=os.path.join(config["reference"]["star_index"], ".done")
    params:
        index=config["reference"]["star_index"],
        overhang=config["star"]["sjdbOverhang"]
    threads: config["threads"]["star_index"]
    shell:
        r"""
        mkdir -p {params.index}
        STAR \
          --runThreadN {threads} \
          --runMode genomeGenerate \
          --genomeDir {params.index} \
          --genomeFastaFiles {input.genome} \
          --sjdbGTFfile {input.gtf} \
          --sjdbOverhang {params.overhang}
        touch {output.done}
        """

# ------------------------------------------------------------
# 3. STAR alignment and BAM indexing
# ------------------------------------------------------------
rule star_align:
    input:
        r1="results/fastp/{sample}_R1.fastq.gz",
        r2="results/fastp/{sample}_R2.fastq.gz",
        index=os.path.join(config["reference"]["star_index"], ".done")
    output:
        bam="results/alignment/{sample}.Aligned.sortedByCoord.out.bam",
        log="results/alignment/{sample}.Log.final.out",
        sj="results/alignment/{sample}.SJ.out.tab",
        reads="results/alignment/{sample}.ReadsPerGene.out.tab"
    params:
        index=config["reference"]["star_index"],
        prefix="results/alignment/{sample}."
    threads: config["threads"]["star_align"]
    shell:
        r"""
        mkdir -p results/alignment
        STAR \
          --runThreadN {threads} \
          --genomeDir {params.index} \
          --readFilesIn {input.r1} {input.r2} \
          --readFilesCommand zcat \
          --outFileNamePrefix {params.prefix} \
          --outSAMtype BAM SortedByCoordinate \
          --outBAMsortingThreadN {threads} \
          --quantMode GeneCounts
        """

rule samtools_index:
    input:
        "results/alignment/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "results/alignment/{sample}.Aligned.sortedByCoord.out.bam.bai"
    threads: 2
    shell:
        "samtools index -@ {threads} {input} {output}"

# ------------------------------------------------------------
# 4. Strandedness check on one sample
# ------------------------------------------------------------
rule strand_test_counts:
    input:
        bam=f"results/alignment/{STRAND_TEST_SAMPLE}.Aligned.sortedByCoord.out.bam",
        gtf=config["reference"]["gtf"]
    output:
        s0="results/qc/strand_test/featureCounts_s0.txt",
        s1="results/qc/strand_test/featureCounts_s1.txt",
        s2="results/qc/strand_test/featureCounts_s2.txt",
        s0_summary="results/qc/strand_test/featureCounts_s0.txt.summary",
        s1_summary="results/qc/strand_test/featureCounts_s1.txt.summary",
        s2_summary="results/qc/strand_test/featureCounts_s2.txt.summary"
    threads: config["threads"]["featurecounts"]
    shell:
        r"""
        mkdir -p results/qc/strand_test
        featureCounts -T {threads} -p --countReadPairs -B -C \
          -t exon -g gene_id -s 0 -a {input.gtf} -o {output.s0} {input.bam}
        featureCounts -T {threads} -p --countReadPairs -B -C \
          -t exon -g gene_id -s 1 -a {input.gtf} -o {output.s1} {input.bam}
        featureCounts -T {threads} -p --countReadPairs -B -C \
          -t exon -g gene_id -s 2 -a {input.gtf} -o {output.s2} {input.bam}
        """

rule strand_test_recommend:
    input:
        s0="results/qc/strand_test/featureCounts_s0.txt.summary",
        s1="results/qc/strand_test/featureCounts_s1.txt.summary",
        s2="results/qc/strand_test/featureCounts_s2.txt.summary"
    output:
        recommendation="results/qc/strand_test/recommendation.txt"
    script:
        "workflow/scripts/recommend_strand.py"

# ------------------------------------------------------------
# 5. Gene-level quantification
# ------------------------------------------------------------
rule featurecounts:
    input:
        bams=expand("results/alignment/{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        gtf=config["reference"]["gtf"]
    output:
        counts="results/counts/gene_counts.txt",
        summary="results/counts/gene_counts.txt.summary"
    params:
        strand=config["featurecounts"]["strand"]
    threads: config["threads"]["featurecounts"]
    shell:
        r"""
        mkdir -p results/counts
        featureCounts \
          -T {threads} \
          -p --countReadPairs \
          -B -C \
          -t exon -g gene_id \
          -s {params.strand} \
          -a {input.gtf} \
          -o {output.counts} \
          {input.bams}
        """

rule expression_matrices:
    input:
        counts="results/counts/gene_counts.txt"
    output:
        raw="results/expression/gene_counts.tsv",
        cpm="results/expression/gene_cpm.tsv",
        tpm="results/expression/gene_tpm.tsv",
        fpkm="results/expression/gene_fpkm.tsv",
        summary="results/expression/expression_summary.txt"
    script:
        "workflow/scripts/make_expression_matrices.py"

# ------------------------------------------------------------
# 6. MultiQC
# ------------------------------------------------------------
rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}.R1_fastqc.zip", sample=SAMPLES),
        expand("results/qc/fastqc/{sample}.R2_fastqc.zip", sample=SAMPLES),
        expand("results/fastp/{sample}.fastp.json", sample=SAMPLES),
        expand("results/alignment/{sample}.Log.final.out", sample=SAMPLES),
        "results/counts/gene_counts.txt.summary"
    output:
        report="results/multiqc/multiqc_report.html",
        data=directory("results/multiqc/multiqc_data")
    threads: 2
    shell:
        "multiqc -f -o results/multiqc results"

# ------------------------------------------------------------
# 7. Exploratory analysis: PCA, correlation, clustering
# ------------------------------------------------------------
rule exploratory:
    input:
        tpm="results/expression/gene_tpm.tsv"
    output:
        pca="results/exploratory/pca.pdf",
        correlation="results/exploratory/sample_correlation.pdf",
        clustering="results/exploratory/sample_clustering.pdf",
        pca_scores="results/exploratory/pca_scores.tsv",
        log_matrix="results/exploratory/log2_tpm_matrix.tsv",
        linkage="results/exploratory/sample_linkage.tsv",
        summary="results/exploratory/exploratory_summary.txt"
    script:
        "workflow/scripts/explore_counts.py"
