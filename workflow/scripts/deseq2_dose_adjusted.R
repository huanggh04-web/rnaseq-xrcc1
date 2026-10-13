suppressPackageStartupMessages(library(DESeq2))

outdir <- "results/deseq2/dose_adjusted"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

read_counts <- function() {
  f_expr <- "results/expression/gene_counts.tsv"
  f_fc   <- "results/counts/gene_counts.txt"

  if (file.exists(f_expr)) {
    x <- read.delim(f_expr, check.names = FALSE, stringsAsFactors = FALSE)
    if (ncol(x) < 2) stop("Expression count file has too few columns")
    colnames(x)[1] <- "gene_id"
    ann <- intersect(c("gene_id", "gene_name", "symbol", "Chr", "Start", "End", "Strand", "Length"), colnames(x))
    sample_cols <- setdiff(colnames(x), ann)
    y <- data.frame(gene_id = x$gene_id, x[, sample_cols, drop = FALSE], check.names = FALSE)
    return(y)
  }

  if (file.exists(f_fc)) {
    x <- read.delim(f_fc, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
    if (!"Geneid" %in% colnames(x)) stop("featureCounts file has no Geneid column")
    ann <- intersect(c("Geneid", "Chr", "Start", "End", "Strand", "Length"), colnames(x))
    sample_cols <- setdiff(colnames(x), ann)
    y <- data.frame(gene_id = x$Geneid, x[, sample_cols, drop = FALSE], check.names = FALSE)
    return(y)
  }

  stop("No count file found: expected results/expression/gene_counts.tsv or results/counts/gene_counts.txt")
}

raw <- read_counts()
sample_cols <- setdiff(colnames(raw), "gene_id")

extract_sample <- function(x) {
  m <- regmatches(x, gregexpr("A(1[0-8]|[1-9])(?![0-9])", x, perl = TRUE))[[1]]
  if (length(m) == 0) return(NA_character_)
  tail(m, 1)
}

samples <- vapply(sample_cols, extract_sample, character(1))

cat("Sample column mapping:\n")
print(data.frame(column = sample_cols, sample = samples, stringsAsFactors = FALSE))

need <- paste0("A", 1:18)

if (any(is.na(samples))) {
  stop("Could not extract A1-A18 from all sample columns")
}
if (!all(need %in% samples)) {
  stop("Not all A1-A18 sample columns were found")
}

idx <- match(need, samples)
counts <- as.matrix(raw[, sample_cols[idx], drop = FALSE])
rownames(counts) <- raw$gene_id
colnames(counts) <- need

if (any(counts < 0 | counts != round(counts))) {
  stop("Count matrix contains non-integer values")
}
storage.mode(counts) <- "integer"

group_vec <- c(
  A1="Ctrl", A2="Ctrl", A3="Ctrl", A4="Ctrl", A5="Ctrl", A6="Ctrl",
  A7="WT", A8="WT", A11="WT", A13="WT", A15="WT", A18="WT",
  A9="Cent", A10="Cent", A12="Cent", A14="Cent", A16="Cent", A17="Cent"
)

dose_vec <- c(
  A1=0, A2=0, A3=0, A4=0, A5=0, A6=0,
  A7=5.25, A8=1.00, A11=66.75, A13=5.75, A15=513.75, A18=3.00,
  A9=353.00, A10=812.00, A12=820.75, A14=1121.25, A16=380.75, A17=79.25
)

coldata_all <- data.frame(
  sample = need,
  group = factor(group_vec[need], levels = c("Ctrl", "WT", "Cent")),
  dose = as.numeric(dose_vec[need]),
  row.names = need
)
coldata_all$log2dose <- as.numeric(scale(log2(coldata_all$dose + 1)))

add_symbol <- function(res) {
  f <- "gene_id2name.tsv"
  if (!file.exists(f)) return(res)

  map <- read.delim(f, header = FALSE, stringsAsFactors = FALSE)
  if (nrow(map) == 0) return(res)

  if (!is.na(map[1, 1]) && grepl("gene", map[1, 1], ignore.case = TRUE)) {
    map2 <- read.delim(f, header = TRUE, stringsAsFactors = FALSE)
    if (ncol(map2) >= 2) map <- map2
  }
  if (ncol(map) < 2) return(res)

  colnames(map)[1:2] <- c("gene_id", "symbol")
  map$gene_id <- sub("\\..*$", "", map$gene_id)
  idx <- match(sub("\\..*$", "", res$gene_id), map$gene_id)
  res$symbol <- map$symbol[idx]
  res
}

summary_df <- data.frame()

run_model <- function(coldata, design_str, contrast, tag) {
  cat("\n=== Running:", tag, "===\n")
  cat("Design:", design_str, "\n")

  smp <- rownames(coldata)
  cnt <- counts[, smp, drop = FALSE]

  dds <- DESeqDataSetFromMatrix(
    countData = cnt,
    colData = coldata,
    design = as.formula(design_str)
  )
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds, quiet = TRUE)

  if (is.null(contrast)) {
    res <- results(dds, name = "log2dose")
  } else {
    res <- results(dds, contrast = contrast)
  }

  res <- as.data.frame(res)
  res$gene_id <- rownames(res)
  res <- res[, c("gene_id", setdiff(colnames(res), "gene_id"))]
  res <- add_symbol(res)

  sig  <- subset(res, !is.na(padj) & padj < 0.05)
  sig1 <- subset(sig, !is.na(log2FoldChange) & abs(log2FoldChange) > 1)

  write.table(res,  file = file.path(outdir, paste0(tag, "_results.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(sig,  file = file.path(outdir, paste0(tag, "_padj005.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(sig1, file = file.path(outdir, paste0(tag, "_padj005_lfc1.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  summary_df <<- rbind(summary_df, data.frame(
    model = tag,
    tested = nrow(res),
    padj005 = nrow(sig),
    padj005_lfc1 = nrow(sig1)
  ))

  cat(sprintf("%s: tested=%d, padj<0.05=%d, padj<0.05 & |LFC|>1=%d\n",
              tag, nrow(res), nrow(sig), nrow(sig1)))
}

# 全部 18 个样本
run_model(coldata_all, "~ group", c("group", "Cent", "WT"),
          "all18_group_Cent_vs_WT")
run_model(coldata_all, "~ log2dose", NULL,
          "all18_dose_only")
run_model(coldata_all, "~ group + log2dose", c("group", "Cent", "WT"),
          "all18_group_plus_dose_Cent_vs_WT")

# 只保留 WT + Cent
coldata_wc <- coldata_all[coldata_all$group %in% c("WT", "Cent"), ]
coldata_wc$group <- droplevels(coldata_wc$group)

run_model(coldata_wc, "~ group", c("group", "Cent", "WT"),
          "WT_Cent_group_Cent_vs_WT")
run_model(coldata_wc, "~ log2dose", NULL,
          "WT_Cent_dose_only")
run_model(coldata_wc, "~ group + log2dose", c("group", "Cent", "WT"),
          "WT_Cent_group_plus_dose_Cent_vs_WT")

write.table(summary_df,
            file = file.path(outdir, "model_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== MODEL SUMMARY ===\n")
print(summary_df)
