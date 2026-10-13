suppressPackageStartupMessages({
  ok_cp <- requireNamespace("clusterProfiler", quietly = TRUE)
  ok_db <- requireNamespace("org.Mm.eg.db", quietly = TRUE)
  ok_gg <- requireNamespace("ggplot2", quietly = TRUE)
})

cat("Package check:\n")
cat("clusterProfiler:", ok_cp, "\n")
cat("org.Mm.eg.db:", ok_db, "\n")
cat("ggplot2:", ok_gg, "\n")

if (!ok_cp || !ok_db || !ok_gg) {
  cat("\nMissing packages. Install with:\n")
  cat("BiocManager::install(c('clusterProfiler', 'org.Mm.eg.db'))\n")
  quit(status = 0)
}

library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)

gene_file <- "results/deseq2/dose_adjusted/WT_Cent_group_plus_dose_lfc1_symbols.txt"
res_file  <- "results/deseq2/dose_adjusted/WT_Cent_group_plus_dose_Cent_vs_WT_results.tsv"
outdir    <- "results/deseq2/dose_adjusted/enrichment"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

genes <- readLines(gene_file)
genes <- unique(genes[!is.na(genes) & genes != "" & genes != "NA"])
cat("Input genes:", length(genes), "\n")

allres <- read.delim(res_file, check.names = FALSE, stringsAsFactors = FALSE)
if ("symbol" %in% colnames(allres)) {
  universe <- unique(allres$symbol[
    !is.na(allres$symbol) & allres$symbol != "" & allres$symbol != "NA"
  ])
} else {
  universe <- NULL
}
cat("Universe genes:", length(universe), "\n")

run_go <- function(ont) {
  cat("\nRunning GO:", ont, "\n")
  ego <- enrichGO(
    gene = genes,
    universe = universe,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    readable = TRUE
  )

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    cat("No significant GO terms for", ont, "\n")
    return(NULL)
  }

  df <- as.data.frame(ego)
  write.table(
    df,
    file = file.path(outdir, paste0("GO_", ont, "_results.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  p <- dotplot(ego, showCategory = 20) +
    ggtitle(paste0("GO ", ont, " - dose-adjusted residual genes"))

  ggsave(
    file.path(outdir, paste0("GO_", ont, "_dotplot.pdf")),
    plot = p, width = 8, height = 10
  )

  cat("Significant GO", ont, "terms:", nrow(df), "\n")
  ego
}

bp <- run_go("BP")
mf <- run_go("MF")
cc <- run_go("CC")

if (!is.null(bp)) {
  cat("\nTop GO BP terms:\n")
  print(head(as.data.frame(bp)[, c("ID", "Description", "GeneRatio", "p.adjust")], 10))
}
