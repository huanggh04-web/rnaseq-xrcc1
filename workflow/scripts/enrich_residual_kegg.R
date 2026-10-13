suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

outdir <- "results/deseq2/dose_adjusted/enrichment"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

genes <- readLines("results/deseq2/dose_adjusted/WT_Cent_group_plus_dose_lfc1_symbols.txt")
genes <- unique(genes[!is.na(genes) & genes != "" & genes != "NA"])

entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
cat("Mapped genes:", nrow(entrez), "\n")

kk <- enrichKEGG(
  gene = entrez$ENTREZID,
  organism = "mmu",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)

if (is.null(kk) || nrow(as.data.frame(kk)) == 0) {
  cat("No significant KEGG pathways\n")
  quit(status = 0)
}

df <- as.data.frame(kk)
write.table(df, file.path(outdir, "KEGG_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

print(head(df[, c("ID","Description","GeneRatio","p.adjust")], 20))
