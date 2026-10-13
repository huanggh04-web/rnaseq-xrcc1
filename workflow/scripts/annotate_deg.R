#!/usr/bin/env Rscript
## 给 DESeq2 结果补上基因名

sym <- read.delim("gene_id2name.tsv", header = FALSE, sep = "\t",
                  col.names = c("gid", "symbol"), stringsAsFactors = FALSE)
cat("映射表:", nrow(sym), "行\n")

## 关键：索引和查询都用「去掉版本号」的 ID
SYM <- setNames(sym$symbol, sub("\\..*$", "", as.character(sym$gid)))

fn <- "results/deseq2/cent_vs_wt_MODEL_D_12samples.tsv"
rD <- read.delim(fn, stringsAsFactors = FALSE)
rD$symbol <- SYM[sub("\\..*$", "", as.character(rD$gene))]

cat("未匹配上的基因数:", sum(is.na(rD$symbol)), "/", nrow(rD), "\n\n")
write.table(rD, fn, sep = "\t", quote = FALSE, row.names = FALSE)

cat("================= Model D (仅 WT+Cent, 12 样本) Top 40 =================\n")
top <- head(rD[rD$padj < 0.05, ], 40)
for (i in seq_len(nrow(top))) {
  cat(sprintf("%-14s %-24s log2FC=%+6.2f  padj=%.1e\n",
              ifelse(is.na(top$symbol[i]), "?", top$symbol[i]),
              top$gene[i], top$log2FoldChange[i], top$padj[i]))
}

cat("\n================= 上调最多的 20 个 =================\n")
up <- head(rD[order(-rD$log2FoldChange), ], 20)
for (i in seq_len(nrow(up)))
  cat(sprintf("%-14s log2FC=%+6.2f  padj=%.1e\n",
              ifelse(is.na(up$symbol[i]), "?", up$symbol[i]), up$log2FoldChange[i], up$padj[i]))

cat("\n================= 下调最多的 20 个 =================\n")
dn <- head(rD[order(rD$log2FoldChange), ], 20)
for (i in seq_len(nrow(dn)))
  cat(sprintf("%-14s log2FC=%+6.2f  padj=%.1e\n",
              ifelse(is.na(dn$symbol[i]), "?", dn$symbol[i]), dn$log2FoldChange[i], dn$padj[i]))

cat(sprintf("\n总共 %d 个差异基因 (padj<0.05)\n", sum(rD$padj < 0.05)))
