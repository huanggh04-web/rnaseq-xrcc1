#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(DESeq2) })

if (!file.exists("results/counts/gene_counts.txt"))
  stop("请在 ~/projects/rnaseq_xrcc1 目录下运行")

cnt <- read.delim("results/counts/gene_counts.txt", comment.char = "#",
                  check.names = FALSE, stringsAsFactors = FALSE)
rownames(cnt) <- cnt$Geneid
cnt <- as.matrix(cnt[, 7:ncol(cnt)])
storage.mode(cnt) <- "integer"
tag <- sub(".*-", "", sub("\\.Aligned.*", "", basename(colnames(cnt))))
colnames(cnt) <- tag

GRP <- c(A1="Ctrl",A2="Ctrl",A3="Ctrl",A4="Ctrl",A5="Ctrl",A6="Ctrl",
         A7="WT",A8="WT",A11="WT",A13="WT",A15="WT",A18="WT",
         A9="Cent",A10="Cent",A12="Cent",A14="Cent",A16="Cent",A17="Cent")

fit <- function(samples, levels) {
  sub <- cnt[, samples, drop = FALSE]
  cd  <- data.frame(row.names = samples, group = factor(GRP[samples], levels = levels))
  dds <- DESeqDataSetFromMatrix(sub, cd, ~ group)
  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
  DESeq(dds, quiet = TRUE)
}
get_res <- function(dds) {
  w <- as.data.frame(results(dds, contrast = c("group","Cent","WT"), alpha = 0.05))
  w$gene <- rownames(w); w <- w[!is.na(w$padj), ]; w[order(w$padj), ]
}

CAT12 <- c("A7","A8","A11","A13","A15","A18",
           "A9","A10","A12","A14","A16","A17")

## 三个模型
dA <- fit(colnames(cnt), c("Ctrl","WT","Cent"));            rA <- get_res(dA)
dB <- fit(setdiff(colnames(cnt), c("A3","A5")),
          c("Ctrl","WT","Cent"));                            rB <- get_res(dB)
dD <- fit(CAT12, c("WT","Cent"));                            rD <- get_res(dD)

cat("\n===== WT vs Cent 在三个模型中的结果 =====\n")
cat(sprintf("  A. 全部 18 样本              : %4d 个 padj<0.05\n", sum(rA$padj < 0.05)))
cat(sprintf("  B. 剔除 A3/A5（16 样本）      : %4d 个 padj<0.05\n", sum(rB$padj < 0.05)))
cat(sprintf("  D. 只用 WT+Cent（12 样本）★   : %4d 个 padj<0.05\n", sum(rD$padj < 0.05)))

sA <- rA$gene[rA$padj < 0.05]; sB <- rB$gene[rB$padj < 0.05]; sD <- rD$gene[rD$padj < 0.05]
cat(sprintf("\n  三者共有的基因: %d\n", length(Reduce(intersect, list(sA, sB, sD)))))
cat(sprintf("  仅模型 D 检出  : %d\n", length(setdiff(sD, union(sA, sB)))))

## 效应量一致性
cm <- Reduce(intersect, list(rA$gene, rB$gene, rD$gene))
dAB <- abs(rA$log2FoldChange[match(cm, rA$gene)] - rB$log2FoldChange[match(cm, rB$gene)])
dAD <- abs(rA$log2FoldChange[match(cm, rA$gene)] - rD$log2FoldChange[match(cm, rD$gene)])
cat(sprintf("\n  效应量 log2FC 最大差异  A-B: %.4f   A-D: %.4f\n",
            max(dAB, na.rm = TRUE), max(dAD, na.rm = TRUE)))

## 带基因名导出模型 D
sym <- read.delim("gene_id2name.tsv", header = FALSE, sep = "\t",
                  col.names = c("gid","symbol"), stringsAsFactors = FALSE)
SYM <- setNames(sym$symbol, sym$gid)
rD$symbol <- SYM[sub("\\..*$", "", rD$gene)]
rD <- rD[, c("gene","symbol","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
write.table(rD, "results/deseq2/cent_vs_wt_MODEL_D_12samples.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n  模型 D 的 Top 30 差异基因:\n")
top <- head(rD, 30)
for (i in seq_len(nrow(top)))
  cat(sprintf("    %-14s %-22s log2FC=%+6.2f  padj=%.2e\n",
              ifelse(is.na(top$symbol[i]), "?", top$symbol[i]),
              top$gene[i], top$log2FoldChange[i], top$padj[i]))

cat("\n完成 → results/deseq2/cent_vs_wt_MODEL_D_12samples.tsv\n")
