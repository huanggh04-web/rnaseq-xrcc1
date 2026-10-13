#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(DESeq2); library(ggplot2) })
dir.create("results/deseq2", showWarnings = FALSE)

## ---- 1. counts ----
cnt <- read.delim("results/counts/gene_counts.txt", comment.char = "#",
                  check.names = FALSE, stringsAsFactors = FALSE)
rownames(cnt) <- cnt$Geneid
cnt <- as.matrix(cnt[, 7:ncol(cnt)])
storage.mode(cnt) <- "integer"

## ---- 2. 从 BAM 文件名解析样本号 ----
tag <- sub(".*-", "", sub("\\.Aligned.*", "", basename(colnames(cnt))))
colnames(cnt) <- tag
cat("解析出的样本号:", paste(tag, collapse = " "), "\n\n")

## ---- 3. 分组 ----
grp <- c(A1="Ctrl",A2="Ctrl",A3="Ctrl",A4="Ctrl",A5="Ctrl",A6="Ctrl",
         A7="WT",  A8="WT",  A11="WT", A13="WT", A15="WT", A18="WT",
         A9="Cent",A10="Cent",A12="Cent",A14="Cent",A16="Cent",A17="Cent")
stopifnot(all(tag %in% names(grp)))

run_dds <- function(samples) {
  sub <- cnt[, samples, drop = FALSE]
  cd  <- data.frame(row.names = samples,
                    group = factor(grp[samples], levels = c("Ctrl","WT","Cent")))
  cd  <- droplevels(cd)
  dds <- DESeqDataSetFromMatrix(sub, cd, ~ group)
  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
  DESeq(dds, quiet = TRUE)
}

save_res <- function(dds, label, contrast) {
  r  <- results(dds, contrast = contrast, alpha = 0.05)
  nm <- paste0(label, "__", paste(contrast[2:3], collapse = "_vs_"))
  out <- as.data.frame(r[order(r$padj), ])
  write.table(out, paste0("results/deseq2/", nm, ".tsv"),
              sep = "\t", quote = FALSE, col.names = NA)
  cat(sprintf("  %-14s  padj<0.05: %5d   padj<0.05 & |LFC|>1: %5d   (检验基因 %d)\n",
              paste(contrast[2:3], collapse=" vs "),
              sum(r$padj < 0.05, na.rm = TRUE),
              sum(r$padj < 0.05 & abs(r$log2FoldChange) > 1, na.rm = TRUE),
              sum(!is.na(r$padj))))
  invisible(r)
}

## ================= A. 全部 18 样本 =================
cat("=========== A. 全部 18 样本 ===========\n")
A <- run_dds(colnames(cnt))
rA_wt_ctrl   <- save_res(A, "all18", c("group","WT","Ctrl"))
rA_cent_ctrl <- save_res(A, "all18", c("group","Cent","Ctrl"))
rA_cent_wt   <- save_res(A, "all18", c("group","Cent","WT"))

vsd <- assay(vst(A, blind = TRUE))
pca <- prcomp(t(vsd), scale. = FALSE)
pv  <- round(100 * summary(pca)$importance[2, 1:2], 1)
df  <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
                  group = factor(grp[colnames(vsd)], levels=c("Ctrl","WT","Cent")),
                  sample = colnames(vsd))
write.table(df, "results/deseq2/pca_vst_all18.tsv", sep="\t", quote=FALSE)
p <- ggplot(df, aes(PC1, PC2, color = group, label = sample)) +
     geom_point(size = 3.5) + geom_text(vjust = -1.2, size = 3.2) +
     labs(title = "PCA on VST-normalized counts (all 18 samples)",
          x = sprintf("PC1 (%.1f%%)", pv[1]),
          y = sprintf("PC2 (%.1f%%)", pv[2])) +
     scale_color_manual(values = c(Ctrl="#D55E00", WT="#0072B2", Cent="#009E73")) +
     theme_bw(base_size = 12)
ggsave("results/deseq2/pca_vst_all18.pdf", p, width = 7.5, height = 5.5)
cat("  → 已出图: results/deseq2/pca_vst_all18.pdf\n\n")

## ================= B. 敏感性分析：剔除 A3 / A5 =================
cat("=========== B. 剔除 A3 / A5（Ctrl 只剩 4 只）===========\n")
keep <- setdiff(colnames(cnt), c("A3","A5"))
B <- run_dds(keep)
rB_wt_ctrl   <- save_res(B, "noA3A5", c("group","WT","Ctrl"))
rB_cent_ctrl <- save_res(B, "noA3A5", c("group","Cent","Ctrl"))
rB_cent_wt   <- save_res(B, "noA3A5", c("group","Cent","WT"))

## ================= C. 一致性检查 =================
cat("\n=========== C. WT vs Cent 的一致性 ===========\n")
common <- intersect(rownames(rA_cent_wt), rownames(rB_cent_wt))
d <- abs(rA_cent_wt[common,"log2FoldChange"] - rB_cent_wt[common,"log2FoldChange"])
cat(sprintf("  log2FC 最大差异: %.4f   中位差异: %.4f\n",
            max(d, na.rm=TRUE), median(d, na.rm=TRUE)))
cat(sprintf("  显著基因数: 全部18样本=%d, 剔除A3A5=%d\n",
            sum(rA_cent_wt$padj < 0.05, na.rm=TRUE),
            sum(rB_cent_wt$padj < 0.05, na.rm=TRUE)))

saveRDS(A, "results/deseq2/dds_all18.rds")
saveRDS(B, "results/deseq2/dds_noA3A5.rds")
cat("\n完成。所有结果在 results/deseq2/\n")
