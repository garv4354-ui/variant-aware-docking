suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(EnhancedVolcano)
  library(org.Hs.eg.db); library(clusterProfiler)
})
set.seed(1)
out <- "../results/01_rnaseq"

## 1. Load data
counts <- read.delim("../data/gse81089/GSE81089_raw_counts.tsv.gz",
                     row.names = 1, check.names = FALSE)
meta <- read.csv("../data/gse81089/sample_info.csv",
                 row.names = 1, check.names = FALSE)
meta <- meta[colnames(counts), ]
stopifnot(identical(rownames(meta), colnames(counts)))

## 2. Labels from the sample title (T = tumor, N = normal)
meta$condition <- factor(ifelse(grepl("N$", meta$title), "normal", "tumor"),
                         levels = c("normal", "tumor"))
print(table(meta$condition))

## 3. DESeq2
dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts),
                              colData = meta[, c("title", "condition")],
                              design = ~ condition)
keep <- rowSums(counts(dds) >= 10) >= 10
dds <- dds[keep, ]
cat("Genes kept after filtering:", nrow(dds), "\n")
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "tumor", "normal"))

## 4. Add gene symbols and save the table
res_df <- as.data.frame(res)
res_df$entrez <- rownames(res_df)
res_df$symbol <- mapIds(org.Hs.eg.db, keys = res_df$entrez,
                        keytype = "ENTREZID", column = "SYMBOL",
                        multiVals = "first")
res_df <- res_df[order(res_df$padj), ]
write.csv(res_df, file.path(out, "deseq2_tumor_vs_normal.csv"), row.names = FALSE)

sig_up   <- which(res_df$padj < 0.05 & res_df$log2FoldChange > 1)
sig_down <- which(res_df$padj < 0.05 & res_df$log2FoldChange < -1)
cat("Significant UP in tumor:", length(sig_up),
    " | DOWN in tumor:", length(sig_down), "\n")

## 5. QC plot (PCA)
vsd <- vst(dds, blind = TRUE)
ggsave(file.path(out, "pca.png"), plotPCA(vsd, intgroup = "condition"),
       width = 6, height = 5, dpi = 200)

## 6. Volcano plot
res_v <- res_df[complete.cases(res_df$padj, res_df$symbol), ]
res_v$padj[res_v$padj == 0] <- min(res_v$padj[res_v$padj > 0]) / 10
png(file.path(out, "volcano.png"), width = 2000, height = 1800, res = 250)
print(EnhancedVolcano(res_v, lab = res_v$symbol,
      x = "log2FoldChange", y = "padj", pCutoff = 0.05, FCcutoff = 1,
      title = "Lung tumor vs normal (GSE81089)", subtitle = NULL))
dev.off()

## 7. Pathway enrichment (GO biological process)
universe <- res_df$entrez[complete.cases(res_df$padj)]
for (dir in c("up", "down")) {
  idx <- if (dir == "up") sig_up else sig_down
  ego <- enrichGO(gene = res_df$entrez[idx], universe = universe,
                  OrgDb = org.Hs.eg.db, ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  write.csv(as.data.frame(ego), file.path(out, paste0("GO_", dir, ".csv")),
            row.names = FALSE)
  if (nrow(as.data.frame(ego)) > 0) {
    ggsave(file.path(out, paste0("GO_", dir, ".png")),
           dotplot(ego, showCategory = 15), width = 8, height = 7, dpi = 200)
  }
}

## 8. Sanity check and target candidates
genes <- c("EGFR", "KRAS", "ALK", "AGER", "SFTPC", "SPP1")
print(res_df[res_df$symbol %in% genes,
             c("symbol", "baseMean", "log2FoldChange", "padj")])

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("DONE\n")
