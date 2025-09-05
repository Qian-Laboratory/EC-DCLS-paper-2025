# Analysis for shLCN2 RNA-seq data
options(timeout = 10000)
options(repos = structure(c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))
# update.packages(ask=F)

# List of packages required
depp <- c(
  "languageserver", "httpgd", "dplyr", "edgeR", "ggplot2", "patchwork", "data.table", "tidyr", "tidyverse", "tidygraph"
)
# Check if the packages were installed if not install
depp_new <- depp[!(depp %in% installed.packages())]
depp_new
if (length(depp_new)) {
  install.packages(depp_new)
}
# remotes::install_github('satijalab/seurat-wrappers')
# remotes::install_github("mojaveazure/seurat-disk")

# List of bioconductor packages required by this analysis
BioDepp <- c(
  "rtracklayer", "org.Hs.eg.db", "biomaRt", "DESeq2", "pheatmap", "airway"
)
# Check if the packages were installed if not install
BioDepp_new <- BioDepp[!(BioDepp %in% installed.packages())]
BioDepp_new
# install.packages("BiocManager")
if (length(BioDepp_new)) {
  BiocManager::install(BioDepp_new)
}

# load required packages
sapply(depp, library, character.only = TRUE)
sapply(BioDepp, library, character.only = TRUE)





# edgeR
# ---- Generate file directory list for counts results of a list of samples
sample_list <- c("IK-1", "IK-2", "IK-4", "IK-shLCN2-1", "IK-shLCN2-2", "IK-shLCN2-3")
file_dir <- paste0("../../data/RNA-seq/IK-shLCN2/NCBIref/", sample_list, "/", sample_list, "_Aligned.sortedByCoord.out.bam.counts")
group <- factor(c(rep("IK", 3), rep("IK-shLCN2", 3)), levels = c("IK", "IK-shLCN2"))


# ---- load data
counts_res <- c()
for (f in file_dir) {
  counts <- read.delim(f, as.is = T, skip = 1, row.names = 1)
  # counts <- counts[which(rownames(counts) %in% gene_list), ]
  # dim(counts)
  counts_res <- cbind(counts_res, counts[, ncol(counts)])
}
rownames(counts_res) <- rownames(counts)
colnames(counts_res) <- sample_list
counts_res <- as.data.frame(counts_res)
head(counts_res)

dim(counts)
gene_anno <- counts[, c(1:5)]
head(gene_anno)

# ---- load data to edgeR
dge <- DGEList(counts = counts_res, group = group, genes = gene_anno)


# # ---- filter low expression genes
# keep <- rowSums(cpm(dge) > 0.5) >= 3 # usually the replication number
keep <- filterByExpr(dge)
table(keep)
dge <- dge[keep, , keep.lib.sizes = F]


# ---- TMM Normalization
dge <- calcNormFactors(dge)


# ---- MDS plot
test <- plotMDS(dge)
xlim <- c(round(min(test$x) - 0.25, 1), round(max(test$x) + 0.25, 1))
ylim <- c(round(min(test$y) - 0.25, 1), round(max(test$y) + 0.25, 1))

pdf("../../figures/MDS_plot.pdf", width = 6, height = 6)
my_cols <- (2:(1 + length(levels(group))))
cols <- my_cols[group]
plotMDS(dge, xlim = xlim, ylim = ylim, main = "MDS Plot", col = cols, cex = 0.7)
legend("topright", legend = levels(group), bty = "n", fill = my_cols, cex = 0.7)
dev.off()



# ---- Exact test
dge <- estimateDisp(dge)
et <- exactTest(dge, pair = c("IK-shLCN2", "IK"))
results <- topTags(et, n = nrow(et$table), adjust.method = "BH")
# results <- data.frame(dge$genes, results[, -c(1:5)], dge$counts, stringsAsFactors = F)
results <- as.data.frame(results)
head(results)
write.table(results, "../result/IK_shLCN2_mRNAseq_STAR_featureCounts_edgeR.txt", row.names = T, col.names = NA, sep = "\t", quote = F)


pdf("../../figures/shLCN2vsNC_smearplot.pdf", width = 6, height = 6)
pvalue <- 0.05
de <- decideTestsDGE(et, p.value = pvalue)
sum_de <- summary(de)
detags <- rownames(dge)[as.logical(de)]
plotSmear(et, de.tags = detags, cex = 0.5, main = "IK-shLCN2 vs NC")
abline(h = c(-1, 1), col = "blue", lty = 2)
legend("bottomright", c(
  paste("FDR < ", pvalue, sep = ""),
  rev(paste(c("Down", "Up"), sum_de[-2, 1], sep = ": ")),
  "---: 2-fold line"
),
text.col = c("black", "red", "blue", "blue"), cex = 0.7, bty = "n"
)
dev.off()





# volcano plot
results$logFC <- -results$logFC
results$group <- case_when(
  results$logFC > 1 & results$PValue < 0.05 ~ "up_regulated",
  results$logFC < -1 & results$PValue < 0.05 ~ "down_regulated",
  TRUE ~ "not_significant"
)
results$group <- factor(results$group, levels = c("up_regulated", "down_regulated", "not_significant"))
head(results, 5)

top10_up <- results %>%
  rownames_to_column(var = "gene") %>%
  filter(PValue < 0.05 & logFC > 1) %>%
  distinct(gene, .keep_all = T) %>%
  top_n(10, logFC)
top10_down <- results %>%
  rownames_to_column(var = "gene") %>%
  filter(PValue < 0.05 & logFC < -1) %>%
  distinct(gene, .keep_all = T) %>%
  top_n(100, -logFC)
sig <- rbind(top10_up, top10_down)
colnames(sig)
sig[, c(1, 7, 9, 10, 11)]
sig <- filter(sig, gene %in% c("LCN2", "S100A10", "S100A9"))

mycol <- c("#8400ff", "#dab3ff", "#b9b9b9")

pdf("./volcano_plot.pdf", width = 6, height = 4)
print(
  ggplot() +
    geom_point(data = results, aes(x = logFC, y = -log10(PValue), color = group), size = 2) +
    scale_colour_manual(name = "", values = alpha(mycol, 0.7)) +
    scale_x_continuous(limits = c(-5, 5), breaks = seq(-5, 5, by = 2)) +
    scale_y_continuous(expand = expansion(add = c(2, 0)), limits = c(0, 40), breaks = seq(0, 40, by = 10)) +
    geom_hline(yintercept = c(-log10(0.05)), size = 0.7, color = "black", lty = "dashed") +
    geom_vline(xintercept = c(-1, 1), size = 0.7, color = "black", lty = "dashed") +
    theme_bw() +
    geom_text_repel(data = sig, aes(x = logFC, y = -log10(PValue), label = gene, color = group))
)
dev.off()
