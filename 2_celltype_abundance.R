# CellType Abundance

# List of packages required
depp <- c(
  "languageserver", "httpgd", "Seurat", "DoubletFinder", "dplyr", "DT", "lme4", "wesanderson",
  "ggplot2", "patchwork", "data.table", "tidyr", "tidyverse", "tidygraph", "circlize",
  "ggraph", "SeuratWrappers", "doParallel", "openxlsx", "viridis", "rstatix", "paletteer"
)
# Check if the packages were installed if not install
depp_new <- depp[!(depp %in% installed.packages())]
if (length(depp_new)) {
  install.packages(depp_new)
}
# remotes::install_github('satijalab/seurat-wrappers')
# remotes::install_github("mojaveazure/seurat-disk")

# List of bioconductor packages required by this analysis
BioDepp <- c(
  "scater", "slingshot", "destiny", "ggrepel", "RColorBrewer", "miloR",
  "GSEABase", "ggpubr", "org.Hs.eg.db", "clusterProfiler", "ComplexHeatmap"
)
# Check if the packages were installed if not install
BioDepp_new <- BioDepp[!(BioDepp %in% installed.packages())]
# install.packages("BiocManager")
if (length(BioDepp_new)) {
  BiocManager::install(BioDepp_new)
}

# load required packages
sapply(depp, library, character.only = TRUE)
sapply(BioDepp, library, character.only = TRUE)




# ---- Create a Milo object
sample.all.original <- readRDS("../../data/sample_integrated_cancer_cb_pc30_filtered.rds")
sample_ca_pca <- subset(sample.all.original, subset = Tissue %in% c("CA", "PCA"))

sample_all_sce <- as.SingleCellExperiment(sample_ca_pca, assay = "RNA")
sample_all_milo <- Milo(sample_all_sce)
sample_all_milo


# ---- Construct KNN graph
sample_all_milo <- buildGraph(sample_all_milo, k = 30, d = 20, reduced.dim = "PCA")


# ---- Defining representative neighbourhoods on the KNN graph
sample_all_milo <- makeNhoods(sample_all_milo, prop = 0.1, k = 30, d = 20, refined = TRUE, reduced_dims = "PCA")
# should see a peak over 20, otherwise you should increase k
pdf("../../figures/milo_all_neigoborhood.pdf", width = 15, height = 5)
plotNhoodSizeHist(sample_all_milo) + geom_vline(xintercept = 20)
dev.off()


# ---- Counting cells in neighbourhoods
head(data.frame(colData(sample_all_milo)))
sample_all_milo <- countCells(sample_all_milo, meta.data = data.frame(colData(sample_all_milo)), sample = "Sample")
head(nhoodCounts(sample_all_milo))


# ---- Defining experimental design
sample_design <- data.frame(colData(sample_all_milo))[, c("Sample", "Tissue")]
# Convert batch info from integer to factor
sample_design$Tissue <- as.factor(sample_design$Tissue)
sample_design <- distinct(sample_design)
rownames(sample_design) <- sample_design$Sample
sample_all_milo <- calcNhoodDistance(sample_all_milo, d = 20, reduced.dim = "PCA")
da_results <- testNhoods(sample_all_milo, design = ~Tissue, design.df = sample_design)
da_results %>%
  arrange(SpatialFDR) %>%
  head()
save(da_results, file = "./da_results.rda")

# add batch effect
sample_design <- data.frame(colData(sample_all_milo))[, c("Sample", "Tissue", "Batch")]
sample_design$Batch <- as.factor(sample_design$Batch)
sample_design <- distinct(sample_design)
rownames(sample_design) <- sample_design$Sample
sample_all_milo_batch <- calcNhoodDistance(sample_all_milo, d = 20, reduced.dim = "PCA")
da_results_batch <- testNhoods(sample_all_milo_batch, design = ~ Batch + Tissue, design.df = sample_design)
da_results_batch %>%
  arrange(SpatialFDR) %>%
  head()
save(da_results_batch, file = "./da_results_batch.rda")


# ---- inspect DA testing results
pdf("../../figures/milo_da_pval_dist.pdf")
ggplot(da_results, aes(PValue)) +
  geom_histogram(bins = 50)
dev.off()

# Then we visualize the test results with a volcano plot (remember that each point here represents a neighbourhood, not a cell).
pdf("../../figures/milo_da_volcano.pdf")
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) +
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)
dev.off()

# Plot single-cell UMAP and neighbourhood graph
sample_all_milo <- buildNhoodGraph(sample_all_milo)
pdf("../../figures/milo_umap.pdf", width = 8, height = 6)
# umap_pl <- plotReducedDim(sample_all_milo, dimred = "UMAP", colour_by = "CellType", text_by = "CellType", text_size = 5) +
#   guides(fill = "none") +
#   scale_color_manual(values = rev(paletteer_d("ggthemes::Classic_20", n = 19))) +
#   NoLegend()
nh_graph_pl <- plotNhoodGraphDA(sample_all_milo, da_results, layout = "UMAP", alpha = 0.1)
# umap_pl + nh_graph_pl + plot_layout(guides = "collect")
nh_graph_pl
dev.off()

# wheather DA is particularly evident in certain cell types
da_results <- annotateNhoods(sample_all_milo, da_results, coldata_col = "CellType")
head(da_results)

pdf("../../figures/milo_ct_fraction.pdf")
ggplot(da_results, aes(CellType_fraction)) +
  geom_histogram(bins = 50)
dev.off()

da_results$celltype <- ifelse(da_results$CellType_fraction < 0.6, "Mixed", da_results$CellType)
pdf("../../figures/milo_ct_fold_change.pdf", width = 6, height = 9)
plotDAbeeswarm(da_results, group.by = "celltype")
dev.off()