# CODEX data analysis
# List of packages required
depp <- c(
  "httpgd", "Seurat", "DoubletFinder", "dplyr", "DT", "lme4","svDialogs",
  "ggplot2", "patchwork", "data.table", "tidyr", "tidyverse", "tidygraph",
  "ggraph", "SeuratWrappers", "doParallel", "openxlsx", "viridis", "rstatix",
  "reticulate", "plyr", "openxlsx", "readxl", "magrittr", "purrr", "cowplot"
)
# Check if the packages were installed if not install
depp_new <- depp[!(depp %in% installed.packages())]
if (length(depp_new)) {
  install.packages(depp_new)
}
# remotes::install_github('satijalab/seurat-wrappers')
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

# List of bioconductor packages required by this analysis
BioDepp <- c(
  "scater", "slingshot", "destiny", "ggrepel", "RColorBrewer",
  "GSEABase", "ggpubr", "org.Hs.eg.db", "clusterProfiler"
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





# Individual sample
## EC4 (IA G1 EC p53wt Ciliated)
inputdata <- read_csv("../../data/EC_codex/EC4CA_1_measurements.csv")
dim(inputdata)
colnames(inputdata)
head(inputdata)
ec4.codex.obj <- LoadAkoya(
  filename = "../../data/EC_codex/EC4CA_1_measurements.csv",
  type = "qupath",
  fov = "EC4"
)
rownames(ec4.codex.obj)
ec4.codex.obj[["B2M"]] <- inputdata$"B2M-BX043: Cell: Mean"
ec4.codex.obj[["cellsize"]] <- ec4.codex.obj$"Cell..Area.µm.2"
ec4.codex.obj[["Sample"]] <- "EC4"

VlnPlot(ec4.codex.obj, features = c("B2M", "nCount_Akoya", "cellsize"), ncol = 3, pt.size = 0, log = TRUE) & geom_hline(yintercept = 5)
ggsave("../../figures/EC4/EC4_vlnplot_low.pdf", width = 8, height = 7)

VlnPlot(ec4.codex.obj, features = c("B2M", "nCount_Akoya", "cellsize"), ncol = 3, pt.size = 0, log = FALSE) & geom_hline(yintercept = 200)
ggsave("../../figures/EC4/EC4_vlnplot_high.pdf", width = 8, height = 7)

ec4.codex.obj <- subset(ec4.codex.obj, subset = B2M > 5 & B2M < 100 & nCount_Akoya > 70 & nCount_Akoya < 1200 & cellsize > 3 & cellsize < 200)

ec4.codex.obj <- NormalizeData(object = ec4.codex.obj, normalization.method = "CLR", margin = 2)
ec4.codex.obj <- ScaleData(ec4.codex.obj)
VariableFeatures(ec4.codex.obj) <- rownames(ec4.codex.obj) # since the panel is small, treat all features as variable.
ec4.codex.obj <- RunPCA(object = ec4.codex.obj)
ec4.codex.obj <- RunUMAP(object = ec4.codex.obj, dims = 1:20)
ec4.codex.obj <- FindNeighbors(ec4.codex.obj, k.param = 20, dims = 1:20)
ec4.codex.obj <- FindClusters(object = ec4.codex.obj, resolution = seq(0.1, 2, 0.1), algorithm = 2, n.start = 1)

FeaturePlot(object = ec4.codex.obj, features = c("B2M", "nCount_Akoya", "cellsize"), raster = FALSE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
ggsave("../../figures/EC4/EC4_umap_qc.pdf", width = 8, height = 7)

ec4.codex.obj <- ScaleData(ec4.codex.obj, vars.to.regress = c("nCount_Akoya"))
VariableFeatures(ec4.codex.obj) <- rownames(ec4.codex.obj) # since the panel is small, treat all features as variable.
ec4.codex.obj <- RunPCA(object = ec4.codex.obj)

ec4.codex.obj <- JackStraw(ec4.codex.obj, num.replicate = 100, prop.freq = 0.01)
ec4.codex.obj <- ScoreJackStraw(ec4.codex.obj, dims = 1:20)
JackStrawPlot(ec4.codex.obj, dims = 1:20) & ElbowPlot(ec4.codex.obj, ndims = 20)
ggsave("../../figures/EC4/EC4_jackstraw.pdf", width = 15, height = 7)

ec4.codex.obj <- RunUMAP(object = ec4.codex.obj, dims = 1:20)
ec4.codex.obj <- FindNeighbors(ec4.codex.obj, k.param = 20, dims = 1:20)
ec4.codex.obj <- FindClusters(object = ec4.codex.obj, resolution = seq(0.1, 2, 0.1), algorithm = 2, n.start = 1)


# ---- seurat cluster
for (i in seq(0.1, 2, 0.1)) {
  print(i)
  DimPlot(ec4.codex.obj, label = TRUE, label.box = TRUE, raster = TRUE, group.by = paste0("Akoya_snn_res.", i))
  ggsave(paste0("../../figures/EC4/EC4_cluster_res_", i, ".pdf"), width = 8, height = 7)
}


# ---- seurat feature
for (gene in rownames(ec4.codex.obj)) {
  print(gene)
  FeaturePlot(object = ec4.codex.obj, features = gene, raster = FALSE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
  ggsave(paste0("../../figures/EC4/EC4_feature_", gene, ".pdf"), width = 8, height = 7)
}


# ---- codex image
ImageDimPlot(ec4.codex.obj, cols = "parade", group.by = "Akoya_snn_res.1", border.size = NA)
ggsave("../../figures/EC4/EC4_neighbor.pdf", width = 8, height = 7)


# ---- feature heatmap
DoHeatmap(
  ec4.codex.obj,
  # feature = c(
  #   rownames(ec4.codex.obj)[18],
  #   rownames(ec4.codex.obj)[26],
  #   rownames(ec4.codex.obj)[35]
  # ),
  group.by = "Akoya_snn_res.1", slot = "scale.data"
) + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
ggsave("../../figures/EC4/EC4_heatmap.pdf", width = 20, height = 8)


# ---- stacked vlnplot
Idents(ec4.codex.obj) <- "Akoya_snn_res.1"
StackedVlnPlot(ec4.codex.obj, sort(rownames(ec4.codex.obj)), pt.size = 0, cols = my36colors)
ggsave("../../figures/EC4/EC4_stackedvlnplot.pdf", width = 40, height = 30)

# ---- cell type annotation
ec4.codex.obj$celltype <- "Undefined"
ec4.codex.obj$celltype[ec4.codex.obj$Akoya_snn_res.1 == "10"] <- "Tumor"
ec4.codex.obj$celltype[ec4.codex.obj$Akoya_snn_res.1 == "16"] <- "Tumor"
ec4.codex.obj$celltype[ec4.codex.obj$Akoya_snn_res.1 == "15"] <- "Tumor"
ec4.codex.obj$celltype[ec4.codex.obj$Akoya_snn_res.1 == "20"] <- "Fibro"
ec4.codex.obj$celltype[ec4.codex.obj$Akoya_snn_res.1 == "18"] <- "Fibro"
ec4.codex.obj$celltype[ec4.codex.obj$Akoya_snn_res.1 == "3"] <- "CD4T"
ec4.codex.obj$celltype[ec4.codex.obj$Akoya_snn_res.1 == "11"] <- "CD4T"
ec4.codex.obj$celltype[ec4.codex.obj$Akoya_snn_res.1 == "6"] <- "CD8T"
ec4.codex.obj$celltype[ec4.codex.obj$Akoya_snn_res.1 == "13"] <- "Plasma"
ec4.codex.obj$celltype[ec4.codex.obj$Akoya_snn_res.1 == "8"] <- "Endo"
ec4.codex.obj$celltype[ec4.codex.obj$Akoya_snn_res.1 == "12"] <- "Myeloid"

ec4.codex.obj$celltype <- factor(ec4.codex.obj$celltype, levels = sort(unique(ec4.codex.obj$celltype)))
ImageDimPlot(
  ec4.codex.obj,
  cols = c(
    "CD4T" = "#990F26", "CD8T" = "#B33E52", "DC/B" = "#CC7A88",
    "Endo" = "#E6B8BF", "Fibro" = "#99600F", "Hypoxia_cells" = "#B3823E",
    "Myeloid" = "#CCAA7A", "Plasma" = "#E6D2B8",
    "Tumor" = "#54990F", "Tumor_BCL6" = "#78B33E", "Tumor_CD24" = "#A3CC7A", "Tumor_Prolif" = "#CFE6B8",
    "Undefined" = "#6e6e6e"
  ),
  group.by = "celltype", border.size = NA
)
ggsave("../../figures/EC4/EC4_neighbor.pdf", width = 8, height = 7)
DimPlot(ec4.codex.obj,
  group.by = "celltype", label = TRUE, label.box = TRUE, raster = TRUE
)
ggsave("../../figures/EC4/EC4_celltype.pdf", width = 8, height = 7)



## EC5 (IA G2 EC p53wt Ciliated)
inputdata <- read_csv("../../data/EC_codex/EC5CA_1_measurements.csv")
head(inputdata)
dim(inputdata)
ec5.codex.obj <- LoadAkoya(
  filename = "../../data/EC_codex/EC5CA_1_measurements.csv",
  type = "qupath",
  fov = "EC5"
)
colnames(ec5.codex.obj@meta.data)[1:10]
ec5.codex.obj[["B2M"]] <- inputdata$"B2M-BX043: Cell: Mean"
ec5.codex.obj[["cellsize"]] <- ec5.codex.obj$"Cell..Area.µm.2"
ec5.codex.obj[["Sample"]] <- "EC5"

VlnPlot(ec5.codex.obj, features = c("B2M", "nCount_Akoya", "cellsize"), ncol = 3, pt.size = 0, log = TRUE) & geom_hline(yintercept = 30)
ggsave("../../figures/EC5/EC5_vlnplot_low.pdf", width = 8, height = 7)
VlnPlot(ec5.codex.obj, features = c("B2M", "nCount_Akoya", "cellsize"), ncol = 3, pt.size = 0, log = FALSE) & geom_hline(yintercept = 1000)
ggsave("../../figures/EC5/EC5_vlnplot_high.pdf", width = 8, height = 7)

ec5.codex.obj <- subset(ec5.codex.obj, subset = B2M > 2 & B2M < 100 & nCount_Akoya > 30 & nCount_Akoya < 1000 & cellsize > 4 & cellsize < 200)

ec5.codex.obj <- NormalizeData(object = ec5.codex.obj, normalization.method = "CLR", margin = 2)
ec5.codex.obj <- ScaleData(ec5.codex.obj)
VariableFeatures(ec5.codex.obj) <- rownames(ec5.codex.obj) # since the panel is small, treat all features as variable.
ec5.codex.obj <- RunPCA(object = ec5.codex.obj)
ec5.codex.obj <- RunUMAP(object = ec5.codex.obj, dims = 1:20)
ec5.codex.obj <- FindNeighbors(ec5.codex.obj, k.param = 20, dims = 1:20)
ec5.codex.obj <- FindClusters(object = ec5.codex.obj, resolution = seq(0.1, 2, 0.1), algorithm = 2, n.start = 1)

FeaturePlot(object = ec5.codex.obj, features = c("B2M", "nCount_Akoya", "cellsize"), raster = FALSE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
ggsave("../../figures/EC5/EC5_umap_qc.pdf", width = 8, height = 7)

ec5.codex.obj <- ScaleData(ec5.codex.obj, vars.to.regress = c("nCount_Akoya"))
VariableFeatures(ec5.codex.obj) <- rownames(ec5.codex.obj) # since the panel is small, treat all features as variable.
ec5.codex.obj <- RunPCA(object = ec5.codex.obj)

ec5.codex.obj <- JackStraw(ec5.codex.obj, num.replicate = 100, prop.freq = 0.1)
ec5.codex.obj <- ScoreJackStraw(ec5.codex.obj, dims = 1:20)
JackStrawPlot(ec5.codex.obj, dims = 1:20) & ElbowPlot(ec5.codex.obj, ndims = 20)
ggsave("../../figures/EC5/EC5_jackstraw.pdf", width = 15, height = 7)

ec5.codex.obj <- RunUMAP(object = ec5.codex.obj, dims = 1:20)
ec5.codex.obj <- FindNeighbors(ec5.codex.obj, k.param = 20, dims = 1:20)
ec5.codex.obj <- FindClusters(object = ec5.codex.obj, resolution = seq(0.1, 2, 0.1), algorithm = 2, n.start = 1)


# ---- seurat cluster
for (i in seq(0.1, 2, 0.1)) {
  print(i)
  DimPlot(ec5.codex.obj, label = TRUE, label.box = TRUE, raster = TRUE, group.by = paste0("Akoya_snn_res.", i))
  ggsave(paste0("../../figures/EC5/EC5_cluster_res_", i, ".pdf"), width = 8, height = 7)
}


# ---- seurat feature
for (gene in rownames(ec5.codex.obj)) {
  print(gene)
  FeaturePlot(object = ec5.codex.obj, features = gene, raster = FALSE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
  ggsave(paste0("../../figures/EC5/EC5_feature_", gene, ".pdf"), width = 8, height = 7)
}


# ---- codex image
ImageDimPlot(ec5.codex.obj, cols = "parade", group.by = "Akoya_snn_res.1", border.size = NA)
ggsave("../../figures/EC5/EC5_neighbor.pdf", width = 8, height = 7)


# ---- feature heatmap
DoHeatmap(
  ec5.codex.obj,
  # feature = c(
  #   rownames(ec5.codex.obj)[18],
  #   rownames(ec5.codex.obj)[26],
  #   rownames(ec5.codex.obj)[35]
  # ),
  group.by = "Akoya_snn_res.1", slot = "scale.data"
) + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
ggsave("../../figures/EC5/EC5_heatmap.pdf", width = 20, height = 8)


# ---- stacked vlnplot
Idents(ec5.codex.obj) <- "Akoya_snn_res.1"
StackedVlnPlot(ec5.codex.obj, sort(rownames(ec5.codex.obj)), pt.size = 0, cols = my36colors)
ggsave("../../figures/EC5/EC5_stackedvlnplot.pdf", width = 40, height = 30)

# ---- cell type annotation
ec5.codex.obj$celltype <- "Undefined"
ec5.codex.obj$celltype[ec5.codex.obj$Akoya_snn_res.1 == "10"] <- "Tumor"
ec5.codex.obj$celltype[ec5.codex.obj$Akoya_snn_res.1 == "16"] <- "Tumor"
ec5.codex.obj$celltype[ec5.codex.obj$Akoya_snn_res.1 == "15"] <- "Tumor"
ec5.codex.obj$celltype[ec5.codex.obj$Akoya_snn_res.1 == "20"] <- "Fibro"
ec5.codex.obj$celltype[ec5.codex.obj$Akoya_snn_res.1 == "18"] <- "Fibro"
ec5.codex.obj$celltype[ec5.codex.obj$Akoya_snn_res.1 == "3"] <- "CD4T"
ec5.codex.obj$celltype[ec5.codex.obj$Akoya_snn_res.1 == "11"] <- "CD4T"
ec5.codex.obj$celltype[ec5.codex.obj$Akoya_snn_res.1 == "6"] <- "CD8T"
ec5.codex.obj$celltype[ec5.codex.obj$Akoya_snn_res.1 == "13"] <- "Plasma"
ec5.codex.obj$celltype[ec5.codex.obj$Akoya_snn_res.1 == "8"] <- "Endo"
ec5.codex.obj$celltype[ec5.codex.obj$Akoya_snn_res.1 == "12"] <- "Myeloid"

ec5.codex.obj$celltype <- factor(ec5.codex.obj$celltype, levels = sort(unique(ec5.codex.obj$celltype)))
ImageDimPlot(
  ec5.codex.obj,
  cols = c(
    "CD4T" = "#990F26", "CD8T" = "#B33E52", "DC/B" = "#CC7A88",
    "Endo" = "#E6B8BF", "Fibro" = "#99600F", "Hypoxia_cells" = "#B3823E",
    "Myeloid" = "#CCAA7A", "Plasma" = "#E6D2B8",
    "Tumor" = "#54990F", "Tumor_BCL6" = "#78B33E", "Tumor_CD24" = "#A3CC7A", "Tumor_Prolif" = "#CFE6B8",
    "Undefined" = "#6e6e6e"
  ),
  group.by = "celltype", border.size = NA
)
ggsave("../../figures/EC5/EC5_neighbor.pdf", width = 8, height = 7)
DimPlot(ec5.codex.obj,
  group.by = "celltype", label = TRUE, label.box = TRUE, raster = TRUE
)
ggsave("../../figures/EC5/EC5_celltype.pdf", width = 8, height = 7)



## EC7 (IB G2 EC p53wt Glandular)
inputdata <- read_csv("../../data/EC_codex/EC7CA_1_measurements.csv")
head(inputdata)
dim(inputdata)
ec7.codex.obj <- LoadAkoya(
  filename = "../../data/EC_codex/EC7CA_1_measurements.csv",
  type = "qupath",
  fov = "EC7"
)
colnames(ec7.codex.obj@meta.data)[1:10]
ec7.codex.obj[["B2M"]] <- inputdata$"B2M-BX043: Cell: Mean"
ec7.codex.obj[["cellsize"]] <- ec7.codex.obj$"Cell..Area.µm.2"
ec7.codex.obj[["Sample"]] <- "EC7"

VlnPlot(ec7.codex.obj, features = c("B2M", "nCount_Akoya", "cellsize"), ncol = 3, pt.size = 0, log = TRUE) & geom_hline(yintercept = 3)
ggsave("../../figures/EC7/EC7_vlnplot_low.pdf", width = 8, height = 7)
VlnPlot(ec7.codex.obj, features = c("B2M", "nCount_Akoya", "cellsize"), ncol = 3, pt.size = 0, log = FALSE) & geom_hline(yintercept = 100)
ggsave("../../figures/EC7/EC7_vlnplot_high.pdf", width = 8, height = 7)

ec7.codex.obj <- subset(ec7.codex.obj, subset = B2M > 3 & B2M < 100 & nCount_Akoya > 50 & nCount_Akoya < 1200 & cellsize > 3 & cellsize < 200)

ec7.codex.obj <- NormalizeData(object = ec7.codex.obj, normalization.method = "CLR", margin = 2)
ec7.codex.obj <- ScaleData(ec7.codex.obj)
VariableFeatures(ec7.codex.obj) <- rownames(ec7.codex.obj) # since the panel is small, treat all features as variable.
ec7.codex.obj <- RunPCA(object = ec7.codex.obj)
ec7.codex.obj <- RunUMAP(object = ec7.codex.obj, dims = 1:20)
ec7.codex.obj <- FindNeighbors(ec7.codex.obj, k.param = 20, dims = 1:20)
ec7.codex.obj <- FindClusters(object = ec7.codex.obj, resolution = seq(0.1, 2, 0.1), algorithm = 2, n.start = 1)

FeaturePlot(object = ec7.codex.obj, features = c("B2M", "nCount_Akoya", "cellsize"), raster = FALSE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
ggsave("../../figures/EC7/EC7_umap_qc.pdf", width = 8, height = 7)

ec7.codex.obj <- ScaleData(ec7.codex.obj, vars.to.regress = c("nCount_Akoya"))
VariableFeatures(ec7.codex.obj) <- rownames(ec7.codex.obj) # since the panel is small, treat all features as variable.
ec7.codex.obj <- RunPCA(object = ec7.codex.obj)

ec7.codex.obj <- JackStraw(ec7.codex.obj, num.replicate = 100, prop.freq = 0.01)
ec7.codex.obj <- ScoreJackStraw(ec7.codex.obj, dims = 1:20)
JackStrawPlot(ec7.codex.obj, dims = 1:20) & ElbowPlot(ec7.codex.obj, ndims = 20)
ggsave("../../figures/EC7/EC7_jackstraw.pdf", width = 15, height = 7)

ec7.codex.obj <- RunUMAP(object = ec7.codex.obj, dims = 1:20)
ec7.codex.obj <- FindNeighbors(ec7.codex.obj, k.param = 20, dims = 1:20)
ec7.codex.obj <- FindClusters(object = ec7.codex.obj, resolution = seq(0.1, 2, 0.1), algorithm = 2, n.start = 1)


# ---- seurat cluster
for (i in seq(0.1, 2, 0.1)) {
  print(i)
  DimPlot(ec7.codex.obj, label = TRUE, label.box = TRUE, raster = TRUE, group.by = paste0("Akoya_snn_res.", i))
  ggsave(paste0("../../figures/EC7/EC7_cluster_res_", i, ".pdf"), width = 8, height = 7)
}


# ---- seurat feature
for (gene in rownames(ec7.codex.obj)) {
  print(gene)
  FeaturePlot(object = ec7.codex.obj, features = gene, raster = FALSE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
  ggsave(paste0("../../figures/EC7/EC7_feature_", gene, ".pdf"), width = 8, height = 7)
}


# ---- codex image
ImageDimPlot(ec7.codex.obj, cols = "parade", group.by = "Akoya_snn_res.1", border.size = NA)
ggsave("../../figures/EC7/EC7_neighbor.pdf", width = 8, height = 7)


# ---- feature heatmap
DoHeatmap(
  ec7.codex.obj,
  # feature = c(
  #   rownames(ec7.codex.obj)[18],
  #   rownames(ec7.codex.obj)[26],
  #   rownames(ec7.codex.obj)[35]
  # ),
  group.by = "Akoya_snn_res.1", slot = "scale.data"
) + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
ggsave("../../figures/EC7/EC7_heatmap.pdf", width = 20, height = 8)


# ---- stacked vlnplot
Idents(ec7.codex.obj) <- "Akoya_snn_res.1"
StackedVlnPlot(ec7.codex.obj, sort(rownames(ec7.codex.obj)), pt.size = 0, cols = my36colors)
ggsave("../../figures/EC7/EC7_stackedvlnplot.pdf", width = 40, height = 30)

# ---- cell type annotation
ec7.codex.obj$celltype <- "Undefined"
ec7.codex.obj$celltype[ec7.codex.obj$Akoya_snn_res.1 == "10"] <- "Tumor"
ec7.codex.obj$celltype[ec7.codex.obj$Akoya_snn_res.1 == "16"] <- "Tumor"
ec7.codex.obj$celltype[ec7.codex.obj$Akoya_snn_res.1 == "15"] <- "Tumor"
ec7.codex.obj$celltype[ec7.codex.obj$Akoya_snn_res.1 == "20"] <- "Fibro"
ec7.codex.obj$celltype[ec7.codex.obj$Akoya_snn_res.1 == "18"] <- "Fibro"
ec7.codex.obj$celltype[ec7.codex.obj$Akoya_snn_res.1 == "3"] <- "CD4T"
ec7.codex.obj$celltype[ec7.codex.obj$Akoya_snn_res.1 == "11"] <- "CD4T"
ec7.codex.obj$celltype[ec7.codex.obj$Akoya_snn_res.1 == "6"] <- "CD8T"
ec7.codex.obj$celltype[ec7.codex.obj$Akoya_snn_res.1 == "13"] <- "Plasma"
ec7.codex.obj$celltype[ec7.codex.obj$Akoya_snn_res.1 == "8"] <- "Endo"
ec7.codex.obj$celltype[ec7.codex.obj$Akoya_snn_res.1 == "12"] <- "Myeloid"

ec7.codex.obj$celltype <- factor(ec7.codex.obj$celltype, levels = sort(unique(ec7.codex.obj$celltype)))
ImageDimPlot(
  ec7.codex.obj,
  cols = c(
    "CD4T" = "#990F26", "CD8T" = "#B33E52", "DC/B" = "#CC7A88",
    "Endo" = "#E6B8BF", "Fibro" = "#99600F", "Hypoxia_cells" = "#B3823E",
    "Myeloid" = "#CCAA7A", "Plasma" = "#E6D2B8",
    "Tumor" = "#54990F", "Tumor_BCL6" = "#78B33E", "Tumor_CD24" = "#A3CC7A", "Tumor_Prolif" = "#CFE6B8",
    "Undefined" = "#6e6e6e"
  ),
  group.by = "celltype", border.size = NA
)
ggsave("../../figures/EC7/EC7_neighbor.pdf", width = 8, height = 7)
DimPlot(ec7.codex.obj,
  group.by = "celltype", label = TRUE, label.box = TRUE, raster = TRUE
)
ggsave("../../figures/EC7/EC7_celltype.pdf", width = 8, height = 7)



## EC10 (II G3 CC p53abn Luminal)
inputdata <- read_csv("../../data/EC_codex/EC10CA_1_measurements.csv")
head(inputdata)
dim(inputdata)
ec10.codex.obj <- LoadAkoya(
  filename = "../../data/EC_codex/EC10CA_1_measurements.csv",
  type = "qupath",
  fov = "EC10"
)
colnames(ec10.codex.obj@meta.data)[1:10]
ec10.codex.obj[["B2M"]] <- inputdata$"B2M-BX043: Cell: Mean"
ec10.codex.obj[["cellsize"]] <- ec10.codex.obj$"Cell..Area.µm.2"
ec10.codex.obj[["Sample"]] <- "EC10"

VlnPlot(ec10.codex.obj, features = c("B2M", "nCount_Akoya", "cellsize"), ncol = 3, pt.size = 0, log = TRUE) & geom_hline(yintercept = 2)
ggsave("../../figures/EC10/EC10_vlnplot_low.pdf", width = 8, height = 7)
VlnPlot(ec10.codex.obj, features = c("B2M", "nCount_Akoya", "cellsize"), ncol = 3, pt.size = 0, log = FALSE) & geom_hline(yintercept = 100)
ggsave("../../figures/EC10/EC10_vlnplot_high.pdf", width = 8, height = 7)

ec10.codex.obj <- subset(ec10.codex.obj, subset = B2M > 2 & B2M < 100 & nCount_Akoya > 30 & nCount_Akoya < 1200 & cellsize > 4 & cellsize < 200)

ec10.codex.obj <- NormalizeData(object = ec10.codex.obj, normalization.method = "CLR", margin = 2)
ec10.codex.obj <- ScaleData(ec10.codex.obj)
VariableFeatures(ec10.codex.obj) <- rownames(ec10.codex.obj) # since the panel is small, treat all features as variable.
ec10.codex.obj <- RunPCA(object = ec10.codex.obj)
ec10.codex.obj <- RunUMAP(object = ec10.codex.obj, dims = 1:20)
ec10.codex.obj <- FindNeighbors(ec10.codex.obj, k.param = 20, dims = 1:20)
ec10.codex.obj <- FindClusters(object = ec10.codex.obj, resolution = seq(0.1, 2, 0.1), algorithm = 2, n.start = 1)

FeaturePlot(object = ec10.codex.obj, features = c("B2M", "nCount_Akoya", "cellsize"), raster = FALSE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
ggsave("../../figures/EC10/EC10_umap_qc.pdf", width = 8, height = 7)

ec10.codex.obj <- ScaleData(ec10.codex.obj, vars.to.regress = c("nCount_Akoya"))
VariableFeatures(ec10.codex.obj) <- rownames(ec10.codex.obj) # since the panel is small, treat all features as variable.
ec10.codex.obj <- RunPCA(object = ec10.codex.obj)

ec10.codex.obj <- JackStraw(ec10.codex.obj, num.replicate = 100, prop.freq = 0.01)
ec10.codex.obj <- ScoreJackStraw(ec10.codex.obj, dims = 1:20)
JackStrawPlot(ec10.codex.obj, dims = 1:20) & ElbowPlot(ec10.codex.obj, ndims = 20)
ggsave("../../figures/EC10/EC10_jackstraw.pdf", width = 15, height = 7)

ec10.codex.obj <- RunUMAP(object = ec10.codex.obj, dims = 1:20)
ec10.codex.obj <- FindNeighbors(ec10.codex.obj, k.param = 20, dims = 1:20)
ec10.codex.obj <- FindClusters(object = ec10.codex.obj, resolution = seq(0.1, 2, 0.1), algorithm = 2, n.start = 1)


# ---- seurat cluster
for (i in seq(0.1, 2, 0.1)) {
  print(i)
  DimPlot(ec10.codex.obj, label = TRUE, label.box = TRUE, raster = TRUE, group.by = paste0("Akoya_snn_res.", i))
  ggsave(paste0("../../figures/EC10/EC10_cluster_res_", i, ".pdf"), width = 8, height = 7)
}


# ---- seurat feature
for (gene in rownames(ec10.codex.obj)) {
  print(gene)
  FeaturePlot(object = ec10.codex.obj, features = gene, raster = FALSE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
  ggsave(paste0("../../figures/EC10/EC10_feature_", gene, ".pdf"), width = 8, height = 7)
}


# ---- codex image
ImageDimPlot(ec10.codex.obj, cols = "parade", group.by = "Akoya_snn_res.1", border.size = NA)
ggsave("../../figures/EC10/EC10_neighbor.pdf", width = 8, height = 7)


# ---- feature heatmap
DoHeatmap(
  ec10.codex.obj,
  # feature = c(
  #   rownames(ec10.codex.obj)[18],
  #   rownames(ec10.codex.obj)[26],
  #   rownames(ec10.codex.obj)[35]
  # ),
  group.by = "Akoya_snn_res.1", slot = "scale.data"
) + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
ggsave("../../figures/EC10/EC10_heatmap.pdf", width = 20, height = 8)


# ---- stacked vlnplot
Idents(ec10.codex.obj) <- "Akoya_snn_res.1"
StackedVlnPlot(ec10.codex.obj, sort(rownames(ec10.codex.obj)), pt.size = 0, cols = my36colors)
ggsave("../../figures/EC10/EC10_stackedvlnplot.pdf", width = 40, height = 30)

# ---- cell type annotation
ec10.codex.obj$celltype <- "Undefined"
ec10.codex.obj$celltype[ec10.codex.obj$Akoya_snn_res.1 == "10"] <- "Tumor"
ec10.codex.obj$celltype[ec10.codex.obj$Akoya_snn_res.1 == "16"] <- "Tumor"
ec10.codex.obj$celltype[ec10.codex.obj$Akoya_snn_res.1 == "15"] <- "Tumor"
ec10.codex.obj$celltype[ec10.codex.obj$Akoya_snn_res.1 == "20"] <- "Fibro"
ec10.codex.obj$celltype[ec10.codex.obj$Akoya_snn_res.1 == "18"] <- "Fibro"
ec10.codex.obj$celltype[ec10.codex.obj$Akoya_snn_res.1 == "3"] <- "CD4T"
ec10.codex.obj$celltype[ec10.codex.obj$Akoya_snn_res.1 == "11"] <- "CD4T"
ec10.codex.obj$celltype[ec10.codex.obj$Akoya_snn_res.1 == "6"] <- "CD8T"
ec10.codex.obj$celltype[ec10.codex.obj$Akoya_snn_res.1 == "13"] <- "Plasma"
ec10.codex.obj$celltype[ec10.codex.obj$Akoya_snn_res.1 == "8"] <- "Endo"
ec10.codex.obj$celltype[ec10.codex.obj$Akoya_snn_res.1 == "12"] <- "Myeloid"

ec10.codex.obj$celltype <- factor(ec10.codex.obj$celltype, levels = sort(unique(ec10.codex.obj$celltype)))
ImageDimPlot(
  ec10.codex.obj,
  cols = c(
    "CD4T" = "#990F26", "CD8T" = "#B33E52", "DC/B" = "#CC7A88",
    "Endo" = "#E6B8BF", "Fibro" = "#99600F", "Hypoxia_cells" = "#B3823E",
    "Myeloid" = "#CCAA7A", "Plasma" = "#E6D2B8",
    "Tumor" = "#54990F", "Tumor_BCL6" = "#78B33E", "Tumor_CD24" = "#A3CC7A", "Tumor_Prolif" = "#CFE6B8",
    "Undefined" = "#6e6e6e"
  ),
  group.by = "celltype", border.size = NA
)
ggsave("../../figures/EC10/EC10_neighbor.pdf", width = 8, height = 7)
DimPlot(ec10.codex.obj,
  group.by = "celltype", label = TRUE, label.box = TRUE, raster = TRUE
)
ggsave("../../figures/EC10/EC10_celltype.pdf", width = 8, height = 7)



## EC11 (IIIC2 G3 EC p53wt EMT-like)
inputdata <- read_csv("../../data/EC_codex/EC11CA_1_measurements.csv")
head(inputdata)
dim(inputdata)
ec11.codex.obj <- LoadAkoya(
  filename = "../../data/EC_codex/EC11CA_1_measurements.csv",
  type = "qupath",
  fov = "EC11"
)
colnames(ec11.codex.obj@meta.data)[1:10]
ec11.codex.obj[["B2M"]] <- inputdata$"B2M-BX043: Cell: Mean"
ec11.codex.obj[["cellsize"]] <- ec11.codex.obj$"Cell..Area.µm.2"
ec11.codex.obj[["Sample"]] <- "EC11"

VlnPlot(ec11.codex.obj, features = c("B2M", "nCount_Akoya", "cellsize"), ncol = 3, pt.size = 0, log = TRUE) & geom_hline(yintercept = 3)
ggsave("../../figures/EC11/EC11_vlnplot_low.pdf", width = 8, height = 7)
VlnPlot(ec11.codex.obj, features = c("B2M", "nCount_Akoya", "cellsize"), ncol = 3, pt.size = 0, log = FALSE) & geom_hline(yintercept = 100)
ggsave("../../figures/EC11/EC11_vlnplot_high.pdf", width = 8, height = 7)

ec11.codex.obj <- subset(ec11.codex.obj, subset = B2M > 3 & B2M < 100 & nCount_Akoya > 70 & nCount_Akoya < 1200 & cellsize > 5 & cellsize < 200)

ec11.codex.obj <- NormalizeData(object = ec11.codex.obj, normalization.method = "CLR", margin = 2)
ec11.codex.obj <- ScaleData(ec11.codex.obj)
VariableFeatures(ec11.codex.obj) <- rownames(ec11.codex.obj) # since the panel is small, treat all features as variable.
ec11.codex.obj <- RunPCA(object = ec11.codex.obj)
ec11.codex.obj <- RunUMAP(object = ec11.codex.obj, dims = 1:20)
ec11.codex.obj <- FindNeighbors(ec11.codex.obj, k.param = 20, dims = 1:20)
ec11.codex.obj <- FindClusters(object = ec11.codex.obj, resolution = seq(0.1, 2, 0.1), algorithm = 2, n.start = 1)

FeaturePlot(object = ec11.codex.obj, features = c("B2M", "nCount_Akoya", "cellsize"), raster = FALSE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
ggsave("../../figures/EC11/EC11_umap_qc.pdf", width = 8, height = 7)

ec11.codex.obj <- ScaleData(ec11.codex.obj, vars.to.regress = c("nCount_Akoya"))
VariableFeatures(ec11.codex.obj) <- rownames(ec11.codex.obj) # since the panel is small, treat all features as variable.
ec11.codex.obj <- RunPCA(object = ec11.codex.obj)

ec11.codex.obj <- JackStraw(ec11.codex.obj, num.replicate = 100, prop.freq = 0.01)
ec11.codex.obj <- ScoreJackStraw(ec11.codex.obj, dims = 1:20)
JackStrawPlot(ec11.codex.obj, dims = 1:20) & ElbowPlot(ec11.codex.obj, ndims = 20)
ggsave("../../figures/EC11/EC11_jackstraw.pdf", width = 15, height = 7)

ec11.codex.obj <- RunUMAP(object = ec11.codex.obj, dims = 1:20)
ec11.codex.obj <- FindNeighbors(ec11.codex.obj, k.param = 20, dims = 1:20)
ec11.codex.obj <- FindClusters(object = ec11.codex.obj, resolution = seq(0.1, 2, 0.1), algorithm = 2, n.start = 1)


# ---- seurat cluster
for (i in seq(0.1, 2, 0.1)) {
  print(i)
  DimPlot(ec11.codex.obj, label = TRUE, label.box = TRUE, raster = TRUE, group.by = paste0("Akoya_snn_res.", i))
  ggsave(paste0("../../figures/EC11/EC11_cluster_res_", i, ".pdf"), width = 8, height = 7)
}


# ---- seurat feature
for (gene in rownames(ec11.codex.obj)) {
  print(gene)
  FeaturePlot(object = ec11.codex.obj, features = gene, raster = FALSE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
  ggsave(paste0("../../figures/EC11/EC11_feature_", gene, ".pdf"), width = 8, height = 7)
}


# ---- codex image
ImageDimPlot(ec11.codex.obj, cols = "parade", group.by = "Akoya_snn_res.1", border.size = NA)
ggsave("../../figures/EC11/EC11_neighbor.pdf", width = 8, height = 7)


# ---- feature heatmap
DoHeatmap(
  ec11.codex.obj,
  # feature = c(
  #   rownames(ec11.codex.obj)[18],
  #   rownames(ec11.codex.obj)[26],
  #   rownames(ec11.codex.obj)[35]
  # ),
  group.by = "Akoya_snn_res.1", slot = "scale.data"
) + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
ggsave("../../figures/EC11/EC11_heatmap.pdf", width = 20, height = 8)


# ---- stacked vlnplot
Idents(ec11.codex.obj) <- "Akoya_snn_res.1"
StackedVlnPlot(ec11.codex.obj, sort(rownames(ec11.codex.obj)), pt.size = 0, cols = my36colors)
ggsave("../../figures/EC11/EC11_stackedvlnplot.pdf", width = 40, height = 30)

# ---- cell type annotation
ec11.codex.obj$celltype <- "Undefined"
ec11.codex.obj$celltype[ec11.codex.obj$Akoya_snn_res.1 == "10"] <- "Tumor"
ec11.codex.obj$celltype[ec11.codex.obj$Akoya_snn_res.1 == "16"] <- "Tumor"
ec11.codex.obj$celltype[ec11.codex.obj$Akoya_snn_res.1 == "15"] <- "Tumor"
ec11.codex.obj$celltype[ec11.codex.obj$Akoya_snn_res.1 == "20"] <- "Fibro"
ec11.codex.obj$celltype[ec11.codex.obj$Akoya_snn_res.1 == "18"] <- "Fibro"
ec11.codex.obj$celltype[ec11.codex.obj$Akoya_snn_res.1 == "3"] <- "CD4T"
ec11.codex.obj$celltype[ec11.codex.obj$Akoya_snn_res.1 == "11"] <- "CD4T"
ec11.codex.obj$celltype[ec11.codex.obj$Akoya_snn_res.1 == "6"] <- "CD8T"
ec11.codex.obj$celltype[ec11.codex.obj$Akoya_snn_res.1 == "13"] <- "Plasma"
ec11.codex.obj$celltype[ec11.codex.obj$Akoya_snn_res.1 == "8"] <- "Endo"
ec11.codex.obj$celltype[ec11.codex.obj$Akoya_snn_res.1 == "12"] <- "Myeloid"

ec11.codex.obj$celltype <- factor(ec11.codex.obj$celltype, levels = sort(unique(ec11.codex.obj$celltype)))
ImageDimPlot(
  ec11.codex.obj,
  cols = c(
    "CD4T" = "#990F26", "CD8T" = "#B33E52", "DC/B" = "#CC7A88",
    "Endo" = "#E6B8BF", "Fibro" = "#99600F", "Hypoxia_cells" = "#B3823E",
    "Myeloid" = "#CCAA7A", "Plasma" = "#E6D2B8",
    "Tumor" = "#54990F", "Tumor_BCL6" = "#78B33E", "Tumor_CD24" = "#A3CC7A", "Tumor_Prolif" = "#CFE6B8",
    "Undefined" = "#6e6e6e"
  ),
  group.by = "celltype", border.size = NA
)
ggsave("../../figures/EC11/EC11_neighbor.pdf", width = 8, height = 7)
DimPlot(ec11.codex.obj,
  group.by = "celltype", label = TRUE, label.box = TRUE, raster = TRUE
)
ggsave("../../figures/EC11/EC11_celltype.pdf", width = 8, height = 7)



## EC12 (IB G2 EC MMRd Glandular)
inputdata <- read_csv("../../data/EC_codex/EC12CA_1_measurements.csv")
head(inputdata)
dim(inputdata)
ec12.codex.obj <- LoadAkoya(
  filename = "../../data/EC_codex/EC12CA_1_measurements.csv",
  type = "qupath",
  fov = "EC12"
)
colnames(ec12.codex.obj@meta.data)[1:10]
ec12.codex.obj[["B2M"]] <- inputdata$"B2M-BX043: Cell: Mean"
ec12.codex.obj[["cellsize"]] <- ec12.codex.obj$"Cell..Area.µm.2"
ec12.codex.obj[["Sample"]] <- "EC12"

VlnPlot(ec12.codex.obj, features = c("B2M", "nCount_Akoya", "cellsize"), ncol = 3, pt.size = 0, log = TRUE) & geom_hline(yintercept = 1)
ggsave("../../figures/EC12/EC12_vlnplot_low.pdf", width = 8, height = 7)
VlnPlot(ec12.codex.obj, features = c("B2M", "nCount_Akoya", "cellsize"), ncol = 3, pt.size = 0, log = FALSE) & geom_hline(yintercept = 100)
ggsave("../../figures/EC12/EC12_vlnplot_high.pdf", width = 8, height = 7)

ec12.codex.obj <- subset(ec12.codex.obj, subset = B2M > 1 & B2M < 100 & nCount_Akoya > 30 & nCount_Akoya < 1200 & cellsize > 3 & cellsize < 200)

ec12.codex.obj <- NormalizeData(object = ec12.codex.obj, normalization.method = "CLR", margin = 2)
ec12.codex.obj <- ScaleData(ec12.codex.obj)
VariableFeatures(ec12.codex.obj) <- rownames(ec12.codex.obj) # since the panel is small, treat all features as variable.
ec12.codex.obj <- RunPCA(object = ec12.codex.obj)
ec12.codex.obj <- RunUMAP(object = ec12.codex.obj, dims = 1:20)
ec12.codex.obj <- FindNeighbors(ec12.codex.obj, k.param = 20, dims = 1:20)
ec12.codex.obj <- FindClusters(object = ec12.codex.obj, resolution = seq(0.1, 2, 0.1), algorithm = 2, n.start = 1)

FeaturePlot(object = ec12.codex.obj, features = c("B2M", "nCount_Akoya", "cellsize"), raster = FALSE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
ggsave("../../figures/EC12/EC12_umap_qc.pdf", width = 8, height = 7)

ec12.codex.obj <- ScaleData(ec12.codex.obj, vars.to.regress = c("nCount_Akoya"))
VariableFeatures(ec12.codex.obj) <- rownames(ec12.codex.obj) # since the panel is small, treat all features as variable.
ec12.codex.obj <- RunPCA(object = ec12.codex.obj)

ec12.codex.obj <- JackStraw(ec12.codex.obj, num.replicate = 100, prop.freq = 0.01)
ec12.codex.obj <- ScoreJackStraw(ec12.codex.obj, dims = 1:20)
JackStrawPlot(ec12.codex.obj, dims = 1:20) & ElbowPlot(ec12.codex.obj, ndims = 20)
ggsave("../../figures/EC12/EC12_jackstraw.pdf", width = 15, height = 7)

ec12.codex.obj <- RunUMAP(object = ec12.codex.obj, dims = 1:20)
ec12.codex.obj <- FindNeighbors(ec12.codex.obj, k.param = 20, dims = 1:20)
ec12.codex.obj <- FindClusters(object = ec12.codex.obj, resolution = seq(0.1, 2, 0.1), algorithm = 2, n.start = 1)


# ---- seurat cluster
for (i in seq(0.1, 2, 0.1)) {
  print(i)
  DimPlot(ec12.codex.obj, label = TRUE, label.box = TRUE, raster = TRUE, group.by = paste0("Akoya_snn_res.", i))
  ggsave(paste0("../../figures/EC12/EC12_cluster_res_", i, ".pdf"), width = 8, height = 7)
}


# ---- seurat feature
for (gene in rownames(ec12.codex.obj)) {
  print(gene)
  FeaturePlot(object = ec12.codex.obj, features = gene, raster = FALSE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
  ggsave(paste0("../../figures/EC12/EC12_feature_", gene, ".pdf"), width = 8, height = 7)
}


# ---- codex image
ImageDimPlot(ec12.codex.obj, cols = "parade", group.by = "Akoya_snn_res.1", border.size = NA)
ggsave("../../figures/EC12/EC12_neighbor.pdf", width = 8, height = 7)


# ---- feature heatmap
DoHeatmap(
  ec12.codex.obj,
  # feature = c(
  #   rownames(ec12.codex.obj)[18],
  #   rownames(ec12.codex.obj)[26],
  #   rownames(ec12.codex.obj)[35]
  # ),
  group.by = "Akoya_snn_res.1", slot = "scale.data"
) + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
ggsave("../../figures/EC12/EC12_heatmap.pdf", width = 20, height = 8)


# ---- stacked vlnplot
Idents(ec12.codex.obj) <- "Akoya_snn_res.1"
StackedVlnPlot(ec12.codex.obj, sort(rownames(ec12.codex.obj)), pt.size = 0, cols = my36colors)
ggsave("../../figures/EC12/EC12_stackedvlnplot.pdf", width = 40, height = 30)

# ---- cell type annotation
ec12.codex.obj$celltype <- "Undefined"
ec12.codex.obj$celltype[ec12.codex.obj$Akoya_snn_res.1 == "10"] <- "Tumor"
ec12.codex.obj$celltype[ec12.codex.obj$Akoya_snn_res.1 == "16"] <- "Tumor"
ec12.codex.obj$celltype[ec12.codex.obj$Akoya_snn_res.1 == "15"] <- "Tumor"
ec12.codex.obj$celltype[ec12.codex.obj$Akoya_snn_res.1 == "20"] <- "Fibro"
ec12.codex.obj$celltype[ec12.codex.obj$Akoya_snn_res.1 == "18"] <- "Fibro"
ec12.codex.obj$celltype[ec12.codex.obj$Akoya_snn_res.1 == "3"] <- "CD4T"
ec12.codex.obj$celltype[ec12.codex.obj$Akoya_snn_res.1 == "11"] <- "CD4T"
ec12.codex.obj$celltype[ec12.codex.obj$Akoya_snn_res.1 == "6"] <- "CD8T"
ec12.codex.obj$celltype[ec12.codex.obj$Akoya_snn_res.1 == "13"] <- "Plasma"
ec12.codex.obj$celltype[ec12.codex.obj$Akoya_snn_res.1 == "8"] <- "Endo"
ec12.codex.obj$celltype[ec12.codex.obj$Akoya_snn_res.1 == "12"] <- "Myeloid"

ec12.codex.obj$celltype <- factor(ec12.codex.obj$celltype, levels = sort(unique(ec12.codex.obj$celltype)))
ImageDimPlot(
  ec12.codex.obj,
  cols = c(
    "CD4T" = "#990F26", "CD8T" = "#B33E52", "DC/B" = "#CC7A88",
    "Endo" = "#E6B8BF", "Fibro" = "#99600F", "Hypoxia_cells" = "#B3823E",
    "Myeloid" = "#CCAA7A", "Plasma" = "#E6D2B8",
    "Tumor" = "#54990F", "Tumor_BCL6" = "#78B33E", "Tumor_CD24" = "#A3CC7A", "Tumor_Prolif" = "#CFE6B8",
    "Undefined" = "#6e6e6e"
  ),
  group.by = "celltype", border.size = NA
)
ggsave("../../figures/EC12/EC12_neighbor.pdf", width = 8, height = 7)
DimPlot(ec12.codex.obj,
  group.by = "celltype", label = TRUE, label.box = TRUE, raster = TRUE
)
ggsave("../../figures/EC12/EC12_celltype.pdf", width = 8, height = 7)





# Integration
codex.obj <- merge(
  x = ec4.codex.obj,
  y = c(ec5.codex.obj, ec7.codex.obj, ec10.codex.obj, ec11.codex.obj, ec12.codex.obj),
  add.cell.ids = c("EC4", "EC5", "EC7", "EC10", "EC11", "EC12"),
  project = "EC_CODEX"
)
saveRDS(codex.obj, "./codex_object.rds")





# ---- split samples
sample.list <- list()
groupings <- unique(codex.obj$Sample)
for (i in groupings) {
  # i <- "EC4"
  cells <- which(x = codex.obj[["Sample", drop = TRUE]] == i)
  cells <- colnames(x = codex.obj)[cells]
  sample.list[[i]] <- subset_opt(object = codex.obj, cells = cells)
}

sample.list <- purrr::map(sample.list, function(x) {
  x <- NormalizeData(object = x, normalization.method = "CLR", margin = 2)
  VariableFeatures(x) <- rownames(x)
  x <- ScaleData(x, features = VariableFeatures(x), vars.to.regress = c("nCount_Akoya"))
  x <- RunPCA(x, features = VariableFeatures(x), npcs = 30, verbose = FALSE)
})
save(sample.list, file = "./sample_list.rda")





# ---- integrate samples
# load("./sample_list.rda")
# codex.obj <- readRDS("./codex_object.rds")
plan("multicore", workers = 20)
plan()
options(future.globals.maxSize = 40 * 1024^3)
gc()
plan()
sample.anchors <- FindIntegrationAnchors(
  sample.list,
  reduction = "rpca",
  dims = 1:20,
  anchor.features = rownames(codex.obj),
  scale = FALSE
)
sample.integrated <- IntegrateData(
  anchorset = sample.anchors,
  dims = 1:20
)
saveRDS(sample.integrated, file = "./sample_integrated.rds")

# sample.integrated <- readRDS("./sample_integrated.rds")
sample.integrated <- ScaleData(sample.integrated)
sample.integrated <- RunPCA(sample.integrated, npcs = 30)
sample.integrated <- FindNeighbors(sample.integrated, k.param = 25, dims = 1:20)
sample.integrated <- FindClusters(object = sample.integrated, resolution = seq(0.1, 3, 0.1), algorithm = 2, n.start = 1)
sample.integrated <- RunUMAP(sample.integrated, dims = 1:20)
saveRDS(sample.integrated, file = "./sample_integrated.rds")

colnames(sample.integrated@meta.data)





# ---- seurat metadata
FeaturePlot(sample.integrated, features = c("B2M", "nCount_Akoya", "cellsize"), order = T, raster = F) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
ggsave("../../figures_metadata.jpg", width = 8, height = 7)
ggsave("../../figures_metadata.jpg", width = 8, height = 7)





# ---- annotation
sample.integrated$celltype <- "Undefined"
sample.integrated$celltype[sample.integrated$integrated_snn_res.3 %in% c(0, 1, 15, 17, 21, 25, 29, 37, 56, 77, 78, 80)] <- "Tumor"
sample.integrated$celltype[sample.integrated$integrated_snn_res.3 %in% c(2, 8, 14, 19, 22, 24, 31, 38, 39, 42, 47, 49, 58, 66, 72, 73, 75, 85, 87)] <- "Fibro"
sample.integrated$celltype[sample.integrated$integrated_snn_res.3 %in% c(3, 48, 54, 63)] <- "CD4T"
sample.integrated$celltype[sample.integrated$integrated_snn_res.3 %in% c(4, 18, 55, 57, 61, 65)] <- "Endo"
sample.integrated$celltype[sample.integrated$integrated_snn_res.3 %in% c(5, 10, 30, 34, 40, 43)] <- "Tumor_Prolif"
sample.integrated$celltype[sample.integrated$integrated_snn_res.3 %in% c(6, 12, 13, 52, 59)] <- "Tumor"
sample.integrated$celltype[sample.integrated$integrated_snn_res.3 %in% c(7, 20, 45, 51, 53, 64, 68, 84)] <- "Myeloid"
sample.integrated$celltype[sample.integrated$integrated_snn_res.3 %in% c(9, 69, 82)] <- "Tumor_MSLN"
sample.integrated$celltype[sample.integrated$integrated_snn_res.3 %in% c(11, 28, 36, 41, 60, 88)] <- "Tumor"
sample.integrated$celltype[sample.integrated$integrated_snn_res.3 %in% c(16)] <- "Tumor"
sample.integrated$celltype[sample.integrated$integrated_snn_res.3 %in% c(23)] <- "Plasma"
sample.integrated$celltype[sample.integrated$integrated_snn_res.3 %in% c(26, 81)] <- "CD4T_ICOS"
sample.integrated$celltype[sample.integrated$integrated_snn_res.3 %in% c(27, 50, 67)] <- "CD8T"
sample.integrated$celltype[sample.integrated$integrated_snn_res.3 %in% c(33, 35, 74, 76)] <- "CD8T_CD103"
sample.integrated$celltype[sample.integrated$integrated_snn_res.3 %in% c(41)] <- "Tumor"
sample.integrated$celltype[sample.integrated$integrated_snn_res.3 %in% c(44)] <- "NK"
sample.integrated$celltype[sample.integrated$integrated_snn_res.3 %in% c(32, 62, 83)] <- "Hypoxia_cells"
sample.integrated$celltype[sample.integrated$integrated_snn_res.3 %in% c(70)] <- "Endo_Prolif"
sample.integrated$celltype[sample.integrated$integrated_snn_res.3 %in% c(71)] <- "Myeloid"
sample.integrated$celltype[sample.integrated$integrated_snn_res.3 %in% c(46, 79, 86)] <- "B"

table(sample.integrated$integrated_snn_res.3, sample.integrated$celltype)
table(sample.integrated$celltype)
table(sample.integrated$integrated_snn_res.3)
table(sample.integrated$integrated_snn_res.0.7) %>% sort()

DimPlot(sample.integrated, label = TRUE, label.box = TRUE, raster = TRUE, group.by = "celltype", label.size = 2) + NoLegend()
ggsave("../../figures_celltype.jpg", width = 8, height = 7)

saveRDS(sample.integrated, file = "./sample_integrated.rds")

sample.integrated$celltype <- factor(sample.integrated$celltype, levels = sort(unique(sample.integrated$celltype)))
sort(unique(sample.integrated$celltype))
for (i in c(4, 5, 7, 10, 11, 12)) {
  print(i)
  # i <- 4
  ImageDimPlot(
    sample.integrated,
    fov = paste0("EC", i),
    # cols = "parade",
    cols = c(
      "B" = "#abfff1",
      "CD4T" = "#FFFFB3",
      "CD4T_ICOS" = "#ddd53e",
      "CD8T" = "#FB8072",
      "CD8T_CD103" = "#fba89f",
      "Endo" = "#ec9ed2",
      "Endo_Prolif" = "#d561dd",
      "Fibro" = "#d0b33d",
      "Hypoxia_cells" = "grey",
      "Myeloid" = "#ff0004",
      "NK" = "#8e3af4",
      "Plasma" = "#77f538",
      "Tumor" = "#373bbf",
      "Tumor_MSLN" = "#5893ba",
      "Tumor_Prolif" = "#b1cdf1"
    ),
    group.by = "celltype", border.size = NA, crop = TRUE, coord.fixed = FALSE
  ) + NoLegend()
  ggsave(paste0("../../figures_neighbor_ec", i, ".pdf"), width = 4.5, height = 6)
}

DimPlot(
  sample.integrated,
  label = TRUE, label.box = TRUE, raster = TRUE, group.by = "celltype", label.size = 2,
  cols = c(
    "B" = "#abfff1",
    "CD4T" = "#FFFFB3",
    "CD4T_ICOS" = "#ddd53e",
    "CD8T" = "#FB8072",
    "CD8T_CD103" = "#fba89f",
    "Endo" = "#ec9ed2",
    "Endo_Prolif" = "#d561dd",
    "Fibro" = "#d0b33d",
    "Hypoxia_cells" = "grey",
    "Myeloid" = "#ff0004",
    "NK" = "#8e3af4",
    "Plasma" = "#77f538",
    "Tumor" = "#373bbf",
    "Tumor_MSLN" = "#5893ba",
    "Tumor_Prolif" = "#b1cdf1"
  )
) + NoLegend()
ggsave("../../figures_celltype.pdf", width = 8, height = 7)

