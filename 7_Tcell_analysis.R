# trajectory analysis for T cells

# List of packages required
depp <- c(
  "languageserver", "httpgd", "Seurat", "DoubletFinder", "dplyr", "DT",
  "ggplot2", "patchwork", "data.table", "tidyr", "tidyverse", "tidygraph",
  "ggraph", "SeuratWrappers", "SeuratDisk", "doParallel", "openxlsx", "viridis",
  "rstatix", "BiocParallel"
)
# Check if the packages were installed if not install
depp_new <- depp[!(depp %in% installed.packages())]
depp_new
if (length(depp_new)) {
  install.packages(depp_new)
}
# remotes::install_github('satijalab/seurat-wrappers')
# remotes::install_github("mojaveazure/seurat-disk")
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
# remotes::install_github("aertslab/SCENIC")
# remotes::install_github("aertslab/SCopeLoomR")

# List of bioconductor packages required by this analysis
BioDepp <- c(
  "scater", "slingshot", "destiny", "ggrepel", "GSEABase",
  "org.Hs.eg.db", "clusterProfiler", "AUCell", "ggpubr", "preprocessCore"
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




 
# ---- Add clonotype information
filter_contig.files <- list.files(
  "/scratch/qlab/project_ec",
  pattern = "filtered_contig_annotations.csv",
  all.files = TRUE, full.names = TRUE, recursive = T
)
head(filter_contig.files)
contig_name <- filter_contig.files[grep("outs", filter_contig.files)]
head(contig_name)

tcr_files <- contig_name[grep("tcr", contig_name)]
bcr_files <- contig_name[grep("bcr", contig_name)]
bcr_files <- bcr_files[-1]

# sample.rdb <- readRDS("../EC_lym.rds")
colnames(sample.rdb@meta.data)
table(sample.rdb$lym_subtype, sample.rdb$orig.ident)

sample.rdb <- readRDS("../sample_all_original_rdb.rds")
sample.rdb <- subset(sample.rdb, subset = CellType_new %in% c("Bcell", "Plasma", "CD4T", "CD8T"))
colnames(sample.rdb@meta.data)
table(sample.rdb$CellType_new, sample.rdb$CellSubType)
table(sample.rdb$CellSubType, sample.rdb$CellType_new)
table(sample.rdb$CellSubType, !is.na(sample.rdb$tcr_clonotype_id))
table(sample.rdb$CellSubType, !is.na(sample.rdb$bcr_clonotype_id))
saveRDS(sample.rdb, "../sample_all_tcr_bcr.rds")



# ---- add function
add_clonotype <- function(file_dir, seurat_obj, type) {
  file_dir <- gsub("/outs/filtered_contig_annotations.csv", "", file_dir)
  mix_name <- gsub("/scratch/qlab/project_ec/", "", file_dir)
  mix_name <- gsub("sc5r", "", mix_name)
  mix_name <- gsub(type, "", mix_name)

  tcr <- read.csv(paste0(file_dir, "/outs/filtered_contig_annotations.csv"), sep = ",")
  tcr$barcode <- paste(mix_name, tcr$barcode, sep = "_")
  tcr <- tcr[!duplicated(tcr$barcode), ]
  names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"

  # Clonotype-centric info
  clono <- read.csv(paste0(file_dir, "/outs/clonotypes.csv"))

  # Slap the AA sequences onto our original table by clonotype_id.
  tcr <- merge(tcr[, c("barcode", "clonotype_id")], clono)

  # Reorder so barcodes are first column and set them as rownames.
  tcr$clonotype_mix <- paste(mix_name, tcr$clonotype_id, sep = "_")
  colnames(tcr) <- paste(type, colnames(tcr), sep = "_")
  return(tcr)
}


# ---- add TCR metadata
tcr_metadata <- lapply(1:length(tcr_files), function(x) {
  add_clonotype(tcr_files[x], sample.rdb, type = "tcr")
})
tcr_metadata <- rbindlist(tcr_metadata)
rownames(tcr_metadata) <- tcr_metadata$tcr_barcode
tcr_metadata[, "tcr_barcode"] <- NULL
dim(tcr_metadata)
head(tcr_metadata)
sample.rdb <- AddMetaData(object = sample.rdb, metadata = tcr_metadata)
head(sample.rdb@meta.data)


# ---- add BCR metadata
bcr_metadata <- lapply(1:length(bcr_files), function(x) {
  add_clonotype(bcr_files[x], sample.rdb, type = "bcr")
})
bcr_metadata <- rbindlist(bcr_metadata)
rownames(bcr_metadata) <- bcr_metadata$bcr_barcode
bcr_metadata[, "bcr_barcode"] <- NULL
dim(bcr_metadata)
head(bcr_metadata)
sample.rdb <- AddMetaData(object = sample.rdb, metadata = bcr_metadata)
head(sample.rdb@meta.data)


# ---- remove cells expressed both T and B clonotype
table(!is.na(sample.rdb$tcr_clonotype_id), !is.na(sample.rdb$bcr_clonotype_id))
sample.rmdup <- subset(
  sample.rdb,
  cells = colnames(sample.rdb)[!(!is.na(sample.rdb$tcr_clonotype_id) & !is.na(sample.rdb$bcr_clonotype_id))]
)
sample.dup <- subset(
  sample.rdb,
  cells = colnames(sample.rdb)[(!is.na(sample.rdb$tcr_clonotype_mix) & !is.na(sample.rdb$bcr_clonotype_mix))]
)
sample.dup@meta.data %>% colnames()
tcr_dup <- sample.dup$tcr_clonotype_mix
bcr_dup <- sample.dup$bcr_clonotype_mix





# CD8T
## add clonotype
cd8T <- readRDS("../../data/sample_integrated_cd8_filtered.rds")
table(cd8T$Patient)
cd8T <- subset(cd8T, subset = Tissue %in% c("CA"), invert = F)


# ---- add tissue and celltype
cd8T$Tissue_merge <- cd8T$Tissue
cd8T$Tissue_merge[cd8T$Tissue %in% c("MBC", "MCP", "MO")] <- "MET"
table(cd8T$Tissue_merge)
cd8T$CellType %>% unique()
cd8T$CellType <- factor(
  cd8T$CellType,
  levels = c("CD8_N", "CD8_EM1", "CD8_EM2", "CD8_EX", "CD8_RM", "CD8_EMRA")
)


# ---- add TCGA molecular class
cd8T$molecular_class <- "unknown"
cd8T$molecular_class[cd8T$Patient %in% c("EC15", "EC16", "EC22")] <- "POLE"
cd8T$molecular_class[cd8T$Patient %in% c("EC1", "EC10", "EC21", "EC25", "EC26")] <- "p53abn"
cd8T$molecular_class[cd8T$Patient %in%
  c("EC2", "EC6", "EC9", "EC12", "EC14", "EC18", "EC24", "EC27")] <- "MMR-d"
cd8T$molecular_class[cd8T$Patient %in%
  c("EC4", "EC5", "EC7", "EC8", "EC11", "EC17", "EC19", "EC23")] <- "NSMP/p53wt"
table(cd8T$molecular_class)


# ---- add TCR clonotype, rename tcr_clonotype_id based on cdr3s_aa
cd8T <- AddMetaData(object = cd8T, metadata = tcr_metadata)
table(!is.na(cd8T$tcr_clonotype_id))
change_id <- data.frame(raw = na.omit(unique(cd8T$tcr_cdr3s_aa)))
change_id$new <- paste0("clonotype", 1:nrow(change_id))
head(change_id)
cd8T$tcr_clonotype_id_new <- change_id$new[match(cd8T$tcr_cdr3s_aa, change_id$raw)]
table(cd8T$CellType, cd8T$tcr_clonotype_id_new > 0)

clonotype_num <- as.data.frame(table(cd8T$tcr_clonotype_id_new))
head(clonotype_num)
clonotype_num %>% arrange(Freq)
table(clonotype_num$Freq > 2) # 410 cells
table(clonotype_num$Freq > 5) # 151 cells

clonotype_num_sample <- as.data.frame(table(cd8T$tcr_clonotype_id_new, cd8T$Sample))
head(clonotype_num_sample)

clonotype_expa_2 <- clonotype_num_sample %>%
  group_by(Var2) %>%
  filter(Freq > 2) %>%
  summarize(freq_2 = n_distinct(Var1)) %>%
  arrange(freq_2)
clonotype_expa_2$cancer_type <- cd8T$cancer_type[match(clonotype_expa_2$Var2, cd8T$Sample)]
clonotype_expa_2

clonotype_expa_5 <- clonotype_num_sample %>%
  group_by(Var2) %>%
  filter(Freq > 5) %>%
  summarize(freq_5 = n_distinct(Var1)) %>%
  arrange(freq_5)
clonotype_expa_5$cancer_type <- cd8T$cancer_type[match(clonotype_expa_5$Var2, cd8T$Sample)]
clonotype_expa_5


## lineage analysis
# ---- lineage analysis
Idents(cd8T) <- "CellType"
sds <- slingshot(
  Embeddings(cd8T, "umap"),
  clusterLabels = cd8T$CellType,
  start.clus = "CD8_N",
  end.clus = c("CD8_EX", "CD8_RM", "CD8_EMRA"),
  stretch = 0
)
lineages <- slingLineages(sds)
lineages
curves <- slingCurves(sds, as.df = TRUE)
curved <- rbind(curves[curves$Lineage == 1, ], curves[curves$Lineage == 2, ], curves[curves$Lineage == 3, ])
curved$Lineage <- paste0("Lineage", curved$Lineage)
table(curved$Lineage)

cd8T$Pseudotime1 <- slingPseudotime(sds)[, 1]
cd8T$Pseudotime1 <- (cd8T$Pseudotime1 - min(cd8T$Pseudotime1, na.rm = T)) / (max(cd8T$Pseudotime1, na.rm = T) - min(cd8T$Pseudotime1, na.rm = T)) * 100
cd8T$Pseudotime2 <- slingPseudotime(sds)[, 2]
cd8T$Pseudotime2 <- (cd8T$Pseudotime2 - min(cd8T$Pseudotime2, na.rm = T)) / (max(cd8T$Pseudotime2, na.rm = T) - min(cd8T$Pseudotime2, na.rm = T)) * 100
cd8T$Pseudotime3 <- slingPseudotime(sds)[, 3]
cd8T$Pseudotime3 <- (cd8T$Pseudotime3 - min(cd8T$Pseudotime3, na.rm = T)) / (max(cd8T$Pseudotime3, na.rm = T) - min(cd8T$Pseudotime3, na.rm = T)) * 100

cd8T$Pseudotime <- rowMeans(cd8T@meta.data[, c("Pseudotime1", "Pseudotime2", "Pseudotime3")], na.rm = T)
head(cd8T@meta.data[, c("Pseudotime1", "Pseudotime2", "Pseudotime3", "Pseudotime")])

data <- FetchData(cd8T,
  vars = c("UMAP_1", "UMAP_2", "Pseudotime", "Pseudotime1", "Pseudotime2", "Pseudotime3")
)
head(data)


ncols <- c("#B22C2CFF", "#85B22CFF", "#E5B17EFF", "#51A3CCFF", "#E57E7EFF", "#BFB2FFFF")
pdf("../../figures/umap_integrated_cd8T_lineage.pdf", width = 5.5, height = 5)
DimPlot(object = cd8T, reduction = "umap", label = F, group.by = "CellType", pt.size = 0.5) +
  geom_path(aes_string("UMAP_1", "UMAP_2", linetype = "Lineage"), curved, size = 1, arrow = arrow(length = unit(0.15, "inches"))) +
  scale_fill_manual(values = ncols) +
  scale_color_manual(values = ncols) +
  NoLegend()
FeaturePlot(cd8T, feature = "Pseudotime") + scale_color_distiller(palette = "YlOrRd", na.value = "grey90") + NoLegend()
dev.off()


pdf("../../figures/cd8T_celltype_density_cancer_subtype_sub.pdf", width = 6, height = 3)
# lineage 1
lineage_plot <- ggplot(cd8T@meta.data[!is.na(cd8T$Pseudotime1), ]) +
  geom_density(aes(x = Pseudotime1, color = cancer_type)) +
  # scale_color_viridis(discrete = TRUE) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(colour = "cancer_type") +
  xlab("Lineage1: CD8_Tex")
celltype_plot <- ggplot(cd8T@meta.data[!is.na(cd8T$Pseudotime1), ]) +
  geom_jitter(aes(x = Pseudotime1, y = 0.00025, color = CellType), height = 0.0025, size = 0.01) +
  scale_color_manual(values = ncols) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank()
  )
(celltype_plot / lineage_plot) + plot_layout(heights = c(0.3, 2))
# lineage 2
lineage_plot <- ggplot(cd8T@meta.data[!is.na(cd8T$Pseudotime2), ]) +
  geom_density(aes(x = Pseudotime2, color = cancer_type)) +
  # scale_color_viridis(discrete = TRUE) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(colour = "cancer_type") +
  xlab("Lineage2: CD8_Trm")
celltype_plot <- ggplot(cd8T@meta.data[!is.na(cd8T$Pseudotime2), ]) +
  geom_jitter(aes(x = Pseudotime2, y = 0.00025, color = CellType), height = 0.0025, size = 0.01) +
  scale_color_manual(values = ncols) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank()
  )
(celltype_plot / lineage_plot) + plot_layout(heights = c(0.3, 2))
# lineage 3
lineage_plot <- ggplot(cd8T@meta.data[!is.na(cd8T$Pseudotime3), ]) +
  geom_density(aes(x = Pseudotime3, color = cancer_type)) +
  # scale_color_viridis(discrete = TRUE) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(colour = "cancer_type") +
  xlab("Lineage3: CD8_Temra")
celltype_plot <- ggplot(cd8T@meta.data[!is.na(cd8T$Pseudotime3), ]) +
  geom_jitter(aes(x = Pseudotime3, y = 0.00025, color = CellType), height = 0.0025, size = 0.01) +
  scale_color_manual(values = ncols) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank()
  )
(celltype_plot / lineage_plot) + plot_layout(heights = c(0.3, 2))
dev.off()

gene <- c(
  "CCR7", "LEF1", "TCF7", "CTLA4", "TOX2", "RUNX3", "GNLY", "ZNF683", "ITGA1",
  "ITGAE", "CX3CR1", "FGFBP2", "FCGR3A", "GZMB", "PDCD1", "PRF1", "IFNG",
  "TOX", "HAVCR2", "CTLA4", "LAG3"
)
pdf("../../figures/cd8T_feature.pdf", width = 10, height = 40)
Idents(cd8T) <- "CellType"
DefaultAssay(cd8T) <- "RNA"
FeaturePlot(cd8T, features = gene, label = T, ncol = 2)
DefaultAssay(cd8T) <- "integrated"
dev.off()

matrix <- GetAssayData(cd8T, assay = "RNA")
expr_matrix <- as.data.frame(t(rbind(
  matrix[gene, , drop = FALSE],
  Pseudotime1 = cd8T$Pseudotime1,
  Pseudotime2 = cd8T$Pseudotime2,
  Pseudotime3 = cd8T$Pseudotime3
)))
table(rownames(expr_matrix) %in% colnames(cd8T))
expr_matrix$Tissue <- cd8T$Tissue
expr_matrix$cancer_type <- cd8T$cancer_type
head(expr_matrix)

df <- pivot_longer(expr_matrix, gene, names_to = "feature", values_to = "expr")
head(df)





# ---- plot
pdf("../../figures/cd8T_exp_trend.pdf", width = 12, height = 4)
ggplot(df) +
  geom_smooth(aes(x = Pseudotime1, y = expr, color = "#DC0000B2"), method = "gam", se = F) +
  geom_smooth(aes(x = Pseudotime2, y = expr, color = "#0073C299"), method = "gam", se = F) +
  geom_smooth(aes(x = Pseudotime3, y = expr, color = "#EFC00099"), method = "gam", se = F) +
  xlab("PseudoTime") +
  ylab("Expression") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_identity(
    name = "Lineages",
    # labels = c("Lineage1: CD8_Tex",d "Lineage2: CD8_Trm", "Lineage3: CD8_Temra"),
    guide = "legend"
  ) +
  facet_wrap(~feature, nrow = 2, scales = "free")
dev.off()


## clonotype
clonotype_num <- as.data.frame(table(cd8T$tcr_clonotype_id_new))
cd8T$Clonotype_num <- clonotype_num$Freq[match(cd8T$tcr_clonotype_id_new, clonotype_num$Var1)]
table(cd8T$Clonotype_num)
cd8T$Clonotype_num_ld <- "Not detected"
cd8T$Clonotype_num_ld[cd8T$Clonotype_num == 1] <- "n = 1"
cd8T$Clonotype_num_ld[cd8T$Clonotype_num > 1 & cd8T$Clonotype_num <= 10] <- "1 < n < 10"
cd8T$Clonotype_num_ld[cd8T$Clonotype_num > 10 & cd8T$Clonotype_num <= 20] <- "10 < n < 20"
cd8T$Clonotype_num_ld[cd8T$Clonotype_num > 20] <- "n > 20"
cd8T$Clonotype_num_ld <- factor(
  cd8T$Clonotype_num_ld,
  levels = c("Not detected", "n = 1", "1 < n < 10", "10 < n < 20", "n > 20")
)
table(cd8T$Clonotype_num_ld)

pdf("../../figures/umap_integrated_cd8T_clonotype.pdf", width = 6.5, height = 5)
print(
  DimPlot(
    object = cd8T, reduction = "umap", label = F, pt.size = 0.5, order = T,
    group.by = "Clonotype_num_ld",
    cols = c("lightgrey", "darkgrey", "#E0C5BEFF", "#DE5C00FF", "#573333FF")
  )
)
print(
  DimPlot(
    object = cd8T, reduction = "umap", label = F, pt.size = 0.5, order = T,
    group.by = "Clonotype_num_ld", split.by = "Tissue_merge",
    cols = c("lightgrey", "darkgrey", "#E0C5BEFF", "#DE5C00FF", "#573333FF")
  )
)
print(
  DimPlot(
    object = cd8T, reduction = "umap", label = F, pt.size = 0.5, order = T,
    group.by = "Clonotype_num_ld", split.by = "Subtype",
    cols = c("lightgrey", "darkgrey", "#E0C5BEFF", "#DE5C00FF", "#573333FF")
  )
)
print(
  DimPlot(
    object = cd8T, reduction = "umap", label = F, pt.size = 0.5, order = T,
    group.by = "Clonotype_num_ld", split.by = "molecular_class",
    cols = c("lightgrey", "darkgrey", "#E0C5BEFF", "#DE5C00FF", "#573333FF")
  )
)
print(
  DimPlot(
    object = cd8T, reduction = "umap", label = F, pt.size = 0.5, order = T,
    group.by = "Clonotype_num_ld", split.by = "cancer_type",
    cols = c("lightgrey", "darkgrey", "#E0C5BEFF", "#DE5C00FF", "#573333FF")
  )
)
dev.off()


## clonotype richness
# ---- L1
summary(cd8T$Pseudotime1)
cd8T$l1_inx <- NA
cd8T$l1_inx[cd8T$Pseudotime1 >= 0 & cd8T$Pseudotime1 < 10] <- 0
cd8T$l1_inx[cd8T$Pseudotime1 >= 10 & cd8T$Pseudotime1 < 20] <- 10
cd8T$l1_inx[cd8T$Pseudotime1 >= 20 & cd8T$Pseudotime1 < 30] <- 20
cd8T$l1_inx[cd8T$Pseudotime1 >= 30 & cd8T$Pseudotime1 < 40] <- 30
cd8T$l1_inx[cd8T$Pseudotime1 >= 40 & cd8T$Pseudotime1 < 50] <- 40
cd8T$l1_inx[cd8T$Pseudotime1 >= 50 & cd8T$Pseudotime1 < 60] <- 50
cd8T$l1_inx[cd8T$Pseudotime1 >= 60 & cd8T$Pseudotime1 < 70] <- 60
cd8T$l1_inx[cd8T$Pseudotime1 >= 70 & cd8T$Pseudotime1 < 80] <- 70
cd8T$l1_inx[cd8T$Pseudotime1 >= 80 & cd8T$Pseudotime1 < 90] <- 80
cd8T$l1_inx[cd8T$Pseudotime1 >= 90 & cd8T$Pseudotime1 < 100] <- 90
table(cd8T$l1_inx)

# unique
richness_unique <- as.data.frame(table(cd8T$tcr_clonotype_id_new, cd8T$l1_inx, cd8T$cancer_type)) %>%
  group_by(Var2, Var3) %>%
  filter(Freq > 0) %>%
  summarize(unique_TCR = n())
richness_unique

# total
richness_total <- as.data.frame(table(cd8T$tcr_clonotype_id_new, cd8T$l1_inx, cd8T$cancer_type)) %>%
  group_by(Var2, Var3) %>%
  mutate(total_TCR = sum(Freq))
richness_total <- unique(richness_total[, c(2, 3, 5)])
richness_total

# richness
richness <- right_join(richness_unique, richness_total)
# richness <- left_join(richness_unique, richness_total)
richness$unique_TCR[is.na(richness$unique_TCR)] <- 0
richness$richness <- richness$unique_TCR / richness$total_TCR
richness <- richness[, c(1, 2, 5)]
names(richness) <- c("time1", "cancer_type", "richness")
unique(richness$cancer_type)

# plot clonotype richness
pdf("../../figures/cd8T_clonotype_density_cancer_subtype_l1_sub.pdf", width = 6, height = 3.2)
ggplot(richness, aes(x = time1, y = richness, group = cancer_type, color = cancer_type)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ylim(0, 1)
dev.off()



# ---- L2
summary(cd8T$Pseudotime2)
cd8T$l2_inx <- NA
cd8T$l2_inx[cd8T$Pseudotime2 >= 0 & cd8T$Pseudotime2 < 10] <- 0
cd8T$l2_inx[cd8T$Pseudotime2 >= 10 & cd8T$Pseudotime2 < 20] <- 10
cd8T$l2_inx[cd8T$Pseudotime2 >= 20 & cd8T$Pseudotime2 < 30] <- 20
cd8T$l2_inx[cd8T$Pseudotime2 >= 30 & cd8T$Pseudotime2 < 40] <- 30
cd8T$l2_inx[cd8T$Pseudotime2 >= 40 & cd8T$Pseudotime2 < 50] <- 40
cd8T$l2_inx[cd8T$Pseudotime2 >= 50 & cd8T$Pseudotime2 < 60] <- 50
cd8T$l2_inx[cd8T$Pseudotime2 >= 60 & cd8T$Pseudotime2 < 70] <- 60
cd8T$l2_inx[cd8T$Pseudotime2 >= 70 & cd8T$Pseudotime2 < 80] <- 70
cd8T$l2_inx[cd8T$Pseudotime2 >= 80 & cd8T$Pseudotime2 < 90] <- 80
cd8T$l2_inx[cd8T$Pseudotime2 >= 90 & cd8T$Pseudotime2 < 100] <- 90
table(cd8T$l2_inx)
# unique
richness_unique <- as.data.frame(table(cd8T$tcr_clonotype_id_new, cd8T$l2_inx, cd8T$cancer_type)) %>%
  group_by(Var2, Var3) %>%
  filter(Freq > 0) %>%
  summarize(unique_TCR = n())
richness_unique

# total
richness_total <- as.data.frame(table(cd8T$tcr_clonotype_id_new, cd8T$l2_inx, cd8T$cancer_type)) %>%
  group_by(Var2, Var3) %>%
  mutate(total_TCR = sum(Freq))
richness_total <- unique(richness_total[, c(2, 3, 5)])
richness_total

# richness
richness <- right_join(richness_unique, richness_total)
# richness <- left_join(richness_unique, richness_total)
richness$unique_TCR[is.na(richness$unique_TCR)] <- 0
richness$richness <- richness$unique_TCR / richness$total_TCR
richness <- richness[, c(1, 2, 5)]
names(richness) <- c("time2", "cancer_type", "richness")
unique(richness$cancer_type)

# plot clonotype richness
pdf("../../figures/cd8T_clonotype_density_cancer_subtype_l2_sub.pdf", width = 6, height = 3.2)
ggplot(richness, aes(x = time2, y = richness, group = cancer_type, color = cancer_type)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ylim(0, 1)
dev.off()



# ---- L3
summary(cd8T$Pseudotime3)
cd8T$l3_inx <- NA
cd8T$l3_inx[cd8T$Pseudotime3 >= 0 & cd8T$Pseudotime3 < 10] <- 0
cd8T$l3_inx[cd8T$Pseudotime3 >= 10 & cd8T$Pseudotime3 < 20] <- 10
cd8T$l3_inx[cd8T$Pseudotime3 >= 20 & cd8T$Pseudotime3 < 30] <- 20
cd8T$l3_inx[cd8T$Pseudotime3 >= 30 & cd8T$Pseudotime3 < 40] <- 30
cd8T$l3_inx[cd8T$Pseudotime3 >= 40 & cd8T$Pseudotime3 < 50] <- 40
cd8T$l3_inx[cd8T$Pseudotime3 >= 50 & cd8T$Pseudotime3 < 60] <- 50
cd8T$l3_inx[cd8T$Pseudotime3 >= 60 & cd8T$Pseudotime3 < 70] <- 60
cd8T$l3_inx[cd8T$Pseudotime3 >= 70 & cd8T$Pseudotime3 < 80] <- 70
cd8T$l3_inx[cd8T$Pseudotime3 >= 80 & cd8T$Pseudotime3 < 90] <- 80
cd8T$l3_inx[cd8T$Pseudotime3 >= 90 & cd8T$Pseudotime3 < 100] <- 90
table(cd8T$l3_inx)
# unique
richness_unique <- as.data.frame(table(cd8T$tcr_clonotype_id_new, cd8T$l3_inx, cd8T$cancer_type)) %>%
  group_by(Var2, Var3) %>%
  filter(Freq > 0) %>%
  summarize(unique_TCR = n())
print(richness_unique, n = 100)

# total
richness_total <- as.data.frame(table(cd8T$tcr_clonotype_id_new, cd8T$l3_inx, cd8T$cancer_type)) %>%
  group_by(Var2, Var3) %>%
  mutate(total_TCR = sum(Freq))
richness_total <- unique(richness_total[, c(2, 3, 5)])
print(richness_total, n = 100)

# richness
richness <- right_join(richness_unique, richness_total)
# richness <- left_join(richness_unique, richness_total)
richness$unique_TCR[is.na(richness$unique_TCR)] <- 0
richness$richness <- richness$unique_TCR / richness$total_TCR
richness <- richness[, c(1, 2, 5)]
names(richness) <- c("time3", "cancer_type", "richness")
unique(richness$cancer_type)

# plot clonotype richness
pdf("../../figures/cd8T_clonotype_density_cancer_subtype_l3_sub.pdf", width = 6, height = 3.2)
ggplot(richness, aes(x = time3, y = richness, group = cancer_type, color = cancer_type)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ylim(0, 1)
dev.off()


## richness
table(cd8T$Sample)
# unique
richness_unique <- as.data.frame(table(cd8T$tcr_clonotype_id_new, cd8T$CellType, cd8T$Sample)) %>%
  group_by(Var2, Var3) %>%
  filter(Freq > 0) %>%
  summarize(unique_TCR = n())


# total
richness_total <- as.data.frame(table(cd8T$tcr_clonotype_id_new, cd8T$CellType, cd8T$Sample)) %>%
  group_by(Var2, Var3) %>%
  mutate(total_TCR = sum(Freq))
richness_total <- unique(richness_total[, c(2, 3, 5)])
richness_total$Var4 <- cd8T$Tissue_merge[match(richness_total$Var3, cd8T$Sample)]
head(richness_total)

# richness
richness <- right_join(richness_unique, richness_total)
# richness <- left_join(richness_unique, richness_total)
richness$unique_TCR[is.na(richness$unique_TCR)] <- 0
richness$Var2 <- factor(
  richness$Var2,
  levels = c("CD8_N", "CD8_EM1", "CD8_EM2", "CD8_EX", "CD8_RM", "CD8_EMRA")
)
richness$richness <- richness$unique_TCR / richness$total_TCR
richness
names(richness) <- c("celltype", "sample", "unique_TCR", "total_TCR", "tissue", "richness")

richness$subtype <- cd8T$Subtype[match(richness$sample, cd8T$Sample)]
richness$molecular_class <- cd8T$molecular_class[match(richness$sample, cd8T$Sample)]
richness$cancer_type <- cd8T$cancer_type[match(richness$sample, cd8T$Sample)]
unique(richness$cancer_type)
richness$cancer_type <- factor(
  richness$cancer_type,
  levels = c("Ciliated", "Glandular", "Luminal_LCN2", "Luminal_ERBB4", "EMT-like", "Glandular_CXCL14")
)
