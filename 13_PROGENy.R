# PROGENy analysis
options(timeout = 10000)

# List of packages required
depp <- c(
  "languageserver", "httpgd", "Seurat", "DoubletFinder", "dplyr", "DT", "ggplot2",
  "patchwork", "data.table", "tidyr", "tidyverse", "tidygraph", "ggraph",
  "SeuratWrappers", "SeuratDisk", "doParallel", "openxlsx", "viridis", "rstatix",
  "KernSmooth", "BiocParallel"
)
# Check if the packages were installed if not install
depp_new <- depp[!(depp %in% installed.packages())]
if (length(depp_new)) {
  install.packages(depp_new)
}
# remotes::install_github('satijalab/seurat-wrappers', force = TRUE)
# remotes::install_github("mojaveazure/seurat-disk", force = TRUE)
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)
# remotes::install_github("aertslab/SCopeLoomR", force = TRUE)
# remotes::install_github("aertslab/SCENIC", force = TRUE)
# devtools::install_github('sunduanchen/Scissor', force = TRUE)

# List of bioconductor packages required by this analysis
BioDepp <- c(
  "scater", "destiny", "ggrepel", "GSEABase", "ggpubr", "org.Hs.eg.db",
  "clusterProfiler", "biomaRt", "preprocessCore", "progeny"
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





# Progeny - Scissor
library(decoupleR)
library(progeny)

sample.cancer <- readRDS("../sample_integrated_Cancer_ca_rdb_rdb_rmCC_regCount.rds")
head(sample.cancer@meta.data)
table(sample.cancer$Cancer_ca_subtype)


# ---- pathway activities
net <- readRDS("~/reference/Progeny/progeny_human_100.rds")
unique(net$source)
head(net)

as_matrix <- function(mat) {
  tmp <- matrix(data = 0L, nrow = mat@Dim[1], ncol = mat@Dim[2])

  row_pos <- mat@i + 1
  col_pos <- findInterval(seq(mat@x) - 1, mat@p[-1]) + 1
  val <- mat@x

  for (i in seq_along(val)) {
    tmp[row_pos[i], col_pos[i]] <- val[i]
  }

  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}
mat <- as_matrix(sample.cancer@assays$RNA@data)

pathway_activity <- progeny(mat, scale = TRUE, organism = "Human", top = 100)
head(pathway_activity)

acts <- run_mlm(
  mat = mat, net = net,
  .source = "source",
  .target = "target",
  .mor = "weight",
  minsize = 5
)
acts
save(acts, file = "../acts.rda")


# ---- Extract mlm and store it in pathwaysmlm in data
sample.cancer[["pathwaysmlm"]] <- acts %>%
  pivot_wider(
    id_cols = "source", names_from = "condition",
    values_from = "score"
  ) %>%
  column_to_rownames("source") %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = sample.cancer) <- "pathwaysmlm"

# Scale the data
sample.cancer <- ScaleData(sample.cancer)
saveRDS(sample.cancer, "../sample_integrated_Cancer_ca_progeny.rds")


# ---- check progeny pathway activities
pdf("./test_500.pdf", width = 8, height = 8)
DimPlot(sample.cancer, reduction = "umap", group.by = "Cancer_ca_subtype", label = TRUE, pt.size = 0.5) +
  NoLegend() + ggtitle("Cell types")
for (path in unique(net$source)) {
  print(
    FeaturePlot(sample.cancer, features = path) &
      scale_colour_gradient2(low = "blue", mid = "white", high = "red")
  ) +
    ggtitle(path)
}
dev.off()



# ---- correlation between progeny pathways and marker genes of scissor_pos cells
mygene <- "LCN2"
expr.mygene <- GetAssayData(sample.cancer, assay = "RNA", slot = "data")[mygene, ]

pathways <- unique(net$source)
df <- foreach(mypath = pathways, .combine = "rbind") %do% {
  # mypath  <- pathways[1]
  score <- GetAssayData(sample.cancer, assay = "pathwaysmlm", slot = "data")[mypath, ]
  mytest <- cor.test(score, expr.mygene)
  mytest$pathway <- mypath
  data.frame(pathway = mypath, gene = mygene, cor = mytest$estimate, pval = mytest$p.value, t = mytest$statistic)
}
df %>% arrange(cor)




# ---- signature gene set: cancer_subtype + scissor_pos genes
markers.LCN2 <- read.xlsx("../../umap_integrated_Cancer_ca_AllMarker_celltype.xlsx") %>%
  dplyr::filter(cluster == "c3_oncogenic", avg_log2FC > 0, p_val_adj < 0.01, pct.1 / pct.2 > 1.25 | pct.1 - pct.2 > 0.1)
# dplyr::pull(gene)
head(markers.LCN2)

load("../../figures/cancer_scissor_pos_c3_oncogenic.RData")
head(marker_ciliated)
scissor.df <- marker_ciliated %>% rownames_to_column(var = "gene")
scissor.df <- scissor.df %>%
  dplyr::filter(gene %in% markers.LCN2$gene, p_val_adj < 0.05) %>%
  arrange(-avg_log2FC)
head(scissor.df)
# scissor.df$gene[scissor.df$gene == "H19"] <- "LCN2"
scissor.df[grep("LCN2", scissor.df$gene), ]


unique(net$source)
mypath <- "TNFa"
scores <- GetAssayData(sample.cancer, assay = "pathwaysmlm", slot = "data")[mypath, ]
registerDoParallel(cores = 10)
NFKB.df <- foreach(gene = markers.LCN2$gene, .combine = "rbind") %dopar% {
  # mygene  <- "KIF26B"
  mygene <- gene
  expr.mygene <- GetAssayData(sample.cancer, assay = "RNA", slot = "data")[mygene, ]
  mytest <- cor.test(scores, expr.mygene, method = "pearson")
  mytest$pathway <- mypath
  data.frame(pathway = mypath, gene = mygene, r = mytest$estimate, pval = mytest$p.value, t = mytest$statistic)
}
NFKB.df %>%
  arrange(-r) %>%
  head(n = 20)
NFKB.df[NFKB.df$gene == "LCN2", ]



mypath <- "JAK-STAT"
scores <- GetAssayData(sample.cancer, assay = "pathwaysmlm", slot = "data")[mypath, ]
registerDoParallel(cores = 10)
EGFR.df <- foreach(gene = markers.LCN2$gene, .combine = "rbind") %dopar% {
  mygene <- gene
  expr.mygene <- GetAssayData(sample.cancer, assay = "RNA", slot = "data")[mygene, ]
  mytest <- cor.test(scores, expr.mygene, method = "pearson")
  mytest$pathway <- mypath
  data.frame(pathway = mypath, gene = mygene, r = mytest$estimate, pval = mytest$p.value, t = mytest$statistic)
}
EGFR.df %>%
  arrange(-r) %>%
  head()
EGFR.df[EGFR.df$gene == "LCN2", ]



NFKB.df <- NFKB.df %>%
  dplyr::rename(NFKB_r = r, NFKB_pval = pval, NFKB_t = t)
EGFR.df <- EGFR.df %>%
  dplyr::rename(EGFR_r = r, EGFR_pval = pval, EGFR_t = t)
combined.df <- NFKB.df %>%
  dplyr::left_join(EGFR.df, by = "gene") %>%
  # dplyr::left_join(scissor.df, by = c("gene"))
  dplyr::left_join(markers.LCN2, by = c("gene"))
head(combined.df)


combined.df %>%
  arrange(-EGFR_t) %>%
  head(20)

combined.df %>%
  filter(NFKB_t > 50 | EGFR_t > 25, avg_log2FC > 0.4) %>%
  pull(gene)




# Progeny - heatmap
pathways <- unique(net$source)
head(sample.cancer@meta.data)

Idents(sample.cancer) <- "Cancer_ca_subtype"
DefaultAssay(sample.cancer) <- "pathwaysmlm"

expr_ca.df <- AverageExpression(sample.cancer, features = pathways, use.scale = F)

markers <- c(
  "TNFa", "JAK-STAT", "EGFR", "VEGF", "NFkB", "Trail", "MAPK", "Hypoxia",
  "p53", "WNT", "PI3K", "Androgen", "TGFb", "Estrogen"
)
markers.df <- tibble(Gene = markers, GeneOrder = 1:length(markers))

expr_combined.long <- expr_ca.df$pathwaysmlm %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  left_join(markers.df, by = "Gene") %>%
  gather(key = celltype, value = expr, -Gene, -GeneOrder)
head(expr_combined.long)

expr_combined.long[expr_combined.long$celltype == "c3-oncogenic", ] %>%
  arrange(-expr) %>%
  pull(Gene)

expr_combined.long <- expr_combined.long %>%
  group_by(Gene) %>%
  mutate(zscore = (expr - mean(expr)) / (sd(expr) + 0.001))
head(expr_combined.long)

plot.order.heatmap <- str_sort(unique(expr_combined.long$celltype), numeric = T)
expr_combined.long$celltype <- factor(expr_combined.long$celltype, levels = plot.order.heatmap)

range(expr_combined.long$expr)
heatmap.markers <- ggplot(data = expr_combined.long) +
  geom_tile(aes(x = celltype, y = reorder(Gene, -GeneOrder), fill = zscore), color = "white", size = 0.1) + # bigger size means bigger spacer
  scale_fill_gradient2(low = "#0571B0", mid = "white", high = "#CA0020")
  # c("#CA0020", "#F4A582", "#F7F7F7", "#92C5DE", "#0571B0") +
  scale_y_discrete(position = "right") +
  xlab("characteristics") +
  theme_grey(base_size = 5) +
  ggtitle("cancer celltype") +
  theme(
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size = 10, colour = "gray50"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave("../heatmap_cancer_Progeny_celltype_markers_scale.pdf", plot = heatmap.markers, width = 2, height = 1.5)