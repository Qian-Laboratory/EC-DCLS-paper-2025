# Cell-Cell Interaction
# List of packages required
depp <- c(
  "languageserver", "httpgd", "Seurat", "DoubletFinder", "dplyr", "DT", "lme4",
  "ggplot2", "patchwork", "data.table", "tidyr", "tidyverse", "tidygraph", "pheatmap",
  "ggraph", "SeuratWrappers", "doParallel", "openxlsx", "viridis", "rstatix"
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
  "scater", "slingshot", "destiny", "ggrepel", "RColorBrewer", "ktplots", "BiocParallel",
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




# prepare input
load("../scEC_data.rda")
sample.all.rdb <- subset(sample.all, subset = CellSubType %in% "Cancer", invert = T)
sample.all.rdb$CellSubType <- factor(
  sample.all.rdb$CellSubType,
  levels = c(
    "CD4T_01_N", "CD4T_02_EM", "CD4T_03_FH", "CD4T_04_TH1", "CD4T_05_TH17", "CD4T_06_Treg",
    "CD8T_01_N", "CD8T_02_EM1", "CD8T_03_EM2", "CD8T_04_EX", "CD8T_05_RM", "CD8T_06_EMRA", "CD8T_07_MAIT",
    "NK_01_NK1", "NK_02_NK2", "NK_03_NK3",
    "dgT", "dgT_V9", "ILC", "ILC_MAIT",
    "B_01_mature_naive", "B_02_memory_IgM+", "B_03_memory_IgM-", "B_04_GC", "B_05_Plasma",
    "DC_01_cDC1", "DC_02_cDC2", "DC_03_migDC", "DC_04_pDC", "DC_05_langDC", "DC_06_ITGAX",
    "Mono_01_CD14", "Mono_02_CD16",
    "Marco_03_CXCL3", "Marco_04_MT1G", "Marco_05_MMP9", "Marco_06_CCL18", "Marco_07_CX3CR1", "Marco_08_FOLR2", "Marco_09_ATG7", "Mast",
    "EC_01_lymEC", "EC_02_TipCell", "EC_03_Artery", "EC_04_Capillary", "EC_05_Activated_PCV", "EC_06_Vein",
    "Fb_01_APOD", "Fb_02_IGFBP3", "CAF_03_MMP11", "CAF_04_ARID1B", "SMC_05_ACTG2", "Pericyte_06_RGS5", "Myfb_07_MYH11",
    "Epi_01_SCGB2A1", "Epi_02_CXCL1", "Epi_03_WNT7A", "Epi_04_LCN2", "Epi_05_FOXJ1", "Epi_06_MMP7",
    "Cancer_01_ciliated_FOXJ1", "Cancer_02_glandular_SCGB2A1", "Cancer_03_luminal_LCN2",
    "Cancer_04_luminal_ERBB4", "Cancer_05_basal_KRT5", "Cancer_06_EMT_like_SULF1",
    "Cancer_07_fibro_like_COL3A1", "Cancer_08_glandular_CXCL14"
  )
)
sample.all.rdb$CellType
sample.all.rdb$CellSubType <- as.character(sample.all.rdb$CellSubType)


# ---- cancer subtype
ciliated_ca <- c("EC4CA", "EC5CA", "EC9CA", "EC17CA", "EC23CA")
glandular_ca <- c("EC1CA", "EC2CA", "EC7CA", "EC12CA", "EC15CA", "EC18CA", "EC24CA", "EC26CA", "EC27CA")
luminal_ca <- c("EC6CA", "EC8CA", "EC10CA", "EC14CA", "EC16CA", "EC21CA", "EC22CA", "EC25CA")
EMT_ca <- c("EC11CA")
CXCL14_ca <- c("EC19CA")

sample.all.rdb$cancer_type <- "unknown"
sample.all.rdb$cancer_type[sample.all.rdb$Patient %in% c("EC19")] <- "Cancer_08_glandular_CXCL14"
sample.all.rdb$cancer_type[sample.all.rdb$Patient %in% c("EC11")] <- "Cancer_06_EMT_like_SULF1"
sample.all.rdb$cancer_type[sample.all.rdb$Patient %in% c("EC4", "EC5", "EC23", "EC17", "EC9")] <- "Ciliated_type"
sample.all.rdb$cancer_type[sample.all.rdb$Patient %in% c("EC7", "EC1", "EC2", "EC12", "EC15", "EC27", "EC24", "EC26", "EC18")] <- "Glandular_type"
sample.all.rdb$cancer_type[sample.all.rdb$Patient %in% c("EC10", "EC6", "EC14", "EC21", "EC8", "EC22", "EC25", "EC16")] <- "Luminal_type"





# compare cpdb - luminal vs others
## CA
# ---- cpdb metadata
cpdb_meta <- data.frame(
  sample_id = sort(unique(sample.all.rdb$Sample[sample.all.rdb$Tissue_merge == "CA"])),
  cellphonedb_folder = paste0("../../22_CCI/h5ad/", sort(unique(sample.all.rdb$Sample[sample.all.rdb$Tissue_merge == "CA"]))),
  sce_file = paste0("../../22_CCI/h5ad/", sort(unique(sample.all.rdb$Sample[sample.all.rdb$Tissue_merge == "CA"])), "/sample_all_RNA.h5ad")
)


# ---- sample metadata
sample_metadata <- as.data.frame(sample.all.rdb@meta.data[sample.all.rdb$Tissue_merge == "CA", c("Sample", "Patient", "cancer_type")] %>% unique())
sample_metadata <- sample_metadata[match(cpdb_meta$sample_id, sample_metadata$Sample), ]
rownames(sample_metadata) <- NULL
sample_metadata$cancer_type_luminal <- "Other_types"
sample_metadata$cancer_type_luminal[sample_metadata$Patient %in% c("EC10", "EC6", "EC14", "EC21", "EC8", "EC22", "EC25", "EC16")] <- "Luminal_type"

head(cpdb_meta)
head(sample_metadata)





# ---- compare cpdb
sample.all.rdb$CellType %>% unique()
sample_metadata$cancer_type_luminal %>% unique()
out <- compare_cpdb(
  cpdb_meta = cpdb_meta,
  sample_metadata = sample_metadata,
  celltypes = unique(sample.all.rdb$CellType[sample.all.rdb$Tissue_merge == "CA"]),
  celltype_col = "CellType",
  method = "wilcox.test",
  groupby = "cancer_type_luminal"
)
save(out, file = "./CellType_cpdb_compare_luminal_ca.rda")





# ---- convert list into dataframe
load("./CellType_cpdb_compare_luminal_ca.rda")
out_table <- out$Other_types_vs_Luminal_type %>% rownames_to_column(var = "LPR")
# split from and to celltype
rows <- strsplit(out_table$LPR, ">@<")
from <- lapply(rows, function(x) x[1])
from <- do.call(rbind, from)
to <- lapply(rows, function(x) x[2])
to <- do.call(rbind, to)
interaction <- lapply(rows, function(x) x[3])
interaction <- do.call(rbind, interaction)
out_table <- as.data.frame(cbind(out_table[, c("LFC", "pval", "padj", "contrast", "celltypes")], interaction = interaction, from = from, to = to))

# rename interaction
out_table$interaction <- gsub(" ", "-", out_table$interaction)
out_table$interaction[grep("^X[1-9]", out_table$interaction)] <- gsub("^X", "", out_table$interaction[grep("^X[1-9]", out_table$interaction)])
sort(unique(out_table$interaction))

# add interaction gene name
interaction_input <- read_csv("~/reference/CellPhoneDB/data/interaction_input_gene.csv")
interaction_input$id_cp_interaction <- gsub("_", "-", interaction_input$id_cp_interaction)
interaction_input$id_cp_interaction <- gsub(":", "-", interaction_input$id_cp_interaction)
out_table$interaction_gene <- interaction_input$interaction_gene[match(out_table$interaction, interaction_input$id_cp_interaction)]
table(is.na(out_table$interaction_gene))

# save result
write_csv(out_table, "./CellType_luminal_ca.csv")
