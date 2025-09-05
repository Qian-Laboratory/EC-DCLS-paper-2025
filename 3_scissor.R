# Scissor analysis
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





library(org.Hs.eg.db)
library(biomaRt)
library(Scissor)

sample.rdb <- readRDS("../sample_integrated_rdb.rds")
Cancer_ca <- readRDS("../sample_integrated_Cancer_ca_rdb_rdb_rmCC_regCount.rds")
Cancer_pca <- readRDS("../sample_integrated_Cancer_pca_rdb.rds")
colnames(sample.rdb@meta.data)
unique(sample.rdb$CellSubType)

sample_cancer <- subset(sample.rdb, subset = CellSubType %in% c(
  "Epi_01_SCGB2A1",
  "Epi_02_CXCL1",
  "Epi_03_WNT7A",
  "Epi_04_LCN2",
  "Epi_05_FOXJ1",
  "Epi_06_MMP7",
  "Cancer_01_ciliated",
  "Cancer_02_glandular",
  "Cancer_03_oncogenic",
  "Cancer_04_ERBB4",
  "Cancer_05_basal_cell",
  "Cancer_06_EMT-like",
  "Cancer_07_fibroblast-like",
  "Cancer_08_CXCL14"
))
DefaultAssay(sample_cancer) <- "RNA"


# ---- scRNA-seq raw counts data:
sc_dataset <- Seurat_preprocessing(
  GetAssayData(Cancer_ca, assay = "RNA", slot = "count"),
  verbose = T
)


# ---- Prepare the bulk data and phenotype
TCGA_all <- read.delim("/share/reference/TCGA/GDC-PANCAN/GDC-PANCAN.htseq_fpkm-uq.tsv.gz", row.names = 1)
phenotype <- read.delim("/share/reference/TCGA/GDC-PANCAN/GDC-PANCAN.basic_phenotype.tsv.gz")
survival <- read.delim("/share/reference/TCGA/GDC-PANCAN/GDC-PANCAN.survival.tsv")

colnames(TCGA_all) <- gsub(r"(\.)", "-", colnames(TCGA_all))
tcga.id <- colnames(TCGA_all)
phenotype.id <- phenotype$sample[grep("COAD", phenotype$project_id)]
survival.id <- survival$sample
ucec.id <- intersect(intersect(phenotype.id, survival.id), tcga.id)


# ---- bulk_dataset
bulk_dataset <- TCGA_all[, ucec.id]
dim(bulk_dataset)
# convert ensemble to gene symbol
library(tinyarray)
bulk_dataset <- trans_exp(TCGA_all[, ucec.id])
bulk_dataset_norm <- normalize.quantiles(x = as.matrix(bulk_dataset))
rownames(bulk_dataset_norm) <- rownames(bulk_dataset)
colnames(bulk_dataset_norm) <- colnames(bulk_dataset)


# ---- bulk_survival
survival_ec <- survival[match(ucec.id, survival$sample), ]
bulk_survival <- data.frame(
  TCGA_patient_barcode = survival_ec$sample,
  OS_time = survival_ec$OS.time,
  Status = survival_ec$OS
)
dim(bulk_survival)

# test consistency
all(colnames(bulk_dataset_norm) == bulk_survival$TCGA_patient_barcode) # 检验一下病人的id顺序是否与bulk表达量矩阵中的顺序相同


# ---- phenotype
phenotype <- bulk_survival[, 2:3]
colnames(phenotype) <- c("time", "status")
head(phenotype)



# ---- Execute Scissor to select the informative cells
# use Scissor to select the phenotype-associated cell subpopulations
infos1 <- Scissor(bulk_dataset_norm, sc_dataset, phenotype,
  alpha = 0.05,
  family = "cox",
  # Save_file = "../Scissor_cancer_survival_test.RData"
)
str(infos1)


# ---- visualize the Scissor selected cells by adding a new annotation in the Seurat object
load("./Scissor_cancer_survival.RData")
Scissor_select <- rep("Background cells", ncol(Cancer_ca))
names(Scissor_select) <- colnames(Cancer_ca)
Scissor_select[infos1$Scissor_pos] <- "Scissor+ cell"
Scissor_select[infos1$Scissor_neg] <- "Scissor- cell"
Cancer_ca <- AddMetaData(Cancer_ca, metadata = Scissor_select, col.name = "scissor")
subtype_ct <- as.data.frame(table(Cancer_ca$Cancer_ca_subtype, Cancer_ca$scissor)) %>%
  group_by(Var1) %>%
  mutate(percent = Freq / sum(Freq)) %>%
  print(n = 300)
subtype_ct <- subtype_ct[c(9:24), ]

pdf("../../figures/cancer_scissor_ucec_os.pdf", width = 10, height = 5)
Idents(Cancer_ca) <- "Cancer_ca_subtype"
print(
  DimPlot(Cancer_ca,
    reduction = "umap", group.by = "scissor", label = F,
    cols = c("grey", "royalblue", "indianred1"),
    pt.size = 1.2, order = c("Scissor+ cell", "Scissor- cell", "Background cells")
  )
)
print(
  ggplot(data = subtype_ct, aes(x = Var1, y = percent, fill = Var2)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(label = Freq),
      vjust = 1.6, color = "white",
      position = position_dodge(0.9), size = 3.5
    ) +
    scale_fill_brewer(palette = "Paired") +
    theme_minimal()
)
dev.off()
