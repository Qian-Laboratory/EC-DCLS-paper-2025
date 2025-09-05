# TCGA deconvolution with BayesPrism

# List of packages required
depp <- c(
  "languageserver", "httpgd", "Seurat", "DoubletFinder", "dplyr", "DT", "lme4",
  "ggplot2", "patchwork", "data.table", "tidyr", "tidyverse", "tidygraph",
  "ggraph", "SeuratWrappers", "doParallel", "openxlsx", "viridis", "rstatix"
)
# Check if the packages were installed if not install
depp_new <- depp[!(depp %in% installed.packages())]
depp_new
if (length(depp_new)) {
  install.packages(depp_new)
}
# remotes::install_github('satijalab/seurat-wrappers')
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
# remotes::install_github("Danko-Lab/BayesPrism/BayesPrism")

# List of bioconductor packages required by this analysis
BioDepp <- c(
  "scater", "slingshot", "destiny", "ggrepel", "ComplexHeatmap", "RColorBrewer",
  "GSEABase", "ggpubr", "org.Hs.eg.db", "clusterProfiler", "BayesPrism"
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




# 529 TCGA data
library(TCGAbiolinks)
library(SummarizedExperiment)

projects <- TCGAbiolinks::getGDCprojects()$project_id
projects <- projects[grepl("^TCGA", projects, perl = TRUE)]

query <- GDCquery(
  project = "TCGA-UCEC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(
  query,
  method = "api",
  directory = "../../../../reference/TCGA/UCEC/ucec_tcga_pan_can_atlas_2018(529)/raw_count/",
  files.per.chunk = 100
)

exprs.sce <- GDCprepare(
  query,
  directory = "../../../../reference/TCGA/UCEC/ucec_tcga_pan_can_atlas_2018(529)/raw_count/"
)
save(exprs.sce, file = "./TCGA-UCEC_mRNA.Rdata")

load("./TCGA-UCEC_mRNA.Rdata")
se_mrna <- exprs.sce[rowData(exprs.sce)$gene_type == "protein_coding", colData(exprs.sce)$sample_type == "Primary Tumor"]
exprs_mrna_count <- assay(se_mrna, "unstranded")
exprs_mrna_tpm <- assay(se_mrna, "tpm_unstrand")
dim(exprs_mrna_count)

# add sample ID
sample_id <- colData(se_mrna)$patient
head(sample_id)
exprs_mrna_count <- cbind(data.frame(sample_id), as.data.frame(t(exprs_mrna_count)))
exprs_mrna_tpm <- cbind(data.frame(sample_id), as.data.frame(t(exprs_mrna_tpm)))

# remove duplicate samples
exprs_mrna_count <- exprs_mrna_count %>%
  as_tibble() %>%
  mutate(meanrow = rowMeans(.[, -1]), .before = 2) %>%
  arrange(desc(meanrow)) %>%
  distinct(sample_id, .keep_all = T) %>%
  select(-meanrow) %>%
  column_to_rownames(var = "sample_id") %>%
  as.data.frame()
dim(exprs_mrna_count)
exprs_mrna_tpm <- exprs_mrna_tpm %>%
  as_tibble() %>% # tibble不支持row name，我竟然才发现！
  mutate(meanrow = rowMeans(.[, -1]), .before = 2) %>%
  arrange(desc(meanrow)) %>%
  distinct(sample_id, .keep_all = T) %>%
  select(-meanrow) %>%
  column_to_rownames(var = "sample_id") %>%
  as.data.frame()
dim(exprs_mrna_tpm)

# add gene symbol
symbol_mrna <- rowData(se_mrna)$gene_name
head(symbol_mrna)
exprs_mrna_count <- cbind(data.frame(symbol_mrna), as.data.frame(t(exprs_mrna_count)))
exprs_mrna_tpm <- cbind(data.frame(symbol_mrna), as.data.frame(t(exprs_mrna_tpm)))

# remove duplicate genes
exprs_mrna_count_rmdup <- exprs_mrna_count %>%
  as_tibble() %>%
  mutate(meanrow = rowMeans(.[, -1]), .before = 2) %>%
  arrange(desc(meanrow)) %>%
  distinct(symbol_mrna, .keep_all = T) %>%
  select(-meanrow) %>%
  column_to_rownames(var = "symbol_mrna") %>%
  as.data.frame()
dim(exprs_mrna_count_rmdup)
exprs_mrna_tpm_rmdup <- exprs_mrna_tpm %>%
  as_tibble() %>%
  mutate(meanrow = rowMeans(.[, -1]), .before = 2) %>%
  arrange(desc(meanrow)) %>%
  distinct(symbol_mrna, .keep_all = T) %>%
  select(-meanrow) %>%
  column_to_rownames(var = "symbol_mrna") %>%
  as.data.frame()
dim(exprs_mrna_tpm_rmdup)
# exprs_mrna_tpm_rmdup[1:5,1:5]
# write.table(exprs_mrna_tpm_rmdup, "../../result/TCGA/TCGA_TPM_exprs.tsv", sep = "\t", quote = FALSE)

# ---- metadata
tcga_529_oncoprint <- read.delim("~/reference/TCGA/UCEC/ucec_tcga_pan_can_atlas_2018(529)/PATIENT_DATA_oncoprint.tsv")
tcga_529_oncoprint[1:10, 1:10]

# ---- phenotype
rownames(tcga_529_oncoprint) <- paste(tcga_529_oncoprint[, 1], tcga_529_oncoprint[, 2], sep = "_")
tcga_529_oncoprint <- tcga_529_oncoprint[, c(3:ncol(tcga_529_oncoprint))]
tcga_529_oncoprint <- tcga_529_oncoprint %>%
  t() %>%
  as.data.frame()
phenotype_ec <- tcga_529_oncoprint[tcga_529_oncoprint$`Subtype_CLINICAL` != "", ]
head(phenotype_ec)
rownames(phenotype_ec) <- gsub(r"(\.)", "-", rownames(phenotype_ec))
write_csv(phenotype_ec[, order(colnames(phenotype_ec))], "../../result/TCGA/phenotype_ec.csv")

# ---- sample ID
ucec.id <- rownames(phenotype_ec)
length(ucec.id)

# ---- expression
head(exprs_mrna_count_rmdup)
bk.dat <- t(exprs_mrna_count_rmdup[, ucec.id])
head(bk.dat)
dim(bk.dat)

save(bk.dat, phenotype_ec, file = "../TCGA_UCEC_529.rda")





# single-cell EC data
# ---- bk.dat
load("../TCGA_UCEC_529.rda")
load("../19_final_analysis/scEC_data.rda")

dim(bk.dat)
dim(phenotype_ec)

# ---- sc.dat
sample_ca <- subset(sample.all, subset = Tissue %in% c("CA"))
sample_ca <- subset(sample_ca, subset = CellSubType %in% c("Cancer"), invert = T)

sample_ca$CellSubType <- factor(
  sample_ca$CellSubType,
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
    "Cancer_01_ciliated_FOXJ1", "Cancer_02_glandular_SCGB2A1", "Cancer_03_luminal_LCN2", "Cancer_04_luminal_ERBB4",
    "Cancer_05_basal_KRT5", "Cancer_06_EMT_like_SULF1", "Cancer_07_fibro_like_COL3A1",
    "Cancer_08_glandular_CXCL14"
  )
)



# BayesPrism Deconvolution
sc.dat <- as.matrix(t(GetAssayData(sample_ca, assay = "RNA", slot = "counts")))
dim(sc.dat)
sc.dat[1:5, 1:5]

# ---- cell type annotation
cell.type.labels <- as.character(sample_ca$CellType)
table(is.na(cell.type.labels))
sort(unique(cell.type.labels))

cell.state.labels <- as.character(sample_ca$CellSubType)
table(is.na(cell.state.labels))
sort(unique(cell.state.labels))

sc.dat.filtered <- cleanup.genes(
  input = sc.dat,
  input.type = "count.matrix",
  species = "hs",
  gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"),
  exp.cells = 5
)
dim(sc.dat.filtered)

# note this function only works for human data. For other species, you are advised to make plots by yourself.
pdf("../../figures/BayesPrism_bk_vs_sc_tcga529.pdf", width = 12, height = 4)
plot.bulk.vs.sc(
  sc.input = sc.dat.filtered,
  bulk.input = bk.dat
)
dev.off()

sc.dat.filtered.pc <- select.gene.type(sc.dat.filtered, gene.type = "protein_coding")
dim(sc.dat.filtered.pc)

# performing pair-wise t test for cell states from different cell types
table(sample_ca$CellSubType) %>% sort()
diff.exp.stat <- get.exp.stat(
  sc.dat = sc.dat[, colSums(sc.dat > 0) > 6],
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  psuedo.count = 0.1, # a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
  cell.count.cutoff = 50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
  n.cores = 50
)
sc.dat.filtered.pc.sig <- select.marker(
  sc.dat = sc.dat.filtered.pc,
  stat = diff.exp.stat,
  pval.max = 0.05,
  lfc.min = 0.1
)
dim(sc.dat.filtered.pc.sig)

# ---- save TPM matrix for Kannsandra deconvolution
exprs_mrna_tpm_rmdup_filtered <- exprs_mrna_tpm_rmdup[colnames(sc.dat.filtered.pc.sig), ]
exprs_mrna_tpm_rmdup_filtered <- exprs_mrna_tpm_rmdup_filtered[!is.na(exprs_mrna_tpm_rmdup_filtered[, 1]), ucec.id]
dim(exprs_mrna_tpm_rmdup_filtered)
exprs_mrna_tpm_rmdup_filtered[1:5, 1:5]
write.table(exprs_mrna_tpm_rmdup_filtered, "../../result/TCGA/TCGA_TPM_exprs_BP_DE.txt", sep = "\t", quote = FALSE)

myPrism <- new.prism(
  reference = sc.dat.filtered.pc.sig,
  mixture = bk.dat,
  input.type = "count.matrix",
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  key = "Cancer",
  outlier.cut = 0.01,
  outlier.fraction = 0.1,
)

plan("multicore", workers = 10)
plan()
options(future.globals.maxSize = 20 * 1000 * 1024^2)
bp.res <- run.prism(prism = myPrism, n.cores = 50)
slotNames(bp.res)
save(bp.res, file = "../bp_res_tcga529_BPGenes.rda")
