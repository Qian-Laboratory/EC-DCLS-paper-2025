# NMF analysis for CCLE

# List of packages required
depp <- c(
  "languageserver", "httpgd", "Seurat", "DoubletFinder", "dplyr", "DT", "lme4",
  "ggplot2", "patchwork", "data.table", "tidyr", "tidyverse", "tidygraph",
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
  "scater", "slingshot", "destiny", "ggrepel",
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





library(NMF)
matrix <- exprs_matrix
matrix <- exprs_matrix[rownames(exprs_matrix) %in% var_genes, ]
colnames(matrix) <- ccle$cell_line_name
table(rowSums(matrix) == 0)
matrix <- matrix[rowSums(matrix) > 0, ]
dim(matrix)

# plan("multicore", workers = 10)
# plan()
# options(future.globals.maxSize = 20 * 1000 * 1024^2)
ranks <- 2:10
estim.coad <- nmf(matrix, ranks, nrun = 50)
duplicated(colnames(matrix))
pdf("../../figures/CCLE_nmf_feature.pdf")
print(plot(estim.coad))
dev.off()

seed <- 20230410
nmf.rank6 <- nmf(matrix,
  rank = 4,
  nrun = 200,
  seed = seed,
  method = "brunet"
)

# set color
jco <- c("#2874C5", "#EABF00", "#C6524A", "#73c630")
index <- extractFeatures(nmf.rank6, "max")
sig.order <- unlist(index)
NMF.Exp.rank6 <- matrix[sig.order, ]
NMF.Exp.rank6 <- na.omit(NMF.Exp.rank6)
group <- predict(nmf.rank6)
table(group)


pdf("../../figures/CCLE_nmf_consensus_rank4_confirm.pdf", width = 8, height = 8.8)
print(
  consensusmap(nmf.rank6,
    labRow = NA,
    annCol = data.frame(
      "cluster" = group[colnames(NMF.Exp.rank6)],
      "seurat_cluster" = ccle$RNA_snn_res.1.2[match(ccle$cell_line_name, colnames(NMF.Exp.rank6))]
    ),
    annColors = list(
      cluster = c("1" = jco[1], "2" = jco[2], "3" = jco[3], "4" = jco[4]),
      seurat_cluster = c("0" = jco[1], "1" = jco[2], "2" = jco[3], "3" = jco[4], "4" = "#e3843b", "5" = "grey")
    )
  )
)
dev.off()



