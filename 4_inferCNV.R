# EC Data Analysis for inferCNVpy
# List of packages required
depp <- c(
  "languageserver", "httpgd", "Seurat", "DoubletFinder", "dplyr", "DT", "SeuratDisk", "ggplot2",
  "patchwork", "data.table", "tidyr", "tidyverse", "tidygraph", "ggraph", "SeuratWrappers",
  "doParallel", "openxlsx", "viridis", "rstatix", "KernSmooth", "BiocParallel", "reticulate"
)
# Check if the packages were installed if not install
depp_new <- depp[!(depp %in% installed.packages())]
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
  "scater", "slingshot", "destiny", "ggrepel", "GSEABase", "ggpubr",
  "org.Hs.eg.db", "clusterProfiler", "AUCell", "GENIE3"
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





# Convert Seurat to loom
sample.rdb <- readRDS("../sample_all_original_rdb.rds")
sample.rdb <- subset(sample.rdb, subset = Tissue == "CA")
table(sample.rdb$Tissue)
table(sample.rdb$CellType)


DefaultAssay(sample.rdb) <- "RNA"
sample.rdb$CellType <- as.character(sample.rdb$CellType)

# first convert to h5Seurat file
SeuratDisk::SaveH5Seurat(sample.rdb, filename = "../sample_integrated_ca_celltype.h5Seurat")

# then convert to h5ad, which can be read by scanpy
SeuratDisk::Convert(
  source = "../sample_integrated_ca_celltype.h5Seurat",
  assay = "RNA",
  dest = "../sample_integrated_ca_celltype.h5ad"
)


DefaultAssay(sample.rdb) <- "RNA"
sample.rdb.sub <- subset(sample.rdb, subset = CellSubType %in% c(
  "Epi_01_SCGB2A1",
  "Epi_02_CXCL1",
  "Epi_03_WNT7A",
  "Epi_04_LCN2",
  "Epi_05_FOXJ1",
  "Epi_06_MMP7",
  "Cancer_01_ciliated_FOXJ1",
  "Cancer_02_glandular_SCGB2A1",
  "Cancer_03_luminal_LCN2",
  "Cancer_04_luminal_ERBB4",
  "Cancer_05_basal_KRT5",
  "Cancer_06_EMT_like_SULF1",
  "Cancer_07_fibro_like_COL3A1",
  "Cancer_08_glandular_CXCL14"
))

# first convert to h5Seurat file
SeuratDisk::SaveH5Seurat(sample.rdb.sub, filename = "../sample_integrated_all_Cancer.h5Seurat")

# then convert to h5ad, which can be read by scanpy
SeuratDisk::Convert(
  source = "../sample_integrated_all_Cancer.h5Seurat",
  assay = "RNA",
  dest = "../sample_integrated_all_Cancer.h5ad"
)




# ---- inferCNVpy(run in python)
import scanpy as sc
import loompy as lp
import numpy as np

adata = sc.read_h5ad("../sample_integrated_all_Cancer.h5ad")

row_attrs = { "Gene": np.array(adata.var_names), }
col_attrs = { "CellID":np.array(adata.obs_names),
   "nGene":np.array(np.sum(adata.X.transpose()>0, axis=0)).flatten(),
   "nUMI":np.array(np.sum(adata.X.transpose(),axis=0)).flatten(), }

adata.var_names

lp.create("../sample_integrated_all_Cancer.loom",
          adata.X.transpose(),
          row_attrs, 
          col_attrs)


## remove genes expressed in less than 20 cells
tmp <- GetAssayData(object = sample.rdb, slot = "data")
sc.obj <- CreateSeuratObject(counts = tmp, min.cells = 20)
genes.kp <- rownames(sc.obj)
sample.rdb.kp <- sample.rdb[genes.kp, , ]

DefaultAssay(sample.rdb.kp) <- "RNA"
sample.rdb.kp$CellSubType <- as.character(sample.rdb.kp$CellSubType)

# first convert to h5Seurat file
SaveH5Seurat(sample.rdb.kp, filename = "./sample_integrated_rdb_r20.h5Seurat")
# then convert to h5ad, which can be read by scanpy
Convert(
  source = "./sample_integrated_rdb_r20.h5Seurat",
  assay = "RNA",
  dest = "./sample_integrated_rdb_r20.h5ad"
)

import scanpy as sc
import loompy as lp
import numpy as np

adata = sc.read_h5ad("../sample_integrated_all_celltype.h5ad")

row_attrs = { "Gene": np.array(adata.var_names), }
col_attrs = { "CellID":np.array(adata.obs_names),
   "nGene":np.array(np.sum(adata.X.transpose()>0, axis=0)).flatten(),
   "nUMI":np.array(np.sum(adata.X.transpose(),axis=0)).flatten(), }

adata.var_names

lp.create("./sample_integrated_rdb_r20.loom",
          adata.X.transpose(),
          row_attrs, 
          col_attrs)





# infer ITH score - python
## loading single cell dataset
import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt

adata = sc.read_h5ad("../sample_integrated_ca_celltype.h5ad")

cnv.io.genomic_position_from_gtf(
  gtf_file = "~/reference/10X/refdata-gex-GRCh38-2020-A/genes/genes.gtf", 
  adata = adata, 
  gtf_gene_id = 'gene_name')

adata.var.loc[:, ["features", "chromosome", "start", "end"]].head()

import pandas as pd
adata.var.loc[:, ["features", "chromosome", "start", "end"]]
data = adata.var.loc[:, ["features", "chromosome", "start", "end"]]
pd.DataFrame(data).to_csv('gene_pos.txt', index=True, mode='a')


## running inferCNV
adata.obs["CellType"]
cnv.tl.infercnv(
    adata = adata,
    reference_key="CellType",
    reference_cat=["Bcell", "CD4T", "CD8T", "NK", "Treg"],
    window_size=100,
)
cnv.pl.chromosome_heatmap(adata, groupby = "Sample")
plt.savefig("../../figures/infercnv_chrom_heatmap.pdf")

cnv.tl.pca(adata)
cnv.pp.neighbors(adata)
cnv.tl.leiden(adata)

adata.obs.columns
cnv.tl.cnv_score(
  adata = adata, 
  # groupby='cnv_leiden',
  groupby='CellType',
  use_rep='cnv', 
  key_added='cnv_score', 
  inplace=True, 
  obs_key=None
)
pd.DataFrame(adata.obs).to_csv('./metadata_cnv_ca_celltype.txt', index=True, mode='w')

import pandas as pd
data1 = adata.obsm["X_cnv"]
data1_dense = data1.todense()
data1df = pd.DataFrame(data1_dense)
pd.DataFrame(data1df).to_csv('./inferCNV.txt', index=True, mode='w')

norm_cells <- rownames(sample.rdb@meta.data[sample.rdb$CellSubType %in% ref_group_names, ])

cnv.tl.copykat(
  adata = adata, 
  gene_ids = "S", 
  segmentation_cut = 0.1, 
  distance = "euclidean", 
  s_name = "copykat_result", 
  min_genes_chr = 5, 
  key_added = "copykat", 
  inplace = True, 
  norm_cell_names = norm_cells)


ca_cells = adata.obs.index[adata.obs.CellSubType.isin([
  "Epi_01_SCGB2A1", "Epi_02_CXCL1", "Epi_03_WNT7A", "Epi_04_LCN2", "Epi_05_FOXJ1", "Epi_06_MMP7", 
  "Cancer_01_ciliated_FOXJ1", "Cancer_02_glandular_SCGB2A1", "Cancer_03_luminal_LCN2",
  "Cancer_04_luminal_ERBB4", "Cancer_05_basal_KRT5", "Cancer_06_EMT_like_SULF1",
  "Cancer_07_fibro_like_COL3A1", "Cancer_08_glandular_CXCL14"])]
adata_ca = adata[ca_cells]
cnv.tl.ithcna(
  adata = adata_ca, 
  groupby = "Sample",
  use_rep='X_cnv', 
  key_added='ithcna_sample', 
  inplace=True
)
cnv.tl.ithgex(
  adata = adata_ca, 
  groupby = "Sample", 
  use_raw=None, 
  layer=None, 
  inplace=True, 
  key_added='ithgex_sample'
)
import pandas as pd
pd.DataFrame(adata_ca.obs).to_csv('../metadata_sample.txt', index=True, mode='w')
