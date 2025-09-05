# WES data analysis
# List of packages required
depp <- c(
  "languageserver", "httpgd", "Seurat", "DoubletFinder", "dplyr", "DT", "lme4",
  "ggplot2", "patchwork", "data.table", "tidyr", "tidyverse", "tidygraph",
  "ggraph", "SeuratWrappers", "doParallel", "openxlsx", "viridis", "rstatix",
  "maftools"
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





# ---- build MAF files
## mutect2 annovar files
samples <- list.files("/staging/qlab/wes/annotation/", pattern = "EC")
data_path <- "/staging/qlab/wes/annotation"
cancer_genes <- read_tsv("~/reference/wes/hotspots_chang2018.tsv") %>%
  dplyr::pull(Gene) %>%
  unique()
oncokb_genes <- read_tsv("~/reference/wes/TCGA_EC_Mutated_Genes.tsv") %>%
  dplyr::filter(OncoKB == "Yes") %>%
  dplyr::pull(Gene)


# ---- merge annovar csv files for different annotation method
assays <- c("mutect2", "freebayes", "manta", "strelka", "tiddit")
for (assay in assays) {
  df_all <- tibble()
  for (mysample in samples) {
    # mysample  <- "EC20CA"
    cat(mysample, "\n")
    if (assay %in% c("mutect2", "freebayes")) {
      df <- suppressMessages(read_csv(paste0(data_path, "/", mysample, "/", mysample, ".", assay, ".filtered.annovar.hg38_multianno.csv.gz")))
    } else if (assay %in% c("manta")) {
      df <- suppressMessages(read_csv(paste0(data_path, "/", mysample, "/", mysample, ".", assay, ".tumor_sv.annovar.hg38_multianno.csv.gz")))
    } else if (assay %in% c("strelka")) {
      df <- suppressMessages(read_csv(paste0(data_path, "/", mysample, "/", mysample, ".", assay, ".variants.annovar.hg38_multianno.csv.gz")))
    } else if (assay %in% c("tiddit")) {
      df <- suppressMessages(read_csv(paste0(data_path, "/", mysample, "/", mysample, ".", assay, ".annovar.hg38_multianno.csv.gz")))
    }
    df_filtered <- df %>%
      dplyr::rename(
        Func.refGene = Func.refGeneWithVer,
        Gene.refGene = Gene.refGeneWithVer,
        GeneDetail.refGene = GeneDetail.refGeneWithVer,
        ExonicFunc.refGene = ExonicFunc.refGeneWithVer,
        AAChange.refGene = AAChange.refGeneWithVer
      ) %>%
      dplyr::select(
        Chr, Start, End, Ref, Alt, Func.refGene, Gene.refGene, GeneDetail.refGene, ExonicFunc.refGene,
        AAChange.refGene, cytoband, AF, AF_popmax, AF_male, AF_female, AF_eas,
        CLNALLELEID, CLNDN, CLNSIG, ICGC_Id
      ) %>%
      dplyr::mutate(Tumor_Sample_Barcode = mysample)

    df_all <- bind_rows(df_all, df_filtered)
  }
  write_csv(df_all, paste0("wes_annotation_", assay, "_merged.csv"))
}

# all anno
# assays <- c("mutect2", "freebayes", "manta", "strelka", "tiddit")
assay <- "mutect2"
df_all <- read_csv(paste0("wes_annotation_", assay, "_merged.csv"))

# exonic anno only
df_all %>%
  dplyr::filter(Func.refGene == "exonic") %>%
  write_csv(paste0("wes_annotation_", assay, "_exonic.csv"))

# exonic and cancer only
df_all %>%
  dplyr::filter(Func.refGene == "exonic" | CLNSIG %in% c("Pathogenic", "Likely_pathogenic")) %>%
  dplyr::filter(Gene.refGene %in% cancer_genes | Gene.refGene %in% oncokb_genes) %>%
  write_csv(paste0("wes_annotation_", assay, "_exonic_cancer.csv"))





# ---- Annovar to MAF
maf_file_all <- annovarToMaf(
  paste0("wes_annotation_", assay, "_merged.csv"),
  Center = "ZJUHW",
  refBuild = "hg38", table = "refGene", sep = ",", tsbCol = "Tumor_Sample_Barcode"
)
maf_file_exonic <- annovarToMaf(
  paste0("wes_annotation_", assay, "_exonic.csv"),
  Center = "ZJUHW",
  refBuild = "hg38", table = "refGene", sep = ",", tsbCol = "Tumor_Sample_Barcode"
)
maf_file_cancer <- annovarToMaf(
  paste0("wes_annotation_", assay, "_exonic_cancer.csv"),
  Center = "ZJUHW",
  refBuild = "hg38", table = "refGene", sep = ",", tsbCol = "Tumor_Sample_Barcode"
)





# ---- read MAF files
maf_all <- read.maf(maf_file_all, useAll = FALSE)
getSampleSummary(maf_all)
maf_exonic <- read.maf(maf_file_exonic, useAll = FALSE)
getSampleSummary(maf_exonic)
maf_cancer <- read.maf(maf_file_cancer, useAll = FALSE)
getSampleSummary(maf_cancer)

# str(maf_all)
maf_all@data %>%
  filter(Hugo_Symbol == "POLQ") %>%
  select(Tumor_Sample_Barcode, Gene.refGene, Func.refGene, Variant_Classification, ExonicFunc.refGene, CLNSIG, ICGC_Id) %>%
  arrange(Tumor_Sample_Barcode)

df_all %>%
  filter(Func.refGene == "exonic") %>%
  filter(ExonicFunc.refGene == "nonsynonymous SNV") %>%
  select(Tumor_Sample_Barcode, Gene.refGene, Func.refGene, ExonicFunc.refGene, CLNSIG, ICGC_Id) %>%
  # filter(Tumor_Sample_Barcode == "EC12CA") %>%
  filter(Gene.refGene %in% c("POLQ"))

df_all %>%
  filter(Func.refGene == "exonic") %>%
  # filter(ExonicFunc.refGene == "nonsynonymous SNV") %>%
  # filter(Tumor_Sample_Barcode == "EC12CA") %>%
  filter(Gene.refGene %in% c("TP53", "POLE", "PMS2", "MLH1", "MSH2", "MSH6")) %>%
  write_csv(paste0("wes_annotation_", assay, "_exonic_cancer_TP53_POLE_MMRd.csv"))
arrange(Tumor_Sample_Barcode) %>%
  # select(-c(GeneDetail.refGene, CLNDN))
  select(-c(GeneDetail.refGene, CLNDN, AAChange.refGene, ExonicFunc.refGene))


## add clinical data
sample.rdb <- readRDS("../09_SCENIC/sample_integrated_rdb.rds")
sample.rdb <- subset(sample.rdb, subset = Tissue %in% c("CA", "PCA"))
sample.rdb <- subset(sample.rdb, subset = CellSubType %in% c("Cancer"), invert = T)





# ---- make clinical data data frame
colnames(sample.rdb@meta.data)
sample.clinical <- sample.rdb@meta.data[, c("Patient", "Tissue", "Age", "Subtype", "Staging", "Grade", "BMI.y", "WeightState", "molecular_class")] %>%
  unique() %>%
  filter(Tissue == "CA")
sample.clinical$Tumor_Sample_Barcode <- paste0(sample.clinical$Patient, sample.clinical$Tissue)
head(sample.clinical)

# ---- add clinical data
maf_all <- read.maf(maf_file_all, useAll = FALSE, clinicalData = sample.clinical)
maf_exonic <- read.maf(maf_file_exonic, useAll = FALSE, clinicalData = sample.clinical)
maf_cancer <- read.maf(maf_file_cancer, useAll = FALSE, clinicalData = sample.clinical)

# ---- maf summary
# Shows sample summry.
getSampleSummary(maf_cancer)
# Shows gene summary.
getGeneSummary(maf_cancer)
# shows clinical data associated with samples
getClinicalData(maf_cancer)
# Shows all fields in MAF
getFields(maf_cancer)
# Writes maf summary to an output file with basename maf.
write.mafSummary(maf = maf_cancer, basename = "maf_cancer")


## add CNV
amp_genes <- read_tsv("~/reference/wes/TCGA_EC_CNA_Genes.tsv") %>%
  dplyr::mutate(Pct = Altered / Profiled) %>%
  dplyr::filter(OncoKB == "Yes", Pct > 0.01, CNA == "AMP") %>%
  dplyr::arrange(-Pct) %>%
  dplyr::pull(Gene)
del_genes <- read_tsv("~/reference/wes/TCGA_EC_CNA_Genes.tsv") %>%
  dplyr::mutate(Pct = Altered / Profiled) %>%
  dplyr::filter(OncoKB == "Yes", Pct > 0.01, CNA == "HOMDEL") %>%
  dplyr::arrange(-Pct) %>%
  dplyr::pull(Gene)

cnv_genes.selected <- c(
  "POLE", "MLH1", "PTEN", "ESR1", "PGR", "PIK3CA", "ARID1A", "PIK3R1", "KRAS", "CTNNB1", "BRCA2", "TP53",
  "PMS2", "MSH2", "MSH6", "AR", "PPP2R1A", "FBXW7", "ERBB2", "CHD4", "CCNE1", "ATM", "BRCA1", "SPOP", "KMT2D"
)

out_all.df <- tibble()
out.df <- tibble()
for (sample in samples) {
  # sample  <- "OC1CA"
  cat(sample, "\n")
  cnv.df <- read_tsv(paste0("/staging/qlab/wes/variant_calling/cnvkit/", sample, "/", sample, ".call.cns"))

  alt.df <- cnv_genes.selected %>%
    purrr::map_df(function(cnv_gene) {
      cnv.df %>%
        rowwise() %>%
        dplyr::mutate(gene_list = list(str_split(gene, ",", simplify = TRUE))) %>%
        dplyr::filter(cnv_gene %in% unlist(gene_list), cn != 2)
    })

  out.df <- 1:nrow(alt.df) %>%
    purrr::map_df(function(i) {
      genes <- unlist(alt.df$gene_list[i])
      cn <- alt.df$cn[i]
      tibble(Gene = genes, Sample_name = sample, CN = cn)
    })

  out_all.df <- bind_rows(out_all.df, out.df)
}

out_all.df <- out_all.df[!is.na(out_all.df$CN), ]
out_maf.df <- out_all.df %>%
  dplyr::mutate(CN = case_when(
    CN == 0 ~ "Del",
    CN == 1 ~ "Shallow_Del",
    CN == 3 ~ "Shallow_Amp",
    CN > 3 ~ "Amp"
  ))
out_maf.df <- out_maf.df[!duplicated(out_maf.df), ]
out_maf.df

out_maf.df %>%
  dplyr::group_by(Gene, Sample_name) %>%
  dplyr::summarise(num = n())

maf_cnv <- read.maf(maf_file_all, cnTable = out_maf.df, clinicalData = sample.clinical, useAll = FALSE)
maf_cnv <- read.maf(maf_file_cancer, cnTable = out_maf.df, clinicalData = sample.clinical, useAll = FALSE)

oncogene_all <- read_tsv("~/reference/wes/TCGA_EC_CNA_Genes.tsv") %>%
  dplyr::mutate(Pct = Altered / Profiled) %>%
  dplyr::filter(Pct > 0.05, OncoKB == "Yes") %>%
  dplyr::pull(Gene) %>%
  union(oncokb_genes)

out_maf.df %>%
  dplyr::filter(CN == "Amp") %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(num = n()) %>%
  dplyr::filter(Gene %in% oncogene_all) %>%
  dplyr::arrange(num)

out_maf.df %>%
  dplyr::filter(CN == "Del") %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(num = n()) %>%
  dplyr::filter(Gene %in% oncogene_all) %>%
  dplyr::arrange(num)

out_maf.df %>%
  dplyr::filter(Gene == "POLE")

out_maf.df %>%
  dplyr::filter(CN %in% c("Amp", "Shallow_Amp")) %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(num = n()) %>%
  dplyr::filter(Gene %in% oncogene_all) %>%
  dplyr::arrange(num)

out_maf.df %>%
  dplyr::filter(CN %in% c("Del", "Shallow_Del")) %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(num = n()) %>%
  dplyr::filter(Gene %in% oncogene_all) %>%
  dplyr::arrange(num)

pdf("../../figures/wes/oncoplot_cnv_select.pdf", height = 10, width = 6)
mycolors <- RColorBrewer::brewer.pal(n = 11, name = "Paired")
# [1] "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00" "#CAB2D6" "#6A3D9A" "#FFFF99"
names(mycolors) <- c(
  "Shallow_Del", "Del", "Frame_Shift_Del", "Missense_Mutation", "Shallow_Amp", "Amp",
  "Nonsense_Mutation", "Multi_Hit", "Frame_Shift_Ins",
  "Splice_Site", "Complex_Event"
)
oncoplot(
  maf = maf_cnv, top = 50, keepGeneOrder = FALSE, colors = mycolors,
  showTumorSampleBarcodes = TRUE, draw_titv = TRUE,
  genes = cnv_genes.selected
)
dev.off()


## MAF visulization
# ---- visulization
pdf(paste0("../../figures/wes/maf_summary_", assay, ".pdf"), width = 12, height = 8)
plotmafSummary(maf = maf_all, rmOutlier = TRUE, addStat = "median", dashboard = TRUE, titvRaw = FALSE)
plotmafSummary(maf = maf_exonic, rmOutlier = TRUE, addStat = "median", dashboard = TRUE, titvRaw = FALSE)
plotmafSummary(maf = maf_cancer, rmOutlier = TRUE, addStat = "median", dashboard = TRUE, titvRaw = FALSE)
dev.off()

pdf(paste0("../../figures/wes/oncoplot_", assay, ".pdf"), width = 10, height = 20)
oncoplot(maf = maf_all, top = 100, showTumorSampleBarcodes = TRUE)
oncoplot(maf = maf_exonic, top = 50, showTumorSampleBarcodes = TRUE)
oncoplot(maf = maf_cancer, top = 50, showTumorSampleBarcodes = TRUE)
oncoplot(maf = maf_cancer, top = 150, showTumorSampleBarcodes = TRUE)
dev.off()
