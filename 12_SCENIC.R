# SCENIC analysis

# List of packages required
depp <- c(
  "languageserver", "httpgd", "Seurat", "DoubletFinder", "dplyr", "DT",
  "ggplot2", "patchwork", "data.table", "tidyr", "tidyverse", "tidygraph",
  "ggraph", "SeuratWrappers", "SeuratDisk", "doParallel", "openxlsx", "viridis",
  "rstatix", "SCENIC", "SCopeLoomR", "KernSmooth", "BiocParallel"
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
  "org.Hs.eg.db", "clusterProfiler", "AUCell", "RcisTarget", "GENIE3"
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
sample.rdb <- readRDS("./sample_integrated_rdb.rds")
DefaultAssay(sample.rdb) <- "RNA"
# first convert to h5Seurat file
SaveH5Seurat(sample.rdb, filename = "../sample_integrated_Cancer.h5Seurat")
# then convert to h5ad, which can be read by scanpy
Convert(
  source = "../sample_integrated_Cancer.h5Seurat",
  assay = "RNA",
  dest = "../sample_integrated_Cancer.h5ad"
)




# ---- run in python
import scanpy as sc
import loompy as lp
import numpy as np

adata = sc.read_h5ad("../sample_integrated_Cancer.h5ad")

row_attrs = { "Gene": np.array(adata.var_names), }
col_attrs = { "CellID":np.array(adata.obs_names),
   "nGene":np.array(np.sum(adata.X.transpose()>0, axis=0)).flatten(),
   "nUMI":np.array(np.sum(adata.X.transpose(),axis=0)).flatten(), }

adata.var_names

lp.create("../sample_integrated_Cancer.loom",
          adata.X.transpose(),
          row_attrs, 
          col_attrs)


# pySCENIC
# 用grnboost2计算adjacencies文件
system("
pyscenic grn \
--num_workers 20 \
--output /share/home/qlab/qlab_cdj/project/project_10x_ec/analysis/adj_Cancer.tsv \
--method grnboost2 \
/share/home/qlab/qlab_cdj/project/project_10x_ec/analysis/sample_integrated_Cancer.loom \
/share/reference/human/hs_hgnc_tfs.txt
")


## pySCENIC downstream analyses
# import dependencies
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import json
import base64
import zlib
import seaborn as sns
from pyscenic.plotting import plot_binarization
from pyscenic.export import add_scenic_metadata
from pyscenic.cli.utils import load_signatures
import matplotlib as mpl
import matplotlib.pyplot as plt
from scanpy.plotting._tools.scatterplots import plot_scatter


adata = sc.read_h5ad("../sample_integrated_Cancer.h5ad")
row_attrs = { "Gene": np.array(adata.var_names), }
col_attrs = { "CellID":np.array(adata.obs_names),
   "nGene":np.array(np.sum(adata.X.transpose()>0, axis=0)).flatten(),
   "nUMI":np.array(np.sum(adata.X.transpose(),axis=0)).flatten(), }
cellAnnot = adata.obs
cellAnnot.head()

# scenic output
f_final_loom = "../sample_integrated_Cancer.loom"
lf = lp.connect( f_final_loom, mode='r', validate=False )
exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID).T

sig = load_signatures("../reg_Cancer.tsv")
adata = add_scenic_metadata(adata, auc_mtx, sig)
adata.head()

AUCELL_MTX_FNAME = "../auc_mtx.csv"
auc_mtx = pd.read_csv(AUCELL_MTX_FNAME, index_col = 0)
auc_mtx.head()

rss_cellType = regulon_specificity_scores(auc_mtx, cellAnnot["CellSubType"])
rss_cellType

from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
import matplotlib.pyplot as plt
from adjustText import adjust_text
import seaborn as sns
from pyscenic.binarization import binarize

rss_cellType = regulon_specificity_scores(auc_mtx, cellAnnot["CellSubType"])

import pandas as pd
pd.DataFrame(rss_cellType).to_csv('../rss_cellType.txt', index=True, mode='w')

## RSS panel plot with all cell types
cats = sorted(list(set(cellAnnot['CellSubType'])))

fig = plt.figure(figsize=(15, 15))
for c,num in zip(cats, range(1,len(cats)+1)):
    x=rss_cellType.T[c]
    ax = fig.add_subplot(3,5,num)
    plot_rss(rss_cellType, c, top_n=10, max_n=None, ax=ax)
    ax.set_ylim( x.min()-(x.max()-x.min())*0.05 , x.max()+(x.max()-x.min())*0.05 )
    for t in ax.texts:
        t.set_fontsize(12)
    ax.set_ylabel('')
    ax.set_xlabel('')
    adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-',color='lightgrey'), precision=0.001 )
 
fig.text(0.5, 0.0, 'Regulon', ha='center', va='center', size='x-large')
fig.text(0.00, 0.5, 'Regulon specificity score (RSS)', ha='center', va='center', rotation='vertical', size='x-large')
plt.tight_layout()
plt.rcParams.update({
    'figure.autolayout': True,
        'figure.titlesize': 'large' ,
        'axes.labelsize': 'medium',
        'axes.titlesize':'large',
        'xtick.labelsize':'medium',
        'ytick.labelsize':'medium'
        })
plt.savefig("../../figures/Cancer_cellType_RSS_top10.pdf", dpi = 600, bbox_inches = "tight")

## Select the top 5 regulons from each cell type
topreg = []
for i,c in enumerate(cats):
    topreg.extend(
        list(rss_cellType.T[c].sort_values(ascending=False)[:5].index)
    )
topreg = list(set(topreg))
topreg

## Generate a Z-score for each regulon to enable comparison between regulons
auc_mtx_Z = pd.DataFrame( index=auc_mtx.index )
for col in list(auc_mtx.columns):
    auc_mtx_Z[ col ] = ( auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof=0)
#auc_mtx_Z.sort_index(inplace=True)
auc_mtx
auc_mtx_Z

## CELL TYPE SPECIFIC REGULATORS - Z-SCORE
To find cell type specific regulators we use a Z score (i.e. the average AUCell score for the cells of a give type are standardized using the overall average AUCell scores and its standard deviation).
df_obs = adata.obs
signature_column_names = list(df_obs.select_dtypes('number').columns)
signature_column_names = list(filter(lambda s: s.startswith('Regulon('), signature_column_names))
df_scores = df_obs[signature_column_names + ['CellSubType']]
df_results = ((df_scores.groupby(by='CellSubType').mean() - df_obs[signature_column_names].mean())/ df_obs[signature_column_names].std()).stack().reset_index().rename(columns={'level_1': 'regulon', 0:'Z'})
df_results['regulon'] = list(map(lambda s: s[8:-1], df_results.regulon))
df_results[(df_results.Z >= 3.0)].sort_values('Z', ascending=False).head()

df_heatmap = pd.pivot_table(data=df_results[df_results.Z >= 1.5].sort_values('Z', ascending=False), index='CellSubType', columns='regulon', values='Z')
#df_heatmap.drop(index='Myocyte', inplace=True) # We leave out Myocyte because many TFs are highly enriched (becuase of small number of cells).
fig, ax1 = plt.subplots(1, 1, figsize=(20, 10))
sns.heatmap(df_heatmap, ax=ax1, annot=True, fmt=".1f", linewidths=.7, cbar=False, square=True, linecolor='gray', 
            cmap="YlGnBu", annot_kws={"size": 10})
ax1.set_ylabel('')
plt.savefig("../../figures/Cancer_cellType_heatmap_regulons.pdf", dpi = 600, bbox_inches = "tight")

## Generate a heatmap
def palplot(pal, names, colors=None, size=1):
    n = len(pal)
    f, ax = plt.subplots(1, 1, figsize=(n * size, size))
    ax.imshow(np.arange(n).reshape(1, n),
              cmap=mpl.colors.ListedColormap(list(pal)),
              interpolation="nearest", aspect="auto")
    ax.set_xticks(np.arange(n) - .5)
    ax.set_yticks([-.5, .5])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    colors = n * ['k'] if colors is None else colors
    for idx, (name, color) in enumerate(zip(names, colors)):
        ax.text(0.0+idx, 0.0, name, color=color, horizontalalignment='center', verticalalignment='center')
    return f

colors = sns.color_palette('bright',n_colors=len(cats) )
colorsd = dict( zip( cats, colors ))
colormap = [ colorsd[x] for x in cellAnnot['CellSubType'] ]

sns.set()
sns.set(font_scale = 0.8)
fig = palplot(colors, cats, size = 1.0)
plt.savefig("../../figures/Cancer_cellType-heatmap-legend-top5.pdf", dpi = 600, bbox_inches = "tight")

sns.set(font_scale=1.2)
g = sns.clustermap(auc_mtx_Z[topreg], annot=False,  square=False,  linecolor='gray',
    yticklabels=False, xticklabels=True, vmin=-2, vmax=6, row_colors=colormap,
    cmap="YlGnBu", figsize=(21,16) )
g.cax.set_visible(True)
g.ax_heatmap.set_ylabel('')
g.ax_heatmap.set_xlabel('')
plt.savefig("../../figures/Cancer_cellType-heatmap-top5.pdf", dpi = 600, bbox_inches = "tight")

## Generate a binary regulon activity matrix:
binary_mtx, auc_thresholds = binarize( auc_mtx, num_workers=25 )
binary_mtx.head()

import pandas as pd
pd.DataFrame(binary_mtx).to_csv('../rss_sc_bin.txt', index=True, mode='w')

## Show the AUC distributions for selected regulons
# select regulons:
r = [ 'NFIB(+)', 'CDX2(+)', 'HOXB13(+)' ]

fig, axs = plt.subplots(1, 3, figsize=(12, 4), dpi=150, sharey=False)
for i,ax in enumerate(axs):
    sns.distplot(auc_mtx[ r[i] ], ax=ax, norm_hist=True, bins=100)
    ax.plot( [ auc_thresholds[ r[i] ] ]*2, ax.get_ylim(), 'r:')
    ax.title.set_text( r[i] )
    ax.set_xlabel('')
    
fig.text(-0.01, 0.5, 'Frequency', ha='center', va='center', rotation='vertical', size='large')
fig.text(0.5, -0.01, 'AUC', ha='center', va='center', rotation='horizontal', size='large')

fig.tight_layout()
fig.savefig('../../figures/Cancer_cellType-binaryPlot2.pdf', dpi=600, bbox_inches='tight')

## Further exploration of modules directly from the network inference output
adjacencies = pd.read_csv("adj_Cancer.tsv", index_col = False, sep = "\t")

from pyscenic.utils import modules_from_adjacencies
modules = list(modules_from_adjacencies(adjacencies, exprMat))

# pick out modules for EBF1:
tf = 'EBF1'
tf_mods = [ x for x in modules if x.transcription_factor==tf ]

for i,mod in enumerate( tf_mods ):
    print( f'{tf} module {str(i)}: {len(mod.genes)} genes' )
print( f'{tf} regulon: {len(regulons[tf+"_(+)"])} genes' )

# write these modules, and the regulon to files:
for i,mod in enumerate( tf_mods ):
    with open( tf+'_module_'+str(i)+'.txt', 'w') as f:
        for item in mod.genes:
            f.write("%s\n" % item)
            
with open( tf+'_regulon.txt', 'w') as f:
    for item in regulons[tf+'_(+)']:
        f.write("%s\n" % item)

from IPython.display import display, Image
display(Image(filename='iRegulon_screenshot_PBMC10k-EBF1.png'))
