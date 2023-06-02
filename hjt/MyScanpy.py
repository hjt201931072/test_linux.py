
import scanpy as sc
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import bbknn



# mpl.rcParams['pdf.fonttype'] = 42
# mpl.rcParams["font.sans-serif"] = "Arial"
# sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')
# results_file='./scrna1.h5ad'
# results_file1='./scrna2.h5ad'
#
#
# adata = sc.datasets.pbmc3k()
adata = sc.read_10x_mtx('hjt\\data1\\filtered_gene_bc_matrices/', var_names='gene_symbols', cache=True)
adata.var_names_make_unique()
adata
obs = pd.read_csv('hjt\\data1\\filtered_gene_bc_matrices\\barcodes.tsv', header=None, index_col=0, sep='\t')
var = pd.read_csv('hjt\\data1\\filtered_gene_bc_matrices\\genes.tsv', header=None, index_col=1, sep='\t')
obs.index.name = ''
var.index.name = ''
var.columns = ['gene_ids']


from scipy.io import mmread
from scipy.sparse import csr_matrix
mtx = mmread('hjt\\data1\\filtered_gene_bc_matrices\\matrix.mtx')
mtx = mtx.T
mtx = csr_matrix(mtx)
adata1 = sc.AnnData(mtx, obs=obs, var=var)

adata1.var_names_make_unique()
sc.pl.highest_expr_genes(adata, n_top=20, )
#
#


#qc
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)#0

sc.settings.set_figure_params(dpi_save=300)
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(8, 3))
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=ax[0], show=False)
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax=ax[1], show=False)
plt.subplots_adjust(wspace=.4)#1
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
#
#normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata


# feature selection
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=2000)
sc.pl.highly_variable_genes(adata)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)


#reduction
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)#3
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)

adata
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

# clustering
sc.tl.leiden(adata)
sc.pl.umap(adata, color='leiden')

sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
#adata.write(results_file)
# adata=sc.read_h5ad(results_file)

markers = ["MS4A1", "TYROBP", "CD14",'FCGR3A', "FCER1A", "CCR7", "IL7R", "PPBP", "CD8A"]
sc.pl.umap(adata, color=markers, ncols=3)
# adata.write(results_file)
new_cluster_names = ['Naive CD4 T', 'Memory CD4', 'CD14 Monocytes','B', 'CD8 T', 'FCGR3A Monocytes','NK', 'DC', 'Platelet']
adata.rename_categories('leiden', new_cluster_names)
sc.settings.set_figure_params(dpi=50, dpi_save=300, figsize=(7, 7))
sc.pl.umap(adata, color='leiden', legend_loc='on data')
# adata.write(results_file1)
# adata=sc.read_h5ad(results_file)
# sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
