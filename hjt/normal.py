import scanpy as sc

# Preprocessing
def preprocess(dataset):

    adata = sc.read_10x_mtx(dataset,var_names='gene_symbols', cache=True)
    # adata=sc.read(dataset)
    adata
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=1000)
    adata = adata[:, adata.var['highly_variable']]
    sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    adata.write('preprocessed.h5ad')

# Clustering
def cluster():
    adata = sc.read('preprocessed.h5ad')
    sc.tl.louvain(adata)
    sc.pl.umap(adata, color=['louvain'])
    adata.write('clustered.h5ad')

# Visualization
def visualize():
    adata = sc.read('clustered.h5ad')
    sc.pl.umap(adata, color=['louvain', 'n_genes'])
preprocess('./data1/filtered_gene_bc_matrices/')
# cluster()
# visualize()