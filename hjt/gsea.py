import scanpy as sc
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

adata = sc.datasets.pbmc3k()
gsea='./gsea.h5ad'
# ==========================
# 质量控制
# ==========================
# 基本过滤
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# 计算质量控制指标
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# 根据基因数量和线粒体百分比进行过滤
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

# ==========================
# 数据标准化
# ==========================
#总计数归一化
sc.pp.normalize_total(adata, target_sum=1e4)
# 对数化
sc.pp.log1p(adata)


# ==========================
# 特征选取
# ==========================
# 识别高度可变的基因
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# 保存原始数据
adata.raw = adata
# 过滤
adata = adata[:, adata.var.highly_variable]

# 将数据缩放到单位方差
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)


# ==========================
# 降维、聚类、可视化
# ==========================
#通过运行umap来降低数据的维数。
sc.pp.neighbors(adata, n_neighbors=10)
sc.tl.umap(adata)

sc.tl.leiden(adata)
sc.pl.umap(adata, color=['leiden'], save='_pbmc3k_leiden.png')

# 通过logreg计算差异表达
sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')

# 每个聚类的前10个高差异表达的基因
rank_genes_df = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
gene_set_list = list(rank_genes_df["0"][:100])

print("gene_set_list: ", gene_set_list)

import gseapy as gp
# simple plotting function
from gseapy.plot import barplot, dotplot

  # default: Human

# 富集查询
enr = gp.enrichr(gene_list=gene_set_list,
                 gene_sets=['KEGG_2021_Human'],
                 organism='Human',
                 outdir='./KEGG_2021_Human',
                 cutoff=0.5,
                 format='png',
                 )
barplot(enr.res2d, title='PBMC 3K', cutoff=0.5, top_term=20, ofname='./figures/pbmc3k_enrichment_barplot.png')
dotplot(enr.res2d, title='PBMC 3K', cmap='viridis_r', size=20,cutoff=0.5, top_term=20,
        ofname='./figures/pbmc3k_enrichment_dotplot.png')