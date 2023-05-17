import scanpy as sc
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd


results_file='./scrna.h5ad'
data=sc.read_h5ad(results_file)
data
sc.tl.draw_graph(data)
sc.pl.draw_graph(data, color='leiden', legend_loc='on data',title = "")
plt.savefig("draw_graph.pdf")

sc.tl.louvain(data) #可以使用resolution调节聚类的簇的数据，如resolution=1.0
sc.tl.paga(data, groups='louvain')
sc.pl.paga(data, color=['louvain', 'MS4A1', 'NKG7', 'PPBP']) ##随便挑选了几个基因
plt.savefig("paga_celltype.pdf")

data.obs['louvain'].cat.categories
data.obs['louvain_anno'] = data.obs['louvain']
data.obs['louvain_anno'] = data.obs['louvain_anno'].cat.rename_categories(['Naive CD4 T', 'Memory CD4', 'CD14 Monocytes', 'B', 'CD8 T', 'FCGR3A Monocytes', 'NK', 'DC', 'Platelet'])
sc.tl.paga(data, groups='louvain_anno')
sc.pl.paga(data, threshold=0.03, show=False) #细胞与细胞之间的关系，距离越近表示关系越接近
plt.savefig("paga_celltype1.pdf")

#利用PAGA重新计算细胞之间的距离
sc.tl.draw_graph(data, init_pos='paga')
sc.pl.draw_graph(data, color=['louvain_anno','MS4A1', 'NKG7', 'PPBP'], legend_loc='on data')
plt.savefig("paga_celltype_graph.pdf")

zeileis_colors = np.array(sc.pl.palettes.zeileis_28)
new_colors = np.array(data.uns['louvain_anno_colors'])
new_colors[[0]] = zeileis_colors[[12]]  # CD4 T colors / green
new_colors[[1]] = zeileis_colors[[5]]  # CD14 Monocytes colors / red
new_colors[[2]] = zeileis_colors[[17]]  # B colors / yellow
new_colors[[3]] = zeileis_colors[[2]]  # CD8 T / grey
new_colors[[4]] = zeileis_colors[[18]]  # NK / turquoise
new_colors[[5]] = zeileis_colors[[6]]  # FCGR3A Monocytes / light blue
new_colors[[6]] = zeileis_colors[[0]]  # Dendritic / dark blue
new_colors[[7]] = zeileis_colors[[25]]  # Megakaryocytes / grey
#new_colors[[10, 17, 5, 3, 15, 6, 18, 13, 7, 12]] = zeileis_colors[[5, 5, 5, 5, 11, 11, 10, 9, 21, 21]]  # CD14 Monocytes colors / red

data.uns['louvain_anno_colors'] = new_colors
sc.pl.paga_compare(
    data, threshold=0.03, title='', right_margin=0.2, size=10, edge_width_scale=0.5,
    legend_fontsize=12, fontsize=12, frameon=False, edges=True)
plt.savefig("paga_compare.pdf")

data.uns['iroot'] = np.flatnonzero(data.obs['louvain_anno']  == 'B')[0] ##假设分化起点为B cells,当然自己分析的时候需要根据数据实际情况选择分化起点
sc.tl.dpt(data)
sc.pl.draw_graph(data, color=['louvain_anno', 'dpt_pseudotime'],  legend_loc='on data',title = ['','pseudotime'], frameon=True)
plt.savefig("paga_peudotime.pdf")
sc.pl.draw_graph(data, color=['louvain', 'dpt_pseudotime'],  legend_loc='on data',title = ['','pseudotime'], frameon=True)
plt.savefig("paga_peudotime1.pdf")