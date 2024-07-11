
import pandas as pd
import anndata
import scanpy as sc

adata = anndata.read('endothelial_atlas.h5ad')
adata.uns['log1p']['base'] = None

#adata = adata[adata.obs['leiden'].isin(['4', '5', '7', '8'])]
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test', pts=True, )
dedf = sc.get.rank_genes_groups_df(adata, None,)
dedf.to_csv('markers_list.leiden.txt', sep='\t')

"""
sc.tl.rank_genes_groups(adata, 'celltype.1124', method='t-test', pts=True, )
dedf = sc.get.rank_genes_groups_df(adata, None,)
dedf.to_csv('markers_list.celltype1124.txt', sep='\t')
"""

