
import pandas as pd
import anndata
import scanpy as sc
sc.set_figure_params(vector_friendly=True, dpi_save=600)

celltype2color = pd.read_csv('../celltype2color.multi.txt', sep='\t', header=0)
celltype2color = dict(celltype2color[['celltype', 'choice1']].to_dict('tight')['data'])

adata = anndata.read('lamprey_atlas.1124.h5ad')
sc.pl.umap(adata, color='celltype.1124', 
        groups=['Endothelial cell'], 
        palette=celltype2color, 
        size=0.5, 
        save='.endothelial_highlight.pdf',
        )

import sys
sys.exit()
adata = anndata.read('endothelial_atlas.h5ad')
sc.tl.dendrogram(adata, groupby='tissue')
sc.pl.dendrogram(adata, groupby='tissue', save='.endothelial_dendrogram.pdf')


