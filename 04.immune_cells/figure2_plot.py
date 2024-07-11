#!/usr/bin/env python3

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import anndata
import scanpy as sc
sc.set_figure_params(vector_friendly=True, dpi_save=600)

adata = anndata.read('./immune_atlas.h5ad')

"""
adata.obs['leiden'] = adata.obs['leiden'].astype(str)
adata.obs.loc[~adata.obs['leiden'].isin(['4', '5', '7', '8']), 'leiden'] = 'other'
adata.obs['leiden'] = adata.obs['leiden'].astype('category')
colors = ['#BB7784', '#8E063B', '#8595E1', '#B5BBE3', '#dddfe6']
sc.pl.umap(adata, color=['leiden'], 
        save='.leiden.highlight.pdf', palette=colors,
        )
sys.exit()
"""

celltype2color = pd.read_csv('../../../celltype2color.multi.txt', sep='\t', header=0)
celltype2color = celltype2color.rename(columns={'celltype': 'celltype.1124'})
celltype2color = celltype2color[['celltype.1124', 'choice1']]

celltype2color_dict = dict(celltype2color.to_dict('tight')['data'])
colors = [celltype2color_dict[x] for x in adata.obs['celltype.1124'].cat.categories.to_list()]
sc.pl.umap(adata, color=['celltype.1124'], 
        save='.celltype.pdf', palette=colors, 
        )
#axes[0].set_aspect('equal')
fig, ax = plt.subplots(dpi=700)
obs = adata.obs.copy()
obs = obs[['celltype.1124']].value_counts(ascending=True).to_frame(name='counts').reset_index()
obs['counts'] = np.log10(obs['counts'])
obs['celltype.1124'] = obs['celltype.1124'].astype(str)
obs = obs.sort_values(by='counts', ascending=True)
obs = obs.merge(celltype2color, how='left', on='celltype.1124')
colors = obs['choice1'].values
obs.plot.barh(x='celltype.1124', y='counts', color=colors, ax=ax, )
fig.savefig('immune.bar.pdf')

tissue2color = pd.read_csv('../../../tissue2color.multi.txt', sep='\t', header=0)
tissue2color = tissue2color[['tissue', 'choice6']]
tissue2color = dict(tissue2color.to_dict('tight')['data'])
colors = [tissue2color[x] for x in adata.obs['tissue'].cat.categories.to_list()]
sc.pl.umap(adata, color=['tissue'], 
        save='.tissue.pdf', palette=colors, 
        )



sys.exit()
sc.pl.umap(adata, color=['leiden'], 
        save='.leiden.pdf', #palette=colors, 
        )

#colors = [celltype2color[x] for x in adata.obs['celltype.1124'].cat.categories.to_list()]
sc.pl.umap(adata, color=['celltype.1124'], 
        save='.celltype1124.pdf', #palette=colors, 
        )


#adata.write('immune_atlas.1124.h5ad')

