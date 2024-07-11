#!/usr/bin/env python3

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import anndata
import scanpy as sc
sc.set_figure_params(vector_friendly=True, dpi_save=600)


def legend_plot(group2color_df, color='hex', ax=None):
    if 'tissue' in group2color_df.columns:
        label = 'tissue'
    else:
        label = 'celltype'
    legend_elements = []
    print(group2color_df)
    #group2color_df = group2color_df.sort_values(by=f'{label}.index', key=lambda x: int(x))
    for index, row in group2color_df.iterrows():
        num = row[f'{label}.index']
        text = row[f'{label}']
        c = row[color]
        element = Line2D(
                [0], [0], marker='o', 
                markeredgecolor='none',
                color=c, 
                label=f'{int(num):02d} {text}'
                )
        legend_elements.append(element)

    ax.legend(handles=legend_elements, loc='center')
    return ax

adata = anndata.read('./lamprey_atlas.1124.h5ad')
obs = adata.obs.copy()

tissue2counts = adata.obs[['tissue']].value_counts().to_frame(name='counts').reset_index()
tissue2counts['tissue'] = tissue2counts['tissue'].astype(str)
tissue2counts = tissue2counts.sort_values(by='counts', ascending=False)
tissue2counts['tissue.index'] = tissue2counts.index.astype(int) + 1
tissue2counts['tissue.index'] = tissue2counts['tissue.index'].astype(str).astype('category')
tissue_obs = obs.reset_index().merge(tissue2counts, how='left', on=['tissue']).set_index('CellID')[['tissue.index']]
adata.obs = adata.obs.merge(tissue_obs, how='left', left_index=True, right_index=True)
tissue2color = pd.read_csv('../tissue2color.multi.txt', sep='\t', header=0)
tissue_meta = tissue2color.merge(tissue2counts, how='left', on='tissue')
tissue2color_dict = dict(tissue_meta[['tissue.index', 'choice6']].to_dict('tight')['data'])
colors = [tissue2color_dict[x] for x in adata.obs['tissue.index'].cat.categories.to_list()]
fig, axes = plt.subplots(1, 3, dpi=600, figsize=(16, 6))
sc.pl.umap(adata, color=['tissue.index'], legend_loc='on data',
        palette=colors, ax=axes[0], 
        )
axes[0].set_aspect('equal')
tissue_meta['counts'] = np.log10(tissue_meta['counts'])
tissue_meta.plot.barh(x='tissue', y='counts', color=tissue_meta['choice6'].values, ax=axes[1], )
axes[2] = legend_plot(tissue_meta, color='choice6', ax=axes[2])
fig.savefig('umap.tissue.index.pdf')

celltype2counts = adata.obs[['celltype.1124']].value_counts().to_frame(name='counts').reset_index()
celltype2counts['celltype.1124'] = celltype2counts['celltype.1124'].astype(str)
celltype2counts = celltype2counts.sort_values(by='counts', ascending=False)
celltype2counts['celltype.index'] = celltype2counts.index.astype(int) + 1
celltype2counts['celltype.index'] = celltype2counts['celltype.index'].astype(str).astype('category')
celltype_obs = obs.reset_index().merge(celltype2counts, how='left', on=['celltype.1124']).set_index('CellID')[['celltype.index']]
adata.obs = adata.obs.merge(celltype_obs, how='left', left_index=True, right_index=True)
celltype2counts = celltype2counts.rename(columns={'celltype.1124': 'celltype'})
celltype2color = pd.read_csv('../celltype2color.multi.txt', sep='\t', header=0)
celltype_meta = celltype2color.merge(celltype2counts, how='left', on='celltype')
celltype2color_dict = dict(celltype_meta[['celltype.index', 'choice1']].to_dict('tight')['data'])
colors = [celltype2color_dict[x] for x in adata.obs['celltype.index'].cat.categories.to_list()]
fig, axes = plt.subplots(1, 2, dpi=600, figsize=(16, 6))
sc.pl.umap(adata, color=['celltype.index'], legend_loc='on data', 
        palette=colors, ax=axes[0], 
        )
axes[0].set_aspect('equal')
axes[1] = legend_plot(celltype_meta, color='choice1', ax=axes[1])
fig.savefig('umap.celltype.index.pdf')


sys.exit()
palette = ['hex', 'choice1', 'choice2', 'choice3', 'choice4', 'choice5', 'choice6']
fig, axes = plt.subplots(1, 7, dpi=700, figsize=(28, 4))
tissue2color = pd.read_csv('../tissue2color.multi.txt', sep='\t', header=0)
tissue2index = adata.obs[['tissue', 'tissue.index']].drop_duplicates()
tissue2color = tissue2color.merge(tissue2index, how='left', on='tissue')
for ax, which in zip(axes, palette):
    t2c = dict(tissue2color[['tissue.index', which]].to_dict('tight')['data'])
    colors = [t2c[x] for x in adata.obs['tissue.index'].cat.categories.to_list()]
    sc.pl.umap(adata, color=['tissue.index'], ax=ax, legend_loc='on data', 
            palette=colors)
    ax.set_title(which)
fig.savefig('umap.tissue_multi.pdf')

fig, axes = plt.subplots(1, 7, dpi=700, figsize=(28, 4))
celltype2color = pd.read_csv('../celltype2color.multi.txt', sep='\t', header=0)
celltype2index = adata.obs[['celltype.1124', 'celltype.index']].drop_duplicates()
celltype2index = celltype2index.rename(columns={'celltype.1124': 'celltype'})
celltype2color = celltype2color.merge(celltype2index, how='left', on='celltype')
for ax, which in zip(axes, palette):
    c2c = dict(celltype2color[['celltype.index', which]].to_dict('tight')['data'])
    colors = [c2c[x] for x in adata.obs['celltype.index'].cat.categories.to_list()]
    sc.pl.umap(adata, color=['celltype.index'], ax=ax, legend_loc='on data', 
            palette=colors)
    ax.set_title(which)
fig.savefig('umap.celltype_multi.pdf')

sys.exit()
adata = anndata.read('./lamprey_atlas.h5ad')

"""
obs = adata.obs.copy()
obs = obs[['tissue']].value_counts(ascending=True).to_frame(name='counts').reset_index()
obs['tissue'] = obs['tissue'].astype(str)
obs = obs.sort_values(by='counts', ascending=True)

tissue2color = pd.read_csv('../tissue2color.txt', sep='\t', header=0)
obs = obs.merge(tissue2color, how='left', on='tissue')
colors = obs['color'].values

fig, ax = plt.subplots(dpi=600)
obs.plot.barh(x='tissue', y='counts', color=colors, ax=ax, )
fig.savefig('barplot.pdf')

tissue2color = dict(tissue2color.to_dict('tight')['data'])
colors = [tissue2color[x] for x in adata.obs['tissue'].cat.categories.to_list()]
sc.pl.umap(adata, color=['tissue'], 
        save='.tissue.pdf', palette=colors, 
        )
"""

obs = pd.read_csv('../lamprey.obs.txt', sep='\t', header=0, index_col=0)
obs = obs[['celltype.1124']]
adata.obs = adata.obs.merge(obs, how='left', left_index=True, right_index=True)
adata.obs['celltype.1124'] = adata.obs['celltype.1124'].astype('category')
adata.write('lamprey_atlas.1124.h5ad')
celltype2color = pd.read_csv('../celltype2color.1124.txt', sep='\t', header=0)
celltype2color = dict(celltype2color.to_dict('tight')['data'])
colors = [celltype2color[x] for x in adata.obs['celltype.1124'].cat.categories.to_list()]
sc.pl.umap(adata, color=['celltype.1124'], 
        save='.celltype1124.replot.pdf', palette=colors, 
        )



