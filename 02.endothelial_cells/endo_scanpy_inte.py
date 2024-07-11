#!/usr/bin/env python3

import sys
from pathlib import Path
import matplotlib.pyplot as plt

import pandas as pd
import anndata
import scanpy as sc
sc.set_figure_params(vector_friendly=True, dpi_save=600)

adata = anndata.read('../07.lamprey_atlas/lamprey.raw.h5ad')

#adata.obs = adata.obs.reset_index()
#adata.obs['CellID'] = adata.obs['CellID'].str.rsplit('_', n=1).str[0]
#adata.obs['CellID'] = adata.obs['CellID'].astype(str) + '.' + adata.obs['tissue'].astype(str)
#adata.obs = adata.obs.set_index('CellID')

obs = pd.read_csv('lamprey.obs.txt', 
        sep='\t', header=0, index_col=0)
obs = obs[['celltype.1124']]
adata.obs = adata.obs.merge(obs, how='left', left_index=True, right_index=True)
adata.obs['celltype.1124'] = adata.obs['celltype.1124'].astype('category')
adata.obs['tissue'] = adata.obs['tissue'].astype('category')

adata = adata[adata.obs['celltype.1124'].isin(["Endothelial cell", "Vascular endothelial cell"])] #, "Mesenchymal cell"])]

sc.pp.filter_cells(adata, min_genes=300)
sc.pp.filter_genes(adata, min_cells=3)

sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata

adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts'])
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack')

sc.external.pp.bbknn(adata, batch_key='batch')
#sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata)
sc.pl.umap(adata, color=['leiden'], save='.leiden.pdf')

celltype2color = pd.read_csv('celltype2color.1124.txt', sep='\t', header=0)
celltype2color = dict(celltype2color.to_dict('tight')['data'])
tissue2color = pd.read_csv('tissue2color.txt', sep='\t', header=0)
tissue2color = dict(tissue2color.to_dict('tight')['data'])

colors = [celltype2color[x] for x in adata.obs['celltype.1124'].cat.categories.to_list()]
sc.pl.umap(adata, color=['celltype.1124'], 
        save='.celltype1124.pdf', palette=colors, 
        )

colors = [tissue2color[x] for x in adata.obs['tissue'].cat.categories.to_list()]
sc.pl.umap(adata, color=['tissue'], 
        save='.tissue.pdf', palette=colors, 
        )

adata.write('endothelial_atlas.h5ad')


