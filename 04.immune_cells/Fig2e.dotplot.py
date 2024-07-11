 #!/usr/bin/env python3

import sys
from pathlib import Path
from collections import defaultdict
import matplotlib.pyplot as plt
import pandas as pd
import anndata
import scanpy as sc

gene_anno = pd.read_csv('../../00.data/genome/gene_annotation.csv', sep=',', header=0)
gene_anno = gene_anno.rename(columns={'gene_name': 'names', 'description': 'gene_name'})
gene_anno = gene_anno.drop_duplicates(subset='names')
gene_anno = gene_anno.drop_duplicates(subset='gene_name')

if not Path('degs.txt').exists():
    degs = pd.read_csv('./markers_list.txt', sep='\t', header=0, index_col=0)
    degs = degs.merge(gene_anno, how='left', on=['names'])
    degs = degs.dropna()
    degs = degs.groupby('group').head(15).reset_index()
    degs = degs.drop_duplicates(subset='names')
    degs = degs.drop_duplicates(subset='gene_name')
    degs.to_csv('degs.txt', sep='\t')
else:
    degs = pd.read_csv('degs.txt', sep='\t', header=0, index_col=0)

adata = anndata.read('immune_atlas.h5ad')
adata = adata.raw.to_adata()
sc.tl.dendrogram(adata, groupby='celltype.1124')
degs = degs[degs['names'].isin(adata.var_names)]
print(degs)
gene_ids = degs['names'].values
gene_names = degs['gene_name'].values
markers = defaultdict(list)
for index, row in degs.iterrows():
    gene = row.names
    if gene not in adata.var_names:
        continue
    markers[row.group].append(row.names)
print(markers)

adata = adata[:, gene_ids]

gene_anno = gene_anno.rename(columns={'names': 'Gene'})
gene_anno = gene_anno.set_index('Gene')
adata.var = adata.var.merge(gene_anno, how='left', left_index=True, right_index=True)

sc.pl.dotplot(adata, markers, 
        groupby='celltype.1124', 
        #gene_symbols='gene_name', 
        dendrogram=True, 
        swap_axes=True, 
        save='.markers.pdf', 
        cmap='YlGnBu'
        )

