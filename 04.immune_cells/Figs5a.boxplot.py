
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import anndata

genes_b = pd.read_csv('B-B.txt', sep=';', header=None, names=['lamprey', 'mouse'])
genes_b['type'] = 'BB'
genes_c = pd.read_csv('C-T.txt', sep=';', header=None, names=['lamprey', 'mouse'])
genes_c['type'] = 'CT'
genes = pd.concat([genes_b, genes_c])
genes['lamprey'] = genes['lamprey'].str.split('_', n=1).str[1]
genes['mouse'] = genes['mouse'].str.split('_', n=1).str[1]
genes['pair'] = genes['lamprey'] + ':' + genes['mouse']

ortholog = pd.read_csv('../../02.blast/lamprey_mouse/ortholog_groups.txt', sep='\t', header=0)
ortholog = ortholog[['rec', 'gene_name']].rename(columns={'rec': 'lamprey', 'gene_name': 'mouse'})
ortholog['pair'] = ortholog['lamprey'] + ':' + ortholog['mouse']

#genes = genes[genes['pair'].isin(ortholog['pair'].values)]

lamprey_markers = pd.read_csv('../01.re-cluster/scanpy_inte/markers_list.celltype1124.txt', sep='\t', header=0, index_col=0)
lamprey_markers = lamprey_markers[lamprey_markers['group'].isin(['VLRB cell', 'VLRC cell'])]
lamprey_markers = lamprey_markers[lamprey_markers['logfoldchanges'] >= 1.2]

mouse_markers = pd.read_csv('markers_list.mouse_immune.txt', sep='\t', header=0, index_col=0)
mouse_markers = mouse_markers[mouse_markers['group'].isin(['B cell', 'T cell'])]
mouse_markers = mouse_markers[mouse_markers['logfoldchanges'] >= 1.2]

def grouped_obs_mean(adata, group_key, layer=None, gene_symbols=None):
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X
    if gene_symbols is not None:
        new_idx = adata.var[idx]
    else:
        new_idx = adata.var_names

    grouped = adata.obs.groupby(group_key)
    out = pd.DataFrame(
            np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
            columns=list(grouped.groups.keys()),
            index=adata.var_names
            )
    for group, idx in grouped.indices.items():
        X = getX(adata[idx])
        out[group] = np.ravel(X.mean(axis=0, dtype=np.float64))
    return out

if not Path('lamprey.exp.txt').exists():
    lamprey = anndata.read('lamprey.h5ad')
    lamprey = lamprey.raw.to_adata()
    lamprey_genes = genes['lamprey'].values
    lamprey = lamprey[:, lamprey_genes]
    lamprey.obs['celltype.1124'] = lamprey.obs['celltype.1124'].astype(str)
    lamprey.obs.loc[~lamprey.obs['celltype.1124'].isin(['VLRB cell', 'VLRC cell']), 'celltype.1124'] = 'other'
    
    lamprey_df = grouped_obs_mean(lamprey, 'celltype.1124')
    lamprey_genes = genes[['lamprey', 'type']].rename(columns={'lamprey': 'Gene'})
    lamprey_df = lamprey_df.merge(lamprey_genes, how='left', on='Gene')
    lamprey_df.to_csv('lamprey.exp.txt', sep='\t')
else:
    lamprey_df = pd.read_csv('lamprey.exp.txt', sep='\t', header=0, index_col=0)
    lamprey_df = lamprey_df.drop_duplicates()
    lamprey_df = lamprey_df[lamprey_df['Gene'].isin(genes['lamprey'].values)]

if not Path('mouse.exp.txt').exists():
    mouse = anndata.read('mouse.h5ad')
    mouse_genes = genes['mouse'].values
    mouse = mouse[:, mouse_genes]
    mouse.obs['name'] = mouse.obs['name'].astype(str)
    mouse.obs.loc[~mouse.obs['name'].isin(['B cell', 'T cell']), 'name'] = 'other'
    
    mouse_df = grouped_obs_mean(mouse, 'name')
    mouse_df = mouse_df.reset_index().rename(columns={'index': 'Gene'})
    mouse_genes = genes[['mouse', 'type']].rename(columns={'mouse': 'Gene'})
    mouse_df = mouse_df.merge(mouse_genes, how='left', on='Gene')
    mouse_df.to_csv('mouse.exp.txt', sep='\t')
else:
    mouse_df = pd.read_csv('mouse.exp.txt', sep='\t', header=0, index_col=0)
    mouse_df = mouse_df.drop_duplicates()
    mouse_df = mouse_df[mouse_df['Gene'].isin(genes['mouse'].values)]

lamprey_df = lamprey_df[lamprey_df['Gene'].isin(lamprey_markers['names'].values)]
lamprey_df = lamprey_df.melt(id_vars=['Gene', 'type'], var_name='group', value_name='exp')

mouse_df = mouse_df[mouse_df['Gene'].isin(mouse_markers['names'].values)]
mouse_df = mouse_df.melt(id_vars=['Gene', 'type'], var_name='group', value_name='exp')

#comm_genes = genes[genes['lamprey'].isin(lamprey_df['Gene'].unique()) & genes['mouse'].isin(mouse_df['Gene'].unique())]
comm_genes = comm_genes[comm_genes['pair'].isin(ortholog['pair'].values)]

fig, axes = plt.subplots(1, 2, figsize=(10, 5), dpi=600, sharey=True)
sns.boxplot(lamprey_df, x='group', y='exp', hue='type', 
        ax=axes[0], color=".8", linecolor="#137", linewidth=.75)
sns.boxplot(mouse_df, x='group', y='exp', hue='type', 
        ax=axes[1], color=".8", linecolor="#137", linewidth=.75)
fig.savefig('test.pdf')

