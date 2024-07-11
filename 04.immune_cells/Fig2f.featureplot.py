
import anndata
import scanpy as sc
sc.set_figure_params(vector_friendly=True, dpi_save=600)

adata = anndata.read('immune_atlas.h5ad')

sc.pl.umap(adata, color=['MSTRG.25883', 'MSTRG.822', 'MSTRG.24553'], save='g.subtype_markers.pdf')
for t in adata.obs['tissue'].unique():
    ad = adata[adata.obs['tissue'] == t]
    sc.pl.umap(ad, color=['MSTRG.25883', 'MSTRG.822'], save=f'g.subtype_markers.{t}.pdf')

