
import anndata

adata = anndata.read('../../00.data/scRNA/merged.h5ad')


human_tissue = ['AdultLiver', 'AdultPancreas', 'FetalLiver', 'FetalIntestine', 'FetalPancreas']
mouse_tissue = ['Intestine', 'Liver', 'Pancreas']
exclude_batch = ['FetalBrain_1', 'FetalFemaleGonad_1', 'FetalHeart_1', 'FetalHeart_2', 'FetalKidney_1', 'FetalKidney_2', 
        'FetalLung_1', 'FetalMaleGonad_1', 'FetalMaleGonad_2', 'FetalStomach_1', 'FetalStomach_2']
lamprey_tissue = ['intestine', 'liver']
adata = adata[
        ((adata.obs['species'] == 'lamprey') & (adata.obs['tissue'].isin(lamprey_tissue))) |
        ((adata.obs['species'] == 'human') & (adata.obs['tissue'].isin(human_tissue))) | 
        ((adata.obs['species'] == 'mouse') & (adata.obs['tissue'].isin(mouse_tissue))) | 
        ((adata.obs['species'] == 'mouse') & (adata.obs['tissue'] == 'Embryo') & (~adata.obs['batch'].isin(exclude_batch))) | 
        ((adata.obs['species'] == 'zebrafish') & (adata.obs['tissue'].isin(zebrafish_tissue)))
        ]
adata.obs.to_csv('pancreas_group.obs.txt', sep='\t', header=True, index=True)
adata.write('pancreas_group.h5ad')


