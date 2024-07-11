# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 08:06:34 2022

@author: Ray
"""

import argparse as arg
from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                                sankey_plot, chord_plot, CellTypeTriangles, 
                                ParalogSubstitutions, FunctionalEnrichment,
                                convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import pickle

parser = arg.ArgumentParser(formatter_class=arg.RawDescriptionHelpFormatter, description = '''This programe is used to perform SAMap.''')
parser.add_argument('-species1', action='store', help='Input expression list.')
parser.add_argument('-species2', action='store', help='Output expression matrix.')
parser.add_argument('-fn1', action='store', help='Species h5ad')
parser.add_argument('-fn2', action='store', help='Species h5ad')
parser.add_argument('-celltype1', action='store', help='Celltype1.')
parser.add_argument('-celltype2', action='store', help='Celltype2.')
parser.add_argument('-NUMITERS', type=int, default=3, action='store', help='NUMITERS.')
parser.add_argument('-cpu', type=int, default=10, help='Number of cpu used.')
parser.add_argument('-mappingtable', action='store', help='Mappingtable path.')
parser.add_argument('-enrichedGene', action='store', help='Enriched genes table')
args = parser.parse_args()

filenames = {args.species1:args.fn1, args.species2:args.fn2}
keys={args.species1:args.celltype1, args.species2:args.celltype2}
sm = SAMAP(
    filenames,
    f_maps = '/dellfsqd1/ST_OCEAN/C_OCEAN/USERS/c-zhangjin/07.project/2.atlas_lampery/3.analysis/1.SAMap/3.mouse/maps/',
    save_processed=False, #if False, do not save the processed results to `*_pr.h5ad`
    keys=keys
   )
#neigh_from_keys={args.species1:True, args.species2:True}
#sm.run(NUMITERS=args.NUMITERS, neigh_from_keys=neigh_from_keys, ncpus=args.cpu)
sm.run(NUMITERS=args.NUMITERS, ncpus=args.cpu)

with open(args.species1+args.species2+'.sm', 'wb') as f:
    pickle.dump(sm, f)

D,MappingTable=get_mapping_scores(sm, keys, n_top=0)
MappingTable.to_csv(args.mappingtable, sep='\t')

plt.figure(figsize=(100,150))
sc.pl.umap(sm.smap.samap.adata, color='species')
plt.savefig('umap_'+ args.species1+ '_'+ args.species2+ '_species.png', bbox_inches='tight', dpi=600)

# plot one specie umap
plt.figure(figsize=(100,150))
sc.pl.umap(sm.sams[args.species1].adata, color=args.celltype1)
plt.savefig('umap_'+ args.species1+ '.png', bbox_inches='tight', dpi=600)

plt.figure(figsize=(100,150))
sc.pl.umap(sm.sams[args.species2].adata, color=args.celltype2)
plt.savefig('umap_'+ args.species2+ '.png', bbox_inches='tight', dpi=600)

# plot one specie umap with samap matrix
sm.sams[args.species1].adata.obsm['X_umap']=sm.sams[args.species1].adata.obsm['X_umap_samap']
sm.sams[args.species2].adata.obsm['X_umap']=sm.sams[args.species2].adata.obsm['X_umap_samap']
plt.figure(figsize=(100,150))
sc.pl.umap(sm.sams[args.species1].adata, color=args.celltype1)
plt.savefig('umap_'+ args.species1+ '_samap.png', bbox_inches='tight', dpi=600)

plt.figure(figsize=(100,150))
sc.pl.umap(sm.sams[args.species2].adata, color=args.celltype2)
plt.savefig('umap_'+ args.species2+ '_samap.png', bbox_inches='tight', dpi=600)

# gene-pair finder
gpf=GenePairFinder(sm,keys=keys)
gene_pairs=gpf.find_all(align_thr=0.1, thr=0.05)
gene_pairs.to_csv(args.enrichedGene, sep='\t', header=True, index=None)

## plot umap for all
#celltype=[args.species1+ '_'+ i for i in sm.smap.samap.adata.obs[args.species1+ '_'+ args.celltype1].to_list() if i !='unassigned'] + [args.species2+ '_'+ i for i in sm.smap.samap.adata.obs[args.species2+ '_'+ args.celltype2].to_list() if i !='unassigned']
#sm.smap.samap.adata.obs['CellType']=celltype
#plt.figure(figsize=(100,150))
#sc.pl.umap(sm.smap.samap.adata, color='CellType')
#plt.savefig('umap_'+ args.species1+ '_'+ args.species2+ '_celltype.png', bbox_inches='tight', dpi=600)


