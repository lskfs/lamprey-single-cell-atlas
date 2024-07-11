### Get the parameters
parser = argparse::ArgumentParser(description = 'Script to integrate regeneration batches')
parser$add_argument('-i', dest = 'input', help = 'input directory contains h5ad files') 
parser$add_argument('-S', dest = 'suffix', default = '.h5ad', help = 'suffix name pattern for find h5ad files to use, default ".h5ad"') 
parser$add_argument('-s', dest = 'sample', help = 'sample id')
parser$add_argument('-d', dest = 'dims', default = 30, help = 'dims for umap, default 30')
parser$add_argument('-r', dest = 'resolution', default = 0.5, help = 'cluster resolution, default 0.5')
parser$add_argument('-n', dest = 'normalize', default = 'sct', choices = c('sct', 'log'), help = 'normalization method for seurat, sct for sctransform, log for LogNormalize, default sct')
parser$add_argument('-p', dest = 'pointSize', default = 1, help = 'point size in the figure')
opts = parser$parse_args()

opts$pointSize <- as.numeric(opts$pointSize)

library(data.table)
library(dplyr)
library(Seurat)
library(anndata)
library(ggplot2)
library(cowplot)
library(patchwork)
library(RColorBrewer)
library(stringr)


obj.integrated <- readRDS('pancreas.evo.atlas.integrated.rds')
obj.integrated <- subset(obj.integrated, subset = species == 'zebrafish', invert = T)
obj.integrated <- subset(obj.integrated, subset = ((species == 'mouse') && (tissue == 'Embryo')), invert = T)

##############################################################
# now let do cluster
##############################################################
DefaultAssay(object = obj.integrated) <- "integrated"
if (opts$normalize == 'log'){
    obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
}
obj.integrated <- RunPCA(object = obj.integrated, verbose = FALSE)
obj.integrated <- FindNeighbors(object = obj.integrated, dims = 1:as.numeric(opts$dims), reduction = "pca")
obj.integrated <- FindClusters(object = obj.integrated, resolution = as.numeric(opts$resolution))
#obj.integrated <- RunUMAP(object = obj.integrated, reduction = "pca", dims = 1:as.numeric(opts$dims))

###################
# save metadata and rds
write.table(obj.integrated@meta.data, file=paste0(opts$sample, '.metadata.txt'), quote=F, sep='\t', row.names=T)
#saveRDS(obj.integrated, file=paste0(opts$sample, '.integrated.rds'))

###################
# prepare palette
cluster_number <- length(unique(obj.integrated@meta.data$seurat_clusters))
if (cluster_number <= 25){
    cluster_Palette <- c('dodgerblue2', '#E31A1C', 'green4', '#6A3D9A', '#FF7F00', 'black', 'gold1', 
                         'skyblue2', '#FB9A99', 'palegreen2', '#CAB2D6', '#FDBF6F', 'gray70', 'khaki2', 
                         'maroon', 'orchid1', 'deeppink1', 'blue1', 'steelblue4', 'darkturquoise', 
                         'green1', 'yellow4', 'yellow3','darkorange4', 'brown')
} else if (cluster_number > 25 && cluster_number <= 70){
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
    cluster_Palette <- unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
}

###################
# draw batch
DimPlot(obj.integrated, reduction = 'umap', group.by = 'seurat_clusters', split.by = 'species',cols = cluster_Palette, ncol = 4, raster=FALSE)
ggsave('umap_species.png')

###################
# draw umap
DimPlot(obj.integrated, reduction = 'umap', label = TRUE, cols = cluster_Palette, repel = TRUE,raster=FALSE)
ggsave('umap.png', width = 12, height = 12)

###################
# find markers
DefaultAssay(obj.integrated) <- 'RNA'
markers <- FindAllMarkers(obj.integrated, min.pct = 0.1, logfc.threshold = 0.25)
markers <- markers[with(markers, order(cluster, -avg_log2FC)), ]
write.table(markers, paste0(opts$sample, '.AllMarkers.xls'), sep = '\t', quote = FALSE, row.names = F)

topn <- markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.table(topn, paste0(opts$sample, '.AllMarkers.top30.xls'), sep='\t', quote = FALSE, row.names = F)



