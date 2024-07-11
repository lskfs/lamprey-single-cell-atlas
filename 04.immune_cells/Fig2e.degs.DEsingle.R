### Get the parameters
parser = argparse::ArgumentParser(description = "find all markers for input rds")
parser$add_argument('-r', '--region', help = 'region')
parser$add_argument('-p', '--part', help = 'part')
parser$add_argument('-o', '--out', help = 'out directory')
opts = parser$parse_args()

library(future)
library(Seurat)
library(SeuratObject)
#library(SeuratDisk)
library(dplyr)
library(data.table)
library(tibble)

#selected_tps <- c('WT', '0hpa1', '12hpa2', '36hpa2', '3dpa2', '5dpa1', '7dpa2', '10dpa1', '14dpa1')
# there is no c31 at 0hpa and 12hpa
#selected_tps <- c('WT', '36hpa2', '3dpa2', '5dpa1', '7dpa2', '10dpa1', '14dpa1') 
#selected_tps <- c('WT', opts$sample)

opts$region <- 'Granulocyte'

rds.file <- paste0(opts$region, '.rds')
if (!file.exists(rds.file)){
    dataset.merge <- readRDS('integrate/lamprey.immune.integrated.rds')
    regions <- unique(dataset.merge@meta.data$celltype)
    other <- setdiff(regions, c(opts$region))
    dataset.merge <- subset(dataset.merge, subset = celltype %in% regions)
    dataset.merge@meta.data <- dataset.merge@meta.data %>% mutate(
                                                              region = case_when(
                                                                                 celltype %in% other ~ 'other',
                                                                                 celltype == opts$region ~ opts$region
                                                                                 )
                                                              )

    DefaultAssay(dataset.merge) <- 'RNA'
    Idents(dataset.merge) <- 'region'
    saveRDS(dataset.merge, rds.file)
}else{
    dataset.merge <- readRDS(rds.file)
}

sce <- as.SingleCellExperiment(dataset.merge, assay = 'RNA')

library(DEsingle)
library(BiocParallel)

# Set the parameters and register the back-end to be used
param <- MulticoreParam(workers = 24, progressbar = TRUE)
register(param)

# Detecting the DE genes in parallelization with 18 cores
group <- factor(dataset.merge@meta.data$region)
results <- DEsingle(counts = sce, group = group, parallel = TRUE, BPPARAM = param)

# Dividing the DE genes into 3 categories at threshold of FDR < 0.05
results.classified <- DEtype(results = results, threshold = 0.05)

# Extract DE genes at threshold of FDR < 0.05
results.sig <- results.classified[results.classified$pvalue.adj.FDR < 0.05, ]

write.table(results.sig, paste0(opts$region, '_degs.txt'), sep='\t', quote = F)




