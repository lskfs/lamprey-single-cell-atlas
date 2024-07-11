
library(Seurat)
library(ggplot2)
library(clustree)
library(tidyverse)
library(cowplot)
library(dplyr)
library(patchwork) 
library(ggrepel)
library(stringr)

features <- read_lines('/dellfsqd2/C_OCEAN/USERS/c-zhangjin/07.project/2.atlas_lampery/3.analysis/3.genelist_score/2.mouse/pancrease.symbol.gmt') %>%
    lapply(str_split, "\\t") %>% 
    unlist(recursive = F) %>% 
    lapply(function(x) setNames(list(x[-c(1:2)]), x[1])) %>% 
    unlist(recursive = F)
score.colnames <- as.vector(names(features))

obj <- readRDS('pancreas.evo.atlas.integrated.rds')
DefaultAssay(obj) <- 'RNA'

name <- 'pancreas_func'
obj.cc <- AddModuleScore(
                      object = obj,
                      features = features,
                      ctrl = 100,
                      nbin = 24,
                      name = name
                     )

cc.columns <- grep(pattern = name, x = colnames(x = obj.cc[[]]), value = TRUE)
cc.scores <- obj.cc[[cc.columns]]
rm(obj.cc)
CheckGC()

obj <- AddMetaData(obj, cc.scores, col.name = score.colnames)

write.table(obj@meta.data, file = 'pancreas_func.metadata.txt', sep = '\t', quote = F)
saveRDS(obj, file = 'pancreas.evo.atlas.integrated.scored.rds')


