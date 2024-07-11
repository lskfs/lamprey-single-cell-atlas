
library(Seurat) 
library(tidyverse)
library(msigdbr) 
library(AUCell)

DefaultAssay(sce2) <- "RNA"
signature=read_lines("lamprey_genelist_pancreas.gmt") %>%
    lapply(str_split, "\\t") %>% 
    unlist(recursive = F) %>% 
    lapply(function(x) setNames(list(x[-c(1:2)]), x[1])) %>% 
    unlist(recursive = F)
for (i in names(signature)){
    print(i)
    print(signature[i])
    sce2 <- AddModuleScore(sce2,
                          features = signature[i],
                          ctrl = 100,
                          name = i)
}
head(sce2@meta.data)

list <- colnames(sce2@meta.data)[11:41]
for (i in list){
pdf(paste0(i,'.pdf'))
p = ggviolin(sce2@meta.data,
         x = "newcell", y = i,
         width = 0.8, color = "black",
         fill = "newcell",
         palette = unique(c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                     '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                     '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                     '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                     '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                     '#968175')),
         xlab = FALSE,
         add = 'mean_sd',
         bxp.errorbar = TRUE,
         bxp.errorbar.width = 0.05,
         size = 0.5,
         legend = "right") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()
}
