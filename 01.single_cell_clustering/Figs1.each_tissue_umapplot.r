
### Get the parameters
parser = argparse::ArgumentParser(description = 'Script to integrate regeneration batches')
parser$add_argument('-i', dest = 'input', help = 'input directory contains h5ad files') 
parser$add_argument('-o', dest = 'outfile', default = '.h5ad', help = 'suffix name pattern for find h5ad files to use, default ".h5ad"') 
opts = parser$parse_args()



allcolour <- c("#B84D64", "#864A68", "#EE7072", "#E32D32", "#998B95", "#5E549A", "#8952A0", "#4552A0", "#384B97", 
            "#2B3B72", "#911310", "#384C99", "#9B8E8C", "#7CA878", "#35A132", "#6B70B0", "#3D6AAA", "#394D9B", "#75ACC3", "#20ACBD", 
            "#38509F", "#959897", "#F4A2A3", "#F69896", "#B6CCD7", "#AF98B5", "#E01516", "#A09C9A", "#F6EDEF")

#library(Seurat)
library(tidyverse)
library(ggrepel)
library(scales)
#show_col(color_cluster) 
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggthemes)
library(ggrastr)
#library(tidydr)
#library(ggunchull)
library(stringr)
#library(ggsci)

umap <- read.csv(opts$input,sep = ',',head=T )

pdf(opts$outfile,width = 10,height = 10)
p <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +
    geom_point(size = 1 , alpha =1 )  + scale_color_manual(values = allcolour)

cell_type_med <- umap %>%
  group_by(cell_type) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

p6 = p + geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med,
                   point.padding=unit(0.5, "lines")) +
  theme_tufte() + 
  theme(legend.position = "none", axis.line = element_line(colour = "black", size = 1, linetype = "solid"))

rasterise(p6, layers='Point', dpi=600)
dev.off()
#ggsave('blood.png',p6,width = 10,height=10)
