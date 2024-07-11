#sc_GSVA
library(ggplot2)
library(dplyr)
library(msigdbr)
library(Seurat)
library(GSVA)
library(pheatmap)
library(patchwork)

setwd('')
list <- read.csv('endothelial.txt')
CD8_signature=read_lines("msigdb.v2023.2.Hs.symbols.gmt") %>%
    lapply(str_split, "\\t") %>% 
    unlist(recursive = F) %>% 
    lapply(function(x) setNames(list(x[-c(1:2)]), x[1])) %>% 
    unlist(recursive = F)
filter_conditions <- list$V1
filtered_names <- names(CD8_signature)[names(CD8_signature) %in% filter_conditions]
filtered_CD8_signature <- CD8_signature[filtered_names]
save(filtered_CD8_signature,file='endo_signature.rdata')
pbmc <- readRDS('../endothelial.rds')

'''
count <- data.frame(lamprey = c("nbisL1-mrna-82", "nbisL1-mrna-14669", "nbisL1-mrna-5358",
                                "MSTRG.1279", "MSTRG.20876", "MSTRG.3730"),
                    gene_name = c("PLIN3", "PLIN4", "PLIN2", "PLIN4", "PLIN4", "PLIN2"))
count=read.csv()

CD8_signature <- list(chr10p13 = c("ACBD7", "BEND7", "BEND7-DT", "BTBD7P1", "C1QL3", "CAMK1D"))

'''
mapping <- setNames(count$lamprey, count$gene_name)


new_CD8_signature <- lapply(CD8_signature, function(gene_list) {

  mapped_genes <- mapping[gene_list]
  return(mapped_genes)
})


output_file <- "endothelial_output.gmt"

file_conn <- file(output_file, "w")


for (name in names(new_CD8_signature)) {

  elements <- new_CD8_signature[[name]]
  

  writeLines(paste(name, paste(elements, collapse = "\t"), sep = "\t"), file_conn)
}

close(file_conn)


setwd("/home/wucheng/jianshu/function/data")
pbmc <-readRDS("pbmc.rds")
expr <- as.data.frame(pbmc@assays$RNA@data) 
meta <- pbmc@meta.data[,c("orig.ident","tissue")] 
m_df = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") 
msigdbr_list = split(x = m_df$gene_symbol, f = m_df$gs_name)
expr=as.matrix(expr) 
kegg <- gsva(expr, msigdbr_list, kcdf="Gaussian",method = "gsva",parallel.sz=10) 
pheatmap(kegg, show_rownames=1, show_colnames=0, annotation_col=meta,fontsize_row=5, filename='gsva_heatmap.png', width=15, height=12)
es <- data.frame(t(kegg),stringsAsFactors=F)  
scRNA <- AddMetaData(pbmc, es)
p1 <- FeaturePlot(scRNA, features = "KEGG_PRIMARY_BILE_ACID_BIOSYNTHESIS", reduction = 'umap')
p2 <- FeaturePlot(scRNA, features = "KEGG_ETHER_LIPID_METABOLISM", reduction = 'umap')
p3 <- FeaturePlot(scRNA, features = "KEGG_RIBOSOME", reduction = 'umap')
p4 <- FeaturePlot(scRNA, features = "KEGG_ASTHMA", reduction = 'umap')
plotc = (p1|p2)/(p3|p4)
ggsave('gsva_featureplot.png', plotc, width = 10, height = 8) 

meta <- meta %>%arrange(meta$seurat_clusters)
data <- kegg[,rownames(meta)]
group <- factor(meta[,"seurat_clusters"],ordered = F)
data1 <-NULL
for(i in 0:(length(unique(group))-1)){
ind <-which(group==i)
dat <- apply(data[,ind], 1, mean)
data1 <-cbind(data1,dat)
}
colnames(data1) <-c("C0","C1","C2","C3","C4","C5","C6","C7","C8")
result<- t(scale(t(data1)))
result[result>2]=2
result[result<-2]=-2
library(pheatmap)
p <- pheatmap(result[1:20,],
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = T,
                color =colorRampPalette(c("blue", "white","red"))(100),
                cellwidth = 10, cellheight = 15,
                fontsize = 10)
pdf(("gsva_celltype.pdf"),width = 7,height = 7)
print(p)
dev.off()



output_file <- "endothelial_output.gmt"


file_conn <- file(output_file, "w")


for (name in names(filtered_CD8_signature)) {

  elements <- filtered_CD8_signature[[name]]
  

  writeLines(paste(name, paste(elements, collapse = "\t"), sep = "\t"), file_conn)
}


close(file_conn)




top5 <- apply(data, 2, function(x) {

  x_sorted <- sort(x, decreasing = TRUE)

  top_5 <- head(x_sorted, 5)

  names(top_5)
})

name <- unique(unlist(top5))

p1 <- pheatmap(df,
              cluster_rows = F,
              cluster_cols = F,
              show_rownames = T,
              show_colnames = T,
              color =colorRampPalette(c("#C5E99B", "#E0E3DA","#A593E0"))(100),
              cellwidth = 10, cellheight = 15,
              fontsize = 10)

