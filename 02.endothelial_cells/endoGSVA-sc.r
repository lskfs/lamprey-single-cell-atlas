#sc_GSVA
library(ggplot2)
library(dplyr)
library(msigdbr)
library(Seurat)
library(GSVA)
library(pheatmap)
library(patchwork)

setwd('f:/0-work/1.项目/4.全组织图谱/9.作图/13.endothelial/血管相关基因集/')
list <- read.csv('endothelial.txt')
CD8_signature=read_lines("f:/0-work/1.项目/4.全组织图谱/8.基因集评分/human/msigdb.v2023.2.Hs.symbols.gmt") %>%
    lapply(str_split, "\\t") %>% 
    unlist(recursive = F) %>% 
    lapply(function(x) setNames(list(x[-c(1:2)]), x[1])) %>% 
    unlist(recursive = F)
filter_conditions <- list$V1
filtered_names <- names(CD8_signature)[names(CD8_signature) %in% filter_conditions]
filtered_CD8_signature <- CD8_signature[filtered_names]
save(filtered_CD8_signature,file='endo_signature.rdata')
pbmc <- readRDS('../endothelial.rds')

# count数据框
count <- data.frame(lamprey = c("nbisL1-mrna-82", "nbisL1-mrna-14669", "nbisL1-mrna-5358",
                                "MSTRG.1279", "MSTRG.20876", "MSTRG.3730"),
                    gene_name = c("PLIN3", "PLIN4", "PLIN2", "PLIN4", "PLIN4", "PLIN2"))#两列对应关系
count=read.csv()#读取成两列的矩阵
# CD8_signature列表
CD8_signature <- list(chr10p13 = c("ACBD7", "BEND7", "BEND7-DT", "BTBD7P1", "C1QL3", "CAMK1D"))#示例，以前边读取的CD8_signature为主

# 映射关系字典
mapping <- setNames(count$lamprey, count$gene_name)

# 使用lapply函数遍历CD8_signature列表，并根据映射关系进行替换
new_CD8_signature <- lapply(CD8_signature, function(gene_list) {
  # 根据映射关系将基因名替换为lamprey值
  mapped_genes <- mapping[gene_list]
  return(mapped_genes)
})

# 指定输出文件路径和文件名
output_file <- "endothelial_output.gmt"
# 打开输出文件
file_conn <- file(output_file, "w")

# 遍历list中的每个集合
for (name in names(new_CD8_signature)) {
  # 获取当前集合的元素列表
  elements <- new_CD8_signature[[name]]
  
  # 将集合名称和元素列表以GMT格式写入文件
  writeLines(paste(name, paste(elements, collapse = "\t"), sep = "\t"), file_conn)
}
# 关闭文件连接
close(file_conn)

#单细胞GSVA
setwd("/home/wucheng/jianshu/function/data")
pbmc <-readRDS("pbmc.rds")
expr <- as.data.frame(pbmc@assays$RNA@data) #表达矩阵
meta <- pbmc@meta.data[,c("orig.ident","tissue")] #类别
m_df = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") #选取物种人类
msigdbr_list = split(x = m_df$gene_symbol, f = m_df$gs_name)
expr=as.matrix(expr) 
kegg <- gsva(expr, msigdbr_list, kcdf="Gaussian",method = "gsva",parallel.sz=10) #gsva
pheatmap(kegg, show_rownames=1, show_colnames=0, annotation_col=meta,fontsize_row=5, filename='gsva_heatmap.png', width=15, height=12)#绘制热图
es <- data.frame(t(kegg),stringsAsFactors=F)  #添加到单细胞矩阵中，可视化相关通路的在umap上聚集情况，可理解为一个通路即一个基因
scRNA <- AddMetaData(pbmc, es)
p1 <- FeaturePlot(scRNA, features = "KEGG_PRIMARY_BILE_ACID_BIOSYNTHESIS", reduction = 'umap')
p2 <- FeaturePlot(scRNA, features = "KEGG_ETHER_LIPID_METABOLISM", reduction = 'umap')
p3 <- FeaturePlot(scRNA, features = "KEGG_RIBOSOME", reduction = 'umap')
p4 <- FeaturePlot(scRNA, features = "KEGG_ASTHMA", reduction = 'umap')
plotc = (p1|p2)/(p3|p4)
ggsave('gsva_featureplot.png', plotc, width = 10, height = 8) #输出图片
##每个细胞类别与功能相关热图
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


# 指定输出文件路径和文件名
output_file <- "endothelial_output.gmt"

# 打开输出文件
file_conn <- file(output_file, "w")

# 遍历list中的每个集合
for (name in names(filtered_CD8_signature)) {
  # 获取当前集合的元素列表
  elements <- filtered_CD8_signature[[name]]
  
  # 将集合名称和元素列表以GMT格式写入文件
  writeLines(paste(name, paste(elements, collapse = "\t"), sep = "\t"), file_conn)
}

# 关闭文件连接
close(file_conn)



# 遍历每一列
top5 <- apply(data, 2, function(x) {
  # 按值的大小排序
  x_sorted <- sort(x, decreasing = TRUE)
  # 取前5个最大值
  top_5 <- head(x_sorted, 5)
  # 返回前5个最大值对应的列名
  names(top_5)
})
# 将列名保存在一个向量中
name <- unique(unlist(top5))

p1 <- pheatmap(df,
              cluster_rows = F,
              cluster_cols = F,
              show_rownames = T,
              show_colnames = T,
              color =colorRampPalette(c("#C5E99B", "#E0E3DA","#A593E0"))(100),
              cellwidth = 10, cellheight = 15,
              fontsize = 10)

