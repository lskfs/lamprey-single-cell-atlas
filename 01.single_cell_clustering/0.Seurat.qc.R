library("argparse")
parser = argparse::ArgumentParser(description = 'scRNA pipeline')
parser$add_argument('-s', '--subCommand', dest = 'subCommand', required = "True", choices = c("runPro","runFilter","runDoubletFinder","runSC","runHarmony") ,help = 'subCommand: [runPro] [runFilter] [runDoubletFinder] [runSC] [runHarmony]')
parser$add_argument('-i', '--input', dest = 'input',required = "True", help = 'input filename')
#parser = parser$add_argument_group('BasicGroup', 'BasicGroup description')
parser$add_argument('--vg', dest = 'vg', default = 3000, type = 'integer', help = 'number of variable genes, default 3000')
parser$add_argument('--dims', dest = 'dims', default = 30, type = 'integer', help = 'number of PC to use, default 30')
parser$add_argument('--topn', dest = 'topn', type = 'integer', default = 100, help = 'save the top N markers')
parser$add_argument('--mtRegEx', dest = 'mtRegEx', default = "^MT", type = 'character', help = 'maximum Feature number, default ^MT')
parser$add_argument('--projectname', dest = 'projectname', type = 'character', help = 'project name of the seurat object')

#subparsers = parser$add_subparsers(help='sub-command help')

#runPro = subparsers$add_parser("runPro",help="get the vlnplot from rawdata")
#runFilter = subparsers$add_parser("runFilter",help="filter")
runProGroup = parser$add_argument_group('runPro', 'get the vlnplot from rawdata')
runFilter = parser$add_argument_group('runFilter', 'filter')
runFilter$add_argument('--minCount', dest = 'minCount', default = 0, type = 'integer', help = 'minimum UMI number')
runFilter$add_argument('--maxCount', dest = 'maxCount', type = 'integer', help = 'maximum UMI number')
runFilter$add_argument('--minFeature', dest = 'minFeature', default = 0, type = 'integer', help = 'minimum Feature number')
runFilter$add_argument('--maxFeature', dest = 'maxFeature', type = 'integer', help = 'maximum Feature number')
runFilter$add_argument('--percentMt', dest = 'percentMt', default = 5, type = 'integer', help = 'maximum Feature number')
runFilter$add_argument('--Foutname', dest = 'Foutname', default = "mysingle", help = 'output name')

#runDoubletFinder = subparsers$add_parser("runDoubletFinder",help="doubletfinder pipeline")
runDoubletFinder = parser$add_argument_group('runDoubletFinder', 'doubletfinder  pipeline')
runDoubletFinder$add_argument('--doubletRate', dest = 'doubletRate', default = 10, type = 'integer', help = 'maximum Feature number')
runDoubletFinder$add_argument('--Doutname', dest = 'Doutname', default = "mysingle", help = 'output name')
#runSC = subparsers$add_parser("runSC",help="seurat normal pipeline")
#runHarmony = subparsers$add_parser("runHarmony",help="harmony pipeline")
runSC = parser$add_argument_group('runSC', 'seurat normal pipeline')
runHarmony = parser$add_argument_group('runHarmony', 'harmony pipeline')
opts = parser$parse_args()

library(SeuratDisk)
library(Seurat)
library(data.table)
library(Matrix)
library(SeuratData)
library(ggplot2)
library(dplyr)
library(harmony)
library(DoubletFinder)
library(tidyverse)
library(patchwork)

ReadFile <- function(x){
        x <- as.vector(x)
        x1 <- x[1]
        x2 <- x[2]
        if(endsWith(x2, ".rds")){
            FileTmp <- readRDS(x2)
        }else{
            FileTmp <- Read10X(data.dir=x2,gene.column = 1)
         #   colnames(FileTmp) <- paste(colnames(FileTmp),x1,sep="_")
            FileTmp <- CreateSeuratObject(counts = FileTmp, min.cells = 3, min.features = 200, project = x1)
        }
        return(FileTmp)
}

RunFilter <- function(myrds){
    myrds <- subset(myrds,subset=nFeature_RNA >= opts$minFeature && nFeature_RNA<= opts$maxFeature && percent.mt <= opts$percentMt && nCount_RNA >= opts$minCount && nCount_RNA <= opts$maxCount)
    return(myrds)
}

RunVlnPlot <- function(myrds){
    p <- VlnPlot(myrds,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3) +  labs(title =unique(myrds$orig.ident[[1]]))
    print(p)
}


RunSC <- function(scobj){
    DefaultAssay(scobj) <- "RNA"
    scobj <- SCTransform(scobj, assay = "RNA", verbose = FALSE, variable.features.n = opts$vg)
    scobj <- RunPCA(scobj)
    scobj <- RunUMAP(scobj, dims = 1:opts$dims)
    scobj <- FindNeighbors(scobj, dims = 1:opts$dims)
    scobj <- FindClusters(
      object = scobj,
      resolution = c(seq(0.1,1,0.1))
    )
    Idents(scobj) <- 'SCT_snn_res.0.5'
    DefaultAssay(scobj) <- 'RNA'
    scobj <- NormalizeData(scobj)
    return(scobj)
}

RunDoubletFinder <- function(rawrds){

	rawrds <- RunSC(rawrds) 
	DefaultAssay(rawrds) <- 'SCT'

	sweep.res.list <- paramSweep_v3(rawrds, PCs = 1:opts$dims, sct = T)
	sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
	bcmvn <- find.pK(sweep.stats)
	pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
	
	DoubletRate = opts$doubletRate
	homotypic.prop <- modelHomotypic(rawrds$SCT_snn_res.0.5)   
	nExp_poi <- round(DoubletRate*ncol(rawrds)) 
	nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
	
	rawrds <- doubletFinder_v3(rawrds, PCs = 1:opts$dims, pN = 0.25, pK = pK_bcmvn,nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
	DFTitle <- grep(pattern="DF.classifications",colnames(rawrds@meta.data),value=T)    
        rawrds$doublet_finder <- rawrds@meta.data[,DFTitle]

        print(head(rawrds@meta.data))
        print(DimPlot(rawrds, reduction = "umap", group.by = "doublet_finder"))
        print(DimPlot(rawrds, reduction = "umap", group.by = "SCT_snn_res.0.5"))
        singletrds <- subset(rawrds, doublet_finder == "Singlet")
        return(singletrds)
}


runHarmony <- function(HarmonyObj){
        DefaultAssay(HarmonyObj) <- 'RNA'
        HarmonyObj <- NormalizeData(HarmonyObj) %>% FindVariableFeatures(nfeatures=opts$vg) %>% ScaleData() %>% RunPCA(verbose=FALSE) %>% RunHarmony("orig.ident")
        harmony_embeddings <- Embeddings(HarmonyObj, 'harmony')
        HarmonyObj <- AddMetaData(HarmonyObj,harmony_embeddings[,1],col.name = 'harmony_1')
        HarmonyObj <- HarmonyObj %>% RunUMAP(reduction = "harmony", dims = 1:opts$dims) %>% FindNeighbors(reduction = "harmony", dims = 1:opts$dims)
        HarmonyObj <- FindClusters(
                object = HarmonyObj,
                resolution = c(seq(0.1,1,0.1))
        )
        Idents(HarmonyObj) <- 'RNA_snn_res.0.5'
   return(HarmonyObj)
}

#rpca

#integration

#---------------Start---------------#
if(file_test("-f",opts$input)){
    if(endsWith(opts$input,".rds")){
        mergedata <- readRDS(opts$input)
    }else{
        FileName <- read.table(opts$input)
        FileList <- apply(FileName,1,ReadFile)
        if(length(FileList) == 1){
            mergedata <- FileList[[1]]
        }else{
            mergedata <- merge(x=FileList[[1]],y=FileList[c(2:length(FileList))])
             }
    }
}else{
    FileTmp <- Read10X(data.dir=opts$input,gene.column = 1)
    mergedata <- CreateSeuratObject(counts = FileTmp, min.cells = 3, min.features = 200)   
}

mergedata[["percent.mt"]]<- PercentageFeatureSet(mergedata,pattern = opts$mtRegEx)
mytime <- format(Sys.time(), "%m%d")

if(opts$subCommand == "runPro"){
    rds.split <- SplitObject(mergedata, split.by = "orig.ident")
    pdf(paste0("RawVlnPlot.",mytime,".pdf"))
    for (i in 1:length(rds.split)) {
       RunVlnPlot(rds.split[[i]]) 
    }
    dev.off()
}

if(opts$subCommand == "runFilter"){
    mysubobj <- RunFilter(mergedata)
    saveRDS(mysubobj,paste0("AfterFilter.",mytime,".",opts$Foutname))
    pdf(paste0("AfterFilter.",mytime,".vlnplot.pdf"))
    print(VlnPlot(mysubobj,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3))
    dev.off()
}

if(opts$subCommand == "runDoubletFinder"){
    pdf(paste0("DoubletFinder.",mytime,".pdf"))
    Singletobj <- RunDoubletFinder(mergedata)
    saveRDS(Singletobj,paste0("Singlet.",mytime,".",opts$Doutname))
    dev.off()
}

if(opts$subCommand == "runSC"){
    mergedata <- RunSC(mergedata) 
    all.genes<- rownames(mergedata)
    mergedata <- ScaleData(mergedata, features=all.genes)
    saveRDS(mergedata,paste0("Seurat.normal.",mytime,".rds"))

    pdf(paste0("Seurat.normal.",mytime,".pdf"))
    print(DimPlot(mergedata,reduction="umap",label=T))
    mergedataMarkers <- FindAllMarkers(mergedata,only.pos=TRUE,min.pct=0.25,logfc.threshold=0.25)
    top10genes <- mergedataMarkers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
    print(DoHeatmap(mergedata,features=top10genes$gene) + NoLegend())
    topgenes <- mergedataMarkers %>% group_by(cluster) %>% top_n(n=opts$topn,wt=avg_log2FC)
    dev.off()

    write.table(topgenes,paste0("Seurat.normal.",mytime,".xls"),sep="\t",quote=F)
}

if(opts$subCommand == "runHarmony"){
    mergedata <- runHarmony(mergedata)
    saveRDS(mergedata,paste0("Harmony.",mytime,".rds"))
    all.genes<- rownames(mergedata)
    mergedata <- ScaleData(mergedata, features=all.genes)

    pdf(paste0("Harmony.",mytime,".pdf"))
    print(DimPlot(mergedata,reduction="umap",label=T))
    mergedataMarkers <- FindAllMarkers(mergedata,only.pos=TRUE,min.pct=0.25,logfc.threshold=0.25)
    top10genes <- mergedataMarkers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
    print(DoHeatmap(mergedata,features=top10genes) + NoLegend())
    topgenes <- mergedataMarkers %>% group_by(cluster) %>% top_n(n=opts$topn,wt=avg_log2FC)
    dev.off()
    write.table(topgenes,paste0("Hamony.",mytime,".xls"),sep="\t",quote=F)
}
