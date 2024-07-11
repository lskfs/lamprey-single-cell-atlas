library(enrichR)
library(ggplot2)
listEnrichrSites()
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2023")

deg = read.csv(file)
    deg <- deg[deg$avg_log2FC>1,]
    file = gsub("\\.csv", "", file)
    dir.create(file)
    cellname <- unique(deg$cluster)
    for(i in cellname){
        #enrichment
        print(i)
        genelist <- deg[deg$cluster==i,]$description
        enriched <- enrichr(c(genelist), dbs)
        gobp <- enriched[[1]]
        kegg <- enriched[[2]]
        gobp$Term <- gsub("\\(.*", "", gobp$Term)
        gobp$logP = log10(gobp$P.value)
        gobp['-logP'] = -gobp$logP
        kegg$logP = log10(kegg$P.value)
        kegg['-logP'] = -kegg$logP
        dir.create(paste0(file,'/',i))
        #setwd(paste0(file,'/',i))
        filename =paste0(file,'/',i,'/')
        write.csv(gobp,paste0(paste0(filename,file,'_',i,"gobp.csv")),quote=F)
        write.csv(kegg,paste0(paste0(filename,file,'_',i,"kegg.csv")),quote=F)
        #gobp
        gobp<- gobp[!duplicated(gobp$Term), ]
        labels = gobp[order(gobp$`-logP`, decreasing = F), 'Term']
        gobp$Term = factor(gobp$Term,levels=labels)
        mixsc <- max(gobp$`-logP`)
        goplot = ggplot(data=gobp[1:20,]) +
        geom_bar(aes(x=Term, y=-logP, fill=-logP), stat='identity',width = 0.75) +
        coord_flip() + 
        scale_fill_gradient(low="#dddfe6", high = "#6C49B8", limits=c(0,mixsc)) + 
        xlab("Term") + 
        ylab("-logP") + 
        theme(axis.text.x=element_text(color="black", size=10), axis.text.y=element_text(color="black", size=10)) + 
        scale_y_continuous(expand=c(0, 0.2)) + 
        scale_x_discrete(expand=c(0.1,0.1)) + 
        theme_bw()
        pdf(paste0(filename,file,'_',i,"gobp.pdf",sep=""),width=7,height=7)
        print(goplot)
        dev.off()
    }