
countData=readRDS("novogene_counts.rds")
library(Rsubread)
options(scipen=999)
library(DESeq2)


design<-data.frame(experiment=colnames(countData[,c(1,6,9,10)]), batch = c("r1","r1","r2","r2"),
                                            condition = c("siC","siK", "siC","siK") )

dLRT <- DESeqDataSetFromMatrix(countData = countData[,c(1,6,9,10)], colData = design, design = ~ batch + condition )
dLRT <- DESeq(dLRT, test="LRT",full= ~ batch + condition , reduced=~ batch )
dLRT_vsd <- varianceStabilizingTransformation(dLRT)

x=assay(dLRT_vsd)
colnames(x)=c("siC_1","siK_1","siC_2","siK_2")
                     
 library(RColorBrewer)
colors <- colorRampPalette(c("white","red"))(100)

ex=x[which(rownames(x) %in% c("WNT1","LRP5","GSK3B","CTNNB1","DNMT1","SNAI2","BIRC5","TP53","BCL2","SIRT1","CDKN1A","MYC","EPCAM","CDK4","CDK6")),]

library(gplots)
heatmap.2(ex,col=colors,scale="row", trace="none")
