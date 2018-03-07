
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
colnames(x)=c("siControl_1","siKAT5_1","siControl_2","siKAT5_2")
                     
 library(RColorBrewer)
colors <- colorRampPalette(c("blue","white","red"))(10)

ex=x[which(rownames(x) %in% c("MDM2","TP53","CDKN1A","PUMA","CDK4","CDK6","SIRT1","BCL2","BIRC5","UHRF1","KAT5",
                              "GSK3B","LRP5","WNT1","CTNNB1","MYC","NICD","NOTCH1","TGFB1","TGFB2","TGFB3","TGFB4","SNAI1",
                              "SNAI2","DNMT1","DNMT3A","DNMT3B","EPCAM","TWIST1","AKT1","AKT2","ZEB1","ZEB2","CDH1","CDH2"
                             )),c(1,3,2,4)]

dataCenter = function(x) { return(t(apply(x, 1, function(y) y - mean(y) ))) }
dataCenter_timeC = function(x) { return(t(
                                  apply(x, 1, funtion(xtest) { unlist(c( xtest[1:4]-mean(as.numeric(xtest[1:4])), 
                                                                xtest[5:8]-mean(as.numeric(xtest[5:8])), 
                                                                xtest[9:12]-mean(as.numeric(xtest[9:12])) ) )
                } )
                                        )) }

                                          
library(gplots)
library(factoextra)

pdf("TIP60KD_model_genes.pdf")
heatmap.2(dataCenter(ex),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="spearman"),srtCol=25,main="HCT116")
dev.off()

illu=read.table("/Users/wone/CSI/Illu-Quant-expression.txt",sep="\t",row.name=1,header=T)
          
exillu=illu[which(rownames(illu) %in% c("MDM2","TP53","CDKN1A","PUMA","CDK4","CDK6","SIRT1","BCL2","BIRC5","UHRF1","KAT5",
                              "GSK3B","LRP5","WNT1","CTNNB1","MYC","NICD","NOTCH1","TGFB1","TGFB2","TGFB3","TGFB4","SNAI1",
                              "SNAI2","DNMT1","DNMT3A","DNMT3B","EPCAM","TWIST1","AKT1","AKT2","ZEB1","ZEB2","CDH1","CDH2"
                             )),]
          
colnames(exillu)=c("WT_0h_1","WT_24h_1","WT_48h_1","WT_72h_1",
                 "DNMT1_0h_1","DNMT1_24h_1","DNMT1_48h_1","DNMT1_72h_1",
                 "TP53_0h_1","TP53_24h_1","TP53_48h_1","TP53_72h_1",
                 "WT_0h_2","WT_24h_2","WT_48h_2","WT_72h_2",
                 "DNMT1_0h_2","DNMT1_24h_2","DNMT1_48h_2","DNMT1_72h_2",
                 "TP53_0h_2","TP53_24h_2","TP53_48h_2","TP53_72h_2")
          
exillu_mean=(exillu[,1:12]+exillu[,13:24])/2
colnames(exillu_mean) = gsub("_1","",colnames(exillu_mean))

colors <- colorRampPalette(c("blue","white","red"))(10)
pdf("hct116_0h.pdf")
heatmap.2(dataCenter(exillu[,c(1,13,9,21,5,17)]),col=colors,scale="none", trace="none",Colv=FALSE,dendrogram="row",
          distfun = function(x) get_dist(x,method="spearman"),srtCol=25,main="HCT116 0h")
dev.off()
pdf("hct116_24h.pdf")
heatmap.2(dataCenter(exillu[,c(1,13,9,21,5,17)+1]),col=colors,scale="none", trace="none",Colv=FALSE,dendrogram="row",
          distfun = function(x) get_dist(x,method="spearman"),srtCol=25,main="HCT116 24h")
dev.off()

pdf("hct116_48h.pdf")
heatmap.2(dataCenter(exillu[,c(1,13,9,21,5,17)+2]),col=colors,scale="none", trace="none",Colv=FALSE,dendrogram="row",
          distfun = function(x) get_dist(x,method="spearman"),srtCol=25,main="HCT116 48h")
dev.off()

pdf("hct116_72h.pdf")
heatmap.2(dataCenter(exillu[,c(1,13,9,21,5,17)+3]),col=colors,scale="none", trace="none",Colv=FALSE,dendrogram="row",
          distfun = function(x) get_dist(x,method="spearman"),srtCol=25,main="HCT116 72h")
dev.off()

######################################################################################################
center_exillu_mean= t(apply(exillu_mean,1,function(xtest){
unlist(c( xtest[1:4]-mean(as.numeric(xtest[1:4])), xtest[5:8]-mean(as.numeric(xtest[5:8])), xtest[9:12]-mean(as.numeric(xtest[9:12])) ) )}))
          

heatmap.2(center_exillu_mean[,c(1,5,9)+1],col=colors,scale="none", trace="none",dendrogram="row",
          distfun = function(x) get_dist(x,method="spearman"),srtCol=45,main="HCT116",Colv=FALSE)
######################################################################################################        
          
