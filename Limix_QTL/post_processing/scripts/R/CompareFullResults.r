#
baseFolder1 <- "./Gene_Mapping1/"
baseFolder2 <- "./GeneLevel2/"
threshold = 0.05
name1 <- "d1"
name2 <- "d2"

##Feature level comparison
fdr_Sig_Full_1 <- read.delim(paste(baseFolder1,"results_global_level_",threshold,".txt",sep=""),as.is=T)
fdr_Sig_Full_2 <- read.delim(paste(baseFolder2,"results_global_level_",threshold,".txt",sep=""),as.is=T)

uniqueFeatures1 <- unique(fdr_Sig_Full_1$feature)
uniqueFeatures2 <- unique(fdr_Sig_Full_2$feature)
library("VennDiagram")
draw.pairwise.venn(area1 = length(uniqueFeatures1) ,area2 = length(uniqueFeatures2),cross.area = length(which(uniqueFeatures1%in%uniqueFeatures2)),category = c(name1,name2))

##Top level QTL comparison based on 1.
fdr_Sig_1 <- read.delim(paste(baseFolder1,"top_results_global_level_",threshold,".txt",sep=""),as.is=T)
#fdr_Sig_Full_2 <- read.delim(paste(baseFolder2,"results_global_level_",threshold,".txt",sep=""),as.is=T)
fdr_Sig_Full_2 <- read.delim(paste(baseFolder2,"results_gene_level_",threshold,".txt",sep=""),as.is=T)
fdr_Sig_Full_2 <- unique(fdr_Sig_Full_2)

fdr_Sig_1["QTL"] <- paste(fdr_Sig_1$snp_id,fdr_Sig_1$feature,sep=":")
fdr_Sig_Full_2["QTL"] <- paste(fdr_Sig_Full_2$snp_id,fdr_Sig_Full_2$feature,sep=":")
org <- nrow(fdr_Sig_1)
fdr_Sig_1 <- fdr_Sig_1[which(fdr_Sig_1$QTL%in%fdr_Sig_Full_2$QTL),]
print(paste("Fraction of top effects from 1 in dataset 2: ",nrow(fdr_Sig_1)/org))
fdr_Sig_Full_2 <- fdr_Sig_Full_2[which(fdr_Sig_Full_2$QTL%in%fdr_Sig_1$QTL),]
dim(fdr_Sig_1)
dim(fdr_Sig_Full_2)

fdr_Sig_1 <- fdr_Sig_1[order(fdr_Sig_1$QTL),]
fdr_Sig_Full_2 <- fdr_Sig_Full_2[order(fdr_Sig_Full_2$QTL),]
all(fdr_Sig_1$QTL==fdr_Sig_Full_2$QTL)

max <- max(c(abs(fdr_Sig_1$beta),abs(fdr_Sig_Full_2$beta)))
plot(fdr_Sig_1$beta,fdr_Sig_Full_2$beta,xlim = c(-max,max),ylim = c(-max,max))

##Top level QTL comparison based on 2.
#fdr_Sig_Full_1 <- read.delim(paste(baseFolder1,"results_gene_level_",threshold,".txt",sep=""),as.is=T)
fdr_Sig_Full_1 <- read.delim(paste(baseFolder1,"results_global_level_",threshold,".txt",sep=""),as.is=T)
fdr_Sig_2 <- read.delim(paste(baseFolder2,"top_results_global_level_",threshold,".txt",sep=""),as.is=T)

fdr_Sig_2["QTL"] <- paste(fdr_Sig_2$snp_id,fdr_Sig_2$feature,sep=":")
fdr_Sig_Full_1["QTL"] <- paste(fdr_Sig_Full_1$snp_id,fdr_Sig_Full_1$feature,sep=":")
org <- nrow(fdr_Sig_2)
fdr_Sig_2 <- fdr_Sig_2[which(fdr_Sig_2$QTL%in%fdr_Sig_Full_1$QTL),]
print(paste("Fraction of top effects from 2 in dataset 1: ",nrow(fdr_Sig_2)/org))
fdr_Sig_Full_1 <- fdr_Sig_Full_1[which(fdr_Sig_Full_1$QTL%in%fdr_Sig_2$QTL),]

fdr_Sig_2 <- fdr_Sig_2[order(fdr_Sig_2$QTL),]
fdr_Sig_Full_1 <- fdr_Sig_Full_1[order(fdr_Sig_Full_1$QTL),]
all(fdr_Sig_Full_1$QTL==fdr_Sig_2$QTL)

max <- max(c(abs(fdr_Sig_Full_1$beta),abs(fdr_Sig_2$beta)))
plot(unique(fdr_Sig_Full_1)$beta,fdr_Sig_2$beta,xlim = c(-max,max),ylim = c(-max,max))
