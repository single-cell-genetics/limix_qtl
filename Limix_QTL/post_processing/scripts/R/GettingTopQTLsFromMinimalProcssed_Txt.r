library(rhdf5)
library(qvalue)
library(dplyr)
library(readr)

filterFile <- "C:/OnlineFolders/BitSync/Current_Work/EBI_HipSci/Trans_iPSC/Expression/FeatureCounts_Genes/Paired/Filter_Min_tpm_0.5_Abundant_in_20perc_Samp.txt"
#filterFile <-  "/nfs/research2/hipsci/processed_data/rna_seq/intersectionProteomics/PEER/TestOptimal_N_Factors/Filter_Abundant_in_50perc_Samp.txt"
filterFeatures <- NULL
if(!(filterFile=="" || is.na(filterFile))){
  filterFeatures = read.delim(filterFile,as.is=T,header=F)[,1]
}

##Settings
baseFolder <- "C:/OnlineFolders/BitSync/Current_Work/EBI_HipSci/QTL_Output/Replication_eQTL_Gen/eQTL_Gen_Cis/"
baseFolder <- "C:/OnlineFolders/BitSync/Current_Work/EBI_HipSci/QTL_Output/Replication_eQTL_Gen/eQTL_Gen_Cis_5p/"

writeGlobalSig = T
writeGlobalSigTop = T
writeFeatureSig = T
writeNominalSig = F
topResultBased = T
threshold = 0.05
threshold2 = 1.0
multipleTestingGlobal = "ST"
# multipleTestingGlobal = "BF"
#################
##Read files.
setwd(baseFolder)
results <- NULL
# filesToRead <- list.files(".",pattern = ".gz",full.names = T)
# 
# for(i in filesToRead){
#   tmp <- readr::read_delim(file = i,delim="\t")
#   if(length(tmp)>0){
#     #if(length(grep(colnames(df),pattern = "n_samples"))>0){
#     #  df <- df[,-grep(colnames(df),pattern = "n_samples")]
#     #}
#     if(!is.null(filterFeatures)){
#       tmp <- tmp[which(tmp$feature_id%in%filterFeatures),]
#     }
#     
#     if(nrow(tmp)>0){
#       results = rbind(results,tmp)  
#     }
#   }
# }
results <- read.delim("./top_qtl_results_all.txt",as.is=T)
results <- results[which(results$feature_id %in% filterFeatures),]
##Remove NaN values and other failed empirical recalibrations.
results <- results[-which(is.na(results$empirical_feature_p_value)),]
results <- results[-which(results$empirical_feature_p_value<0),]

results <- results[order(results$empirical_feature_p_value, results$p_value,decreasing = F),]

if(topResultBased){
  resultsFull <- results
  results <- resultsFull[which(!duplicated(resultsFull$feature_id)),]
}


##Multiple testing
if(multipleTestingGlobal=="ST"){
  results["global_corr_p_value"] <- qvalue(results$empirical_feature_p_value)$qvalues
} else if (multipleTestingGlobal=="BF"){
  observedFeatures <- length(unique(results$feature_id))
  results["global_corr_p_value"] <- results$empirical_feature_p_value*observedFeatures
  results$global_corr_p_value[results$global_corr_p_value>1]<-1
}

results <- results[order(results$global_corr_p_value,results$empirical_feature_p_value, results$p_value,decreasing = F),]
length(which(results$global_corr_p_value<threshold))

if(writeGlobalSigTop){
  write.table(paste(baseFolder,"top_results_global_level_",threshold,".txt",sep=""),x = results[intersect(which(results$global_corr_p_value<threshold),which(!duplicated(results$feature))),],sep="\t",row.names=F,quote=F)
}


if(writeGlobalSigTop & topResultBased){
  write.table(paste(baseFolder,"top_results_global_level_",threshold,".txt",sep=""),x = results[intersect(which(results$global_corr_p_value<threshold),which(!duplicated(results$feature))),],sep="\t",row.names=F,quote=F)
}

if(writeGlobalSig & !topResultBased){
  write.table(paste(baseFolder,"results_global_level_",threshold,".txt",sep=""),x = results[results$global_corr_p_value<threshold,],sep="\t",row.names=F,quote=F)
}

if(writeFeatureSig & !topResultBased){
  write.table(paste(baseFolder,"results_gene_level_",threshold,".txt",sep=""),x = results[results$feature_corr_p_value<threshold,],sep="\t",row.names=F,quote=F)
}

if(writeNominalSig){
  write.table(paste(baseFolder,"results_nominal_level_",threshold2,".txt",sep=""),x = resultsFull[resultsFull$p_value<threshold2,],sep="\t",row.names=F,quote=F)
}

if(writeGlobalSig & topResultBased){
  if(length(which(results$global_corr_p_value<threshold))>0){
    minFeatureCorrected <- max(results$empirical_feature_p_value[which(results$global_corr_p_value<threshold)])
    write.table(paste(baseFolder,"results_global_level_",threshold,".txt",sep=""),x = resultsFull[resultsFull$empirical_feature_p_value<minFeatureCorrected,],sep="\t",row.names=F,quote=F)    
  }
}
