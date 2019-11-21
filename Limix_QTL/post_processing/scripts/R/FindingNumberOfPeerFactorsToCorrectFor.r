library(rhdf5)
library(qvalue)
library(dplyr)
#################

baseFolder <- "C:/OnlineFolders/BitSync/Current_Work/EBI_HipSci/Trans_iPSC/Expression/FeatureCounts_Genes/PEER/TestOptimal_N_Factors/"
subFolderBase <- "OutMappingPeerFactors."
chrSpecific = F
range <- 2:101
writeGlobalSig = T
writeGeneSig = F
threshold = 0.05
#multipleTestingGlobal = "BF"
multipleTestingGlobal = "ST"
topFeaturesOnly = FALSE


####
setwd(baseFolder)
results <- NULL
for(i in range){
  if(chrSpecific){
    tmp <- h5dump(file = paste(subFolderBase,i,"/qtl_results_",i,".h5",sep=""),)
  } else {
    #tmp <- h5dump(file = paste(subFolderBase,i,"/qtl_results_all.h5",sep=""),)
    tmp <- h5dump(file = paste(subFolderBase,i,"/qtl_results_2.h5",sep=""),)
  }
  for (j in names(tmp)) tmp[[j]][["feature"]] <- j
  df <- bind_rows(tmp)
  colnames(df)[which(colnames(df)=="corr_p_value")] <- "feature_corr_p_value"
  if(topFeaturesOnly){
    df <- df[match(unique(df$feature), df$feature),]
  }
  if(multipleTestingGlobal=="ST"){
    df["global_corr_p_value"] <- qvalue(df$feature_corr_p_value)$qvalues
  } else if (multipleTestingGlobal=="BF"){
    df["global_corr_p_value"] <- df$feature_corr_p_value*length(unique(df$feature))
    df$global_corr_p_value[df$global_corr_p_value>1]<-1
  }
  
  df <- df[order(df$global_corr_p_value,decreasing = F),]
  
  if(writeGlobalSig){
    print(paste(length(which(df$global_corr_p_value<threshold))," SNPs have a genetic effect on Factor: ",i,sep=""))
    write.table(paste(baseFolder,"results_global_level_",i,"_",threshold,".txt",sep=""),x = df[df$global_corr_p_value<threshold,],sep="\t",row.names=F,quote=F)
  }
  
  if(writeGeneSig){
    print(paste(length(which(df$gene_corr_p_value<threshold))," SNPs have a genetic effect on Factor: ",i," after correction",sep=""))
    write.table(paste(baseFolder,"results_gene_level_",i,"_",threshold,".txt",sep=""),x = df[df$gene_corr_p_value<threshold,],sep="\t",row.names=F,quote=F)
  } 

}

#Count number of affected genes.
affectedFeatures <- NULL
posMainEffects <- NULL
negMainEffects <- NULL
oppositeEffectsVsMain <- NULL
sharedEffectsVsMain <- NULL
for(i in range){
  tmp <- read.delim(paste(baseFolder,"results_global_level_",i,"_",threshold,".txt",sep=""),as.is=T)
  tmp["sign"] <- "+"
  tmp["key"] <- paste(tmp$feature,tmp$snp_id,sep="-")
  tmp$sign[tmp$beta<0] <- "-"
  tmp$sign[tmp$beta<0] <- "-"
  if(is.null(posMainEffects) || is.null(negMainEffects)){
    tmp.first <- tmp[match(unique(tmp$feature), tmp$feature),]
    posMainEffects = tmp.first$key[which(tmp.first$sign=="+")]
    negMainEffects = tmp.first$key[which(tmp.first$sign=="-")]
  }
  
  affectedFeatures <- c(affectedFeatures,length(unique(tmp$feature)))
  sharedEffectsVsMain <- c(sharedEffectsVsMain,sum(tmp$key[which(tmp$sign=="-")]%in%negMainEffects) + sum(tmp$key[which(tmp$sign=="+")]%in%posMainEffects))
  oppositeEffectsVsMain <- c(oppositeEffectsVsMain,sum(tmp$key[which(tmp$sign=="+")]%in%negMainEffects) + sum(tmp$key[which(tmp$sign=="-")]%in%posMainEffects))
}
plot(affectedFeatures)
#points(sharedEffectsVsMain)
plot(sharedEffectsVsMain)
plot(oppositeEffectsVsMain)
plot(sharedEffectsVsMain,oppositeEffectsVsMain)
plot(sharedEffectsVsMain,affectedFeatures)

interestMax = which(affectedFeatures == max(affectedFeatures))
#Picked 40 PEER Factors, no factors skipped.

