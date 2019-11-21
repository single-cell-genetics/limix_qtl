library(rhdf5)
library(qvalue)
library(dplyr)

##Settings
baseFolder <- "./"

snpName = "12_28156081_C_G"
featureName = ""

# multipleTestingGlobal = "BF"
#################
##Read files.
setwd(baseFolder)
observedFeatures <- 0
results <- NULL
filesToRead <- list.files("./",pattern = ".h5",full.names = T)
if(length(filesToRead)==1){
  tmp <- h5dump(file = "./qtl_results_all.h5")
  if(length(tmp)>0){
    for (j in names(tmp)) tmp[[j]][["feature"]] <- j
    observedFeatures = observedFeatures+length(tmp)
    df <- bind_rows(tmp)
    if(multipleTestingGlobal=="BF"){
      df <- df[df$corr_p_value<threshold,]
    }
    if(nrow(df)>0){
      results = rbind(results,df)  
    }
  }
} else {
  for(i in filesToRead){
    baseName <- gsub(gsub(i,pattern = ".h5",replacement = ""),pattern = "./qtl_results",replacement = "")
    if(length(list.files("./",pattern = paste(baseName,"\\.",sep="")))==3 | length(list.files("./",pattern = paste(baseName,"\\.",sep="")))==4){
      tmp <- h5dump(file = i)
    } else {
      print(paste("Skipping",i,"because not necessary data is available or to many files with the same name."))
    }
    if(length(tmp)>0){
      for (j in names(tmp)) tmp[[j]][["feature"]] <- j
      observedFeatures = observedFeatures+length(tmp)
      df <- bind_rows(tmp)
      if(length(snpName)>0){
        if(length(which(df$snp_id==snpName))>0){
          df = df[which(df$snp_id==snpName),]
        } else {
          next
        }
      }
      if(featureName!=""){
        if(length(which(df$feature==featureName))>0){
          df = df[which(df$feature==featureName),]
        } else {
          next
        }
      }
      if(nrow(df)>0){
        results = rbind(results,df)  
      }
    }
  }
}
colnames(results)[which(colnames(results)=="corr_p_value")] <- "feature_corr_p_value"
if(length(which(is.na(results$feature_corr_p_value)))!=0){
  results <- results[-which(is.na(results$feature_corr_p_value)),]
}

##Multiple testing
if(multipleTestingGlobal=="ST"){
  results["global_corr_p_value"] <- qvalue(results$feature_corr_p_value)$qvalues
} else if (multipleTestingGlobal=="BF"){
  results["global_corr_p_value"] <- results$feature_corr_p_value*observedFeatures
  results$global_corr_p_value[results$global_corr_p_value>1]<-1
}


results <- results[order(results$global_corr_p_value,decreasing = F),]

if(writeGlobalSigTop){
  write.table(paste(baseFolder,"top_results_global_level_",threshold,".txt",sep=""),x = results[intersect(which(results$global_corr_p_value<threshold),which(!duplicated(results$feature))),],sep="\t",row.names=F,quote=F)
}

if(writeGlobalSig){
  write.table(paste(baseFolder,"results_global_level_",threshold,".txt",sep=""),x = results[results$global_corr_p_value<threshold,],sep="\t",row.names=F,quote=F)
}

if(writeFeatureSig){
  write.table(paste(baseFolder,"results_gene_level_",threshold,".txt",sep=""),x = results[results$feature_corr_p_value<threshold,],sep="\t",row.names=F,quote=F)
}

if(writeNominalSig){
  write.table(paste(baseFolder,"results_nominal_level_",threshold2,".txt",sep=""),x = results[results$p_value<threshold2,],sep="\t",row.names=F,quote=F)
}
