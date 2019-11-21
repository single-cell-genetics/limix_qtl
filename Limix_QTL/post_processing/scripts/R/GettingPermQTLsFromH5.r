library(rhdf5)
library(qvalue)
library(dplyr)

##Settings
baseFolder <- "C:/OnlineFolders/BitSync/Current_Work/EBI_HipSci/T/"

#################
##Read files.
setwd(baseFolder)
observedFeatures <- 0
results <- NULL
snpAnnotation <- NULL
featureAnnotation <- NULL
filesToRead <- list.files(".",pattern = ".h5",full.names = T)
for(i in filesToRead){
  tmp=h5dump(file = i)
  if(length(tmp)>0){
    for (j in names(tmp)) tmp[[j]][["feature"]] <- j
    observedFeatures = observedFeatures+length(tmp)
    df <- bind_rows(tmp)
    if(nrow(df)>0){
      results = rbind(results,df)
    }
  }
}
rm(df,tmp)
results["QTL"] <- paste(results$snp_id, results$feature,sep="-")
if(length(which(duplicated(results$QTL)))>0){
  results <- results[-which(duplicated(results$QTL)),]
}

write.table(paste(baseFolder,"permutationInformation.txt",sep=""),x = results,sep="\t",row.names=F,quote=F)

