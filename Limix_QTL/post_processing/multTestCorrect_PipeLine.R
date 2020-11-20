library("qvalue")
library("multtest")
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<3 || length(args)>4) {
  stop("At least three argument must be supplied (input folder, output folder, FDR level). (Optionally a filtering file to filter genes).n", call.=FALSE)
}

##Parameters
folder = args[1]
folderOut = args[2]
fdrLevel = as.numeric(args[3])

filterFile= NULL
if(length(args)==4){
	filterFile= args[4]
}

qtlResults <- read.delim(paste(folder,"/top_qtl_results_all.txt",sep=""), as.is=T)

qtlResults <- qtlResults[order(qtlResults$empirical_feature_p_value, qtlResults$p_value),]
qtlResults["QTL"] <- paste(qtlResults$feature_id,qtlResults$snp_id,sep = "-")

#Check for double QTLs.
if(length(unique(qtlResults$QTL)) != nrow(qtlResults)){
	qtlResults <- qtlResults[-which(duplicated(qtlResults$QTL)),]
}
#Check for double genes.
if(length(unique(qtlResults$feature_id)) != nrow(qtlResults)){
	qtlResults <- qtlResults[-which(duplicated(qtlResults$feature_id)),]
}

#Filtering to a subset of the genes.
if(!is.null(filterFile)){
	relevantGenes <- read.delim(filterFile, as.is=T, header=F)[,1]
	qtlResults <- qtlResults[which(qtlResults$feature_id %in% relevantGenes),]
}

qtlResults$global_corrected_pValue <- qvalue(qtlResults$empirical_feature_p_value)$qvalues
qtlResults$global_corrected_pValue_BH <- multtest::mt.rawp2adjp(qtlResults$empirical_feature_p_value,proc = "BH")$adjp[,2]
qtlResults$global_corrected_pValue_BF <- multtest::mt.rawp2adjp(qtlResults$empirical_feature_p_value,proc = "Bonferroni")$adjp[,2]

##To extract data from the full file.
minimal_p <- max(qtlResults$empirical_feature_p_value[which(qtlResults$global_corrected_pValue<fdrLevel)])

write.table(qtlResults, paste(folderOut,"/top_qtl_results_all_FDR.txt",sep=""), quote = F, row.names=F, sep="\t" )
write.table(qtlResults[qtlResults$global_corrected_pValue<fdrLevel,], paste(folderOut,"/top_qtl_results_all_FDR",fdrLevel,".txt",sep=""), quote = F, row.names=F, sep="\t" )

qtlResultsAll <- read.delim(paste(folder,"/qtl_results_all.txt",sep=""), as.is=T)
qtlResultsAll["QTL"] <- paste(qtlResultsAll$feature_id,qtlResultsAll$snp_id,sep = "-")

qtlResultsAll <- qtlResultsAll[order(qtlResultsAll$empirical_feature_p_value, qtlResultsAll$p_value),]
qtlResultsAll["QTL"] <- paste(qtlResultsAll$feature_id,qtlResultsAll$snp_id,sep = "-")
#Check for double QTLs.
if(length(unique(qtlResultsAll$QTL)) != nrow(qtlResultsAll)){
	qtlResultsAll <- qtlResultsAll[-which(duplicated(qtlResultsAll$QTL)),]
}

if(!is.null(filterFile)){
	qtlResultsAll <- qtlResultsAll[which(qtlResultsAll$feature_id %in% relevantGenes),]
}

write.table(qtlResultsAll[qtlResultsAll$empirical_feature_p_value<=minimal_p,], paste(folderOut,"/qtl_results_all_FDR",fdrLevel,".txt",sep=""), quote = F, row.names=F, sep="\t" )

