library(matrixStats)
library(readr)
#### gene annotation
gene_anno = read.delim("/vol/projects/wli/projects/genetics/sw/limix/data/anno.nochr.txt",as.is=T)
testCombinations = NULL
#
nGenes = 50
startPos = 0
endOffset = 1000000000
sink("/vol/projects/wli/projects/genetics/sw/limix/data/chunks.txt")
for(chr in unique(gene_anno$chromosome)){
  #print(chr)
  annotationRel = gene_anno[which(gene_anno$chromosome==chr),]
  annotationRel = annotationRel[order(annotationRel$start,annotationRel$end),]
  ##First go through the list to fix genes to they all touch.
  annotationRel$start[1] = startPos
  for(i in 2:nrow(annotationRel)){
    if(i == nrow(annotationRel)){
      annotationRel$end[i] = annotationRel$end[i]+endOffset
    }
    #If "overlapping" than we don't need to do anything.
    if((annotationRel$start[i]>annotationRel$end[i-1])){
      #print(i)
      distance = (annotationRel$start[i]-annotationRel$end[i-1])/2
      annotationRel$start[i] = ceiling(annotationRel$start[i]-distance)
      annotationRel$end[i-1] = floor(annotationRel$end[i-1]+distance)
    }
  }
  chunks = seq(1, nrow(annotationRel),nGenes)
  #Need to add the last one as a separate entry.
  if(chunks[length(chunks)] < nrow(annotationRel)){
    chunks = c(chunks,(nrow(annotationRel)+1))
  }
  for(i in 1:(length(chunks)-1)){
    print(paste(chr,":",annotationRel$start[chunks[i]],"-",annotationRel$end[(chunks[i+1]-1)],sep=""))
  }
}
sink()