HeuristicMutExSorting<-function(mutPatterns){

  mutPatterns<-sign(mutPatterns)

  if(dim(mutPatterns)[1]==1){
    mutPatterns<-matrix(c(mutPatterns[,order(mutPatterns,decreasing=TRUE)]),
                        1,ncol(mutPatterns),
                        dimnames = list(rownames(mutPatterns),colnames(mutPatterns)))

    return(mutPatterns)
  }

  if(dim(mutPatterns)[2]==1){
    mutPatterns<-matrix(c(mutPatterns[order(mutPatterns,decreasing=TRUE),]),
                        nrow(mutPatterns),1,
                        dimnames = list(rownames(mutPatters),colnames(mutPatterns)))
    return(mutPatterns)
  }

  nsamples<-ncol(mutPatterns)

  coveredGenes<-NA
  uncoveredGenes<-rownames(mutPatterns)

  if (length(uncoveredGenes)>1){

    idNull<-which(colSums(mutPatterns)==0)
    nullCol<-matrix(c(mutPatterns[,idNull]),nrow(mutPatterns),length(idNull),dimnames = list(rownames(mutPatterns),colnames(mutPatterns)[idNull]))

    idNonNull<-which(colSums(mutPatterns)>0)
    mutPatterns<-matrix(c(mutPatterns[,idNonNull]),nrow(mutPatterns),length(idNonNull),dimnames=list(rownames(mutPatterns),colnames(mutPatterns)[idNonNull]))

    coveredSamples<-NA
    uncoveredSamples<-colnames(mutPatterns)
    BS<-NA

    while(length(uncoveredGenes)>0 & length(uncoveredSamples)>0){

      patterns<-matrix(c(mutPatterns[uncoveredGenes,uncoveredSamples]),
                       nrow = length(uncoveredGenes),
                       ncol = length(uncoveredSamples),
                       dimnames = list(uncoveredGenes,uncoveredSamples))

      if(length(uncoveredGenes)>1){
        bestInClass<-findBestInClass(patterns)
      }else{
        bestInClass<-uncoveredGenes
      }

      if(is.na(BS[1])){
        BS<-bestInClass
      }else{
        BS<-c(BS,bestInClass)
      }

      if(is.na(coveredGenes[1])){
        coveredGenes<-bestInClass
      }else{
        coveredGenes<-c(coveredGenes,bestInClass)
      }

      uncoveredGenes<-setdiff(uncoveredGenes,coveredGenes)
      toCheck<-matrix(c(patterns[bestInClass,uncoveredSamples]),nrow = 1,ncol=ncol(patterns),dimnames = list(bestInClass,uncoveredSamples))

      if (length(coveredGenes)==1){
        coveredSamples<-names(which(colSums(toCheck)>0))
      }else{
        coveredSamples<-c(coveredSamples,names(which(colSums(toCheck)>0)))
      }

      uncoveredSamples<-setdiff(uncoveredSamples,coveredSamples)

    }

    BS<-c(BS,uncoveredGenes)

    CID<-rearrangeMatrix(mutPatterns,BS)

    FINALMAT<-mutPatterns[BS,CID]

    FINALMAT<-cbind(FINALMAT,nullCol[rownames(FINALMAT),])

    return(FINALMAT)
  }

}
findBestInClass<-function(patterns){

  if(nrow(patterns)==1){
    return(rownames(patterns))
  }

  if(ncol(patterns)==1){
    return(rownames(patterns)[1])
  }

  genes<-rownames(patterns)

  exclCov<-rep(NA,length(genes))
  names(exclCov)<-genes
  for (g in genes){
    residGenes<-setdiff(genes,g)
    if (length(residGenes)>1){
      exclCov[g]<-sum(patterns[g,]-colSums(patterns[residGenes,]))
    }else{
      exclCov[g]<-sum(patterns[g,]-patterns[residGenes,])
    }
  }

  return(names(sort(exclCov,decreasing=TRUE))[1])
}
rearrangeMatrix<-function(patterns,GENES){

  remainingSamples<-colnames(patterns)

  toAdd<-NULL

  for (g in GENES){
    remainingGenes<-setdiff(GENES,g)

    P1<-matrix(c(patterns[g,remainingSamples]),length(g),length(remainingSamples),dimnames = list(g,remainingSamples))
    P2<-matrix(c(patterns[remainingGenes,remainingSamples]),length(remainingGenes),length(remainingSamples),
               dimnames=list(remainingGenes,remainingSamples))

    if(length(remainingGenes)>1){
      DD<-colnames(P1)[order(P1-colSums(P2),decreasing=TRUE)]
    }else{
      DD<-colnames(P1)[order(P1-P2,decreasing=TRUE)]
    }

    toAdd<-c(toAdd,names(which(patterns[g,DD]>0)))
    remainingSamples<-setdiff(remainingSamples,toAdd)
    if(length(remainingSamples)==0){
      break
    }
  }

  toAdd<-c(toAdd,remainingSamples)

  return(toAdd)
}
