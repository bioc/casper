## METHODS

setMethod("getNreads",signature(pc='pathCounts'),function(pc) {
  if(pc@stranded) {
    ans <- c(sapply(pc@counts[['plus']], sum), sapply(pc@counts[['minus']], sum))
  } else {
    ans <- sapply(pc@counts[[1]], sum)
  }
  ans
}
          )

mergePC <- function(pc, DB){
  tmp <- vector(mode='list', length=length(names(DB@transcripts)))
  names(tmp) <- names(DB@transcripts)
  for(x in 1:length(pc)){
    y <- pc[[x]]@counts[[1]]
    sel <- names(y)[unlist(lapply(y, length)>0)]
    tmp[sel] <- y[sel]
  }  
  pc <- new("pathCounts", counts=list(tmp), stranded=pc[[1]]@stranded, denovo=pc[[1]]@denovo)
  pc
}

procPaths <- function(reads, DB, mc.cores, verbose){
    if(verbose) cat("Finding overlaps between reads and exons\n")
    chrs <- as.character(unique(seqnames(reads)@values))
    DB <- subsetGenome(genomeDB=DB, chr=chrs)    
    over<-findOverlaps(reads, DB@exonsNI)    
    if("id" %in% colnames(values(reads))) {
      readid<-values(reads)$id[queryHits(over)]
    } else if("names" %in% colnames(values(reads))) {
      readid<-values(reads)$names[queryHits(over)]
    } else readid <- as.integer(names(reads)[queryHits(over)])
    readside<-values(reads)$rid[queryHits(over)]
    exid <- as.integer(names(DB@exonsNI)[subjectHits(over)])
    exst<-start(DB@exonsNI)[subjectHits(over)]
    if(verbose) cat("Counting paths\n")
    rm(over)
    gc()
    pCounts<-.Call("pathCounts", readid, readside, exst, exid)
    rm(readid, readside, exst, exid)
    pCounts<-lapply(pCounts[1:2], function(x) x[1:pCounts[[3]]])
    names(pCounts[[2]])<-pCounts[[1]]
    pCounts<-pCounts[[2]]
    pCounts <- pCounts[grepl("-[0-9]", names(pCounts))]
    
    sel <- strsplit(names(pCounts), split='-|\\.')
    sel1 <- lapply(sel, "[", 2)
    sel1 <- unlist(sel1)
    
    nislEx <- elementLengths(DB@islands)
    nislEx <- rep(names(DB@islands), unlist(nislEx))
    islEx <- names(unlist(DB@islands))
    islEx <- sub(".*\\.","", islEx)
    names(islEx) <- nislEx
    isl <- match(sel1, islEx)
    isl <- names(islEx)[isl]
    splCounts <- split(pCounts, isl)
    splCounts <- lapply(splCounts, function(x) x[grepl("-", names(x))])
    
    if(DB@denovo){
      sel <- sapply(sel, "[", -1)
      tmp <- split(sel, isl)
      if(mc.cores>1) {
        require(parallel)
        tmp1 <- parallel::mclapply(names(tmp), function(x){
          n <- sapply(tmp[[x]], length)
          nn <- unlist(tmp[[x]])
          names(nn) <- rep(names(splCounts[[x]]), n)
          nnn <- nn %in% names(DB@islands[[x]])
          nnnn <- tapply(nnn, names(nn), all)
          splCounts[[x]][names(splCounts[[x]]) %in% names(nnnn)[nnnn]]
        }, mc.cores=mc.cores)
      } else {
        tmp1 <- lapply(names(tmp), function(x){
          n <- sapply(tmp[[x]], length)
          nn <- unlist(tmp[[x]])
          names(nn) <- rep(names(splCounts[[x]]), n)
          nnn <- nn %in% names(DB@islands[[x]])
          nnnn <- tapply(nnn, names(nn), all)
          splCounts[[x]][names(splCounts[[x]]) %in% names(nnnn)[nnnn]]
        })
      }
      names(tmp1) <- names(tmp)
      splCounts <- tmp1
    }
    
    ans <- vector(length(DB@islands), mode='list')
    names(ans) <- names(DB@islands)
    ans[names(splCounts)] <- splCounts
    ans
  }

setGeneric("pathCounts", function(reads, DB, mc.cores=1, verbose=FALSE) standardGeneric("pathCounts"))
setMethod("pathCounts", signature(reads='procBam'), function(reads, DB, mc.cores, verbose) {
  if (class(reads) != 'procBam') stop('reads must be an object of class procBam')
  if(!reads@stranded) {
    counts <- procPaths(reads=reads@pbam, DB=DB, mc.cores=mc.cores, verbose=verbose)
    ans <- new("pathCounts", counts=list(counts), denovo=DB@denovo, stranded=reads@stranded)
  }
  else {
    plusDB <- genomeBystrand(DB, strand="+")
    plus <- procPaths(reads@plus, plusDB, mc.cores)
    minusDB <- genomeBystrand(DB, strand="-")
    minus <- procPaths(reads@minus, minusDB, mc.cores)
    ans <- new("pathCounts", counts=list(plus=plus, minus=minus), denovo=DB@denovo, stranded=reads@stranded)
  }
  ans
})
setMethod("pathCounts", signature(reads='list'), function(reads, DB, mc.cores, verbose) {
  ans <- mclapply(reads, function(x) pathCounts(x, DB=DB, mc.cores=1, verbose=verbose), mc.cores=mc.cores, mc.preschedulle=TRUE)
  mergePC(pc=ans, DB=DB)
})
          
