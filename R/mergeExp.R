mergeExp <- function(..., sampleNames, keep=c('transcript','gene_id','island_id')) {
  esets <- list(...)
  if (is.list(esets[[1]])) esets <- unlist(esets)
  if (length(unique(sapply(esets,nrow))) != 1) stop('Number of rows do not match')
  if (any(sapply(esets,ncol) != 1)) {
    sampleNames <- unlist(lapply(esets, Biobase::sampleNames))
  } else {   #all ExpressionSets have 1 column
    if (missing(sampleNames)) {
      sampleNames <- paste('Sample',1:length(esets),sep='')
    } else {
      if (length(sampleNames)!=length(esets)) stop('sampleNames length does not match the input')
    }
  }
  if (length(esets)>1) {
    n <- featureNames(esets[[1]])
    readCount <- rowSums(do.call("cbind",lapply(esets, function(z) fData(z)[n,grep('explCnts',fvarLabels(z))])))
    sel <- !(keep %in% c('transcript','gene_id','island_id'))
    cursel <- c('transcript','gene_id','island_id')
    if(sum(sel)>0) cursel <- c('transcript','gene_id','island_id',fvarLabels(esets[[1]])[grep(keep[sel],fvarLabels(esets[[1]]))])
    fData(esets[[1]]) <- fData(esets[[1]])[,cursel]
    fvarLabels(esets[[1]])[-1:-3] <- paste(sampleNames[1],fvarLabels(esets[[1]])[-1:-3],sep='.')
    for (i in 2:length(esets)) {
      #if (any(!(keep %in% fvarLabels(esets[[i]])))) stop("Some variables in argument 'keep' were not found")
        cursel <- c('transcript','gene_id','island_id')   
        if(sum(sel)>0) cursel <- c('transcript','gene_id','island_id',fvarLabels(esets[[i]])[grep(keep[sel],fvarLabels(esets[[i]]))])
        fData(esets[[i]]) <- fData(esets[[i]])[,cursel]
        fvarLabels(esets[[i]])[-1:-3] <- paste(sampleNames[i],fvarLabels(esets[[i]])[-1:-3],sep='.')
        esets[[i]] <- esets[[i]][n,]
#      sampleNames(esets[[i]]) <- sampleNames[i]
    }
    ans <- do.call("combine",esets)
    if (ncol(ans) != length(sampleNames)) stop("Length of sampleNames does not match output. Check that all input sampleNames are different")
    sampleNames(ans) <- sampleNames
  } else {
    readCount <- fData(esets[[1]])[,'explCnts']
    ans <- esets[[1]]
  }
  fData(ans)$readCount <- readCount
  return(ans)
}

