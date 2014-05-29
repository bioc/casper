mergeExp <- function(..., sampleNames,  keep=c('transcript','gene_id','island_id')) {
  esets <- list(...)
  if(class(esets[[1]])=='list') esets <- unlist(esets)
  if (length(unique(sapply(esets,nrow))) != 1) stop('Number of rows do not match')
  if (any(sapply(esets,ncol) != 1)) stop("All ExpressionSet objects in '...' should have 1 column")
  if (missing(sampleNames)) sampleNames <- paste('Sample',1:length(esets),sep='')
  if (length(esets)>1) {
    n <- featureNames(esets[[1]])
    readCount <- rowSums(do.call("cbind",lapply(esets, function(z) fData(z)[n,'explCnts'])))
    sel <- !(keep %in% c('transcript','gene_id','island_id'))
    fData(esets[[1]]) <- fData(esets[[1]])[,keep]
    fvarLabels(esets[[1]])[sel] <- paste(sampleNames[1],fvarLabels(esets[[1]])[sel],sep='.')
    sampleNames(esets[[1]]) <- sampleNames[1]
    for (i in 2:length(esets)) {
      if (any(!(keep %in% fvarLabels(esets[[i]])))) stop("Some variables in argument 'keep' were not found")
      fData(esets[[i]]) <- fData(esets[[i]])[,keep]
      sel <- !(keep %in% c('transcript','gene_id','island_id'))
      fvarLabels(esets[[i]])[sel] <- paste(sampleNames[i],fvarLabels(esets[[i]])[sel],sep='.')
      esets[[i]] <- esets[[i]][n,]
      sampleNames(esets[[i]]) <- sampleNames[i]
    }
    ans <- do.call("combine",esets)
  } else {
    readCount <- fData(esets[[1]])[,'explCnts']
    ans <- esets[[1]]
  }
  fData(ans)$readCount <- readCount
  return(ans)
}
