relexprByGene <- function(x) {
  #Input
  # - x: ExpressionSet with gene_id in fData(x) containing relative expressions (typically, adding up to 1 for each island_id
  #Output
  # ExpressionSet where expression adds up to 1 per gene_id
  if (!is(x,'ExpressionSet')) stop('x must be of class ExpressionSet')
  if (any((exprs(x)>1) | (exprs(x)<0))) stop('Input x should contain relative expressions. Values outside [0,1] detected')
  sumpi <- aggregate(exprs(x), by=list(fData(x)$gene_id), FUN=sum)
  names(sumpi)[1] <- 'gene_id'
  exprsx <- data.frame(fData(x)[,c('transcript','gene_id')],exprs(x))
  exprsx <- merge(exprsx, sumpi, by='gene_id')
  pi <- exprsx[,3:(2+ncol(x))]/exprsx[,(3+ncol(x)):(2+2*ncol(x)),drop=FALSE]
  rownames(pi) <- as.character(exprsx$transcript)
  pi <- pi[featureNames(x),,drop=FALSE]
  names(pi) <- sampleNames(x)
  exprs(x) <- as.matrix(pi)
  return(x)
}
