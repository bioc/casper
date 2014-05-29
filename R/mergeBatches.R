setMethod("mergeBatches", signature(x='ExpressionSet',y='simulatedSamples'), function(x, y, mc.cores=1) {
  if (mc.cores>1) {
    if ('parallel' %in% loadedNamespaces()) {
      ans <- parallel::mclapply(y, function(z) mergeBatches(x,z$simExpr), mc.cores=mc.cores)
    } else {
      stop("parallel package has not been loaded!")
    }
  } else {
    ans <- lapply(y, function(z) mergeBatches(x,z$simExpr))
  }
  return(ans)
}
)

setMethod("mergeBatches", signature(x='ExpressionSet',y='ExpressionSet'), function(x, y, mc.cores=1) {
  #Merge ExpressionSets x and y, apply quantile norm and adjust for batch effects via partial residuals
  n <- names(fData(x))[names(fData(x)) %in% names(fData(y))]
  nx <- names(fData(x))
  nx[nx %in% n] <- paste(nx[nx %in% n],'x',sep='.')
  names(fData(x)) <- nx
  ny <- names(fData(y))
  ny[ny %in% n] <- paste(ny[ny %in% n],'y',sep='.')
  names(fData(y)) <- ny
  xnew <- combine(x,y)
  xnew$batch <- factor(rep(c('batch1','batch2'),c(ncol(x),ncol(y))))
  xnewnorm <- quantileNorm(xnew)
  xnewadj <- anovaAdjustment(xnewnorm, adjustmentVariable='batch')
  if (all(c('readCount.x','readCount.y') %in% fvarLabels(xnewadj))) fData(xnewadj)$readCount <- rowSums(fData(xnewadj)[,c('readCount.x','readCount.y')])
  cat('.')
  return(xnewadj)
}
)

### AUXILIZARY FUNCTIONS: QUANTILE NORMALIZATION & ANOVA ADJUSTMENT
setMethod("quantileNorm", signature(x="ExpressionSet"), function(x) {
    exprs(x) <- quantileNorm(exprs(x))
    return(x)
}
)

setMethod("quantileNorm", signature(x="matrix"), function(x) {
  probsx <- seq(0,1,length=nrow(x))
  qx <- quantile(x,probs=probsx,na.rm=TRUE)
  f <- approxfun(probsx,qx)
  for (i in 1:ncol(x)) {
    probsy <- rank(x[,i])/(nrow(x)+1)
    x[,i] <- f(probsy)
  }
  return(x)
}
)


anovaAdjustment <- function(x,adjustmentVariable,na.rm=FALSE) {
  #Perform ANOVA adjustment, i.e. remove effect of nuisance covariates via lmFit. This can be used to remove experimentation date or batch processing effects on gene expression.
  # - x: ExpressionSet with expression values to be adjusted
  # - adjustmentVariable: name of adjustment variable, which must be contained in pData(x).
  # Output: ExpressionSet with adjusted expression values
  if (missing(adjustmentVariable)) stop('adjustmenVariable must be specified')
  UseMethod("anovaAdjustment")
}


anovaAdjustment.ExpressionSet <- function(x,adjustmentVariable,na.rm=FALSE) {
  if (!(adjustmentVariable %in% names(pData(x)))) { stop('Specified variable not found in pData') }
  design <- model.matrix(~ -1 + as.factor(pData(x)[,adjustmentVariable]))
  beta <- exprs(x) %*% design %*% diag(1/colSums(design))
  linpred <- beta %*% t(design)
  exprs(x) <- exprs(x) - linpred + rowMeans(exprs(x),na.rm=na.rm)
  return(x)
}
