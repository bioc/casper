
variantMargExpr <- function(x,minProbExpr=0.5, minExpr=0.05) {
  #Marginal expression for each variant (obtained via model averaging) and marginal post prob of being expressed
  # - minProbExpr: variants with marginal post prob < minProbExpr are not reported
  # - minExpr: variants with expression < minExpr are not reported
  # Note: at least one variant is always reported. If no variants satisfy minProbExpr and minExpr, the variant with largest expression is reported
  if(!(all(is.na(x@posprob$posprob))))
  {
      pospr <- x@posprob$posprob/sum(x@posprob$posprob)
      names(pospr) <- x@posprob$model
      pospr <- pospr[as.character(x@expression$model)]
      ans <- by(data.frame(pospr*x@expression$expr,pospr),INDICES=list(var=x@expression$varName),FUN=colSums,simplify=FALSE)
      n <- names(ans)
      ans <- matrix(unlist(ans),ncol=2,byrow=TRUE)
      colnames(ans) <- c('expr','probExpressed')
      rownames(ans) <- n
      sel <- ans[,'probExpressed']>minProbExpr & ans[,'expr']>minExpr
      if (any(sel)) ans <- ans[sel,,drop=FALSE] else ans <- ans[which.max(ans[,'expr']),,drop=FALSE]
      ans[,'expr'] <- ans[,'expr']/sum(ans[,'expr'])
  }
  else
  {
      ans <- NULL
  }
  return(ans)
}

relativeExpr <- function(expr, summarize='modelAvg', minProbExpr=0.5, minExpr=0.05){
  if (class(expr)!='denovoGenomeExpr') stop("expr must be of class 'denovoGenomeExpr'")
  if (!(summarize %in% c("bestModel", "modelAvg"))) stop("summarize must be one of 'bestModel' or 'modelAvg'")
  if (summarize=='bestModel'){
    ans <- lapply(as.list(expr), function(x){
      
      if(!(all(is.na(x@posprob$posprob))))
      {
        pospr <- x@posprob$posprob
        names(pospr) <- x@posprob$model
        
        bestModel <- x@posprob$model[which.max(x@posprob$posprob)]
        pospr <- pospr[as.character(bestModel)]
        
        expr2var <- x@expression[x@expression$model==bestModel,]
        
        ans <- matrix(nrow=length(expr2var$varName),ncol=2,byrow=TRUE)
        colnames(ans) <- c('expr','probExpressed')
        rownames(ans) <- expr2var$varName
        ans[,"expr"] <- expr2var$expr
        ans[,"probExpressed"] <- rep(pospr, length(expr2var$varName))
        return(ans)    
      }
      else
      {
        ans <- NULL
      }
    })
    ans <- do.call("rbind", unname(ans))
  } else {
    ans <- lapply(as.list(expr), variantMargExpr, minProbExpr=minProbExpr, minExpr=minExpr)
    ans <- do.call("rbind", unname(ans))
  }
  ans
}




denovoExpr <- function(x, pc, rpkm=TRUE, summarize='modelAvg', minProbExpr=0.5, minExpr=0.05) {
  sel <- sapply(as.list(x), function(z) !is.na(posprob(z)[1,'posprob']))
  if (any(sel)) {
      pis <- relativeExpr(x[sel], summarize=summarize, minProbExpr=minProbExpr, minExpr=minExpr)
  } else {
      pis <- matrix(nrow=0,ncol=2); colnames(pis) <- c('expr','probExpressed')
  }
  if (any(!sel)) {
    pis.noreads <- do.call(rbind,lapply(as.list(x[!sel]), function(z) z@expression))
    rownames(pis.noreads) <- pis.noreads[,'varName']
    pis.noreads <- cbind(pis.noreads[,c('expr'),drop=FALSE],rep(NA,nrow(pis.noreads)))
    colnames(pis.noreads) <- colnames(pis)
    pis <- rbind(pis, pis.noreads)
  }
  fdata <- variants(x)
  names(fdata)[1:2] <- c('island_id','transcript')
  rownames(fdata) <- as.character(fdata$transcript)
  if (pc@stranded) {       
    nreads <- c(unlist(lapply(pc@counts$plus,sum)), unlist(lapply(pc@counts$minus,sum)))
    nreads <- nreads[names(x)]
  } else {
    nreads <- unlist(lapply(pc@counts[[1]][names(x)], sum))
  }
  if (rpkm) {
    len <- x@txLength
    apost <- nreads + (x@priorq-1)
    thest <- apost/sum(apost)
    theta <- data.frame(island_id=names(nreads), thest=thest)
    #
    exprsx <- as.matrix(logrpkm(fdata[,c('island_id','transcript')], th=theta, pi=pis[,'expr',drop=FALSE], len=len))
    colnames(exprsx) <- NULL
    fdata <- fdata[rownames(exprsx),]
    ans <- new("ExpressionSet", exprs=exprsx, featureData=new("AnnotatedDataFrame",fdata))
  } else {
    exprsx <- pis[,1,drop=FALSE]
    fdata <- fdata[rownames(exprsx),]
    ans <- new("ExpressionSet", exprs=exprsx, featureData=new("AnnotatedDataFrame",fdata))
  }
  fData(ans)$explCnts <- nreads[as.character(fData(ans)$island_id)]
  fData(ans)$probExpressed <- pis[featureNames(ans),'probExpressed']
  return(ans)  
}
