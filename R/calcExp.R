procExp <- function(distrs, genomeDB, pc, readLength, islandid, rpkm=TRUE, priorq=2, priorqGeneExpr=2, citype='none', niter=10^3, burnin=100, mc.cores=1, verbose=FALSE, totReads) {
  #Format input
  startcdf <- distrs@stDis(seq(0,1,.001))
  lendis <- as.double(distrs@lenDis/sum(distrs@lenDis))
  lenvals <- as.integer(names(distrs@lenDis))
  readLength <- as.integer(readLength)
  priorq <- as.double(priorq)
  citype <- as.integer(switch(citype, none=0, asymp=1, exact=2))
  if (length(citype)==0) stop("citype should take the value 'none', 'asymp', or 'exact'")
  niter <- as.integer(niter)
  burnin <- as.integer(burnin)
  verbose <- as.integer(verbose)
  if ((niter<=burnin) & citype==2) stop("Too many burnin iterations specified. Decrease burnin or increase niter")
  #if (missing(islandid)) islandid <- names(genomeDB@islands)[elementLengths(genomeDB@islands)>1]
  if (missing(islandid)) islandid <- names(genomeDB@islands)

  exons <- as.integer(names(genomeDB@islands@unlistData))
  names(exons) <- rep(names(genomeDB@islands), elementLengths(genomeDB@islands))
  exons <- split(unname(exons), names(exons))
  exonwidth <- width(genomeDB@islands@unlistData)
  names(exonwidth) <- rep(names(genomeDB@islands), elementLengths(genomeDB@islands))
  exonwidth <- split(unname(exonwidth), names(exonwidth))
  strand <- as.character(strand(genomeDB@islands@unlistData))[cumsum(c(1, elementLengths(genomeDB@islands)[-length(genomeDB@islands)]))]
  names(strand) <- names(genomeDB@islands)
  
  if (!all(islandid %in% names(exons))) stop('islandid not found in genomeDB@islands')
  if (!all(islandid %in% names(pc))) stop('islandid not found in pc')
  if (!all(islandid %in% names(genomeDB@transcripts))) stop('islandid not found in genomeDB@transcripts')
  
      #Define basic function
  
  f <- function(z) {
    islandid <- as.integer(z)
    exons <- exons[z]
    exonwidth <- exonwidth[z]
    transcripts <- genomeDB@transcripts[z]
    tmp <- strand[z]
    strand <- vector(mode='integer', length=length(tmp))
    sel <- tmp=='+'
    strand[sel] <- 1
    sel <- tmp=='-'
    strand[sel] <- -1
    sel <- tmp=='*'
    strand[sel] <- 0
    strand <- as.list(as.integer(strand))
    pc <- pc[z]
    ans <- calcKnownMultiple(exons=exons,exonwidth=exonwidth,transcripts=transcripts,islandid=as.list(islandid),pc=pc,startcdf=startcdf,lendis=lendis,lenvals=lenvals,readLength=readLength,priorq=priorq, strand=strand, citype=citype, niter=niter, burnin=burnin, verbose=verbose)
    if (citype==0) {
      ans <- lapply(ans, function(z) { res=vector("list",2); res[[1]]= z[[1]]; names(res[[1]])= z[[2]]; res[[2]]=z[[5]]; res })
    } else if (citype==1) {
      ans <- lapply(ans, function(z) { res=vector("list",3); res[[1]]= z[[1]]; res[[2]]= z[[3]]; names(res[[1]])= names(res[[2]])= z[[2]]; res[[3]] = z[[5]]; res })
    } else if (citype==2) {
      ans <- lapply(ans, function(z) { res=vector("list",4); res[[1]]= z[[1]]; res[[2]]= z[[3]]; res[[3]]= matrix(z[[4]],nrow=niter-burnin); names(res[[1]])= names(res[[2]])= z[[2]]; res[[4]] = z[[5]]; res })
    }
    ans
  }

  #Run
  sel <- !sapply(pc[islandid], is.null)
  all <- islandid
  islandid <- islandid[sel]
  sel <- elementLengths(genomeDB@islands[islandid])>1
  islandid <- islandid[sel]
  if(length(islandid)>0){
    if (verbose) cat("Obtaining expression estimates...\n")
    if (mc.cores>1 && length(islandid)>mc.cores) {
      if ('parallel' %in% loadedNamespaces()) {
                    #split into smaller jobs
        islandid <- split(islandid, cut(1:length(islandid), mc.cores))
        ans <- parallel::mclapply(islandid,f,mc.cores=mc.cores, mc.preschedule=FALSE)
        ans <- unlist(ans, recursive=F)
        names(ans) <- unlist(islandid)
      } else stop('parallel library has not been loaded!')
    } else {
      ans <- f(islandid)
      names(ans) <- islandid
    }
  } else ans <- NULL

  nreads <- unlist(lapply(pc[all], sum))
  nreads[names(ans)] <- unlist(sapply(ans, '[[', length(ans[[1]])))  #count only reads with positive probability under provided variants
  
  sel <- sapply(ans, function(z) sum(z[[1]])==0) #genes with a single variant
  ans[sel] <- NULL

  miss <- lapply(genomeDB@transcripts[all[!(all %in% names(ans))]], names)  #careful: some genes repeated in miss though they appear only once in transcripts, e.g. '20306' sum(names(miss)=='20306')

  #Format as ExpressionSet

  if (verbose) cat("Formatting output...\n")
  transcript <- unlist(lapply(ans, function(z) names(z[[1]])))
  gene <- rep(names(ans),sapply(ans,function(z) length(z[[1]])))
  fdata <- data.frame(transcript=transcript, gene=gene)
  if(length(miss)>0){
    transcript <- c(transcript, unname(unlist(miss)))
    gene <- c(gene, rep(names(miss), sapply(miss, length)))
    fdata <- data.frame(transcript=transcript, gene=gene)
    ntx.miss <- sapply(miss,length) 
    misse <- rep(1/ntx.miss,ntx.miss) 
    a0 <- rep(priorq*ntx.miss,ntx.miss) 
    missSE <- sqrt(priorq * (a0-priorq) / (a0^2 * (a0+1))) 
    names(misse) <- unlist(miss)
    exprsx <- matrix(c(unlist(lapply(ans,'[[',1)), misse) ,ncol=1)
    rownames(exprsx) <- c(unlist(lapply(ans, function(z) names(z[[1]]))), names(misse))
  } else {
    exprsx <- matrix(unlist(lapply(ans, '[[', 1)), ncol=1)
    rownames(exprsx) <- unlist(lapply(ans, function(z) names(z[[1]])))
  }
    exprsx <- round(exprsx,10)
    if (citype==1) {
      se <- unlist(lapply(ans,'[[',2))
      if(length(miss)>0) se <- c(se, missSE)
      if (!rpkm) {
        alpha <- .05
        fdata$ci95.low <- qnorm(alpha/2,mean=exprsx,sd=se)
        fdata$ci95.high <- qnorm(1-alpha/2,mean=exprsx,sd=se)
        sel <- round(fdata$ci95.low,10)<0 & !is.na(se); fdata$ci95.low[sel] <- 0
        se[sel][se[sel]==0] <- 1e-20
        fdata$ci95.high[sel] <- unlist(mapply(function(m,s) casper:::qtnorm(1-alpha/2,mean=m,sd=s,lower=0,upper=1), as.list(exprsx[sel,1]), as.list(se[sel])))
        sel <- round(fdata$ci95.high,10)>1 & !is.na(se); fdata$ci95.high[sel] <- 1
        fdata$ci95.low[sel] <- unlist(mapply(function(m,s) casper:::qtnorm(alpha/2,mean=m,sd=s,lower=0,upper=1), as.list(exprsx[sel,1]), as.list(se[sel])))
      }
      fdata$se <- se
    }
    if (citype==2) {
      if(sum(unlist(lapply(ans, function(x) is.na(x[[3]]))))==0){
        if (mc.cores>1) {
          if ('parallel' %in% loadedNamespaces()) {
           se <- unlist(parallel::mclapply(ans, function(z) (colMeans(z[[3]]^2) - colMeans(z[[3]])^2) * nrow(z[[3]])/(nrow(z[[3]]) - 1), mc.cores=mc.cores))
         } else stop('parallel library has not been loaded!')
        } else se <- unlist(lapply(ans, function(z) (colMeans(z[[3]]^2) - colMeans(z[[3]])^2) * nrow(z[[3]])/(nrow(z[[3]]) - 1)))
        
        if(length(miss)>0) se <- c(se, missSE)
        fdata$se <- se
        if (!rpkm) {
          q <- lapply(ans,function(z) apply(z[[3]],2,quantile,probs=c(.025,.975)))
          q <- t(do.call(cbind,q))
          if(length(miss)>0) q <- rbind(q, matrix(NA, nrow=length(misse), ncol=2))
          fdata$ci95.low <- q[,1]; fdata$ci95.high <- q[,2]
          tmp <- cbind(exprsx, fdata[,c('ci95.low','ci95.high')])
          tmp <- as.data.frame(t(apply(tmp, 1, function(x){
            if(!any(is.na(x))){
              y <- x
              if(x[2]>x[1]) y[2]=y[1]
              if(x[3]<x[1]) y[3]=y[1]
            } else y <- x
            y
            })))
          fdata$ci95.low <- tmp$ci95.low
          fdata$ci95.high <- tmp$ci95.high
        }
      }
    }
  rownames(exprsx) <- rownames(fdata) <- fdata$transcript

  fdata <- cbind(fdata, nreads[as.character(fdata$gene)])
  colnames(fdata)[ncol(fdata)] <- "explCnts"

  #Compute log(RPKM)
  if (rpkm) {
    if (citype != 0) se.logpi <- fdata$se / exprsx
    nreads <- nreads[unique(as.character(fdata$gene))]
    totReads <- sum(nreads) + priorqGeneExpr*length(nreads) 
    txLength <- txLength(genomeDB=genomeDB)
    apost <- nreads + (priorqGeneExpr-1)
    thest <- apost/sum(apost)
    theta <- data.frame(gene=names(nreads), thest=thest)
    exprsx <- as.matrix(logrpkm(fdata, th=theta, pi=exprsx, genomeDB=genomeDB))
    colnames(exprsx) <- NULL
    if (citype != 0) {
      se.theta <- sqrt(apost * (1-apost/totReads) / (totReads * (totReads+1))) 
      se.logtheta <- se.theta / thest
      fdata$se <- sqrt(se.logpi^2 + se.logtheta[as.character(fdata$gene)]^2)
      alpha <- 0.05
      err <- qnorm(alpha/2)*fdata$se; fdata$ci95.low <- exprsx + err; fdata$ci95.high <- exprsx - err #Added
    }
  }
  #Return ExpressionSet

  fdata <- new("AnnotatedDataFrame",fdata)
  ans <- new("ExpressionSet",exprs=exprsx,featureData=fdata)
  ans
}

mergeExp <- function(minus, plus){
  exp <- rbind(exprs(plus), exprs(minus))
  fdata <- new("AnnotatedDataFrame", rbind(fData(plus), fData(minus)))
  ans <- new("ExpressionSet",exprs=exp,featureData=fdata)
  ans
}

calcExp <- function(distrs, genomeDB, pc, readLength, islandid, rpkm=TRUE, priorq=2, priorqGeneExpr=2, citype='none', niter=10^3, burnin=100, mc.cores=1, verbose=FALSE) {
  totReads <- getNreads(pc); totReads <- sum(totReads) + priorqGeneExpr*length(totReads) #Modified: totReads includes prior sample size
  
  if (missing(readLength)) stop("readLength must be specified")
  if (genomeDB@denovo) stop("genomeDB must be a known genome")
  if (pc@denovo) stop("pc must be a pathCounts object from known genome")
  if(pc@stranded){
    if(missing(islandid)) islandid <- c(names(pc@counts$plus), names(pc@counts$minus))
    is <- as.character(strand(genomeDB@islands@unlistData))[cumsum(c(1, elementLengths(genomeDB@islands)[-length(genomeDB@islands)]))]
    names(is) <- names(genomeDB@islands)
    plusGI <- islandid[is[islandid]=="+"]
    plus <- NULL
    if(length(plusGI)>0) {
      plusDB <- genomeBystrand(genomeDB, "+")
      plus <- procExp(distrs, plusDB, pc=pc@counts$plus, readLength=readLength, islandid=plusGI, rpkm=rpkm, priorq=priorq, priorqGeneExpr=priorqGeneExpr, niter=niter, burnin=burnin, mc.cores=mc.cores, citype=citype, totReads=totReads)
    }
    minusGI <- islandid[is[islandid]=="-"]
    minus <- NULL
    if(length(minusGI)>0) {
      minusDB <- genomeBystrand(genomeDB, "-")
      minus <- procExp(distrs, minusDB, pc=pc@counts$minus, readLength=readLength, islandid=minusGI, rpkm=rpkm, priorq=priorq, priorqGeneExpr=priorqGeneExpr, niter=niter, burnin=burnin, mc.cores=mc.cores, citype=citype, totReads=totReads)
    } 
    if(is.null(plus) & is.null(minus)) stop("No counts in islandid genes")
    if(!(is.null(plus) | is.null(minus))) { ans <- mergeExp(plus, minus) }
    else { if(is.null(plus)) {ans <- minus } else ans <- plus }
  } else {
    if(missing(islandid)) islandid <- names(pc@counts[[1]])
    if(sum(!unlist(lapply(pc@counts[islandid], is.null)))) stop("No counts in islandid genes")
    ans <- procExp(distrs, genomeDB, pc=pc@counts[[1]], readLength=readLength, islandid=islandid, rpkm=rpkm, priorq=priorq, priorqGeneExpr=priorqGeneExpr, niter=niter, burnin=burnin, mc.cores=mc.cores, citype=citype, totReads=totReads)
  }

  ans
}
    
calcKnownMultiple <- function(exons, exonwidth, transcripts, islandid, pc, startcdf, lendis, lenvals, readLength, priorq, strand, citype, niter, burnin, verbose) {
  ans <- .Call("calcKnownMultiple",exons,exonwidth,transcripts,islandid,pc,startcdf, lendis, lenvals, readLength, priorq, strand, citype, niter, burnin, verbose, PACKAGE="casper")
  return(ans)
}


lhoodGrid <- function(pc, distrs, genomeDB, readLength, islandid, grid, priorq=2) {
  if (length(islandid)>1) stop("Only 1 gene is allowed, please specify a single gene in argument 'islandid'")
  if (missing(readLength)) stop("readLength must be specified")
  if (genomeDB@denovo) stop("genomeDB must be a known genome")
  if (pc@denovo) stop("pc must be a pathCounts object from known genome")

  #Format input
  startcdf <- distrs@stDis(seq(0,1,.001))
  lendis <- as.double(distrs@lenDis/sum(distrs@lenDis))
  lenvals <- as.integer(names(distrs@lenDis))
  readLength <- as.integer(readLength)
  priorq <- as.double(priorq)

  exons <- as.integer(names(genomeDB@islands[[islandid]]))
  exonwidth <- width(genomeDB@islands[[islandid]])
  strand <- as.character(strand(genomeDB@islands@unlistData))[cumsum(c(1, elementLengths(genomeDB@islands)[-length(genomeDB@islands)]))]
  names(strand) <- names(genomeDB@islands)
  strand   <- strand[islandid]
  transcripts <- genomeDB@transcripts[[islandid]]
  if (length(transcripts)==1) stop("Single transcript specified, estimation not run")
  strand <- as.integer(ifelse(strand=='+', 1, -1))
  if (pc@stranded) {
     if (strand==1) pc <- pc@counts[['plus']][[islandid]] else pc <- pc@counts[['minus']][[islandid]]
  } else {
     pc <- pc@counts[[1]][[islandid]]
  }
  islandid <- as.integer(islandid)
  
  if (missing(grid)) {
    grid <- lapply(1:(length(transcripts)-1), function(z) seq(-10,10,length=10^(5/(length(transcripts)-1))))
    #grid <- lapply(1:length(transcripts), function(z) seq(0.0001,.9999,length=10^(5/length(transcripts))))
  } else {
    if (!is.list(grid)) stop("grid must be of class list")
    if (length(grid) != length(transcripts)-1) stop("grid must have one row less than the number of known transcripts")
  }
  milogit <- function(theta) cbind(rep(1,nrow(theta)),exp(theta))/(1+rowSums(exp(theta))) #matrix conversion
  gridmat <- t(milogit(expand.grid(grid)))
  
  ans <- .Call("lhoodGrid",gridmat,exons,exonwidth,transcripts,pc,startcdf,lendis,lenvals,readLength,priorq,strand)
  names(ans[[2]]) <- ans[[3]]; ans[[3]] <- NULL
  names(ans) <- c('logpos','piest','logpos.piest','S','pathprob','marginalLhood')
  rownames(ans$pathprob) <- names(ans$piest)
  ans$pathprob <- t(ans$pathprob)
  names(ans$marginalLhood) <- c('laplace','IS')
  ans$grid <- grid
  return(ans)  
}


