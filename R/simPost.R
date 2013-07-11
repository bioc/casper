mergePCs <- function(pcs, genomeDB, mc.cores=1){
    allisl <- unique(unlist(lapply(pcs, function(x) names(x@counts[[1]]))))
      alltxs <- lapply(pcs, function(x) unlist(unname(x@counts[[1]])))
      alltxsu <- unlist(alltxs, recursive=F)
      alltxsn <- unique(names(alltxsu))
      ans <- matrix(ncol=length(pcs), nrow=length(alltxsn))
      rownames(ans) <- alltxsn
      for(i in 1:length(alltxs)) ans[rownames(ans) %in% names(alltxs[[i]]),i] <- alltxs[[i]]
      ans <- rowSums(ans, na.rm=T)
    mode(ans) <- "integer"
    #nans <- names(ans)
    #ans <- as.integer(ans)
    #names(ans) <- nans
      counts <- casper:::splitPaths(paths=ans, DB=genomeDB, mc.cores=mc.cores, stranded=pcs[[1]]@stranded, geneid=allisl)
      new("pathCounts", counts=counts, denovo=pcs[[1]]@denovo, stranded=pcs[[1]]@stranded)
  }

simMAE <- function(nsim, islandid=NULL, n, r, f, burnin=1000, pc, distrs, usePilot=FALSE, retTxsError=FALSE, genomeDB, mc.cores=1, mc.cores.int=1, verbose=FALSE) {
   if (length(r) != length(n)) stop("length(n) not equal to length(r)")
   if(is.null(islandid)) islandid <- names(genomeDB@transcripts)
   U <- NULL
   txe <- list()
   pcs <- list()
   for (j in 1:length(n)) {
     d <- distrs
     nmean <-  mean(rep(as.numeric(names(d@lenDis)), d@lenDis))
     nmean <- f[j] - nmean
     names(d@lenDis) <- as.numeric(names(d@lenDis)) + round(nmean)
     if(verbose) cat(paste("Generating posterior samples j =",j, "\n"))
      pis <- simPost(islandid=islandid, nsim=nsim, distrs=d, genomeDB=genomeDB, pc=pc, readLength=r[j], mc.cores=mc.cores.int*mc.cores, verbose=verbose)
     if(verbose) cat(paste("Running simulations for j =",j, "\n"))
     if(mc.cores.int>1) {
       require(parallel)
       res <- parallel::mclapply(1:nsim, function(i){
         sim.pc <-  simPostPred(islandid=islandid, nreads=n[j], pis=pis[i,], pc=pc, distrs=d, rl=r[j], genomeDB=genomeDB, verbose=verbose)
         if(usePilot) sim.pc$pc <- mergePCs(pcs=list(sim.pc$pc,pc), genomeDB=genomeDB)
         sim.exp <- exprs(calcExp(islandid=islandid, distrs=d, genomeDB=genomeDB, pc=sim.pc$pc, readLength=r[j], rpkm=FALSE, mc.cores=mc.cores))
         list(maes=abs(sim.exp[colnames(pis),]-pis[i,]), pc=sim.pc$pc)
       }, mc.cores=mc.cores.int, mc.preschedule=FALSE)
     } else {
       res <- lapply(1:nsim, function(i){
         sim.pc <-  simPostPred(islandid=islandid, nreads=n[j], pis=pis[i,], pc=pc, distrs=d, rl=r[j], genomeDB=genomeDB, verbose=verbose)
         if(usePilot) sim.pc$pc <- mergePCs(pcs=list(sim.pc$pc,pc), genomeDB=genomeDB)
         sim.exp <- exprs(calcExp(islandid=islandid, distrs=d, genomeDB=genomeDB, pc=sim.pc$pc, readLength=r[j], rpkm=FALSE))
         list(maes=abs(sim.exp[colnames(pis),]-pis[i,]), pc=sim.pc$pc)
       })
     }
     Un <- do.call(cbind, lapply(res, "[[", "maes"))
     df <- data.frame(MAE= colMeans(Un), Nreads=rep(n[j], nsim), ReadLength=rep(r[j], nsim), frLength=rep(f[j], nsim))
     U <- rbind(U, df)
     if(retTxsError) {
       rownames(Un) <- colnames(pis)
       txe[[j]] <- Un
       pcs[[j]] <- lapply(res, "[[", "pc")
     }
   }
   if(retTxsError){
     names(txe) <- paste(n, r, f, sep='-')
     names(pcs) <- paste(n, r, f, sep='-')
     return(list(txe=txe, U=U, pcs=pcs))
   }
   return(U)
}


procsimPost <- function(nsim, distrs, genomeDB, pc, readLength, islandid, initvals, useinit, relativeExpr=TRUE, priorq=2, burnin=1000, mc.cores=1, verbose=FALSE) {
  startcdf <- distrs@stDis(seq(0,1,.001))
  lendis <- as.double(distrs@lenDis/sum(distrs@lenDis))
  lenvals <- as.integer(names(distrs@lenDis))
  readLength <- as.integer(readLength)
  priorq <- as.double(priorq)
  citype <- as.integer(2)
  niter <- as.integer(nsim+burnin)
  burnin <- as.integer(burnin)
  verbose <- as.integer(verbose)
  useinit <- as.integer(useinit)
  if (niter<=burnin) stop("Too many burnin iterations specified. Decrease burnin or increase niter")
  if (missing(islandid)) islandid <- names(genomeDB@islands)[sapply(genomeDB@islands,length)>1]

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
    ans <- casper:::calcKnownMultiple(exons=exons,exonwidth=exonwidth,transcripts=transcripts,islandid=as.list(islandid),pc=pc,startcdf=startcdf,lendis=lendis,lenvals=lenvals,readLength=readLength,priorq=priorq, strand=strand, citype=citype, niter=niter, burnin=burnin, verbose=verbose)
    
    if(length(ans)==1) {
      trans <-  ans[[1]][[2]]
      ntrans <- length(trans)
    }
    else {
      trans <- sapply(ans, function(z) z[[2]])
      ntrans <- unlist(lapply(trans, length))
      trans <- unlist(trans)
    }
    l <- unlist(lapply(ans, function(z) length(z[[4]])))
    def <- which(l!=nsim*ntrans)
    for(i in def) ans[[i]][[4]] <- rep(0, length(ans[[i]][[4]])*nsim)
    ans <- lapply(ans, function(z) { res = matrix(z[[4]],nrow=nsim); res })
    ans <- do.call(cbind,ans)
    colnames(ans) <- trans
    ans
}

    #Run
    sel <- !sapply(pc[islandid], is.null)
    all <- islandid
    islandid <- islandid[sel]
  if (verbose) cat("Obtaining expression estimates...\n")
    if (mc.cores>1 && length(islandid)>mc.cores) {
      if ('parallel' %in% loadedNamespaces()) {
                    #split into smaller jobs
        islandid <- split(islandid, cut(1:length(islandid), mc.cores))
        ans <- parallel::mclapply(islandid,f,mc.cores=mc.cores)
        ans <- do.call(cbind,ans)
      } else stop('parallel library has not been loaded!')
    } else {
      ans <- f(islandid)
    }
    if (verbose) cat("Formatting output...\n")
  ans
}

simPost <- function(nsim, distrs, genomeDB, pc, readLength, islandid, initvals, relativeExpr=TRUE, priorq=2, burnin=1000, mc.cores=1, verbose=FALSE) {
  if (missing(initvals)) {
    useinit <- 0
    initvals <- NULL
  }
  else useinit <-  1
  if (missing(readLength)) stop("readLength must be specified")
  if (genomeDB@denovo) stop("genomeDB must be a known genome")
  if (pc@denovo) stop("pc must be a pathCounts object from known genome")
  if(pc@stranded){
    if(missing(islandid)) islandid <- c(names(pc@counts$plus), names(pc@counts$minus))
    plusGI <- islandid[genomeDB@islandStrand[islandid]=="+"]
    plus <- NULL
    if(length(plusGI)>0) {
      plusDB <- genomeBystrand(genomeDB, "+")
      plus <- procsimPost(nsim, distrs, plusDB, pc=pc@counts$plus, readLength=readLength, islandid=plusGI, initvals=initvals, useinit=useinit, relativeExpr=relativeExpr, priorq=priorq, burnin=burnin, mc.cores=mc.cores, verbose=verbose)
    }
    minusGI <- islandid[genomeDB@islandStrand[islandid]=="-"]
    minus <- NULL
    if(length(minusGI)>0) {
      minusDB <- genomeBystrand(genomeDB, "-")
      minus <- procsimPost(nsim, distrs, minusDB, pc=pc@counts$minus, readLength=readLength, islandid=minusGI, initvals=initvals, useinit=useinit, relativeExpr=relativeExpr, priorq=priorq, burnin=burnin, mc.cores=mc.cores, verbose=verbose)
    }
    if(is.null(plus) & is.null(minus)) stop("No counts in islandid genes")
    if(!(is.null(plus) | is.null(minus))) { ans <- mergeExp(plus, minus) }
    else { if(is.null(plus)) {ans <- minus } else ans <- plus }
}
else {
    if(missing(islandid)) islandid <- names(pc@counts[[1]])
    if(sum(!unlist(lapply(pc@counts[islandid], is.null)))) stop("No counts in islandid genes")
    ans <- procsimPost(nsim, distrs, genomeDB, pc=pc@counts[[1]], readLength=readLength, islandid=islandid, initvals=initvals, useinit=useinit, relativeExpr=relativeExpr, priorq=priorq, burnin=burnin, mc.cores=mc.cores, verbose=verbose)
}
  ans
}

#calcKnownMultiple <- function(exons, exonwidth, transcripts, islandid, pc, startcdf, lendis, lenvals, readLength, initvals, useinit, priorq, strand, citype, niter, burnin, verbose) {
#  ans <- .Call("calcKnownMultiple",exons,exonwidth,transcripts,islandid,pc,startcdf, lendis, lenvals, readLength, initvals, useinit, priorq, strand, citype, niter, burnin, verbose)
#  return(ans)
#}

simPostPred <- function(nreads, islandid=NULL, pis, pc, distrs, rl, genomeDB, seed=1, mc.cores=1, verbose=FALSE) {
  a <- tapply(pis, getIsland(txid=names(pis), genomeDB=genomeDB), sum)
  def <- names(a)[a==0]
  pc@counts[[1]][def] <- NULL
  tdef <- unlist(lapply(genomeDB@transcripts[def],names))
  pis <- pis[!(names(pis)%in%tdef)]
  geneExpr <- NULL
  nreadsPerGene <- getNreads(pc)
  if(!is.null(islandid)) nreadsPerGene <- nreadsPerGene[names(nreadsPerGene) %in% islandid]
  thetas <- as.vector(rdirichlet(1, nreadsPerGene+1))
  names(thetas) <- names(nreadsPerGene)
  nreadsPerGeneSim <- rmultinom(1, nreads, thetas)
  nreadsPerGeneSim <- nreadsPerGeneSim[nreadsPerGeneSim>0,]
  nonzero <- names(nreadsPerGeneSim)
  txs <- lapply(genomeDB@transcripts[nonzero], names)
  ntxs <- lapply(txs, length)
  isl <- rep(names(txs), ntxs)
  txs <- unname(unlist(txs))
  names(txs) <- isl
  if (!all(txs %in% names(pis))) {
    miss <-  subset(txs, !txs %in% names(pis))
    tab <- table(names(miss))
    vars <- as.numeric(tab)
    names(vars) <- names(tab)
    one <- subset(vars, vars == 1)
    names(one) <- miss[names(miss) %in% names(one)]
    pis <- append(one, pis)
    nonone <- subset(vars, vars > 1)
    pi.miss <- NULL
    if(length(nonone)>0){
      pi.miss <- lapply(1:length(nonone), function(i){
        pi <- as.numeric(rdirichlet(1, rep(2, nonone[i])))
        names(pi) <- miss[names(miss) %in% names(nonone[i])]
        pi
      })
    pis <- append(unlist(pi.miss), pis)
    }
  }
  pc <- simReads(nonzero, nSimReads=nreadsPerGeneSim, pis=pis, rl=rl, seed=seed, distrs=distrs, genomeDB=genomeDB, mc.cores=mc.cores, repSims=FALSE, writeBam=FALSE, verbose=verbose)
  gc()
  if(verbose) cat("Finished simulations\n")
  notinpc <- names(genomeDB@transcripts)[!(names(genomeDB@transcripts) %in% names(pc@counts[[1]]))]
  pcs <- vector(length=length(pc@counts[[1]]) + length(notinpc), mode='list')
  names(pcs) <- c(names(pc@counts[[1]]), notinpc)
  pcs[names(pc@counts[[1]])] <- pc@counts[[1]]
  pc@counts[[1]] <- pcs
  list(pc=pc, pis=pis, distrsim=distrs)
}
