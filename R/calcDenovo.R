#########################################################################
## DEFINITION AND METHODS FOR CLASS denovoGeneExpr AND denovoGenomeExpr
#########################################################################

require(methods)

setClass("denovoGeneExpr", representation(posprob = "data.frame", expression = "data.frame", variants = "matrix", integralSum= "numeric", npathDeleted= "numeric"))

valid_denovoGeneExpr <- function(object) {
  msg <- NULL
  if (any(!(c('model','posprob') %in% colnames(object@posprob)))) msg <- "Incorrect colnames in 'posprob'"
  if (any(!(c('model','expr','varName') %in% colnames(object@expression)))) msg <- "Incorrect colnames in 'expression'"
  if (is.null(msg)) { TRUE } else { msg }
}

setValidity("denovoGeneExpr", valid_denovoGeneExpr)

setMethod("show", signature(object="denovoGeneExpr"), function(object) {
  cat("denovoGeneExpr object\n")
  cat("\nPosterior model probabilities\n")
  show(head(object@posprob))
  cat("...\nEstimated expression (conditional on each model)\n")
  show(head(object@expression))
  cat("...\nUse posprob() to access posterior probabilities; variants() to get exons in each variant\n")
}
)

setGeneric("posprob", function(object) standardGeneric("posprob"))
setMethod("posprob", signature(object="denovoGeneExpr"), function(object) {
  object@posprob
}
)

setMethod("variants", signature(object="denovoGeneExpr"), function(object) {
  object@variants
}
)

setGeneric("variants<-", function(object,value) standardGeneric("variants<-"))
setReplaceMethod("variants", "denovoGeneExpr", function(object, value) { object@variants <- value; object })



#Class denovoGenomeExpr
setClass("denovoGenomeExpr", representation(islands = "list", txLength = "numeric", priorq = "numeric"))

valid_denovoGenomeExpr <- function(object) {
  msg <- NULL
  if (any(sapply(object@islands,class)!='denovoGeneExpr')) msg <- "All elements must be of class denovoGeneExpr"
  if (is.null(msg)) { TRUE } else { msg }
}

setValidity("denovoGenomeExpr", valid_denovoGenomeExpr)

setMethod("show", signature(object="denovoGenomeExpr"), function(object) {
  cat("denovoGenomeExpr object with",length(object@islands),"gene islands\n\n")
  cat("Individual islands can be accessed via '[[' and '['\n")
  cat("Use relativeExpr() to obtain expression estimates; variants() to see exons in each variant\n")
}
)

setMethod("[", signature(x="denovoGenomeExpr"), function(x, i, ...) { new("denovoGenomeExpr", islands=x@islands[i]) })
setMethod("[[", signature(x="denovoGenomeExpr"), function(x, i, j, ...) { x@islands[[i]] } )
setMethod("as.list", signature(x="denovoGenomeExpr"), function(x) {x@islands})
setMethod("names", signature(x="denovoGenomeExpr"), function(x) names(as.list(x)))

setMethod("variants", signature(object="denovoGenomeExpr"), function(object) {
  ans <- lapply(as.list(object), 'variants')
  data.frame(islandid= rep(names(object), sapply(ans,nrow)), do.call(rbind,ans))
}
)





#########################################################################
## Function calcDenovo
#########################################################################

calcDenovo <- function(distrs, targetGenomeDB, knownGenomeDB=targetGenomeDB, pc, readLength, islandid, priorq=3, mprior, minpp=0.001, selectBest=FALSE, method='submodels', niter, exactMarginal=TRUE, integrateMethod='plugin', verbose=TRUE, mc.cores=1) {
  if (integrateMethod=='plugin') {
      integrateMethod <- as.integer(0)
  } else if (integrateMethod=='Laplace') {
      integrateMethod <- as.integer(1)
  } else { integrateMethod <- as.integer(2) }
  if (missing(readLength)) stop("readLength must be specified")
  if (class(targetGenomeDB)!='annotatedGenome') stop("targetGenomeDB must be of class 'annotatedGenome'")
  if (class(knownGenomeDB)!='annotatedGenome') stop("knownGenomeDB must be of class 'annotatedGenome'")
  if (class(pc)!="pathCounts") stop("pc must be of class 'pathCounts'")
  if (missing(mprior)) { mprior <- modelPrior(knownGenomeDB, verbose=verbose) }
  if (!all(c('nvarPrior','nexonPrior') %in% slotNames(mprior))) stop("Incorrect mprior. Please use modelPrior to generate it.")
  modelUnifPrior <- as.integer(0)
  nvarPrior <- as.list(data.frame(t(mprior@nvarPrior$nbpar)))
  nexonPrior <- as.list(data.frame(t(mprior@nexonPrior$bbpar)))
  #Set Uniform prior on the model space (deprecated)
  #nvarPrior <- list(nbpar=matrix(c(0,0),nrow=1),obs=NA,pred=NA)
  #nexonPrior <- list(bbpar=matrix(c(0,0),nrow=1),obs=NA,pred=NA)
  #modelUnifPrior <- as.integer(1)
  if (!(method %in% c('auto','rwmcmc','priormcmc','allmodels','submodels'))) stop("method must be auto, rwmcmc, priormcmc, allmodels or submodels")

  #Format input
  if (verbose) cat("Formatting input...\n")
  sseq <- seq(0,1,.001)
  startcdf <- as.double(distrs@stDis(sseq))

  lenvals <- as.integer(names(distrs@lenDis))
  lenvals <- as.integer(seq(min(lenvals),max(lenvals),1))
  lendis <- rep(0,length(lenvals)); names(lendis) <- as.character(lenvals)
  lendis[names(distrs@lenDis)] <- as.double(distrs@lenDis/sum(distrs@lenDis))

  readLength <- as.integer(readLength)
  priorq <- as.double(priorq)
  nvarPrior <- lapply(nvarPrior,as.double)
  nexonPrior <- lapply(nexonPrior,as.double)
  minpp <- as.double(minpp)
  selectBest <- as.integer(selectBest)
  method <- as.integer(switch(method, auto=0, allmodels=1, rwmcmc=2, priormcmc=3, submodels=4))
  verbose <- as.integer(verbose)
  exactMarginal <- as.integer(exactMarginal)
  if (missing(islandid)) {
    islandid <- names(targetGenomeDB@islands)[elementLengths(targetGenomeDB@islands)>1]
    islandid <- islandid[islandid %in% names(pc@counts[[1]])]
  }
  else if(is.null(islandid))
  {
    islandid <- names(targetGenomeDB@islands)[elementLengths(targetGenomeDB@islands)>1]
    islandid <- islandid[islandid %in% names(pc@counts[[1]])]
  }
  if (!all(islandid %in% names(targetGenomeDB@islands))) stop('islandid not found in targetGenomeDB@islands')
  if (!all(islandid %in% unlist(lapply(pc@counts, names)))) stop('islandid not found in pc')
  if (!all(islandid %in% names(targetGenomeDB@transcripts))) stop('islandid not found in targetGenomeDB@transcripts')
  exons <- as.integer(names(targetGenomeDB@islands@unlistData))
  names(exons) <- rep(names(targetGenomeDB@islands), elementLengths(targetGenomeDB@islands))
  exons <- split(unname(exons), names(exons))
  exonwidth <- width(targetGenomeDB@islands@unlistData)
  names(exonwidth) <- rep(names(targetGenomeDB@islands), elementLengths(targetGenomeDB@islands))
  exonwidth <- split(unname(exonwidth), names(exonwidth))
  strand <- as.character(strand(targetGenomeDB@islands@unlistData))[cumsum(c(1, elementLengths(targetGenomeDB@islands)[-length(targetGenomeDB@islands)]))]
  names(strand) <- names(targetGenomeDB@islands)
  
  if (missing(niter)) {
     niter <- as.list(as.integer(ifelse(sapply(exons,length)>20,10^3,10^4)))
  } else {
     niter <- as.list(as.integer(rep(niter[1],length(islandid))))
  }
  names(niter) <- islandid
  if (length(targetGenomeDB@knownVars)>0) warning('knownVars in targetGenomeDB are not used by calcDenovo')
  
  #Define basic function
  f <- function(z) {
    islandid <- as.integer(z)
    exons <- exons[z]
    exonwidth <- exonwidth[z]
    transcripts <- targetGenomeDB@transcripts[z]
    prioradj <- lapply(knownGenomeDB@transcripts[z], function(zz) as.double(c(length(zz),mean(sapply(zz,length)))))
    knownVars <- lapply(1:length(z), function(z) character(0))
#    if (length(genomeDB@knownVars)>0) { knownVars <- genomeDB@knownVars[z] } else { knownVars <- lapply(1:length(z), function(z) character(0)) }
    tmp <- strand[z]
    strand <- vector(mode='integer', length=length(tmp))
    sel <- tmp=='+'
    strand[sel] <- 1
    sel <- tmp=='-'
    strand[sel] <- -1
    sel <- tmp=='*'
    strand[sel] <- 0
    strand <- as.list(as.integer(strand))
    #pc <- pc[z] pc's have to be subset from previous step to deal with strandedness !!!!!!
    ans <- calcDenovoMultiple(exons=exons,exonwidth=exonwidth,transcripts=transcripts,knownVars=knownVars,islandid=as.list(islandid),pc=pc@counts[[1]][z],startcdf=startcdf,lendis=lendis,lenvals=lenvals,readLength=readLength,modelUnifPrior=modelUnifPrior,nvarPrior=nvarPrior,nexonPrior=nexonPrior,prioradj=prioradj,priorq=priorq,minpp=minpp,selectBest=selectBest,method=method,niter=niter[z],exactMarginal=exactMarginal,verbose=verbose, integrateMethod=integrateMethod, strand=strand)
    lapply(1:length(ans), function(y) formatDenovoOut(ans[[y]], targetGenomeDB@islands[z][[y]]))
  }

  formatZeroExpr <- function(ids){
    isl <- targetGenomeDB@islands[ids]
    txs <- targetGenomeDB@transcripts[ids]
    exo <- lapply(txs, function(z) { ans <- cbind(names(z), sapply(z, paste, collapse=',')); colnames(ans) <- c('varName','exons'); rownames(ans) <- NULL; ans })
    expr <- lapply(exo, function(x) data.frame(model=rep(0, nrow(x)), expr=rep(1/nrow(x), nrow(x)), varName=x[,'varName']))
    posprob <- lapply(ids, function(x) data.frame(model=0, posprob=NA, priorprob=NA))
    names(posprob) <- ids
    res <- lapply(ids, function(x) new("denovoGeneExpr", variants=exo[[x]], expression=expr[[x]], posprob=posprob[[x]]))
    names(res) <- ids
    res
  }
    
  runCalc <- function(islandid) {
    sel <- !(sapply(pc@counts[[1]][islandid], is.null))
    #sel <- !(sapply(pc@counts[[1]][islandid], is.null) | strand[islandid]=='*')
    #if (verbose==1) cat("Note: Islands with transcripts from both strands will not be processed at the moment\n")
    all <- islandid
    islandid <- islandid[sel]
    if(sum(sel)==0) stop("No islands left to process due to strand or lack of path counts\n")
    if (mc.cores>1 && length(islandid)>mc.cores) {
      require(parallel)
      if ('parallel' %in% loadedNamespaces()) {
        #ans <- mclapply(islandid, f, mc.cores=mc.cores)
        #split into smaller jobs
        nsplit <- ceiling(max(length(islandid), mc.cores)/mc.cores)
        islandidList <- lapply(1:min(length(islandid), mc.cores), function(z) islandid[seq(z,length(islandid),by=mc.cores)])
        ans <- parallel::mclapply(islandidList,f,mc.cores=min(length(islandidList), mc.cores))
        ans <- do.call(c,ans); names(ans) <- unlist(islandidList); ans <- ans[islandid]
      } else stop('parallel library has not been loaded!')
    } else {
      ans <- f(islandid)
      names(ans) <- islandid
    }
    res <- vector(mode="list", length=length(all))
    names(res) <- all
    z <- formatZeroExpr(all[!sel])
    res[names(z)] <- z
    res[names(ans)] <- ans
    res
  }

  #Initialize transcripts for new islands with known orientation
  sel <- names(targetGenomeDB@transcripts)[sapply(targetGenomeDB@transcripts,is.null) & strand=='*']
  if (length(sel)>0) targetGenomeDB@transcripts[sel] <- tapply(as.integer(names(targetGenomeDB@islands[sel]@unlistData)), rep(names(targetGenomeDB@islands[sel]), elementLengths(targetGenomeDB@islands[sel])), function(x) list(as.numeric(x))) 
  islandidUnknown <- islandid[islandid %in% names(targetGenomeDB@transcripts)[sapply(targetGenomeDB@transcripts,is.null)]]
  if (length(islandidUnknown)>0) { islandidini <- islandid; islandid <- islandid[!(islandid %in% islandidUnknown)] }


  #Run
  if (verbose==1) cat("Performing model search (this may take a while)\n")
  if (length(islandidUnknown)==0) {
    ans <- runCalc(islandid)
  } else {
    #Islands with known strand
    ans <- vector("list",length(islandidini)); names(ans) <- islandidini
    if (length(islandid)>0) ans[islandid] <- runCalc(islandid)
    
    #Islands with unknown strand. Run 2 strands and select the one with largest post prob
    targetGenomeDB@transcripts[islandidUnknown] <- lapply(targetGenomeDB@islands[islandidUnknown],function(z) list(var1=as.integer(names(z))))
    strand[islandidUnknown] <- '+'
    ansforw <- runCalc(islandidUnknown)
    strand[islandidUnknown] <- '-'
    targetGenomeDB@transcripts[islandidUnknown] <- lapply(targetGenomeDB@transcripts[islandidUnknown],rev)
    targetGenomeDB@islands[islandidUnknown] <- lapply(targetGenomeDB@islands[islandidUnknown], rev)
    exons[islandidUnknown] <- lapply(exons[islandidUnknown], rev)
    exonwidth[islandidUnknown] <- lapply(exonwidth[islandidUnknown], rev)
    ansrev <- runCalc(islandidUnknown)
    difndel <- sapply(ansforw,function(z) z@npathDeleted) - sapply(ansrev,function(z) z@npathDeleted)
    difmax <- sapply(ansforw,function(z) z@integralSum['logmax']) - sapply(ansrev,function(z) z@integralSum['logmax'])
    difsum <- log(sapply(ansforw,function(z) z@integralSum['sum'])) - log(sapply(ansrev,function(z) z@integralSum['sum']))
    sel <- ifelse(difndel<0 || (difndel==0 && (difsum+difmax)>=0), TRUE, FALSE)
    ans[islandidUnknown[sel]] <- ansforw[sel]; ans[islandidUnknown[!sel]] <- ansrev[!sel]
  }
  if (verbose==1) cat("\n")
  ans <- new("denovoGenomeExpr", islands=ans, priorq=priorq)
  #compute tx length
  fdata <- variants(ans)
  names(fdata)[1:2] <- c('gene','transcript')
  rownames(fdata) <- as.character(fdata$transcript)
  tx <- strsplit(as.character(fdata$exons),split=',')
  tx <- data.frame(tx=rep(fdata$transcript, sapply(tx,length)) , exon=unlist(tx))
  txLength <- txLengthBase(tx=tx, genomeDB=targetGenomeDB)
  ans@txLength <- txLength
  return(ans)
}


formatDenovoOut <- function(ans, genesel) {
  ans[[1]] <- data.frame(ans[[1]])
  colnames(ans[[1]]) <- c('model','posprob','priorprob')
  ans[[1]] <- ans[[1]][order(ans[[1]][,'posprob'],decreasing=TRUE),]
  ans[[2]] <- data.frame(ans[[2]],ans[[3]])
  colnames(ans[[2]]) <- c('model','expr','varName')
  ans[[3]] <- NULL
  colnames(ans[[3]]) <- c('varName','exons')
  names(ans[[4]]) <- c('sum','logmax')
  names(ans) <- c('posprob','expression','variants','integralSum','npathDeleted')
  new("denovoGeneExpr",posprob=ans$posprob,expression=ans$expression,variants=ans$variants,integralSum=ans$integralSum,npathDeleted=ans$npathDeleted)
}

calcDenovoMultiple <- function(exons, exonwidth, transcripts, knownVars, islandid, pc, startcdf, lendis, lenvals, readLength, modelUnifPrior, nvarPrior, nexonPrior, prioradj, priorq, minpp, selectBest, method, niter, exactMarginal, verbose, integrateMethod, strand) {
  ans <- .Call("calcDenovoMultiple",exons,exonwidth,transcripts,knownVars,islandid,pc,startcdf,lendis,lenvals,readLength,modelUnifPrior,nvarPrior,nexonPrior,prioradj,priorq,minpp,selectBest,method,niter,exactMarginal,verbose,integrateMethod,strand)
  return(ans)
}



