
mergePCWr <- function(x, genomeDB){
  pc <- lapply(x, '[[', 'pc')
  tmp <- vector(mode='list', length=length(genomeDB@transcripts))
  names(tmp) <- names(genomeDB@transcripts)
  pcs <- unlist(lapply(pc, function(y) y@counts[[1]]), recursive=F)
  tmp[names(pcs)] <- pcs
  new("pathCounts", counts=list(tmp), stranded=x[[1]]$pc@stranded, denovo=FALSE)
}

mergeDisWr <- function(distrs, pcs){
  lenDis <- lapply(distrs, function(x) x@lenDis)
  lenDis <- lenDis[unlist(lapply(lenDis, function(x) !any(names(x)==0)))]
  lenDis <- lenDis[sapply(lenDis,sum)>0]
  if(length(lenDis)>1){
    minlen <- min(unlist(lapply(lenDis, function(x) min(as.numeric(names(x))))))
    maxlen <- max(unlist(lapply(lenDis, function(x) max(as.numeric(names(x))))))
    tmp <- lapply(lenDis, function(x){
      tmp <- vector(mode='numeric', length=maxlen-minlen+1)
      names(tmp) <- minlen:maxlen
      tmp[names(x)] <- x
      tmp
    } )
    tmp <- do.call(cbind, tmp)
    tmp <- as.array(rowSums(tmp))
  } else tmp <- lenDis[[1]]
  distr <- new('readDistrs', lenDis=tmp)
  th <- seq(0,1,length=10000)
  if (missing(pcs)) w <- sapply(distrs, function(z) sum(z@lenDis)) else w <- sapply(1:length(distrs), function(x) sum(getNreads(pcs[[x]])))
  tmp <- lapply(1:length(distrs), function(x) { all <- distrs[[x]]@stDis(th)*w[x] } )
  tmp <- do.call(cbind, tmp)
  tmp <- rowMeans(tmp)
  tmp <- tmp/tmp[length(tmp)]
  tmp <- approxfun(th, tmp)
  distr@stDis <- tmp
  distr
}

wrapKnown <- function(bamFile, verbose=FALSE, seed=1, mc.cores.int=1, mc.cores=1, genomeDB, readLength, rpkm=TRUE, priorq=2, priorqGeneExpr=2, citype='none', niter=10^3, burnin=100, keep.pbam=FALSE, keep.multihits=TRUE, chroms=NULL) {
  
  if(!exists(as.character(substitute(genomeDB)))) stop("No genomeDB found")
  if (missing(readLength)) stop("readLength must be specified")
  if (genomeDB@denovo) stop("genomeDB must be a known genome")
  if(!(citype %in% c('none', 'asymp', 'exact'))) stop("citype should take the value 'none', 'asymp', or 'exact'") 
  
  if (length(bamFile)==1) {
    ans <- wrapKnownSingle(bamFile=bamFile,verbose=verbose,seed=seed,mc.cores.int=mc.cores.int,mc.cores=mc.cores,genomeDB=genomeDB,readLength=readLength,rpkm=rpkm,priorq=priorq,priorqGeneExpr=priorqGeneExpr,citype=citype,niter=niter,burnin=burnin,keep.pbam=keep.pbam,keep.multihits=keep.multihits,chroms=chroms)
  } else if (length(bamFile)>1) {
    x <- vector("list",length(bamFile))
    for (i in 1:length(bamFile)) {
      cat(paste("\n PROCESSING",bamFile[i],"\n"))
      x[[i]] <- wrapKnownSingle(bamFile=bamFile[i],verbose=verbose,seed=seed,mc.cores.int=mc.cores.int,mc.cores=mc.cores,genomeDB=genomeDB,readLength=readLength,rpkm=rpkm,priorq=priorq,priorqGeneExpr=priorqGeneExpr,citype=citype,niter=niter,burnin=burnin,keep.pbam=keep.pbam,keep.multihits=keep.multihits,chroms=chroms)
    }
    cat("\n MERGING ALL FILES...\n")
    ans <- vector("list",3); names(ans) <- c('pc','distr','exp')
    ans$exp <- mergeExp(lapply(x,'[[','exp'), sampleNames=sub('.bam$','',bamFile), keep=c('transcript','island_id','gene_id','explCnts'))
    ans$distr <- lapply(x,'[[','distr'); names(ans$distr) <- sub('.bam$','',bamFile)
    ans$pc <- lapply(x,'[[','pc'); names(ans$distr) <- sub('.bam$','',bamFile)
  } else {
    stop("Invalid length(bamFile)")
  }
  return(ans)
}

wrapKnownSingle <- function(bamFile, verbose=FALSE, seed=1, mc.cores.int=1, mc.cores=1, genomeDB, readLength, rpkm=TRUE, priorq=2, priorqGeneExpr=2, citype='none', niter=10^3, burnin=100, keep.pbam=FALSE, keep.multihits=TRUE, chroms=NULL) {


  what <- c('qname','strand','pos','mpos','cigar')
  #what <- scanBamWhat(); what <- what[!(what %in% c('seq','qual','qwidth','flag','mapq','mrnm','mpos','isize'))]
  if(!keep.multihits) what <- c(what, 'mapq')
  t <- scanBamHeader(bamFile)[[1]][["targets"]]
  which <- GRanges(names(t), IRanges(1, unname(t)))
  if(!is.null(chroms)) {
    which <- which[as.character(seqnames(which)) %in% chroms]
  } else {
    which <- which[!grepl("_",as.character(seqnames(which)))]
    which <- which[!as.character(seqnames(which))=='chrM']
  }
  sel <- as.vector(seqnames(which)) %in% unique(genomeDB@exon2island$seqnames)
  if (all(!sel)) stop("Did not find any of the specified chromosomes")
  if (any(!sel)) warning(paste("Did not find in genomeDB chromosomes",paste(as.vector(seqnames(which))[!sel],collapse=' '),'. Skipping them'))
  which <- which[sel,]
  if(sum(grepl("_", as.character(seqnames(which))))>0 | sum(grepl("M", as.character(seqnames(which))))>0) cat("Warning, non standard chromosomes included in bam (it is not recommended to include mitochondrial chromosome nor random or unstable chromosomes)") 
  flag <- scanBamFlag(isPaired=TRUE,hasUnmappedMate=FALSE)

 
   
## Define function for one chromosome

  procChr <- function(i) {
    param <- ScanBamParam(flag=flag,what=what, which=which[i], tag='XS')
    cat("Processing chromosome: ", as.character(seqnames(which[i])), "\n")
    bam <- scanBam(file=bamFile,param=param)
    if (length(bam[[1]]$qname)>0) {  #if some reads satisfied the criteria
      if(!keep.multihits) {
        single.hit <- which(bam[[1]][['mapq']]>0)
        bam[[1]] <- lapply(bam[[1]], '[', single.hit)
      }
    
      bam[[1]]$qname <- as.integer(as.factor(bam[[1]]$qname))
      if(verbose) cat(paste("Replaced qname for chr", as.character(seqnames(which[i])), "\n"))
      if(keep.pbam) {
        ans <- vector("list",3); names(ans) <- c("pbam","distr","pc")
        ans$pbam <- procBam(bam=bam[[1]], stranded=FALSE, seed=as.integer(seed), verbose=FALSE, keep.junx=FALSE, rname=as.character(seqnames(which)[i]))
        #ans$distr <- getDistrs(DB=genomeDB, bam=bam[[1]], verbose=verbose, readLength=readLength)
        #cat("Removing bam object\n")
        rm(bam); gc()
        ans$distr <- getDistrs(DB=genomeDB, pbam=ans$pbam, verbose=FALSE)
        ans$pc <- pathCounts(reads=ans$pbam, DB=genomeDB, mc.cores=mc.cores, verbose=FALSE)
      } else {
        ans <- vector("list",2); names(ans) <- c("distr","pc")
        pbam <- procBam(bam=bam[[1]], stranded=FALSE, seed=as.integer(seed), verbose=FALSE, keep.junx=FALSE, rname=as.character(seqnames(which)[i]))
        #ans$distr <- getDistrs(DB=genomeDB, bam=bam[[1]], verbose=verbose, readLength=readLength)
        #cat("Removing bam object\n")
        rm(bam); gc()
        ans$pc <- pathCounts(reads=pbam, DB=genomeDB, mc.cores=mc.cores, verbose=FALSE)
        ans$distr <- getDistrs(DB=genomeDB, pbam=pbam, verbose=FALSE)
        #rm(pbam); gc()
      }
    } else {  #if no reads satisfied the criteria, return integer instead of list
      ans <- as.integer(-1)
    }
    cat("Finished chromosome ", as.character(seqnames(which[i])), "\n")
    return(ans)
  }

  ## Run for all chromosomes, mclapply or for loop
  if (mc.cores.int>1 ){
    if ('parallel' %in% loadedNamespaces()) {
      ans <- parallel::mclapply(1:length(which), procChr, mc.cores=mc.cores.int, mc.preschedule=FALSE)
    } else stop('parallel library has not been loaded!')
  } else {
    ans <- list()
    for(i in 1:length(which)) ans[[i]] <- procChr(i)
  }
  gc()
  ans <- ans[sapply(ans,is.list)] #discard chromosomes with no reads
  if (length(ans)>0) {
    if(keep.pbam) allpbam <- lapply(ans, "[[", "pbam")
    allpc <- mergePCWr(ans, genomeDB)
    alldistr <- suppressWarnings(mergeDisWr(lapply(ans, '[[', 'distr'), lapply(ans, '[[', 'pc')))
    exp <- calcExp(distrs=alldistr, genomeDB=genomeDB, pc=allpc, readLength=readLength, rpkm=rpkm, priorq=priorq, priorqGeneExpr=priorqGeneExpr, citype=citype, niter=niter, burnin=burnin, mc.cores=mc.cores, verbose=verbose)
    if(keep.pbam) {
      ans <- list(pc=allpc, distr=alldistr, exp=exp, pbam=allpbam)
    } else ans <- list(pc=allpc, distr=alldistr, exp=exp)
  } else {
    ans <- list(pc=NA, distr=NA, exp=NA)
  }
  return(ans)
}

