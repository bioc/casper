require(methods)

#setGeneric("wrapKnown", function(bamFile, mc.cores) standardGeneric("wrapKnown"))
#setMethod("wrapKnown", signature(bamFile='character', mc.cores='missing'), function(bamFile, mc.cores){
  

 # setMethod("wrapKnown", signature(bamFile='character', mc.cores='numeric'), function(bamFile, mc.cores){
  
          


setGeneric("subsetPbam", function(pbam, strand, chr) standardGeneric("subsetPbam"))
setMethod("subsetPbam", signature(pbam="procBam", strand="character", chr='missing'),
          function(pbam, strand) {
            if(strand=="+") {
              pbam@pbam <- pbam@plus
              pbam@junx <- pbam@pjunx
            } else {
              pbam@pbam <- pbam@minus
              pbam@junx <- pbam@pminus
            }
            pbam
          }
          )
setMethod("subsetPbam", signature(pbam='procBam', strand='missing', chr='character'),
          function(pbam, chr) {
            pans <- pbam@pbam[seqnames(pbam@pbam) %in% chr]
            #names(pans) <- names(pbam@pbam)[as.character(seqnames(pbam@pbam)) %in% chr]
            jans <- pbam@junx[seqnames(pbam@junx) %in% chr]
            #names(jans) <- names(pbam@junx)[as.character(seqnames(pbam@junx)) %in% chr]
            pbam@pbam <- pans
            pbam@junx <- jans
            pbam
          }
          )


buildRD<-function(reads){
  sel <- reads$end<reads$start
  st <- en <- integer(length(reads$start))
  st[sel] <- reads$end[sel]; st[!sel] <- reads$start[!sel]
  en[!sel] <- reads$end[!sel]; en[sel] <- reads$start[sel]
  if(sum(grepl("chr", unique(reads$chrom)))>0) {
    ans<-GRanges(ranges=IRanges(start=st, end=en), seqnames=reads$chrom, rid=reads$rid, strand=reads$rstrand, XS=Rle(reads$strand), names=reads$key)
  } else {
    ans<-GRanges(ranges=IRanges(start=st, end=en), seqnames=paste("chr", reads$chrom, sep=""), rid=reads$rid, strand=reads$rstrand, XS=Rle(reads$strand), names=reads$key)
    #names(ans) <- reads$key
  }
  if("flag" %in% names(reads)) mcols(ans)$flag <- reads$flag
  return(ans)
}

uniquifyQname<-function(bam, seed=1){   	
  qname <- vector(mode='character', length=length(bam$qname))
  dups <- .Call("uniqQname", bam$qname, length(bam$qname), bam$pos, bam$mpos, qname)
  probPos <- dups[[1]] %in% dups[[2]]
  bam <- lapply(bam, function(x) x[!probPos])    
  qdup <- bam$qname[bam$qname %in% bam$qname[duplicated(bam$qname)]] 
  uni <- bam$qname %in% qdup
  res <- lapply(bam, "[", uni)
  names(res)<-names(bam)
  res
}

nbReads <- function(bam0) {
  tab <- table(bam0$cigar)
  count <- sapply(gregexpr("M",names(tab)),length)
  tjunx <- sum(tab*(count-1))
  treads <- sum(tab*count)
  c(tjunx=tjunx, treads=treads)
}

procB <- function(bam, strnd, seed=1, verbose=FALSE, rname='null', keep.junx=FALSE, keep.flag=FALSE){
  echrom=TRUE
  if(!rname=='null') echrom=FALSE
  lev=NULL
  if(class(bam$strand)=='factor') lev <- levels(bam$strand)
  if(!echrom) bam$rname<-as.character(bam$rname)
  if(verbose) cat("Calculating total number of reads...\n")
  if(class(bam$qname)!='integer')   bam$qname <- as.integer(as.factor(bam$qname))
  nreads <- nbReads(bam)
  njunx=1
  if(keep.junx) njunx <- nreads['tjunx']
  nreads <- nreads['treads']+10
  
  if(verbose) cat("done.\nProcessing cigars and building read's object...\n")
  len=vector(mode="integer", length=nreads)
  strs=vector(mode="integer", length=nreads)
  rid=vector(mode="integer", length=nreads)
  key=vector(mode="character", length=nreads)
  chrom=vector("character", length=1)
  if(echrom) chrom=vector(mode="character", length=nreads)
  flag=vector(mode='integer', length=1)
  if(keep.flag) flag=vector(mode='integer', length=nreads)
  if(!("flag" %in% names(bam))) bam$flag <- 1
  strand=vector(mode="integer", length=nreads)
  jchrom=vector("character", length=1)
  if(echrom) jchrom=vector(mode="character", length=njunx)
  jstrs=vector(mode="integer", length=njunx)
  jlen=vector(mode="integer", length=njunx)
  if(class(bam$rname)=='factor') bam$rname <- as.character(bam$rname)
  if(echrom) {
    data<-.Call("procBam", bam$qname, bam$rname, bam$pos, bam$mpos, bam$cigar, as.integer(bam$strand), length(bam$pos), as.integer(nreads), as.integer(bam$flag), as.integer(njunx), len, strs, key, chrom, rid, strand, jchrom, jstrs, jlen, flag)
  } else data<-.Call("procBam", bam$qname, rname, bam$pos, bam$mpos, bam$cigar, as.integer(bam$strand), length(bam$pos), as.integer(nreads), as.integer(bam$flag), as.integer(njunx), len, strs, key, chrom, rid, strand, jchrom, jstrs, jlen, flag)
  names(data) <- c('end', 'start', 'key', 'chrom', 'rid', 'rstrand', 'jchrom', 'jstart', 'jend', 'flag')
  if(!echrom) {
    data$chrom <- rep(rname, length(data$start))
    data$jchrom <- rep(rname, length(data$jstart))
  }
  if(!keep.flag) data$flag <- NULL
  if(keep.junx) {
    junx <- GRanges(ranges=IRanges(start=data$jstart[data$jstart!=0], end=data$jend[data$jstart!=0]), seqnames=data$jchrom[data$jstart!=0], XS=Rle(rep(strnd, sum(data$jstart!=0))))
    junx <- lapply(unique(as.character(seqnames(junx))), function(x){
      if(sum(seqnames(junx)==x)){
        y <- junx[seqnames(junx) == x]
        z <- paste(start(y), end(y), sep='.')
        zz <- table(z)
        jst <- strsplit(names(zz), split='.', fixed=T)
        jend <- as.numeric(sapply(jst, '[', 2))
        jst <- as.numeric(sapply(jst, '[', 1))
        GRanges(IRanges(jst, jend), seqnames=x, counts=as.numeric(zz))
      }
    })
    junx <- mergeGRanges(junx)
  }
  else junx <- GRanges(IRanges(0,0), seqnames=1)
  if(!is.null(lev)) data[[6]] <- lev[data[[6]]]
  sel<-data[[1]]!=0 & !is.na(data[[3]])
  data<-lapply(data, "[", sel)
  data[[length(data)+1]] <- rep(strnd, length(data[[1]]))
  names(data)[length(data)] <- 'strand'
  if(verbose) cat("done.\n")
  id <- data$key
  idnum <- rep(0,length(data$key))
  idnum[which(data$key[-length(data$key)] != data$key[-1])+1] <- 1
  data$key <- as.integer(cumsum(idnum))
  ans<-buildRD(data)
  list(pbam=ans, junx=junx)
}

setMethod("procBam", signature(bam='list',stranded='missing',seed='missing', verbose='missing', rname='missing') ,
          function(bam, ...) procBam(bam=bam, stranded=FALSE, seed=as.integer(1), verbose=FALSE, rname='null', keep.junx, keep.flag)
          )
setMethod("procBam", signature(bam='list',stranded='missing',seed='missing', verbose='logical', rname='missing') ,
          function(bam, stranded, seed, verbose, rname, ...) procBam(bam=bam, stranded=FALSE, seed=as.integer(1), verbose=verbose, rname='null', keep.junx, keep.flag)
          )
setMethod("procBam", signature(bam='list',stranded='missing',seed='integer', verbose='missing', rname='missing') ,
          function(bam, seed, ...) procBam(bam=bam, stranded=FALSE, seed=seed, verbose=FALSE, rname='null', keep.junx, keep.flag)
          )
setMethod("procBam", signature(bam='list',stranded='logical',seed='missing', verbose='logical', rname='missing') ,
          function(bam, stranded, seed, verbose, ...) procBam(bam=bam, stranded=stranded, seed=as.integer(1), verbose=verbose, rname='null', keep.junx, keep.flag)
          )
setMethod("procBam", signature(bam='list',stranded='logical',seed='missing', verbose='missing', rname='missing') ,
          function(bam, stranded, seed, verbose, ...) procBam(bam=bam, stranded=stranded, seed=as.integer(1), verbose, rname='null', keep.junx, keep.flag)
          )
setMethod("procBam", signature(bam='list',stranded='logical',seed='integer', verbose='logical', rname='missing') ,
          function(bam, stranded, seed, verbose, rname, ...) procBam(bam=bam, stranded=stranded, seed=seed, verbose=verbose, rname="null", keep.junx, keep.flag)
          )
setMethod("procBam", signature(bam='list',stranded='logical',seed='integer', verbose='logical', rname='character') ,
          function(bam, stranded, seed, verbose, rname, keep.junx, keep.flag) {
            byList <- ifelse(class(bam[[1]])=='list', TRUE, FALSE)
            if(!byList) {
              ans <- procBamF(bam=bam, stranded=stranded, seed=seed, verbose=verbose, rname=rname, keep.junx=keep.junx, keep.flag=keep.flag)
            } else {
              ans <- lapply(bam, function(x) procBamF(bam=x, stranded=stranded, seed=seed, verbose=verbose, rname=rname, keep.junx=keep.junx, keep.flag=keep.flag))
              if(length(ans)>1) {
                ans <- mergePbam(ans, fixNames=TRUE)
                lens <- seqnames(ans@pbam)@lengths
                fix <- rep(c(0, (as.numeric(names(ans@pbam)[cumsum(lens)])+1)), c(seqnames(ans@pbam)@lengths, 0))
                names(ans@pbam) <- as.numeric(names(ans@pbam))+fix
              }
            }
            ans
          }
            )

mergeGRanges <- function(gr, fixNames=FALSE) {
  values <- do.call('rbind',unname(lapply(gr, values)))
  colnames(values) <- sub("values.", "", colnames(values))
  seqn <- unlist(RleList(unname(lapply(gr, seqnames))))
  ranges <- unlist(IRangesList(lapply(gr, ranges)))
  ans <- GRanges(ranges=ranges, seqnames=seqn)
  names(ans) <- unlist(sapply(gr, names))
  if(fixNames){
    lens <- sapply(gr, length)
    lens <- rep(c(0, lens[-length(lens)]), lens)
    n <- as.numeric(names(ans))+lens
    names(ans) <- n
  }
  values(ans) <- values
  ans
}

mergePbam <- function(ans, fixNames=FALSE){
  if(!all(sapply(ans, function(x) x@stranded)) & sum(sapply(ans, function(x) x@stranded))!=0) stop("all pbams must be either stranded or non-stranded")
  if(ans[[1]]@stranded){
    plus <- casper:::mergeGRanges(lapply(ans, function(x) slot(x, 'plus')), fixNames)
    minus <- casper:::mergeGRanges(lapply(ans, function(x) slot(x, 'minus')), fixNames)
    pjunx <- casper:::mergeGRanges(lapply(ans, function(x) slot(x, 'pjunx')), fixNames)
    mjunx <- casper:::mergeGRanges(lapply(ans, function(x) slot(x, 'mjunx')), fixNames)
    new("procBam", pbam=NULL, plus=plus, minus=minus,junx=NULL, pjunx=pjunx, mjunx=mjunx, stranded=ans[[1]]@stranded)
  } else{
    pbam <- mergeGRanges(lapply(ans, function(x) slot(x, 'pbam')), fixNames)
    junx <- mergeGRanges(lapply(ans, function(x) slot(x, 'junx')), fixNames)
    ans <- new("procBam",pbam=pbam, junx=junx,stranded=ans[[1]]@stranded)
  }
  ans
}

    
procBamF<-function(bam, stranded=FALSE, seed=1,  verbose=FALSE, rname='null', keep.junx, keep.flag){
  require(IRanges)
  if(stranded) {
    minus <- bam$tag$XS=='-'
    if(keep.junx) mjunx <- GRanges(IRanges(0,0), seqnames=1)
    if(sum(minus)>0){
      minus <- lapply(bam, '[', minus)
      minus <- procB(minus, "-", seed=seed, verbose=verbose, rname=rname, keep.junx=keep.junx, keep.flag=keep.flag)
      if(keep.junx) mjunx <- minus$junx
      minus <- minus$pbam
    } else minus <- GRanges(IRanges(0,0), seqnames=1)
    plus <- bam$tag$XS=='+'
    mplus <- GRanges(IRanges(0,0), seqnames=1)
    if(sum(plus)>0){
      plus <- lapply(bam, '[', plus)
      plus <- procB(plus, "+", seed=seed, verbose=verbose, rname=rname, keep.junx=keep.junx, keep.flag=keep.flag)
      if(keep.junx) pjunx <- plus$junx
      plus <- plus$pbam
    } else plus <- GRanges(IRanges(0,0), seqnames=1)
    if(keep.junx) {
      ans <- list(plus=plus, pjunx=pjunx, minus=minus, mjunx=mjunx, stranded=TRUE)
    } else ans <- list(plus=plus, minus=minus, stranded=TRUE)
  } else {
    pbam <- procB(bam, "*", seed=seed, verbose=verbose, rname=rname, keep.junx=keep.junx, keep.flag=keep.flag)
    if(keep.junx) {
      ans <- list(pbam=pbam$pbam, junx=pbam$junx, stranded=FALSE)
    } else ans <- list(pbam=pbam$pbam, stranded=FALSE)    
  }
  if(!ans$stranded) {
    if(!keep.junx) ans$junx <- GRanges(IRanges(0,0), seqnames=1)
    ans <- new("procBam",pbam=ans$pbam, junx=ans$junx,stranded=ans$stranded)
  } else {
    if(!keep.junx) { ans$pjunx <- GRanges(IRanges(0,0), space=NULL); ans$mjunx <- GRanges(IRanges(0,0), seqnames=1)}
    ans <- new("procBam",pbam=NULL, plus=ans$plus, minus=ans$minus,junx=NULL, pjunx=ans$pjunx, mjunx=ans$mjunx, stranded=ans$stranded)
  }
}
