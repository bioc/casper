startDist <- function(st,fragLength,txLength, nreads=NULL) {
                                        # Estimate relative start distribution under left­truncation (st < 1 ­ fragLength/txLength)
                                        # ­ st: relative start (i.e. start/txLength)
                                        # ­ fragLength: fragment length
                                        # ­ txLength: transcript length
                                        # Output: cumulative probability function (actually, a linear interpolation)
  if(!is.null(nreads)){
    if(nreads<length(st)){
      ran <- sample(1:length(st), size=nreads, replace=T)
      st <- st[ran]
      fragLength <- fragLength[ran]
      txLength <- txLength[ran]
    }
  }
  
  trunc <- 1-fragLength/txLength
  sel <- trunc>(st+1e-10)
  trunc <- 1-trunc[sel]
  st <- 1-st[sel]

  fit <- summary(survfit(Surv(time=trunc,time2=st,event=rep( TRUE,length(st))) ~ 1))
  s <- 1-fit$time
  pcum <- fit$surv
  f <- approxfun(s,pcum)
  sseq <- seq(0,1,.001)
  startcdf <- f(sseq)
  startcdf[1] <- 0; startcdf[length(startcdf)] <- 1
  f <- approxfun(sseq[!is.na(startcdf)], startcdf[!is.na(startcdf)])
  return(f)
}

firstBamReads <- function(bam, nreads) {
  if (nreads < length(bam$qname)) {
    id <- unique(bam$qname[1:nreads])
    sel <- bam$qname %in% id
    bam <- lapply(bam,function(z) z[sel])
  }
  return(bam)
}

getDistrs<-function(DB, bam, islandid=NULL, verbose=FALSE, nreads=4*10^6){
  
  #Format exons as RangedData
  if(verbose) cat("Calculating fragment length distribution\n")

  exonsRD <- DB@exonsNI  
  #Select a sample of reads
  bam <- firstBamReads(bam, nreads=nreads)

  #Find fragment length distribution for fragments aligning to exons larger than 1000 bases

  
  if (!all(c('qname','rname','qwidth','pos','mpos') %in% names(bam))) stop('bam must contain elements qname, rname, qwidth, pos, mpos')
#  tab <- table(bam$qname)
#  sel <- bam$qname %in% names(tab[tab==2])
#  bam <- bam[!(names(bam) %in% 'tag')]
#  bam <- lapply(bam, '[', sel)
  
  d <- bam$mpos - bam$pos
  sel <- d<0; n <- bam$qname[sel]; sp <- bam$rname[sel]; names(sp) <- n; en <- bam$pos[sel]+bam$qwidth[sel]-1; names(en) <- n
  sel <- d>0; st <- bam$pos[sel]; names(st) <- bam$qname[sel];
  sel <- match(n, names(st))
  st <- st[sel[!is.na(sel)]]
  en <- en[names(st)]; sp <- sp[names(st)]
  sel <- st<en; st <- st[sel]; en <- en[sel]; sp <- sp[sel]

  if(!any(unique(sp) %in% names(exonsRD))) {
    if(any(grepl('chr', names(exonsRD)))) {
      sp <- as.factor(paste('chr', as.character(sp), sep=''))
    } else if(any(grepl('chr', unique(sp)))) sp <- sub('chr', '', sp)
  }

  if(!any(unique(sp) %in% levels(seqnames((exonsRD))))) {
    if(any(grepl('chr', levels(seqnames(exonsRD))))) {
      sp <- as.factor(paste('chr', as.character(sp), sep=''))
    } else if(any(grepl('chr', unique(sp)))) sp <- sub('chr', '', sp)
  }
  
  if(!any(unique(sp) %in% levels(seqnames(exonsRD)))) stop('Different chromosome names in bam and genome')
  frags <- GRanges(IRanges(start=st,end=en),seqnames=sp)
  n <- levels(seqnames(frags))[levels(seqnames(frags)) %in% levels(seqnames(exonsRD))]
  fragsL<-frags[levels(seqnames(frags)) %in% n]
  over <- suppressWarnings(findOverlaps(fragsL, subset(exonsRD, width(exonsRD)>1000), type="within"))
  ld<-table(width(fragsL)[queryHits(over)])
  ld <- ld[ld/sum(ld) > 0.0001]
 
    
  

  #Find fragment start distribution for fragments aligning to transcripts in genes with only one annotated tx

  if(verbose) cat("Calculating fragment start distribution\n")
  if(is.null(islandid)){
    sel <- unlist(lapply(DB@transcripts, length))==1
  } else sel= islandid
  oneTx <- DB@transcripts[sel]
  oneTx <- unlist(oneTx, recursive=F)
  names(oneTx) <- sub("[0-9]+\\.", "", names(oneTx))
  oneExons <- exonsRD[names(exonsRD) %in% unlist(oneTx)]
  oneExons <- oneExons[as.character(unlist(oneTx))]
  
  islandStrand <- as.character(strand(DB@islands@unlistData))[cumsum(c(1, elementLengths(DB@islands)[-length(DB@islands)]))]
  names(islandStrand) <- names(DB@islands)
  

  exon_rank <- unlist(lapply(oneTx, function(x) 1:length(x)))
  exon_strand <- islandStrand[as.character(DB@exon2island$island[match(names(oneExons), rownames(DB@exon2island))])]
  strand(oneExons) <- exon_strand
  values(oneExons)$exon_rank <- exon_rank
  txid <- cumsum( values(oneExons)$exon_rank==1 )
  wid <- end(oneExons) - start(oneExons) + 1
  values(oneExons)$exstnogap <- unlist(tapply(wid, txid, function(z) c(0, cumsum(z[-length(z)]))))
  values(oneExons)$txlength <- unlist(tapply((end(oneExons) - start(oneExons) + 1),INDEX=txid,function(z) rep(sum(z),length(z))))
  
  frags <- frags[width(frags)<max(width(oneExons))]
  over <- suppressWarnings(findOverlaps(frags, oneExons, type="within"))

  exstnogap <- values(oneExons)$exstnogap[subjectHits(over)]
  txlength <- values(oneExons)$txlength[subjectHits(over)]
  exst <- start(oneExons)[subjectHits(over)]
  exen <- end(oneExons)[subjectHits(over)]
  readst <- start(frags)[queryHits(over)]
  readen <- end(frags)[queryHits(over)]
  str <- as.character(strand(oneExons))[subjectHits(over)]
  frlen <- width(frags)[queryHits(over)]
  
  stDis <- double(length(readst))
  sel <- str=='+';  stDis[sel] <- (exstnogap[sel]+readst[sel]-exst[sel])/txlength[sel]
  sel <- str=='-'; stDis[sel] <- (exstnogap[sel]+exen[sel]-readen[sel])/txlength[sel]

   stDis <- startDist(stDis, frlen, txlength)
  new("readDistrs",lenDis=ld,stDis=stDis)
}

