require(methods)

## METHODS

setMethod("getIsland",signature(entrezid='character',txid='missing',genomeDB='annotatedGenome'),function(entrezid, txid, genomeDB) {
  as.character(genomeDB@aliases[genomeDB@aliases$gene_id==entrezid,'island_id'][1])
}
)

setMethod("getIsland",signature(entrezid='missing',txid='character',genomeDB='annotatedGenome'),function(entrezid, txid, genomeDB) {
  if (all(txid %in% rownames(genomeDB@aliases))) {
    ans <- as.character(genomeDB@aliases[txid,'island_id'])
  } else {
    warning("Exact match for txid not found, using 'grep' instead. Multiple matches may be returned")
    ans <- as.character(genomeDB@aliases[grep(txid,rownames(genomeDB@aliases)),'island_id'])
  }
  return(ans)
}
)


setMethod("getChr",signature(entrezid='missing',txid='missing',islandid='character', genomeDB='annotatedGenome'), function(entrezid, txid, islandid, genomeDB) {
  as.character(genomeDB@exon2island$seqnames[match(TRUE,genomeDB@exon2island$island == islandid)])
}
)

setMethod("getChr",signature(entrezid='character',txid='missing',islandid='missing', genomeDB='annotatedGenome'), function(entrezid, txid, islandid, genomeDB) {
  islandid <- getIsland(entrezid=entrezid, genomeDB=genomeDB)
  getChr(islandid=islandid, genomeDB=genomeDB)
}
)

setMethod("getChr",signature(entrezid='missing',txid='character',islandid='missing', genomeDB='annotatedGenome'), function(entrezid, txid, islandid, genomeDB) {
  islandid <- getIsland(txid=txid, genomeDB=genomeDB)
  getChr(islandid=islandid, genomeDB=genomeDB)
}
)

setMethod("transcripts", signature(genomeDB='annotatedGenome',txid='missing',islandid='character'), function(genomeDB, txid, islandid) {
  IRangesList(lapply(genomeDB@transcripts[[islandid]],function(z) ranges(genomeDB@islands[[islandid]][as.character(z)])))
}
)

setMethod("transcripts", signature(genomeDB='annotatedGenome',txid='character',islandid='missing'), function(genomeDB, txid, islandid) {
  islandid <- getIsland(txid=txid,genomeDB=genomeDB)
  transcripts(genomeDB=genomeDB, islandid=islandid)[txid]
}
)

setMethod("transcripts", signature(genomeDB='annotatedGenome', txid='missing',islandid='missing'), function(genomeDB, txid, islandid) {
  tx <- unlist(genomeDB@transcripts,recursive=FALSE)
  names(tx) <- unlist(sapply(genomeDB@transcripts,names))
  tx <- data.frame(tx=rep(names(tx),sapply(tx,length)), exon=unlist(tx))
  txranges <- genomeDB@exonsNI[as.character(tx$exon)]
  names(txranges) <- NULL
  ans <- GRanges(ranges=IRanges(start(txranges),end(txranges)), space=as.character(tx$tx), seqnames= seqnames(txranges))
  #ans <- RangedData(IRanges(start(txranges),end(txranges)), space=as.character(tx$tx), chr=as.character(seqnames(txranges)))  #old version, superceeded
  return(ans)
}
)


setGeneric("txLength", function(islandid, txid, genomeDB) standardGeneric("txLength"))
setMethod("txLength", signature(islandid='missing', txid='missing', genomeDB='annotatedGenome'), function(islandid, txid, genomeDB) {
  txid <- rownames(genomeDB@aliases)
  txLength(txid=txid, genomeDB=genomeDB)
}
)

setMethod("txLength", signature(islandid='character', txid='missing', genomeDB='annotatedGenome'), function(islandid, txid, genomeDB) {
  txid <- rownames(genomeDB@aliases)[genomeDB@aliases$island_id==islandid]
  txLength(txid=txid, genomeDB=genomeDB)
}
)

setMethod("txLength", signature(islandid='missing', txid='character', genomeDB='annotatedGenome'), function(islandid, txid, genomeDB) {
  if(length(genomeDB@txLength)==0){
    tx <- unlist(genomeDB@transcripts,recursive=FALSE)
    nn <- strsplit(names(tx),'\\.')
    nn <- sapply(nn, function(z) paste(z[-1],collapse='.'))
    names(tx) <- nn
    #names(tx) <- sapply(strsplit(names(tx),'\\.'),'[[',2)
    tx <- tx[txid]
    tx <- data.frame(tx=rep(names(tx),sapply(tx,length)), exon=unlist(tx))
    #Get exon id & width into data.frame
    w <- unlist(width(genomeDB@islands))
    r <- unlist(ranges(genomeDB@islands))
    names(w) <- sapply(strsplit(names(r),'\\.'),'[[',2)
    w <- data.frame(exon=as.integer(names(w)), width=w)
    #Merge and by
    exonw <- merge(tx,w,by='exon',all.x=TRUE)
    txw <- by(exonw$width, INDICES=factor(exonw$tx), FUN=sum)
    ans <- as.integer(txw); names(ans) <- names(txw)
  } else ans <- genomeDB@txLength[txid]
  return(ans)
}
)


setMethod("txLength", signature(islandid='missing', txid='data.frame', genomeDB='annotatedGenome'), function(islandid, txid, genomeDB) {
  if(length(genomeDB@txLength)==0){
    tx <- unlist(genomeDB@transcripts,recursive=FALSE)
    names(tx) <- sapply(strsplit(names(tx),'\\.'),'[[',2)
    tx <- tx[txid]
    tx <- data.frame(tx=rep(names(tx),sapply(tx,length)), exon=unlist(tx))
    ans <- txLengthBase(tx=tx, genomeDB=genomeDB)
  } else ans <- genomeDB@txLength[txid]
  return(ans)
}
)

txLengthBase <- function(tx, genomeDB) {
  #Get exon id & width into data.frame
  w <- unlist(width(genomeDB@islands))
  r <- unlist(ranges(genomeDB@islands))
  names(w) <- sapply(strsplit(names(r),'\\.'),'[[',2)
  w <- data.frame(exon=as.integer(names(w)), width=w)
  #Merge and by
  exonw <- merge(tx,w,by='exon',all.x=TRUE)
  txw <- by(exonw$width, INDICES=factor(exonw$tx), FUN=sum)
  ans <- as.integer(txw); names(ans) <- names(txw)
  return(ans)
}



matchTranscripts <- function(queryDB, subjectDB, maxbp=10) {
  #Find transcripts in queryDB matching a transcript in subjectDB. The best match for each transcript in newDB is returned, unless difference is >maxbp
  # Ouput is data.frame with tx1, len1, tx2, len2, bpintersect. (tx2 minimizes len2-bpintersect)
  #Format as GRanges objects
  txnew <- transcripts(queryDB)
  chrnew <- txnew[['chr']]; txname <- space(txnew)
  txnew <- GRanges(ranges=IRanges(start(txnew),end(txnew)), seqnames=chrnew, tx_id=txname) 
  #txnew <- RangedData(IRanges(start(txnew),end(txnew)), space=chrnew, tx_id=txname) #old version, superceeded
  #
  txknown <- transcripts(subjectDB)
  chrknown <- txknown[['chr']]; txname <- space(txknown)
  txknown <- GRanges(ranges=IRanges(start(txknown),end(txknown)), seqnames=chrknown, tx_id=txname)
  #txknown <- RangedData(IRanges(start(txknown),end(txknown)), space=chrknown, tx_id=txname) #old version, superceeded
  #Find overlaps
  o <- as.matrix(findOverlaps(txnew, txknown))
  match1 <- txnew[o[,'queryHits'],]
  match2 <- txknown[o[,'subjectHits'],]
  #Count bp in common
  maxst <- ifelse(start(match1) <= start(match2),start(match2),start(match1))
  minen <- ifelse(end(match1) <= end(match2),end(match1),end(match2))
  bpintersect <- minen - maxst + 1
  matches <- data.frame(id=paste(match1[['tx_id']], match2[['tx_id']], sep='//'), bpintersect)
  bpintersect <- sqldf::sqldf("select id, sum(bpintersect) from matches group by id")
  #Format output and filter based on maxbp
  n <- strsplit(as.character(bpintersect$id), split='//')
  n1 <- sapply(n,'[[',1); n2 <- sapply(n,'[[',2)
  len1 <- txLength(genomeDB=queryDB)[n1]
  len2 <- txLength(genomeDB=subjectDB)[n2]
  ans <- data.frame(tx1=n1,len1=len1,tx2=n2,len2=len2,bpintersect=bpintersect[,'sum(bpintersect)'])
  maxlen <- ifelse(ans$len1>ans$len2, ans$len1, ans$len2)
  ans <- ans[ans$bpintersect >= maxlen - maxbp,]
  #If >1 matches, return best match
  tab <- table(as.character(ans$tx1))
  onematch <- ans[ans$tx1 %in% names(tab)[tab==1],]
  multimatch <- ans[ans$tx1 %in% names(tab)[tab>1],]
  bestmatch <- do.call(rbind,by(multimatch, INDICES=as.character(multimatch$tx1), FUN= function(z) z[which.min(z$len2-z$bpintersect),]))
  ans <- rbind(onematch,bestmatch)
  return(ans)
}


setGeneric("subsetGenome", function(islands, chr, genomeDB) standardGeneric("subsetGenome"))
setMethod("subsetGenome", signature(islands='character', chr='missing', genomeDB='annotatedGenome'), function(islands, chr, genomeDB) {
  islands <- unique(islands)
  isl <- subset(genomeDB@islands, names(genomeDB@islands) %in% islands)
  txs <- unlist(sapply(genomeDB@transcripts[islands], names))
  exs <- unique(names(genomeDB@islands[islands]@unlistData))
  alia <- genomeDB@aliases[txs,]
  ex2is <- genomeDB@exon2island[as.character(genomeDB@exon2island$island) %in% islands,]
  txs <- genomeDB@transcripts[islands]
  exs <- subset(genomeDB@exonsNI, names(genomeDB@exonsNI) %in% exs)
  new("annotatedGenome", islands=isl, transcripts=txs, exonsNI=exs, aliases=alia, exon2island=ex2is, dateCreated=genomeDB@dateCreated, denovo=genomeDB@denovo, genomeVersion=genomeDB@genomeVersion)
})
setMethod("subsetGenome", signature(islands='missing', chr='character', genomeDB='annotatedGenome'), function(islands, chr, genomeDB) {
  islands <- unique(as.character(genomeDB@exon2island[genomeDB@exon2island$seqnames %in% chr,]$island))
  subsetGenome(islands=islands, genomeDB=genomeDB)
})

### FUNCTIONS

makeIslands <- function(exons){
  txs <- as.integer(as.factor(names(exons)))
  totEx <- length(exons)
  uniex <- unique(exons)
  nexR <- length(uniex)
  islands <- rep(0, nexR)
  nexons <- names(exons)
  exons <- as.character(sprintf("%d",exons))
  names(exons) <- nexons
  tabex <- table(exons)
  tabex <- tabex[exons]
  tabtx <- table(names(exons))
  tabtx <- tabtx[names(exons)]
  ans<-.Call("makeGeneIslands", as.integer(exons), islands, uniex, txs, totEx, nexR, as.integer(tabex), as.integer(tabtx))
  names(ans) <- as.character(sprintf("%d",uniex))
  ans
}

generateNOexons<-function(exByTx, startId=1, mc.cores){
  exByTx<-unlist(exByTx)
  strand(exByTx) <- "*"
  exonsNI <- disjoin(exByTx)
  names(exonsNI) <- startId:(length(exonsNI)+startId-1)
  overEx <- findOverlaps(exonsNI, exByTx)
  exid <- names(exonsNI)[queryHits(overEx)]
  exkey <- split(exid, values(exByTx)$exon_id[subjectHits(overEx)])
  if(mc.cores>1 && requireNamespace("parallel", quietly=TRUE)) {
    exkey<-parallel::mclapply(exkey, unique, mc.cores=mc.cores)
  }
  else {
    exkey<-lapply(exkey, unique, mc.cores=mc.cores)
  }
  res<-list(exkey=exkey, exons=exonsNI)
  res 
}

genomeBystrand <- function(DB, strand){
  is <- as.character(strand(DB@islands@unlistData))[cumsum(c(1, elementNROWS(DB@islands)[-length(DB@islands)]))]
  sel <- names(DB@islands)[is==strand]
  islands <- DB@islands[sel]
  transcripts <- DB@transcripts[sel]
  exonsNI <- DB@exonsNI[names(DB@exonsNI) %in% as.character(unlist(transcripts)),]
  exon2island <- DB@exon2island[rownames(DB@exon2island) %in% as.character(unlist(transcripts)),]
  txid <- unlist(lapply(transcripts, names))
  aliases <- DB@aliases[DB@aliases$tx %in% txid,]
  ans <- new("annotatedGenome", aliases=aliases, denovo=TRUE, exonsNI=exonsNI, transcripts=transcripts, exon2island=exon2island, dateCreated=Sys.Date(), genomeVersion=DB@genomeVersion, islands=islands)
  ans
}


rmDuplicateTxs <- function(txs, txname, exid) {
  txs <- txs[match(unique(mcols(txs)[,txname]), mcols(txs)[,txname]),]
  exid <- exid[mcols(txs)[,txname]]
  txs <- txs[match(unique(exid), exid),] #keep single tx per identical exon sequence
  aliases <- values(txs)[c('tx_id',txname,'gene_id')] 
  aliases[,'gene_id'] <- as.character(aliases[,'gene_id']) 
  aliases <- as.data.frame(aliases)
  aliases <- cbind(aliases, exid=exid[as.character(aliases[,txname])]) #add sequence of exons to aliases
  reptxs <- tapply(names(exid), INDEX=exid, FUN=paste, collapse=',')
  reptxs <- reptxs[as.character(aliases$exid)]
  names(reptxs) <- aliases$exid
  aliases <- cbind(aliases, tx=reptxs)
  #old code
  #alitx <- split(as.character(aliases[,2]), aliases[,4])
  #names(alitx) <- unlist(lapply(alitx, function(x) ifelse(length(x)>1, x[x %in% txs@elementMetadata[,txname]], x)))
  #nalitx <- rep(names(alitx), unlist(lapply(alitx, length)))
  #alitx <- unlist(alitx)
  #names(alitx) <- nalitx
  #aliases <- cbind(aliases, tx=names(alitx)[match(aliases[,2], alitx)])
  return(list(txs=txs,aliases=aliases))
}


createGenome <- function(txs, Exons, genome, mc.cores, verbose=TRUE) {
  exid <- sapply(txs@elementMetadata$exon_id, function(x) paste(unlist(x), collapse="."))
  names(exid) <- values(txs)[,'tx_name']
  txs <- rmDuplicateTxs(txs, txname='tx_name', exid=exid)
  aliases <- txs$aliases; txs <- txs$txs
  Exons<-Exons[names(Exons) %in% txs@elementMetadata$tx_id,]
  txnames<-match(names(Exons), txs@elementMetadata$tx_id)
  names(Exons)<-txs@elementMetadata$tx_name[txnames]
  if (verbose) cat("Finding non-overlapping exons\n")
  txStrand <- as.character(strand(txs))
  names(txStrand) <- mcols(txs)[,"tx_name"] #values(txs)$tx_name
  exonsNI <- generateNOexons(Exons, startId=1, mc.cores=mc.cores)  
  exkey <- exonsNI$exkey
  exonsNI <- exonsNI$exons
  #Find transcript structure for new exons
  if (verbose) cat("Remapping transcript structure to new exons\n")
  s<-ifelse(as.character(Exons@unlistData@strand)=="+", 1, -1)
  exids <- values(Exons@unlistData)$exon_id
  idx <- values(Exons@unlistData)$exon_rank
  idx <- cumsum(idx==1)
  exids <- s*exids
  exon_ids <- tapply(exids, idx, function(x) if(any(x<0)) {rev(-x)} else x)
  names(exon_ids) <- names(Exons)
  exlen <- sapply(exkey, length)
  newTxs <- as.integer(unlist(exkey[as.character(sprintf("%d",unlist(exon_ids)))]))
  lenkey <- exlen[as.character(sprintf("%d",unlist(exon_ids)))]
  lenkey <- tapply(lenkey, rep(names(exon_ids), sapply(exon_ids, length)), sum)
  lenkey <- lenkey[names(exon_ids)]
  newTxs <- split(newTxs, rep(names(lenkey), lenkey))
  newTxs <- newTxs[names(Exons)]
  geneids <- mcols(txs)[,"gene_id"] #values(txs)$gene_id
  #geneids <- unlist(geneids)
  names(geneids)<- mcols(txs)[,"tx_id"] #unlist(values(txs)$tx_id)
  # Make islands
  ex2tx <- unlist(newTxs)
  names(ex2tx) <- rep(names(newTxs), unlist(lapply(newTxs, length)))
  islands <- makeIslands(ex2tx)
  if (verbose) cat("Splitting transcripts\n")
  extxs <- unlist(lapply(newTxs, "[", 1))    
  sel <- match(extxs, names(exonsNI))
  tx2island <- islands[as.character(sel)]
  names(tx2island) <- names(newTxs)
  transcripts <- newTxs[names(tx2island)]
  sel <- txStrand[names(transcripts)]=='-'
  transcripts[sel] <- sapply(transcripts[sel], rev)
  transcripts <- split(transcripts, tx2island)
  islandStrand <- tapply(txStrand, tx2island[names(txStrand)], unique)
  ss <- sapply(islandStrand, length)
  islandStrand[ss>1] <- "*"
  id2tx <- data.frame(island_id=rep(names(transcripts),sapply(transcripts,length)) , txname=unlist(sapply(transcripts,names)))
  rownames(id2tx) <- id2tx$txname
  aliases$island_id <- id2tx[rownames(aliases),'island_id']
  exon2island <- as.data.frame(ranges(exonsNI))
  exon2island <- exon2island[,-4]
  exon2island$seqnames <- as.character(seqnames(exonsNI))
  rownames(exon2island) <- names(exonsNI)
  exon2island$island <- islands[rownames(exon2island)]
  exp <- exonsNI[unlist(islandStrand[as.character(exon2island$island)]) == '+']
  strand(exp) <- "+"
  tmp <- split(exp, islands[names(exp)])
  exu <- exonsNI[unlist(islandStrand[as.character(exon2island$island)]) == '*']
  strand(exu) <- "*"
  tmpu <- split(exu, islands[names(exu)])
  exm <- exonsNI[unlist(islandStrand[as.character(exon2island$island)]) == '-']
  strand(exm) <- "-"
  tmpm <- split(rev(exm), islands[names(rev(exm))])
  tmp <- c(tmp, tmpm, tmpu)
  tmp <- tmp[names(transcripts)]
  ans <- new("annotatedGenome", islands=tmp, transcripts=transcripts, exon2island=exon2island, aliases=aliases, exonsNI=exonsNI, dateCreated=Sys.Date(), genomeVersion=genome, denovo=FALSE)
  txL <- txLength(genomeDB=ans)
  ans@txLength <- txL
  return(ans)
}


setMethod("procGenome", signature(genDB="TxDb"), function(genDB, genome, mc.cores=1, verbose=TRUE) {
  #  genDB<-makeTranscriptDbFromUCSC(genome=genome, tablename="refGene")
  if (verbose) cat("Processing Exons and Transcrips\n")
  txs <- GenomicFeatures::transcripts(genDB,columns=c("tx_id","tx_name","gene_id","exon_id"))
  Exons <- exonsBy(genDB, by="tx")
  createGenome(txs=txs, Exons=Exons, genome=genome, mc.cores=mc.cores, verbose=verbose)
} 
)

setMethod("procGenome", signature(genDB="GRanges"), function(genDB, genome, mc.cores=1, verbose=TRUE) {
 ## Store gene type to add to aliases

    if('gene_type' %in% colnames(values(genDB))){
        gt <- as.character(genDB$gene_type)
        names(gt) <- genDB$transcript_id
    }

      #keep original strand information to revert strand information
    genDB_org <- genDB
    
  #store the tx ids of those transcripts for which the strand information is missing
    tx_unknownStrand <- unique(as.character(as.data.frame(genDB[strand(genDB) == "*" & genDB$type=='transcript'])[,"transcript_id"]))
    
      # set all unknown strands to '+'
    strand(genDB[which(strand(genDB) == "*")]) <- "+"
    
    
  if (verbose) cat("Formatting GTF table with GenomicFeatures tools...\n")
  if (all(c("gene_id", "transcript_id") %in% colnames(mcols(genDB)))) {
    tables <- .prepareGTFTables(genDB, exonRankAttributeName=NA) #function copied from GenomicFeatures
  } else {
    stop("Columns named 'gene_id' and 'transcript_id' not found")
  }
  if (any(is.na(genDB$gene_id))) stop("Missing values in geneDB$gene_id are not allowed")
  if (any(is.na(genDB$transcript_id))) stop("Missing values in geneDB$transcript_id are not allowed")
  chroms <- unique(tables[["transcripts"]][["tx_chrom"]])
  chrominfo <- data.frame(chrom = chroms, length = rep(NA,length(chroms)))

  
  #Eliminate transcripts with unknown strands (usually transcripts with 1 exon)
  #sel <- (tables$transcripts$tx_strand %in% c('+','-'))
  #tables$transcripts <- tables$transcripts[sel,]
  #sel <- (tables$splicings$exon_strand %in% c('+','-'))  
  #tables$splicings <- tables$splicings[sel,]
  
  tables$genes <- tables$genes[tables$genes$tx_name %in% tables$transcripts$tx_name,]
  if (verbose) cat("Making TranscriptDb object...\n")
  txdb <- makeTxDb(transcripts=tables[["transcripts"]], splicings=tables[["splicings"]], genes=tables[["genes"]], chrominfo=chrominfo, reassign.ids=TRUE, )
  #txdb <- GenomicFeatures::makeTxDb(transcripts=tables[["transcripts"]], splicings=tables[["splicings"]], genes=tables[["genes"]], chrominfo=chrominfo, reassign.ids=TRUE, )
  txs <- GenomicFeatures::transcripts(txdb,columns=c("tx_id","tx_name","gene_id","exon_id"))
  Exons <- exonsBy(txdb, by="tx")
  #Remove transcripts with no exons
  noexons <- sapply(txs$exon_id,length)==0
  if (any(noexons)) {  
    warning('Some of the specified transcripts had no exons. They were removed')
    txs <- txs[!noexons,]
}
  ans <- createGenome(txs=txs, Exons=Exons, genome=genome, mc.cores=mc.cores, verbose=verbose)
  if(length(tx_unknownStrand)>0) ans <- restoreUnknownStrand(ans, tx_unknownStrand, genDB_org, mc.cores)
  if('gene_type' %in% colnames(values(genDB)))  ans@aliases$gene_type <- gt[rownames(ans@aliases)]
  ans
})

restoreUnknownStrand <- function(ans, tx_unknownStrand, genDB, mc.cores)
{

  # get the island ids for every transcript with unknown strand   
    islWithTxUnknownStrand <- unique(as.character(ans@aliases[tx_unknownStrand,]$island_id))
  
  # if the transcript is in an island on its own then change the strand of the island to *

    txPerIsl <- elementNROWS(ans@islands[islWithTxUnknownStrand])
    
    islSingleTx <- islWithTxUnknownStrand[txPerIsl==1]
  
    if(length(islSingleTx)>0) strand(ans@islands[islSingleTx]@unlistData) <- "*"
  
  islWithTxUnknownStrand <- islWithTxUnknownStrand[!islWithTxUnknownStrand %in% islSingleTx]
  
  # if there are only transcripts with an unknown strand in an island according to the original annotation, then change back to "*"

    if(length(islWithTxUnknownStrand)>0){
        txs <- unlist(lapply(ans@transcripts[islWithTxUnknownStrand], names))
        str <- as.character(strand(genDB[match(txs, genDB$transcript_id),]))
        isl <- rep(islWithTxUnknownStrand, unlist(lapply(ans@transcripts[islWithTxUnknownStrand], length)))
        
        isall <- tapply(str, isl, function(x) all(x == "*"))
        islAllStarStrand <- names(isall)[isall]

        strand(ans@islands[islAllStarStrand]@unlistData) <- "*"
        islWithTxUnknownStrand <- islWithTxUnknownStrand[!islWithTxUnknownStrand %in% islAllStarStrand] 
    }


#    if ('parallel' %in% loadedNamespaces()) {  
#      islAllStarStrand <- unlist(parallel::mclapply(islWithTxUnknownStrand, function(x){
#          if(all(strand(genDB[which(genDB$transcript_id %in% names(ans@transcripts[[x]]))]) == "*"))
#              return(x)
#      }, mc.cores=mc.cores))
#    } else {
#      islAllStarStrand <- unlist(lapply(islWithTxUnknownStrand, function(x){
#          if(all(strand(genDB[which(genDB$transcript_id %in% names(ans@transcripts[[x]]))]) == "*"))
#              return(x)
#      }))
#    }
  
#    strand(ans@islands[islAllStarStrand]@unlistData) <- "*"
#    islWithTxUnknownStrand <- islWithTxUnknownStrand[!islWithTxUnknownStrand %in% islAllStarStrand]
#}
    

  # if the strand of the island is +, then leave it like that

    if(length(islWithTxUnknownStrand)>0){        
#        islAllPlusStrand <- tapply(islAllStarStrand, rep(islWithTxUnknownStrand, unlist(lapply(ans@islands[islWithTxUnknownStrand], length))), unique)
#        islAllPlusStrand <- names(islAllPlusStrand[islAllPlusStrand=='+'])

         txs <- unlist(lapply(ans@transcripts[islWithTxUnknownStrand], names))
         str <- as.character(strand(genDB[match(txs, genDB$transcript_id),]))
         isl <- rep(islWithTxUnknownStrand, unlist(lapply(ans@transcripts[islWithTxUnknownStrand], length)))
         ispl <- tapply(str, isl, function(x) all(x %in% c("+", "*")))
         islAllPlusStrand <- names(ispl)[ispl]

        islWithTxUnknownStrand <- islWithTxUnknownStrand[!islWithTxUnknownStrand %in% islAllPlusStrand]
    }
  
  # if the strand of the island is *, then check if all other transcripts in this island are (1) - or (2) not                  
  
  
  # (1) if all others are on the - strand, then change to - 
  # or in other words if none of the transcripts in the original annotation is on the '+' strand, then they were all either '-' or '*'
    if(length(islWithTxUnknownStrand)>0){
#        islAllMinusStrand <- tapply(islAllStarStrand, rep(islWithTxUnknownStrand, unlist(lapply(ans@islands[islWithTxUnknownStrand], length))), unique)
        #islAllMinusStrand <- names(islAllMinusStrand[islAllMinusStrand=='-'])

        txs <- unlist(lapply(ans@transcripts[islWithTxUnknownStrand], names))
        str <- as.character(strand(genDB[match(txs, genDB$transcript_id),]))
        isl <- rep(islWithTxUnknownStrand, unlist(lapply(ans@transcripts[islWithTxUnknownStrand], length)))
        ismi <- tapply(str, isl, function(x) all(x %in%  c("-", "*")))
        islAllMinusStrand <- names(ismi)[ismi]
        
        strand(ans@islands[islAllMinusStrand]@unlistData) <- "-"
        islWithTxUnknownStrand <- islWithTxUnknownStrand[!islWithTxUnknownStrand %in% islAllMinusStrand]
    }
       
  # (2) all remaining cases should be: the island has transcripts on both strands, even when excluding the transcript with unknown strand
  # --> change to *
    if(length(islWithTxUnknownStrand)>0)  strand(ans@islands[islWithTxUnknownStrand]@unlistData) <- "*"
  
  return (ans)
}

