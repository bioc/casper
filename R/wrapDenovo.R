# #temp 
# load("/Volumes/biostats/mstobbe/Kinases/RNASeq/alttranscripts/pG_genDB_denovoCuff_knownCompl_ren_2014.RData")
# load("/Volumes/640 GB/RNASeq/deNovo/pg_mm10known.RData")
# 
# #test data
# bamFiles <- "/Volumes/ga/mstobbe/RNAseq_nebreda_stracker/tophat_lane1152013/sorted_hits.bam"
# genomeDB <- pG_genDB_denovoCuff_knownCompl_ren
# readLength <- 101
# citype <- "asymp"
# runWrapKnownFirst <- TRUE
# genomeDB_mPrior <- pg_mm10known
# maxExons <- 40
# smooth=TRUE
# verbose <- TRUE
# method <- "auto"
# priorq_denovo <- 3
# exactMarginal <- TRUE
# mc.cores_denovo <- 2
# keep.multihits <- TRUE
# 
# load("/Volumes/640 GB/RNASeq/deNovo/sampleForSim_knownNew2013_relExpr.RData")
# 
# results <- wrapKnown(bamFile=bamFiles, verbose=verbose, seed=1, mc.cores.int=1, 
#                                         mc.cores=1, genomeDB=pG_genDB_denovoCuff_knownCompl_ren,readLength=readLength, rpkm=TRUE, 
#                                         priorq=2, priorqGeneExpr=2, citype="none", niter=10^3, 
#                                         burnin=100, keep.pbam=TRUE, keep.multihits=keep.multihits, chroms="chr13")  
# 
# o <- wrapDenovo(bamFiles=bamFiles, genomeDB=pG_genDB_denovoCuff_knownCompl_ren,mc.cores.int=2, readLength=101, genomeDB_mPrior=pg_mm10known)
# o_pbam <- wrapDenovo(bamFiles=bamFiles, genomeDB=pG_genDB_denovoCuff_knownCompl_ren,mc.cores.int=2, readLength=101, genomeDB_mPrior=pg_mm10known, output_wrapKnown=results)

wrapDenovo <- function(bamFiles, verbose=TRUE, seed=1, mc.cores.int=1, mc.cores=1, targetGenomeDB, 
                       readLength, rpkm=TRUE, priorq=2, priorqGeneExpr=2, citype="none", 
                       niter=10^3, burnin=100, returnPbam = FALSE, keep.multihits=TRUE, chroms=NULL,
                       output_wrapKnown, keepPbamInMemory=FALSE, knownGenomeDB,
                       islandid, priorq_denovo=3,method="auto", niter_denovo,
                       exactMarginal=TRUE,  integrateMethod = "plugin", mc.cores_denovo=1,maxExons=40, smooth=TRUE)
{
    # parameter checks to prevent that mistakes are noticed half way
    if(missing(bamFiles))
        stop("Please provide bam files.")
    if(missing(knownGenomeDB))
        stop("No knownGenomeDB found with current annotation")
    if(missing(targetGenomeDB))
        stop("No targetGenomeDB found with predicted annotation")
    if(missing(readLength))
        stop("readLength must be specified")
    if(!(citype %in% c('none', 'asymp', 'exact')))
        stop("citype should take the value 'none', 'asymp', or 'exact'")    
    if (class(knownGenomeDB)!='annotatedGenome') 
        stop("knownGenomeDB must be of class 'annotatedGenome'")
    if (class(targetGenomeDB)!='annotatedGenome') 
      stop("targetGenomeDB must be of class 'annotatedGenome'")
    if (!(method %in% c('auto','rwmcmc','priormcmc','allmodels','submodels'))) 
        stop("method must be auto, rwmcmc, priormcmc, allmodels or submodels")
    if(missing(niter_denovo))
        niter_denovo <- NULL
    if(missing(islandid))
        islandid <- NULL
    
    pcs_all_samples <- c()
    distr_all_samples <- c()
    pbams_all_samples <- c()
    
    # Compute (or collect) distributions and pathcounts based provided genome with only the current/known annotation
    if(missing(output_wrapKnown))
    {
        if(verbose)
          cat("Computing distributions and pathcounts\n")
        
        for(b in 1:length(bamFiles))
        {          
            cur_distrsAndPCs <- getDistrsAndPCs(bamFiles[b], verbose=verbose, seed=seed, mc.cores.int=mc.cores.int, 
                                                 mc.cores=mc.cores, genomeDB=knownGenomeDB, keep.pbam=keepPbamInMemory, 
                                                 keep.multihits=keep.multihits, chroms=chroms)
            
            pcs_all_samples <- c(pcs_all_samples, cur_distrsAndPCs$pc)
            distr_all_samples <- c(distr_all_samples, cur_distrsAndPCs$distr)
            
            if(keepPbamInMemory)
                pbams_all_samples <- c(pbams_all_samples, cur_distrsAndPCs$pbam)
            
            rm(cur_distrsAndPCs)
        }
    }
    else
    {          
        for(i in 1:length(output_wrapKnown))  
        {
            pcs_all_samples <- c(pcs_all_samples, output_wrapKnown[i]$pc)
            distr_all_samples <- c(distr_all_samples, output_wrapKnown[i]$distr)
            
            if("pbam" %in% names(output_wrapKnown))
            {
                pbams_all_samples <- c(pbams_all_samples, output_wrapKnown$pbam)
            }
            else if(keepPbamInMemory)
            {
                keepPbamInMemory <- FALSE
                warning("keepPbamInMemory set to FALSE as pBams were not provided")
            }
        }
                
        rm(output_wrapKnown)
    }
      
    # merge distributions and pathcounts
    if(length(pcs_all_samples) > 1)
    {
        if(verbose)
          cat("Merging distributions and pathcounts\n")
      
        mergedPCs_allSamples <- mergePCs(pcs=pcs_all_samples, knownGenomeDB, mc.cores)    
        mergedDistr_allSamples <- suppressWarnings(mergeDisWr(distr_all_samples, pcs_all_samples))  
    }
    else
    {
        mergedPCs_allSamples <- pcs_all_samples[[1]]
        mergedDistr_allSamples <- distr_all_samples[[1]]
    }
  
    rm(pcs_all_samples)
    rm(distr_all_samples)
    gc()
    
    # run calcDenovo based on the merged distributions and pathcounts of all samples and the provided genome
    # with only the current/known annotation
    mprior <- modelPrior(genomeDB=knownGenomeDB, maxExons=maxExons, smooth=smooth, verbose=verbose)            

    if(verbose)
      cat("Running calcDenovo\n")
    
    out_calcDenovo <- calcDenovo(distrs=mergedDistr_allSamples, targetGenomeDB=targetGenomeDB, knownGenomeDB=knownGenomeDB, pc=mergedPCs_allSamples, 
                                   readLength=readLength, priorq=priorq_denovo, mprior=mprior, minpp=0, selectBest=FALSE, 
                                   method=method, exactMarginal=exactMarginal, verbose=verbose, integrateMethod, niter=niter_denovo,
                                   mc.cores=mc.cores_denovo, islandid= islandid)
    
     # Build denovo genome based on the output of calcDenovo
    
    if(verbose)
      cat("Constructing denovo genome object\n")
    
    genomeDB_denovo <- constructDenovoGenomeObj(vars_info=variants(out_calcDenovo), genomeDB=targetGenomeDB, mc.cores=max(c(mc.cores_denovo, mc.cores.int, mc.cores)))
    
    # Run wrapKnown to get expression estimates, distributions and pathcounts based on the denovo genome
    
    if(verbose)
      cat("Running wrapKnown with denovo Genome\n")
    
    if(!keepPbamInMemory)
    {
        out_wrapKnown_denovoGenome <- wrapKnown(bamFile=bamFiles, verbose=verbose, seed=seed, mc.cores.int=mc.cores.int, 
                                                   mc.cores=mc.cores, genomeDB=genomeDB_denovo,readLength=readLength, rpkm=rpkm, 
                                                   priorq=priorq, priorqGeneExpr=priorqGeneExpr, citype=citype, niter=niter, 
                                                   burnin=burnin, keep.pbam=returnPbam, keep.multihits=keep.multihits, chroms="chr13")  
     }
     else
     {
         out_wrapKnown_denovoGenome <- wrapKnownStartFromPBam(pbams=pbams_all_samples, verbose=verbose, seed=seed, mc.cores.int=mc.cores.int, 
                                                              mc.cores=mc.cores, genomeDB=genomeDB_denovo,readLength=readLength, rpkm=rpkm, 
                                                              priorq=priorq, priorqGeneExpr=priorqGeneExpr, citype=citype, niter=niter, 
                                                              burnin=burnin, keep.pbam=returnPbam, keep.multihits=keep.multihits, chroms="chr13")  
     }
  
    # Prepare output
    pc_denovo <- out_wrapKnown_denovoGenome$pc
    distr_denovo <- out_wrapKnown_denovoGenome$distr
    
    if(returnPbam)
      pbam_denovo <- out_wrapKnown_denovoGenome$pbam
    else
      pbam_denovo <- NULL
    
    eset_denovo <- out_wrapKnown_denovoGenome$exp    
    fData(eset_denovo)[,"newTx"] <- substr(fData(eset_denovo)$transcript, start = 1, stop = nchar("CASP")) == "CASP"

    # get probabilities of each variants being expressed
    out_relExpr <- relativeExpr(out_calcDenovo, summarize="modelAvg", minProbExpr=0, minExpr=0)
    fData(eset_denovo) <- merge(fData(eset_denovo), out_relExpr, by=0, all=TRUE, sort=FALSE)
    rownames(fData(eset_denovo)) <- fData(eset_denovo)[,1]
    fData(eset_denovo) <- fData(eset_denovo)[,-1]
    fData(eset_denovo) <- fData(eset_denovo)[,-grep("expr", colnames(fData(eset_denovo)))]

    return(list(genomeDB_denovo=genomeDB_denovo,eset_denovo=eset_denovo, pc_denovo=pc_denovo, distr_denovo=distr_denovo,pbam_denovo=pbam_denovo))
}

# merge pathcounts from multiple samples
mergePCs <- function(pcs, genomeDB, mc.cores=1)
{
  allisl <- unique(unlist(lapply(pcs, function(x) names(x@counts[[1]]))))
  alltxs <- parallel::mclapply(pcs, function(x) unlist(unname(x@counts[[1]])), mc.cores=mc.cores)
  alltxsu <- unlist(unname(alltxs), recursive=F)
  tmp <- by(alltxsu, names(alltxsu), sum)
  tmpn <- names(tmp)
  tmp1 <- as.numeric(tmp)
  names(tmp1) <- tmpn
  counts <- splitPaths(paths=tmp1, DB=genomeDB, mc.cores=mc.cores, stranded=pcs[[1]]@stranded, geneid=allisl)
  new("pathCounts", counts=counts, denovo=pcs[[1]]@denovo, stranded=pcs[[1]]@stranded)
}

# Construct genomeDB object based on the denovo predictions
constructDenovoGenomeObj <- function(vars_info, genomeDB, mc.cores)
{
    # build genome GRanges object to be used by procGenome    
    genomeDB_denovo_GRanges <- do.call(c, parallel::mclapply(1:nrow(vars_info),function(x)
                      {
                            cur_island <- as.character(vars_info[x,"islandid"])
                            
                            # get exons
                            cur_exons <- unlist(strsplit(x=as.character(vars_info[x,"exons"]), split=","))
                            
                            # get coordinates of each of the exons
                            exonCoor <- genomeDB@islands[[cur_island]][which(names(genomeDB@islands[[cur_island]]) %in% cur_exons),]
                              
                            # get strand information
                            strand(exonCoor) <- determineStrand(cur_exons, unique(as.character(strand(genomeDB@islands[[cur_island]]))))
                            
                            # get gene id
                            if(vars_info[x,"varName"] %in% genomeDB@aliases[,"tx_name"])
                            {
                                gene_id <- unique(as.character(genomeDB@aliases[which(genomeDB@aliases[,"tx_name"] == vars_info[x,"varName"]),"gene_id"]))
                                sourceTx <- "orgAnnotation" 
                            }
                            else
                            {
                                gene_id <- paste(unique(as.character(genomeDB@aliases[which(genomeDB@aliases[,"island_id"] == cur_island),"gene_id"])), collapse="__")
                                sourceTx <- "denovo_transcript" 
                            }
                                              
                            mcols(exonCoor)[,"gene_id"] <- gene_id
                            mcols(exonCoor)[,"transcript_id"] <- vars_info[x,"varName"] 
                            mcols(exonCoor)[,"source"] <- as.factor(sourceTx)
                            mcols(exonCoor)[,"type"] <- as.factor("exon")
                            
                            return(exonCoor)
                        }, mc.cores=mc.cores))
             
    genomeDB_denovo_GRanges <- deduceStrandSingleExonTx(genomeDB_denovo_GRanges, vars_info)
    
    names(genomeDB_denovo_GRanges) <- NULL
    
    save(genomeDB_denovo_GRanges, file="genomeDB_denovo_GRanges.RData")
    
    genomeDB_denovo<- procGenome(genDB=genomeDB_denovo_GRanges, genome="casper_denovo", mc.cores=mc.cores)
    
    return(genomeDB_denovo)
}

# Based on the order in which exons are mentioned determine the strand 
# if this is not possible, then use the strand of the whole island
determineStrand <- function(exons, strand_island)
{
    if(as.numeric(exons[1]) > as.numeric(exons[length(exons)]))
    {
        strand <- "-"
    }
    else if(as.numeric(exons[1]) < as.numeric(exons[length(exons)]))
    {
        strand <- "+"
    }
    else
    {
        strand <- strand_island 
    }
}

# For single exons we cannot determine the strand by the order of the exons
# If all transcripts of an island are on the same strand we already assigned this to the single exon transcript as well.
# If it is a mixed island, then the only case in which we can deduce the strand is if all other transcripts, i.e. excluding the
# single exon transcrip without strand information, are on the same strand.
deduceStrandSingleExonTx <- function(genomeDB_denovo, vars_info)
{
   tx_unknownStrand <- as.character(genomeDB_denovo[which(genomeDB_denovo@strand == "*"),]$transcript_id)
   
   for(i in 1:length(tx_unknownStrand))
   {
       curIsl <- as.character(vars_info[which(vars_info$varName == tx_unknownStrand[i]), "islandid"])
       txInIsl <- as.character(vars_info[which(vars_info$islandid == curIsl), "varName"])
       
       txInIsl <- txInIsl[which(!(txInIsl %in% tx_unknownStrand[i]))]
       
       strand_otherTx <- unique(as.character(strand(genomeDB_denovo[which(genomeDB_denovo$transcript_id %in% txInIsl),])))
       
       if(length(strand_otherTx) == 1)
       {
          if(strand_otherTx == "+")
          {
              strand(genomeDB_denovo[which(genomeDB_denovo$transcript_id %in% tx_unknownStrand[i]),]) <- "-"
          }
          else if(strand_otherTx == "-")
          {
              strand(genomeDB_denovo[which(genomeDB_denovo$transcript_id %in% tx_unknownStrand[i]),]) <- "+"
          }
       }       
   }
   
   return(genomeDB_denovo)
}

# Only the getDistr and pathCounts part of the wrapKnown function
getDistrsAndPCs <- function(bamFile, verbose, seed, mc.cores.int, mc.cores, genomeDB, keep.pbam, keep.multihits, chroms)
{
  if (is.null(genomeDB)) 
    stop("No genomeDB found")

  what <- c('qname','strand','pos','mpos','cigar')

  if(!keep.multihits) what <- c(what, 'mapq')
  t <- scanBamHeader(bamFile)[[1]][["targets"]]
  which <- GRanges(names(t), IRanges(1, unname(t)))
  if(!is.null(chroms)) {
    which <- which[as.character(seqnames(which)) %in% chroms]
  } else {
    which <- which[!grepl("_",as.character(seqnames(which)))]
    which <- which[!as.character(seqnames(which))=='chrM']
    sel <- as.vector(seqnames(which)) %in% names(seqlengths(genomeDB@islands))
    if (any(!sel)) warning(paste("Did not find in genomeDB chromosomes",paste(as.vector(seqnames(which))[!sel],collapse=' '),'. Skipping them'))
    which <- which[sel,]
  }
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
        rm(bam); gc()
        ans$distr <- getDistrs(DB=genomeDB, pbam=ans$pbam, verbose=FALSE)
        ans$pc <- pathCounts(reads=ans$pbam, DB=genomeDB, mc.cores=mc.cores, verbose=FALSE)
      } else {
        ans <- vector("list",2); names(ans) <- c("distr","pc")
        pbam <- procBam(bam=bam[[1]], stranded=FALSE, seed=as.integer(seed), verbose=FALSE, keep.junx=FALSE, rname=as.character(seqnames(which)[i]))
        rm(bam); gc()
        ans$pc <- pathCounts(reads=pbam, DB=genomeDB, mc.cores=mc.cores, verbose=FALSE)
        ans$distr <- getDistrs(DB=genomeDB, pbam=pbam, verbose=FALSE)
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
    if(keep.pbam)
    {
       allpbam <- lapply(ans, "[[", "pbam")
       names(allpbam) <-  seqnames(which)
    }
    allpc <- mergePCWr(ans, genomeDB)
    alldistr <- suppressWarnings(mergeDisWr(lapply(ans, '[[', 'distr'), lapply(ans, '[[', 'pc')))
    if(keep.pbam) {
      ans <- list(pc=allpc, distr=alldistr, pbam=allpbam)
    } else ans <- list(pc=allpc, distr=alldistr)
  } else {
    ans <- list(pc=NA, distr=NA)
  }
  return(ans)
}

# If the pBams were already generated we can skip this stap in wrapKnown
wrapKnownStartFromPBam <- function(pbams, verbose=FALSE, seed=1, mc.cores.int=1, mc.cores=1, genomeDB, readLength, rpkm=TRUE, priorq=2, priorqGeneExpr=2, citype='none', niter=10^3, burnin=100, keep.pbam=FALSE, keep.multihits=TRUE, chroms=NULL) 
{
  if (length(pbams)==1) {
    ans <- wrapKnownStartFromPBamSingle(pbam=pbams,verbose=verbose,seed=seed,mc.cores.int=mc.cores.int,mc.cores=mc.cores,genomeDB=genomeDB,readLength=readLength,rpkm=rpkm,priorq=priorq,priorqGeneExpr=priorqGeneExpr,citype=citype,niter=niter,burnin=burnin,keep.pbam=keep.pbam,keep.multihits=keep.multihits,chroms=chroms)
  } else if (length(pbams)>1) {
    x <- vector("list",length(pbams))
    for (i in 1:length(pbams)) {
      cat(paste("\n PROCESSING",pbams[i],"\n"))
      x[[i]] <- wrapKnownStartFromPBamSingle(pbam=pbams[i],verbose=verbose,seed=seed,mc.cores.int=mc.cores.int,mc.cores=mc.cores,genomeDB=genomeDB,readLength=readLength,rpkm=rpkm,priorq=priorq,priorqGeneExpr=priorqGeneExpr,citype=citype,niter=niter,burnin=burnin,keep.pbam=keep.pbam,keep.multihits=keep.multihits,chroms=chroms)
    }
    cat("\n MERGING ALL FILES...\n")
    if(keep.pbam)
    {
      ans <- vector("list",4); names(ans) <- c('pc','distr','exp', 'pbam')
      ans$exp <- mergeExp(lapply(x,'[[','exp'), sampleNames=names(pbams), keep=c('transcript','island_id','gene_id','explCnts'))
      ans$distr <- lapply(x,'[[','distr'); names(ans$distr) <- names(pbams)
      ans$pc <- lapply(x,'[[','pc'); names(ans$distr) <- names(pbams)
      ans$pbam <- lapply(x,'[[','pbam'); names(ans$pbam) <- names(pbams)
    }
    else
    {
      ans <- vector("list",3); names(ans) <- c('pc','distr','exp')
      ans$exp <- mergeExp(lapply(x,'[[','exp'), sampleNames=names(pbams), keep=c('transcript','island_id','gene_id','explCnts'))
      ans$distr <- lapply(x,'[[','distr'); names(ans$distr) <- names(pbams)
      ans$pc <- lapply(x,'[[','pc'); names(ans$distr) <- names(pbams)
    }
  } else {
    stop("Invalid length(bamFile)")
  }
  return(ans)
}

wrapKnownStartFromPBamSingle <- function(pbam, verbose=FALSE, seed=1, mc.cores.int=1, mc.cores=1, genomeDB, readLength, rpkm=TRUE, priorq=2, priorqGeneExpr=2, citype='none', niter=10^3, burnin=100, keep.pbam=FALSE, keep.multihits=TRUE, chroms=NULL) 
{
    
  ## Get distributions and pathcounts for one chromosome
  
  procChr <- function(i) {
    
    cat("Processing chromosome: ", names(pbam)[i], "\n")
    if (length(pbam[[i]]$pbam)>0) {  #if some reads satisfied the criteria
     
      if(keep.pbam) {
        ans <- vector("list",3); names(ans) <- c("pbam","distr","pc")
        ans$pbam <- pbam[i]
        ans$distr <- getDistrs(DB=genomeDB, pbam=ans$pbam, verbose=FALSE)
        ans$pc <- pathCounts(reads=ans$pbam, DB=genomeDB, mc.cores=mc.cores, verbose=FALSE)
      } else {
        ans <- vector("list",2); names(ans) <- c("distr","pc")
        ans$distr <- getDistrs(DB=genomeDB, pbam=pbam[i], verbose=FALSE)
        ans$pc <- pathCounts(reads=pbam[i], DB=genomeDB, mc.cores=mc.cores, verbose=FALSE)
      }
    } else {  #if no reads satisfied the criteria, return integer instead of list
      ans <- as.integer(-1)
    }
    cat("Finished chromosome ", names(pbam)[i], "\n")
    return(ans)
  }
  
  ## Run for all chromosomes, mclapply or for loop
  if (mc.cores.int>1 ){
    if ('parallel' %in% loadedNamespaces()) {
      ans <- parallel::mclapply(1:length(pbam), procChr, mc.cores=mc.cores.int, mc.preschedule=FALSE)
    } else stop('parallel library has not been loaded!')
  } else {
    ans <- list()
    for(i in 1:length(pbam)) ans[[i]] <- procChr(i)
  }
  
  gc()
  ans <- ans[sapply(ans,is.list)] #discard chromosomes with no reads
  
  if (length(ans)>0) 
  {
    if(keep.pbam) 
    {
       allpbam <- lapply(ans, "[[", "pbam")
       names(allpbam) <- seqnames(which)
    }
    
    allpc <- mergePCWr(ans, genomeDB)
    alldistr <- suppressWarnings(mergeDisWr(lapply(ans, '[[', 'distr'), lapply(ans, '[[', 'pc')))
    
    exp <- calcExp(distrs=alldistr, genomeDB=genomeDB, pc=allpc, readLength=readLength, rpkm=rpkm, priorq=priorq, priorqGeneExpr=priorqGeneExpr, citype=citype, niter=niter, burnin=burnin, mc.cores=mc.cores, verbose=verbose)
    
    if(keep.pbam) 
    {
      ans <- list(pc=allpc, distr=alldistr, exp=exp, pbam=allpbam)
    } 
    else 
      ans <- list(pc=allpc, distr=alldistr, exp=exp)
  } 
  else 
  {
    ans <- list(pc=NA, distr=NA, exp=NA)
  }
  return(ans)
}
