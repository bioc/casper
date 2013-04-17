\name{simReads}
\alias{simReads}
\title{
Function to simulate paired end reads following given read start and
fragment length distributions and gene and variant expressions.
}
\description{
This function generates path counts and bam files with simulated paired
end reads according to given read start distribution, fragment length
distribution and gene and variant expressions.
}
\usage{
simReads(islandid, nSimReads, pis, rl, seed, writeBam, distrs, genomeDB,
repSims=FALSE, bamFile=NULL, stranded=FALSE, verbose=TRUE, chr=NULL, mc.cores=1)
}
\arguments{
  \item{islandid}{
    Island ID's from the genomeDB object to simulate reads
  }
  \item{nSimReads}{
    Named numeric vector with number of fragments to
    simulate in each island.}
  \item{pis}{
    Named numeric vector with relative expression of
    transcripts. Expressions add up to one for each island to simulate.
    }
  \item{rl}{
    Read length
  }
  \item{seed}{
    Seed of the random numbers generator
}
  \item{writeBam}{
    Set to 1 to generate bam files with the simulated reads
}
  \item{distrs}{
    Object of class 'readDistrs' with read start and fragment length
    distributions
}
  \item{genomeDB}{
    Object of class 'annotatedGenome' with the genome to genererate
    reads from
}
\item{repSims}{
    Set to TRUE to return relative read starts and fragment lengths from the simulation
}
  \item{bamFile}{
    Name of the bam file to write reads to. Must end with '.bam' 
}
  \item{stranded}{
    Set to TRUE to preserve gene strand when generating reads. The 'XS'
    tag will be added to reads in the bam file and the returned 'pc'
    object will be stranded
}
\item{verbose}{Set to \code{TRUE} to print progress}
\item{chr}{
  Characters vector with chromosomes to simulate. Defaults to
  whole genome simulations.
  }
\item{mc.cores}{Number of cores to use in function}
}

\value{
  \item{Nsim}{Numerical vector with the number of reads simulated for
    each island.}
  \item{pc}{Object of class 'pathCounts' with simulated path counts}
  \item{sims}{
    Only if 'repSims' is set to TRUE. List with vectors of length 'n'
    with the following elements:
    -'varl': Length of variant for corresponding read
    -'st' Start of fragment relative to variant start (not in genomic
      coordinates)
    -len:Fragment length
    -'strand':Strand of gene for simulated read
  }
}
\author{
Camille Stephan-Otto Attolini
}

\examples{
data(hg19DB)
data(K562.r1l1)
distrs <- getDistrs(hg19DB,bam=K562.r1l1)

islandid <- c('10319','463')
txs <- unlist(lapply(hg19DB@transcripts[islandid], names))
pis <- vector(mode='numeric', length=length(txs))
npis <- sapply(hg19DB@transcripts[islandid],length)
pis[1:npis[1]] <- rep(1/npis[1],npis[1])
pis[-1:-npis[1]] <- rep(1/npis[2],npis[2])
names(pis) <- txs
nSimReads <- c(100, 100)
names(nSimReads) <- islandid

simpc <- simReads(islandid=islandid, nSimReads=nSimReads, pis=pis,
rl=75, repSims=TRUE, seed=1, writeBam=FALSE, distrs=distrs,genomeDB=hg19DB)

}
\keyword{ datagen }
