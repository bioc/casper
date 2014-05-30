\name{relexprByGene}
\alias{relexprByGene}
\title{
  Compute relative expressions within each gene
}
\description{
  Transforms relative expressions that add up to 1 within each gene
  island (the default output of casper) to relative expressions that add
  up to 1 per gene.
}
\usage{
relexprByGene(x)
}
\arguments{
  \item{x}{\code{ExpressionSet} containing relative expressions.
    (typically, adding up to 1 for each island_id)
    Column \code{gene_id} in \code{fData(x)} should contain a unique gene identifier.}
}
\value{
  \code{ExpressionSet} with relative expressions adding up to one for
  each \code{gene_id}.
}
\author{
  David Rossell
}
\examples{
#See help(calcExp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }% __ONLY ONE__ keyword per line
