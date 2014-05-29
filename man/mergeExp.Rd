\name{mergeExp}
\alias{mergeExp}
\title{
  Merge splicing variant expression from multiple samples
}
\description{
  \code{mergeExp} combines the output of \code{calcExp} from multiple
  samples, i.e. multiple \code{ExpressionSet} objects, into a single \code{ExpressionSet}
}
\usage{
mergeExp(..., sampleNames,  keep=c('transcript','gene_id','island_id'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{\code{ExpressionSet} objects to be combined.}
  \item{sampleNames}{Character vector indicating the name of each
    sample. Defaults to 'Sample1', 'Sample2', etc.}
  \item{keep}{Variables in the featureData of each individual
    \code{ExpressionSet} to keep in the merged output.}
}
\details{
  \code{mergeExp} runs some checks to ensure that object can be combined
  (e.g. making sure that measurements are obtained on same set of
  genes), then sorts and formats each input
  \code{ExpressionSet}.

  A label with the sample name is appended to
  variables in the featureData that appear in multiple samples,
  e.g. variable 'se' reporting standard errors (obtained by setting
  citype='asymp' in calcExp).
}
\value{
  Object of class \code{ExpressionSet} combining the input
  ExpressionSets. Its \code{featureData} contains the columns indicated
  in the keep argument, plus a column \code{readCount} with the total
  number of reads mapped to each gene (or gene island, when multiple
  genes have overlapping exons).
}
\author{
  David Rossell
}
\seealso{
  \code{calcExp} to obtain an \code{ExpressionSet} for an individual sample.
}
\examples{
#See example in calcExp
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
