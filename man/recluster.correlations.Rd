\name{recluster.correlations}
\alias{recluster.correlations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Groups correlated vectors in clusters
}
\description{
Uses in-memory matrix of all pairs correlations previously computed with \code{cluster.correlations}
and returns new clustering for a given new threshold
}
\usage{
recluster.correlations(m, cutoff)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
     a set of equal-sized vectors in rows of matrix \code{m}
}
  \item{cutoff}{
     correlation threshold at which clustering is performed (analogous to \code{h} in \code{cutree})
}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Returns an input-sized vector labeling the input vectors with cluster numbers.
  These labels index input vectors that seeded the clustering (they have no deep meaning).
}
%\references{
%}
\author{
Khristian Kotov
}
%\note{
%}

\seealso{
 \code{\link{setcores}}, \code{\link{cleanup}}, \code{\link{cluster.correlations}}, \code{\link{dist}}, \code{\link{hclust}}
}
\examples{
sample <- matrix( rnorm(1000), ncol=100 )
setcores(4)
cl <- cluster.correlations(sample,0.5)
unique(cl)
cl <- recluster.correlations(sample,0.1)
unique(cl)
cleanup()
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
