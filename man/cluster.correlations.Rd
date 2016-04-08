\name{cluster.correlations}
\alias{cluster.correlations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Groups correlated vectors in clusters
}
\description{
Compute all pairs correlations for the input vectors and cluster them for
a given threshold.
}
\usage{
cluster.correlations(m, cutoff)
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
  Returns a list of clusters ordered by size in decreasing order. Each cluster is a list of
  indices of the input vectors, merged into this cluster.
}
%\references{
%}
\author{
Khristian Kotov
}
%\note{
%}

\seealso{
 \code{\link{setcores}}, \code{\link{cleanup}}, \code{\link{recluster.correlations}}, \code{\link{dist}}, \code{\link{hclust}}
}
\examples{
sample <- matrix( rnorm(1000), ncol=100 )
setcores(4)

# let's see how many clusters we get at 0.5 correlation threshold
cl <- cluster.correlations(sample,0.5)
length(cl)

# recluster with another threshold
cl <- recluster.correlations(sample,0.1)
length(cl)

# don't forget to release the memory if you don't plan to recluster the data again
cleanup()
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
