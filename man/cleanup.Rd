\name{cleanup}
\alias{cleanup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Release memory allocated in cluster.correlations
}
\description{
\code{cleanup} deletes matrix of all pairs correlations calculated
in \code{cluster.correlations} that is quadratic in size of its' input.
While this matrix is still stored in the memory, one can cheaply
recluster the input vectors with a different cutoff threshold.
Calling \code{cluster.correlations} on a different input will
automatically delete the old matrix.
}
\usage{
cleanup()
}
%- maybe also 'usage' for other objects documented here.
%\details{
%}
\value{
Always returns NULL
}
%\references{
%}
\author{
Khristian Kotov
}
%\note{
%  ~~further notes~~
%}

\seealso{
cluster.correlations, recluster.correlations, setcores, cross.correlations
}
\examples{
sample <- matrix( rnorm(1000), ncol=100 )
cluster.correlations(sample,0.5)
cleanup()
}
