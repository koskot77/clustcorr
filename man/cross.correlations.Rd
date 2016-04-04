\name{cross.correlations}
\alias{cross.correlations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find most correlated vectors between two matricies 
}
\description{
Compute all pairs correlations for the input vectors from two matricies and return the best matches
}
\usage{
cross.correlations(m1, m2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m1,m2}{
     two sets of equal-sized vectors in rows of the matricies \code{m1} and \code{m2}
}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Returns a list of two data frames containing the best matches from \code{m1} to \code{m2} and \code{m2} to \code{m1}
}
%\references{
%}
\author{
Khristian Kotov
}
%\note{
%}

\seealso{
 \code{\link{setcores}}, \code{\link{cleanup}}, \code{\link{dist}}, \code{\link{hclust}}
}
\examples{
sample1 <- matrix( rnorm(1000), ncol=100 )
sample2 <- matrix( rnorm(1000), ncol=100 )
setcores(4)
cl <- cross.correlations(sample1,sample2)
cl[[1]]
cl[[2]]
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
