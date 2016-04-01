\name{clustcorr-package}
\alias{clustcorr-package}
\alias{clustcorr}
\docType{package}
\title{
\packageTitle{clustcorr}
}
\description{
\packageDescription{clustcorr}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{clustcorr}
\packageIndices{clustcorr}
}
\author{
\packageAuthor{clustcorr}

Maintainer: \packageMaintainer{clustcorr}
}
%\references{
%}
\keyword{ package }
\seealso{
\code{\link[stats:stats-package]{stats}}
}
\examples{
sample <- matrix( rnorm(1000), ncol=100 )
setcores(4)
cl <- cluster.correlations(sample,0.5) 
unique(cl)
cl <- cluster.correlations(sample,0.1) 
unique(cl)
cleanup()
}
