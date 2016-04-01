\name{setcores}
\alias{setcores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Enables multithreading
}
\description{
\code{setcores} defines maximum number of threads to be run in parallel by \code{cluster.correlations}
(the default is one).
}
\usage{
setcores(4)
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
cluster.correlations, setcores
}
\examples{
sample <- matrix( rnorm(1000), ncol=100 )
setcores(4)
cluster.correlations(sample,0.5)
}
