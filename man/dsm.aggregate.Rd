\name{dsm.aggregate}
\alias{dsm.aggregate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Sum over prediction grid cells to produce overall estimate of abundance }
\description{
  Produce aggregated abundance estimates for each specified layer.
}
\usage{
dsm.aggregate(predicted.values = predicted.values, prediction.layers = prediction.layers, path)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{predicted.values}{ Estimated abundance in each cell of the prediction grid }
  \item{prediction.layers}{ Matrix indicating membership of each cell into groups for each specified layer }
  \item{path}{ Temporary directory in which labels for classes in each layer are stored }
}
%\details{
%}
\value{
A table is printed containing the aggregated abundances for each classification in the specified layers.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Eric Rexstad, RUWPA, \email{ericr@mcs.st-and.ac.uk} }
%\note{ ~~further notes~~ 
%{
\seealso{\code{\link{dsm.predict}}}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

