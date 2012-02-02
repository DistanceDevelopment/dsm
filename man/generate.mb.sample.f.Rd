\name{generate.mb.sample.f}
\alias{generate.mb.sample.f}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generates vector of residulas to be mapped onto the data }
\description{
  For parametric bootstrapping purposes, residuals are needed, which are then mapped back onto data
  to form pseudo-values for resampling and refitting
  
}
\usage{
generate.mb.sample.f(k = required.num.blocks, len = len.block, 
                     which.blocks = bs.blocks, x = data, unit.info = unit.info, 
			   n.units = num.sampling.units)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{ Number of blocks }
  \item{len}{ Length of each block }
  \item{which.blocks}{ Blocks over which bootstrap will be performed }
  \item{x}{ Observed data }
  \item{unit.info}{ Details regardinng each sampling unit }
  \item{n.units}{ Number of segments into which survey was subdivided }
}
%\details{
%}
\value{
  Vector of residuals is returned.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Louise Burt, RUWPA, \email{louise@mcs.st-and.ac.uk} }
%\note{ ~~further notes~~ 
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }
