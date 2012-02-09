\name{get.dsm.fit.command.f}
\alias{get.dsm.fit.command.f}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Construct portion of calling sequence to gam() or glm()}
\description{
  Convenience function to assist in constructing R command to fit gam or glm to response variable.
}
\usage{
get.dsm.fit.command.f(model = model.details)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ String identifying model requested for error terms }
}
%\details{
%}
\value{
  Returns a string that is concatenated with other text to form call to \code{gam()} or \code{glm()}
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Eric Rexstad, RUWPA, \email{ericr@mcs.st-and.ac.uk} }
%\note{ ~~further notes~~ 
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{  }

