\name{dsm.predict}
\alias{dsm.predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Produce estimated abundance at each cell in the prediction grid }
\description{
  Function that extends the fitted density surface model from the covered region
  to the entire study region courtesy of predicting with covariates available in
  each of the cells of the prediction grid.
}
\usage{
dsm.predict(gam.model = gam.model, newdata = NULL, field = FALSE, off = NULL, silent = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gam.model}{ Fitted response model object }
  \item{newdata}{ Values of covariates at each of the prediction grid cells }
  \item{field}{ Boolean that defines whether offset (area) is a field or a constant }
  \item{off}{ Offset }
  \item{silent}{ Boolean indicating whether computed sum should be printed, or just returned unprinted }
}
%\details{
%}
\value{
  List including
  \item{result }{predict.gam object}
  \item{call }{Call that led to the generation of the prediction object}
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Eric Rexstad, RUWPA, \email{ericr@mcs.st-and.ac.uk} }
\note{ 
  Presently, only density surface models fitted with \code{mgcv()} can be used for prediction
  with this incarnation of \code{dsm.predict()}.
}
\seealso{ \code{\link[mgcv]{predict.gam}} }
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{  }
