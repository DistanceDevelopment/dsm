\name{dsm.fit}
\alias{dsm.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit a density surface model to line transect data accounting for detectability }
\description{
  Function constructs and then invokes a call to \code{gam()} and returns the result of
  the fitting of the density surface model.
}
\usage{
dsm.fit(ddfobject = ddfobject, phat = NULL, response = response,
        formula = formula, model.defn = list(fn = "gam", family = "quasipoisson"),
	  obsdata = obsdata, segdata = segdata, wghts = NULL,
	  link = "log", off = NULL, convert.units = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ddfobject}{ Detection function object containing detection probability for each object detected in survey }
  \item{phat}{ Vector of estimated detection probabilities for each object (if not derived from fitted detection function) }
  \item{response}{ Nature of response variable (abundance, density, or count and groups or individuals) }
  \item{formula}{ Density surface model to be fitted }
  \item{model.defn}{ GAM or GLM along with error structure to be modelled }
  \item{obsdata}{ Data from the observation layer of the Distance project (sightings and covariates) }
  \item{segdata}{ Data from the subsample1 layer of Distance project (covariates for the segments) }
  \item{wghts}{ Optional weights to apply to observations (pertinent if density is response, and effective area is weight) }
  \item{link}{ Link function applied to the response variable }
  \item{off}{ Offset }
  \item{convert.units}{ Conversion factor to convert distances for the transect to distances for perpendicular distance measurement  }
}
%\details{
%}
\value{
  List consisting of:

  \item{result }{Result of the fitted density surface model}
  \item{call.dsm }{Match call of the request to \code{gam()}}
}
\references{ Hedley, S. and S. T. Buckland. 2004. Spatial models for line transect sampling. JABES 9:181-199. \cr
Wood, S.N. 2006. Generalized Additive Models: An Introduction with R. CRC/Chapman & Hall.}
\author{ Eric Rexstad, RUWPA, \email{ericr@mcs.st-and.ac.uk} }
\note{ 
Note that the gamma parameter to \code{gam()} is hardwired here; set to a value of 1.4 (from advice
in Wood (2006)) such that the \code{gam()} is inclined to not 'overfit.'
}
\seealso{ \code{\link[mgcv]{gam}} }
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{  }

