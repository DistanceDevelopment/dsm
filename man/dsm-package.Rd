\name{DSM-package}
\alias{DSM-package}
\alias{DSM}
\docType{package}
\title{
Density surface modelling of distance sampling data
}
\description{
Performs point and interval estimation of abundance, density or
counts from line transect data.  This is done by first segmenting transects into
sections (outside of this package).  The package then invokes GAMs and perhaps
GLMs for estimation purposes.  Variance estimation is via parametric moving block
bootstrapping.  This package is the analytical engine embedded in Distance 6 
with user interaction with the program via a VB interface.

}
\details{
\tabular{ll}{
Package: \tab DSM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2007-07-17\cr
License: \tab GPL\cr
}
}
\author{
Eric Rexstad, Research Unit for Wildlife Population Assessment
with functions written by Louise Burt, Research Unit for Wildlife Population Assessment

Maintainer: Eric Rexstad \email{ericr@mcs.st-and.ac.uk}
}
\references{
Hedley, S.L., and S.T. Buckland. 2004. Spatial models for line transect sampling. Journal of Agricultural, Biological, and Envirnmental Statistics 9:181-199.
}

\keyword{ package }
\seealso{
\url{http://www.ruwpa.st-and.ac.uk/distance/}
}
\examples{
Not a stand-alone R package, rather to be used in conjunction with Distance 6 interface
}
