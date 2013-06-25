\name{dsm.cor}
\alias{dsm.cor}
\title{Check for autocorrelation in residuals}
\usage{
  dsm.cor(dsm.obj, Transect.Label = "Transect.Label",
    max.lag = 10,
    resid.type = c("deviance", "pearson", "scaled.pearson", "working", "response"),
    fun = cor, ylim = c(0, 1))
}
\arguments{
  \item{dsm.obj}{a fitted dsm object.}

  \item{Transect.Label}{label for the transect (default:
  \code{Transect.Label}). Using different labels can be
  useful when transects are split over geographical
  features.}

  \item{resid.type}{the type of residuals used, see
  \code{\link{residuals.gam}}.}

  \item{fun}{the function to use, by default
  \code{\link{cor}}, must take two column vectors as
  arguments.}

  \item{max.lag}{maximum lag to calulate at.}

  \item{ylim}{user defined limits in y direction.}
}
\value{
  a plot or a vector of \code{fun} applied at the lags.
}
\description{
  Once a DSM has been fitted to data, this function can be
  used to check for autocorrelation in the residuals.
}
\section{Details}{
  Note we assume that \code{Sample.Label} is in order
  within each transect, so sort data beforehand. The
  current iteration of this function will only plot
  correlations nicely, other things are up to you but you
  can get the function to return the data (by assigning the
  result to an object).
}
\examples{
library(Distance)
 library(dsm)

 data(mexdolphins)

 hr.model <- ds(mexdolphins$distdata, max(mexdolphins$distdata$distance), key = "hr", adjustment = NULL)
 mod1<-dsm(N~s(x,y), hr.model, mexdolphins$segdata, mexdolphins$obsdata)

 dsm.cor(mod1,resid.type="d",max.lag=9)
}
\author{
  David L. Miller
}

