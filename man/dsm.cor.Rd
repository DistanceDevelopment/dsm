\name{dsm.cor}
\alias{dsm.cor}
\title{Check for autocorrelation in residuals}
\usage{
dsm.cor(dsm.obj, Transect.Label = "Transect.Label",
  Segment.Label = "Segment.Label", max.lag = 10,
  resid.type = "scaled.pearson", fun = cor, ylim = c(0, 1),
  subset = "all", ...)
}
\arguments{
  \item{dsm.obj}{a fitted dsm object.}

  \item{Transect.Label}{label for the transect (default:
  \code{Transect.Label}). Using different labels can be
  useful when transects are split over geographical
  features or when transects are surveyed multiple times.}

  \item{Segment.Label}{label for the segments (default:
  \code{Segment.Label}).The result of calling
  \code{order()} must make sense.}

  \item{resid.type}{the type of residuals used, see
  \code{\link{residuals.gam}} and
  \code{\link{residuals.gam}}. Defaults to
  \code{"scaled.pearson"} in the GAM case and
  \code{"normalized"} in the GAMM case (which are
  equivalent).}

  \item{fun}{the function to use, by default
  \code{\link{cor}}, must take two column vectors as
  arguments.}

  \item{max.lag}{maximum lag to calulate at.}

  \item{ylim}{user defined limits in y direction.}

  \item{subset}{which subset of the data should the
  correlation function be calculated on?}

  \item{...}{other options to pass to \code{plot}.}
}
\value{
a plot or a vector of \code{fun} applied at the lags.
}
\description{
Once a DSM has been fitted to data, this function can be
used to check for autocorrelation in the residuals.
}
\section{Details}{
  Within each \code{Transect.Label}, segments will be
  sorted according to their \code{Segment.Labels}. This may
  require some time to get right for your particular data.
  If one has multiple surveys where transects are
  revisited, for example, one may want to make
  \code{Transect.Label} a unique transect-survey id.
  Neither label need to be included in the model, they must
  just be present in the \code{$data} field in the model.
  This usually means that they have to be in the segment
  data passed to \code{dsm}.

  The current iteration of this function will only plot
  correlations nicely, other things are up to you but you
  can get the function to return the data (by assigning the
  result to an object).

  If there are NA values in the residuals then the
  correlogram will not be calculated. This usually occurs
  due to NA values in the covariates (so the smoother will
  not have fitted values there). Code like
  `any(is.na(dsm.obj$data))` might be helpful.
}
\examples{
library(Distance)
 library(dsm)

 data(mexdolphins)

 hr.model <- ds(mexdolphins$distdata, max(mexdolphins$distdata$distance),
                key = "hr", adjustment = NULL)
 mod1<-dsm(N~s(x,y), hr.model, mexdolphins$segdata, mexdolphins$obsdata)

 dsm.cor(mod1,resid.type="d",max.lag=9,Segment.Label="Sample.Label")
}
\author{
David L. Miller
}

