\name{dsm.var.movblk}
\alias{dsm.var.movblk}
\title{Variance estimation via parametric moving block bootstrap}
\usage{
  dsm.var.movblk(dsm.object, pred.data, n.boot, block.size,
    off.set, ds.uncertainty = FALSE,
    samp.unit.name = "Transect.Label",
    progress.file = NULL, bs.file = NULL, bar = TRUE)
}
\arguments{
  \item{dsm.object}{object returned from
  \code{\link{dsm}}.}

  \item{pred.data}{a \code{data.frame} that holds
  prediction points, must have the correct columns for
  other environmental covariates. It also MUST have a
  column called \code{cell.area} which gives the area for
  each prediction cell}

  \item{n.boot}{number of bootstrap resamples.}

  \item{block.size}{number of segments in each block.}

  \item{off.set}{offset to be used in the model, see
  \code{\link{offsets}}.  Note that this should NOT be
  \code{log()}'d.}

  \item{ds.uncertainty}{incorporate uncertainty in the
  detection function? See Details, below. Note that this
  feature is EXPERIMENTAL at the moment.}

  \item{samp.unit.name}{name sampling unit to resample
  (default 'Transect.Label').}

  \item{progress.file}{path to a file to be used (usually
  by Distance) to generate a progress bar (default
  \code{NULL} -- no file written).}

  \item{bs.file}{path to a file to store each boostrap
  round. This stores all of the bootstrap results rather
  than just the summaries, enabling outliers to be detected
  and removed. (Default \code{NULL}).}

  \item{bar}{should a progress bar be printed to screen?
  (Default \code{TRUE}).}
}
\description{
  Estimate the variance in abundance over an area using a
  moving block bootstrap. Two procedures are implemented,
  one incorporating detection function uncertainty, one
  not.
}
\section{Details}{
  Setting \code{ds.uncertainty=TRUE} will incorporate
  detection function uncertainty directly into the
  bootstrap. This is done by generating observations from
  the fitted detection function and then re-fitting a new
  detection function (of the same form), then calculating a
  new effective strip width. Rejection sampling is used to
  generate the observations (except in the half-normal
  case) so the procedure can be rather slow. Note that this
  is currently not supported with covariates in the
  detection function.

  Setting \code{ds.uncertainty=FALSE} will incorporate
  detection function uncertainty using the delta method.
  This assumes that the detection function and the spatial
  model are INDEPENDENT. This is probably not reasonable.
}
\examples{
\dontrun{
library(Distance)
library(dsm)

# load the Gulf of Mexico dolphin data (see ?mexdolphins)
data(mexdolphins)

# fit a detection function and look at the summary
hr.model <- ds(mexdolphins$distdata, max(mexdolphins$distdata$distance),
               key = "hr", adjustment = NULL)
summary(hr.model)

# fit a simple smooth of x and y
mod1<-dsm(N~s(x,y), hr.model, mexdolphins$segdata, mexdolphins$obsdata)
summary(mod1)

# create an offset (in metres)
# each prediction cell is 444km2
off.set <- 444*1000*1000

# calculate the variance by 500 moving block bootstraps
mod1.movblk <- dsm.var.movblk(mod1, mexdolphins$preddata, n.boot = 500,
   block.size = 3, samp.unit.name = "Transect.Label", off.set = off.set,
   bar = TRUE, bs.file = "mexico-bs.csv", ds.uncertainty = TRUE)
}
}

