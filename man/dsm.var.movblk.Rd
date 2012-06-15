\name{dsm.var.movblk}
\alias{dsm.var.movblk}
\title{Variance estimation via parametric moving block bootstrap}
\usage{
  dsm.var.movblk(dsm.object, pred.data, n.boot, block.size,
    off.set, ds.uncertainty = FALSE,
    samp.unit.name = "Transect.Label",
    progress.file = NULL, bs.file = NULL)
}
\arguments{
  \item{dsm.object}{object returned from
  \code{\link{dsm.fit()}}.}

  \item{pred.data}{a \code{data.frame} that holds
  prediction points, must have the correct columns for
  other environmental covariates. It also MUST have a
  column called \code{cell.area} which gives the area for
  each prediction cell}

  \item{n.boot}{number of bootstrap resamples.}

  \item{block.size}{number of segments in each block.}

  \item{off.set}{offset to be used in the model, see
  \code{\link{offsets}}.  Note that this is NOT logged.}

  \item{ds.uncertainty}{incorporate uncertainty in the
  detection function?  Note that this doesn't work for
  covariate models for the detection function at the
  moment.}

  \item{samp.unit.name}{name sampling unit to resample
  (default 'Transect.Label').}

  \item{progress.file}{path to a file to be used (usually
  by Distance) to generate a progress bar (default
  \code{NULL} -- no file written).}

  \item{bs.file}{path to a file to store each boostrap
  round (this stores all of the bootstrap results rather
  than just the summaries, enabling outliers to be detected
  and removed.}
}
\description{
  Estimate the variance in abundance over an area using a
  moving block bootstrap. Two procedures are implemented,
  one incorporating detection function uncertainty, one
  not.
}

