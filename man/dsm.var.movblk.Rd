\name{dsm.var.movblk}
\alias{dsm.var.movblk}
\title{Variance estimation via parametric moving block bootstrap}
\usage{
  dsm.var.movblk(n.boot, dsm.object, pred.data,
    ds.uncertainty = FALSE,
    samp.unit.name = "Transect.Label", block.size,
    cell.size.name = NULL, cell.size = NULL, bpfile = NULL)
}
\arguments{
  \item{n.boot}{number of bootstrap resamples.}

  \item{dsm.object}{object returned from
  \code{\link{dsm.fit()}}.}

  \item{pred.grid}{a \code{data.frame} that holds
  prediction points, must have the correct columns for
  other environmental covariates.}

  \item{ds.uncertainty}{incorporate uncertainty in the
  detection function?}

  \item{samp.unit.name}{name sampling unit to resample
  (default 'Transect.Label').}

  \item{block.size}{number of segments in each block.}

  \item{cell.size.name}{name of the column in pred.data
  corresponding to the size of the prediction cells
  (default NULL, see below).}

  \item{cell.size}{scalar (or vector of length
  \code{nrow(pred.data)}) of the sizes of the prediction
  cells (default \code{NULL}, either \code{cell.size} or
  \code{cell.size.name} MUST be specified).}

  \item{bpfile}{path to a file to be used (usually by
  Distance) to get generate a progress bar (default
  \code{NULL} -- no file written).}
}
\description{
  Estimate the variance in abundance over an area using a
  moving block bootstrap. Two procedures are implemented,
  one incorporating detection function uncertainty, one
  not.
}

