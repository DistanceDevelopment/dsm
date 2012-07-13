\name{plot.dsm.var}
\alias{plot.dsm.var}
\title{Create plots of abundance uncertainty}
\usage{
  \method{plot}{dsm.var} (object, poly = NULL, ...,
    limits = NULL, xlab = "x", ylab = "y")
}
\arguments{
  \item{object}{a \code{dsm.var} object}

  \item{poly}{a \code{list} or \code{data.frame} with
  columns \code{x} and \code{y}, which gives the
  coordinates of a polygon to draw.}

  \item{\dots}{any arguments that can usually be passed to
  a}

  \item{limits}{limits for the fill colours
  \code{\link{plot}} method}
}
\value{
  a plot
}
\description{
  Note that the prediction data set must have \code{x} and
  \code{y} columns even if these were not used in the
  model.
}
\author{
  David L. Miller
}

