\name{plot.dsm.var}
\alias{plot.dsm.var}
\title{Create plots of abundance uncertainty}
\usage{
  \method{plot}{dsm.var} (object, poly = NULL,
    limits = NULL, breaks = NULL, legend.breaks = NULL,
    xlab = "x", ylab = "y", observations = TRUE,
    plot = TRUE, ...)
}
\arguments{
  \item{object}{a \code{dsm.var} object}

  \item{poly}{a \code{list} or \code{data.frame} with
  columns \code{x} and \code{y}, which gives the
  coordinates of a polygon to draw.}

  \item{limits}{limits for the fill colours}

  \item{breaks}{breaks for the colour fill}

  \item{legend.breaks}{breaks as they should be displayed}

  \item{xlab}{label for the \code{x} axis}

  \item{ylab}{label for the \code{y} axis}

  \item{observations}{should observations be plotted?}

  \item{plot}{actually plot the map, or just return a
  \code{ggplot2} object?}

  \item{\dots}{any arguments that can usually be passed to
  a}
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

