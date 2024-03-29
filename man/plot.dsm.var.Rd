% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.dsm.var.R
\name{plot.dsm.var}
\alias{plot.dsm.var}
\title{Create plots of abundance uncertainty}
\usage{
\method{plot}{dsm.var}(
  x,
  poly = NULL,
  limits = NULL,
  breaks = NULL,
  legend.breaks = NULL,
  xlab = "x",
  ylab = "y",
  observations = TRUE,
  plot = TRUE,
  boxplot.coef = 1.5,
  x.name = "x",
  y.name = "y",
  gg.grad = NULL,
  ...
)
}
\arguments{
\item{x}{a \code{dsm.var} object}

\item{poly}{a \code{list} or \code{data.frame} with columns \code{x} and \code{y}, which gives
the coordinates of a polygon to draw. It may also optionally have a column
\code{group}, if there are many polygons.}

\item{limits}{limits for the fill colours}

\item{breaks}{breaks for the colour fill}

\item{legend.breaks}{breaks as they should be displayed}

\item{xlab}{label for the \code{x} axis}

\item{ylab}{label for the \code{y} axis}

\item{observations}{should observations be plotted?}

\item{plot}{actually plot the map, or just return a \code{ggplot2} object?}

\item{boxplot.coef}{control trimming (as in
\code{\link{summary.dsm.var}}), only has an effect if the bootstrap
file was saved.}

\item{x.name}{name of the variable to plot as the \code{x} axis.}

\item{y.name}{name of the variable to plot as the \code{y} axis.}

\item{gg.grad}{optional \code{\link[ggplot2:ggplot]{ggplot}} gradient object.}

\item{\dots}{any other arguments}
}
\value{
a plot
}
\description{
Note that the prediction data set must have \code{x} and \code{y} columns even if
these were not used in the model.
}
\section{Details}{


In order to get plotting to work with \code{\link{dsm_var_prop}} and
\code{\link{dsm_var_gam}}, one must first format the data correctly
since these functions are designed to compute very general summaries. One
summary is calculated for each element of the list \code{pred} supplied to
\code{\link{dsm_var_prop}} and \code{\link{dsm_var_gam}}.

For a plot of uncertainty over a prediction grid, \code{pred} (a \code{data.frame}),
say, we can create the correct format by simply using \code{pred.new <- split(pred,1:nrow(pred))}.
}

\seealso{
\code{\link{dsm_var_prop}}, \code{\link{dsm_var_gam}},
\code{\link{dsm_var_movblk}}
}
\author{
David L. Miller
}
