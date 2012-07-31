\name{summary.dsm.var}
\alias{summary.dsm.var}
\title{Summarize the variance of a density surface model}
\usage{
  \method{summary}{dsm.var} (object, alpha = 0.05,
    boxplot.coef = 1.5, bootstrap.subregions, ...)
}
\arguments{
  \item{object}{a \code{dsm.var} object}

  \item{alpha}{alpha level for confidence intervals}

  \item{boxplot.coef}{the value of \code{coef} used to
  calculate the outliers see \code{\link{boxplot}}.}

  \item{bootstrap.subregions}{list of vectors of logicals
  or indices for subregions for which variances need to be
  calculated (only for bootstraps (see
  \code{\link{dsm.var.prop}} for how to use subregions with
  variance propagation).}

  \item{\dots}{unused arguments for S3 compatibility}
}
\value{
  a summary object
}
\description{
  Gives a brief summary of a fitted \code{dsm} variance
  object.
}
\author{
  David L. Miller
}
\seealso{
  dsm.var.movblk dsm.var.prop
}

