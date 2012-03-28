\name{dsm.var.prop}
\alias{dsm.var.prop}
\title{Variance propogation for DSM models}
\usage{
  dsm.var.prop(dsm.obj, pred.grids,
    pred.area = lapply(pred.grids, function(x) rep(1/nrow(x), nrow(x))),
    seglen.varname = "Effort", type.pred = "response")
}
\arguments{
  \item{dsm.obj}{an object returned from running
  \code{\link{dsm.fit}}.}

  \item{pred.grids}{list of prediction grids. Each element
  should be a \code{data.frame} with the same columns as
  the original data.}

  \item{pred.area}{a list with as many elements as there
  are in \code{pred.grids}. Each element is a vector as
  long as the number of rows in the corresponding element
  of \code{pred.grids}. These give the area associated with
  each prediction point. Default gives equal weighting to
  each point such that the "area" sums to one.}

  \item{seglen.varname}{name for the column which holds the
  segment length (default value "Effort").}

  \item{type.pred}{should the predictions be on the
  "response" or "link" scale?  (default "response").}
}
\value{
  a list with elements \tabular{ll}{\code{model} \tab the
  fitted model object\cr \code{pred.var} \tab covariances
  of the regions given in \code{pred.grids}. Diagonal
  elements are the variances in order.  }
}
\description{
  Rather than use a bootstrap to calculate the variance in
  a \code{dsm} model, use the clever variance propogation
  trick from Williams et al. (2011).
}
\details{
  The idea is to refit the spatial model but including the
  Hessian of the offset as an extra term. Variance
  estimates using this new model can then be used to
  calculate the variance of abundance estimates which
  incorporate detection function uncertainty. Further
  mathematical details are given in the paper in the
  references below.

  Many prediction grids can be supplied by supplying a list
  of \code{data.frame}s to the function.

  Based on (much more general) code from Mark Bravington
  and Sharon Hedley.
}
\author{
  Mark V. Bravington, Sharon L. Hedley. Bugs added by David
  L. Miller.
}
\references{
  Williams, R., Hedley, S.L., Branch, T.A., Bravington,
  M.V., Zerbini, A.N. and Findlay, K.P. (2011). Chilean
  Blue Whales as a Case Study to Illustrate Methods to
  Estimate Abundance and Evaluate Conservation Status of
  Rare Species. Conservation Biology 25(3), 526–535.
}

