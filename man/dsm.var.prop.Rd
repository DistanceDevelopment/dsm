\name{dsm.var.prop}
\alias{dsm.var.prop}
\title{Variance propogation for DSM models}
\usage{
  dsm.var.prop(dsm.obj, pred.data, off.set,
    seglen.varname = "Effort", type.pred = "response")
}
\arguments{
  \item{dsm.obj}{an object returned from running
  \code{\link{dsm}}.}

  \item{pred.data}{either: a single prediction grid or list
  of prediction grids. Each grid should be a
  \code{data.frame} with the same columns as the original
  data.}

  \item{off.set}{a a vector or list of vectors with as many
  elements as there are in \code{pred.data}. Each vector is
  as long as the number of rows in the corresponding
  element of \code{pred.data}. These give the area
  associated with each prediction point.}

  \item{seglen.varname}{name for the column which holds the
  segment length (default value "Effort").}

  \item{type.pred}{should the predictions be on the
  "response" or "link" scale?  (default "response").}
}
\value{
  a list with elements \tabular{ll}{\code{model} \tab the
  fitted model object\cr \code{pred.var} \tab covariances
  of the regions given in \code{pred.data}. Diagonal
  elements are the variances in order\cr \code{bootstrap}
  \tab logical, always \code{FALSE}\cr \code{pred.data}
  \tab as above\cr \code{off.set} \tab as above\cr
  \code{model}\tab the fitted model with the extra term\cr
  \code{dsm.object} \tab the original model, as above\cr
  \code{model.check} \tab simple check of subtracting the
  coefficients of the two models to see if there is a large
  difference\cr \code{deriv} \tab numerically calculated
  Hessian of the offset\cr.  }
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
  Rare Species. Conservation Biology 25(3), 526-535.
}

