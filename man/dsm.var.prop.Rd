\name{dsm.var.prop}
\alias{dsm.var.prop}
\title{Variance propogation for DSM models}
\usage{
  dsm.var.prop(dsm.obj, pred.data, off.set,
    seglen.varname = "Effort", type.pred = "response")
}
\arguments{
  \item{dsm.obj}{an object returned from running
  \code{\link{dsm.fit}}.}

  \item{pred.data}{the prediction grid. A \code{data.frame}
  with the same columns as the original data. If there is a
  column named \code{off.set}, it will be ignored in favour
  of the argument to this function.}

  \item{off.set}{a vector as long as the number of rows in
  \code{pred.data}.  These give the offset for each
  prediction point. Can also just be a scalar, if all
  entries are the same.}

  \item{seglen.varname}{name for the column which holds the
  segment length (default value "Effort").}

  \item{type.pred}{should the predictions be on the
  "response" or "link" scale?  (default "response").}
}
\value{
  a list with elements \tabular{ll}{\code{model} \tab the
  fitted model object\cr \code{pred.var} \tab variances in
  abundance of the the prediction region.  }
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

