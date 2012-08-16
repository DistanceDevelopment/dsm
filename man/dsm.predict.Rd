\name{dsm.predict}
\alias{dsm.predict}
\title{Produce estimated abundance at each cell in the prediction grid}
\usage{
  dsm.predict(model, newdata = NULL, off.set = NULL)
}
\arguments{
  \item{model}{fitted model object (created by a call to
  \code{\link{gam}} or \code{\link{dsm.fit}})}

  \item{newdata}{spatially referenced covariates e.g.  sea
  temperature, depth, distance to shore, etc.  Note
  covariates in this dataframe must have names *identical*
  to variable names used in fitting the GAM.}

  \item{off.set}{offset for the cells in the prediction
  grid (\code{newdata}).  NB this should NOT be logged, it
  gets logged below.}
}
\value{
  array of predicted values.
}
\description{
  Function that extends the fitted density surface model
  from the covered region to the entire study region
  courtesy of predicting with covariates available in each
  of the cells of the prediction grid.
}
\note{
  Presently, only density surface models fitted with
  \code{\link{mgcv}} can be used for prediction with this
  incarnation of \code{dsm.predict}.
}
\author{
  Eric Rexstad, David L. Miller
}
\references{
  Hedley, S. and S. T. Buckland. 2004. Spatial models for
  line transect sampling. JABES 9:181-199.

  Wood, S.N. 2006. Generalized Additive Models: An
  Introduction with R. CRC/Chapman & Hall.
}
\seealso{
  predict.gam
}

