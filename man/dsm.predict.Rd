\name{dsm.predict}
\alias{dsm.predict}
\title{Produce estimated abundance at each cell in the prediction grid}
\usage{
  dsm.predict(gam.model, newdata = NULL, field = FALSE,
    off = NULL, silent = FALSE)
}
\arguments{
  \item{gam.model}{fitted model object (created by previous
  call to \code{gam()}/ \code{dsm.fit()})}

  \item{newdata}{spatially referenced (latitude and
  longitude) covariates e.g., sea temperature, depth,
  distance to shore, etc.  Known as
  'prediction.covariates.dat.r' in D6 DSM statement of
  direction prediction grid data frame.}

  \item{field}{logical indicating whether offset is
  provides as constant (\code{FALSE}) or field in grid
  dataframe (\code{TRUE}).}

  \item{off.set}{area of cells in prediction grid
  (\code{newdata}).  Provided by Distance Note covariates
  in this dataframe must have names *identical* to variable
  names used in fitting the \code{gam()} If entered between
  double quotes, will be a field named in prediction grid,
  otherwise a numeric constant Assume if a constant comes
  in, it will already have had the link transform applied,
  whereas if field, apply link fn.}

  \item{silent}{Set to true when used in conjunction with
  variance bootstraps; in this case we don't wish to see
  the abundance estimate printed This identical nature of
  the names in the model object and dataframe is not
  currently enforced; however, some code along the line of
  \code{for (j in 1:length(names(newd))){
  print(any(i<-grep(names(newd)[[j]],
  dimnames(attr(b$terms,"factors"))[[1]]))) }} would
  perform some checking of this

  Failure to have all covariates in \code{newdata} results
  in message of this sort: Error in eval(expr, envir,
  enclos) : object "x2" not found In addition: Warning
  message: not all required variables have been supplied in
  newdata!  in: predict.gam(b, newd, type = "response")}
}
\value{
  predicted list consisting of: \tabular{ll}{\code{result}
  \tab one-dimensional array of predicted values (without
  se's)\cr \code{call} \tab the call that generated this
  instance of the function}
}
\description{
  Function that extends the fitted density surface model
  from the covered region to the entire study region
  courtesy of predicting with covariates available in each
  of the cells of the prediction grid.
}
\note{
  Presently, only density surface models fitted with
  \code{mgcv()} can be used for prediction with this
  incarnation of \code{dsm.predict()}. Note the default
  arguments in \code{predict.gam()}
  \code{newdata.guaranteed=FALSE} and
  \code{na.action=na.pass} (creation of NA as predicted
  value where a covariate is missing) seem believable;
  leave them at their defaults, *but*
  \code{type="response"} is not the default argument, and
  needs to be specified.
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

