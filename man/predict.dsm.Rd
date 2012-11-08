\name{predict.dsm}
\alias{predict.dsm}
\title{Predict from a fitted density surface model}
\usage{
  \method{predict}{dsm} (object, newdata = NULL,
    off.set = NULL, type = "response", ...)
}
\arguments{
  \item{object:}{a fitted \code{\link{dsm}} object as
  produced by \code{gam()}.}

  \item{newdata}{spatially referenced covariates e.g.  sea
  temperature, depth, distance to shore, etc.  Note
  covariates in this dataframe must have names *identical*·
  to variable names used in fitting the GAM.}

  \item{off.set}{area of each of the cells in the
  prediction grid. Ignored if if their is already a column
  in \code{newdata} called \code{off.set}.  The link
  function is applied to the offset inside this function.}

  \item{type}{what scale should the results be on. The
  default is \code{"response"} and is almost always what
  you want.}

  \item{\dots}{any other arguments passed to
  \code{\link{predict.gam}} or \code{\link{predict.glm}}.}
}
\value{
  predicted values on the response scale
  (density/abundance).
}
\description{
  Make predictions outside (or inside) the covered area.
}
\author{
  David L. Miller
}

