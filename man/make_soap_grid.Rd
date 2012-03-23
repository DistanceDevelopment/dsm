\name{make_soap_grid}
\alias{make_soap_grid}
\title{Create a knot grid for the internal part of a soap film smoother.}
\usage{
  make_soap_grid(bnd, n.grid)
}
\arguments{
  \item{bnd}{list with elements \code{x} and \code{y} which
  give the locations of the boundary vertices. The first
  and last elements should be the same.}

  \item{n.grid}{either one number giving the number of
  points along the \code{x} and \code{y} axes that should
  be used to create the grid, or a vector giving the number
  in the \code{x} direction, then \code{y} direction.}
}
\value{
  a list with elements \code{x} and \code{y}, containing
  the knot locations.
}
\description{
  This routine simply creates a grid of knots (in the
  correct format) to be used as in the "internal" part of
  the soap film smoother
}
\author{
  David L Miller
}

