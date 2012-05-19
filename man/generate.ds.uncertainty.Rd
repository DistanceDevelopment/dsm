\name{generate.ds.uncertainty}
\alias{generate.ds.uncertainty}
\title{Generate data from a fitted detection function}
\usage{
  generate.ds.uncertainty(ds.object)
}
\arguments{
  \item{ds.object}{a fitted detection function object (as
  returned by a call to \code{ddf.ds()}.}
}
\description{
  When \code{ds.uncertainty} is \code{TRUE}, this procedure
  generates data from the fitted detection function
  (assuming that it is correct).
}
\note{
  This function changes the random number generator seed.
  To avoid any potential side-effects, use something like:
  \code{seed <- get(".Random.seed",envir=.GlobalEnv)}
  \code{### Run some code}
  \code{assign(".Random.seed",seed,envir=.GlobalEnv)}
}
\author{
  David L. Miller
}

