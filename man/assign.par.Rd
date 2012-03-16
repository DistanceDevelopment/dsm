\name{assign.par}
\alias{assign.par}
\title{Extraction and assignment of parameters to vector}
\usage{
  assign.par(ddfobj, fpar)
}
\arguments{
  \item{ddfobj}{distance sampling object (see
  \code{\link{create.ddfobj}})}

  \item{fpar}{parameter vector}
}
\value{
  index==FALSE, vector of parameters that were requested or
  index==TRUE, vector of 3 indices for scale, shape,
  adjustment
}
\description{
  Assigns parameters of a particular type (scale, shape,
  adjustments or g0 (p(0))) from the vector of parameters
  in \code{ddfobj}. All of the parameters are kept in a
  single vector for optimization even though they have very
  different uses.  \code{assign.par} parses them from the
  vector based on a known structure and assigns them into
  \code{ddfobj}.  \code{getpar} extracts the requested
  types to be extracted from \code{ddfobj}.
}
\note{
  Internal functions not intended to be called by user.
}
\author{
  Jeff Laake
}
\seealso{
  getpar
}
\keyword{utility}

