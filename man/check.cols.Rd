\name{check.cols}
\alias{check.cols}
\title{Check column names exist}
\usage{
  check.cols(ddf.obj, segment.data, observation.data,
    strip.width, segment.area)
}
\arguments{
  \item{ddf.obj}{a ddf object from `mrds`}

  \item{segment.data}{segment data as defined in
  \code{\link{dsm}}}

  \item{observation.data}{observation data as defined in
  \code{\link{dsm}}}

  \item{strip.width}{strip width if strip transects are
  being used}

  \item{segment.area}{area of segments}
}
\value{
  nothing, but throws an error if something went wrong
}
\description{
  Internal function to check that supplied `data.frames`
  have the correct columns and checks that sample labels
  are all unique.
}
\author{
  David Lawrence Miller
}

