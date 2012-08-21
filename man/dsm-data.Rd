\name{dsm-data}
\alias{dsm-data}
\title{Data format for DSM}
\description{
  Two data frame must be provided to \code{\link{dsm.fit}}.
  They are referred to as \code{obsdata} and \code{segdata}
  (for observation and segment data, respectively.
}
\details{
  \code{obsdata} - the observation data frame must have the
  following columns: \tabular{ll}{\code{object} \tab unique
  object identifier \cr \code{Sample.Label} \tab the
  identifier for the segment that the observation occurred
  in \cr \code{size} \tab the size of each observed group
  (i.e. 1 for individuals) \cr \code{distance} \tab
  perpendicular/radial distance to observation}

  \code{segdata} - the segment data frame must have the
  following columns: \tabular{ll}{ \code{x} \tab centreline
  of the transect (i.e. "across the transect")\cr \code{y}
  \tab centre in the direction of the transect\cr
  \code{Effort} \tab the effort (in terms of length of the
  segment)\cr \code{Transect.Label} \tab identifier for
  this transect\cr \code{Sample.Label} \tab identifier for
  the segment (unique!)\cr \code{???} \tab any other
  environmental covariates that are required}
}

