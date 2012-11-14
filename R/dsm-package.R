#' Density surface modelling
#'
#' Some blurb will eventually go here.
#'
#'
#' @name dsm-package
#'
NULL


#' Offsets
#'
#' This will be documentation on calculating offsets.
#'
#' @name offsets
#'
NULL


#' Data format for DSM
#'
#' Two data frame must be provided to \code{\link{dsm}}. They are referred
#' to as \code{observation.data} and \code{segment.data} (for observation 
#' and segment data, respectively.
#'
#' \code{observation.data} - the observation data frame must have the 
#' following columns:
#' \tabular{ll}{\code{object} \tab unique object identifier \cr
#'    \code{Sample.Label} \tab the identifier for the segment that the 
#'      observation occurred in \cr
#'    \code{size} \tab the size of each observed group (i.e. 1 for 
#'      individuals) \cr
#'    \code{distance} \tab perpendicular/radial distance to observation}
#'
#' \code{segment.data} - the segment data frame must have the following columns:
#'  \tabular{ll}{
#'    \code{Effort} \tab the effort (in terms of length of the segment)\cr
#'    \code{Transect.Label} \tab identifier for this transect\cr
#'    \code{Sample.Label} \tab identifier for the segment (unique!)\cr
#'    \code{???} \tab environmental covariates, for example \code{x} and 
#'          \code{y}}
#'
#' @name dsm-data
NULL
