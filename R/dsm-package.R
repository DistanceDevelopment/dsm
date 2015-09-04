#' Density surface modelling
#'
#' \code{dsm} implements spatial models for distance sampling data.
#'
#' Further information on distance sampling methods and example code is available at \url{http://distancesampling.org/R/}.
#'
#' For help with distance sampling and this package, there is a Google Group \url{https://groups.google.com/forum/#!forum/distance-sampling}.
#'
#' A example analysis is available at \url{http://distancesampling.org/R/vignettes/mexico-analysis.html}.
#'
#' @name dsm-package
#' @import mgcv nlme mrds
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
#' Two \code{data.frame}s must be provided to \code{\link{dsm}}. They are referred to as \code{observation.data} and \code{segment.data} (for observation and segment data, respectively).
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
#'    \code{Sample.Label} \tab identifier for the segment (unique!)\cr
#'    \code{???} \tab environmental covariates, for example \code{x} and 
#'          \code{y}}
#'
#' @name dsm-data
NULL

#' Pan-tropical spotted dolphins in the Gulf of Mexico
#'
#' Data from a combination of several NOAA shipboard surveys conducted on pan-tropical spotted dolphins in the Gulf of Mexico. 47 observations of groups of dolphins The group size was recorded, as well as the Beaufort sea state at the time of the observation. Coordinates for each observation and bathymetry data were also available as covariates for the analysis. A complete example analysis (and description of the data) is provided at \url{http://distancesampling.org/R/vignettes/mexico-analysis.html}.
#'
#' @references Halpin, P.N., A.J. Read, E. Fujioka, B.D. Best, B. Donnelly, L.J. Hazen, C. Kot, K. Urian, E. LaBrecque, A. Dimatteo, J. Cleary, C. Good, L.B. Crowder, and K.D. Hyrenbach. 2009. OBIS-SEAMAP: The world data center for marine mammal, sea bird, and sea turtle distributions. Oceanography 22(2):104-115
#'
#' NOAA Southeast Fisheries Science Center. 1996. Report of a Cetacean Survey of Oceanic and Selected Continental Shelf Waters of the Northern Gulf of Mexico aboard NOAA Ship Oregon II (Cruise 220)
#'
#' @name mexdolphins
#' @docType data
#' @keywords datasets
NULL
