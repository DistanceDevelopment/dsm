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
#' Two \code{data.frame}s must be provided to \code{\link{dsm}}. They are referred to as \code{observation.data} and \code{segment.data}.
#'
#' The \code{segment.data} table has the sample identifiers which define the segments, the corresponding effort (line length) expended and the environmental covariates that will be used to model abundance/density. \code{observation.data} provides a link table between the observations used in the detection function and the samples (segments), so that we can aggregate the observations to the segments (i.e. \code{observation.data} is a "look-up table" between the observations and the segments).
#'
#' \code{observation.data} - the observation \code{data.frame} must have (at least) the following columns:
#'
#' \tabular{ll}{\code{object} \tab unique object identifier \cr
#'    \code{Sample.Label} \tab the identifier for the segment that the observation occurred in \cr
#'    \code{size} \tab the size of each observed group (e.g 1 if all animals occurred individually)\cr
#'    \code{distance} \tab distance to observation}
#'
#' One can often also use \code{observation.data} to fit a detection function (so additional columns for detection function covariates are allowed in this table).
#'
#' \code{segment.data}: the segment \code{data.frame} must have (at least) the following columns:
#'
#'  \tabular{ll}{
#'    \code{Effort} \tab the effort (in terms of length of the segment)\cr
#'    \code{Sample.Label} \tab identifier for the segment (unique!)\cr
#'    ??? \tab environmental covariates, for example location (projected latitude and longitude), and other relevant covariates (sea surface temperature, foliage type, altitude, bathymetry etc).}
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
#' @aliases obsdata segdata preddata pred.polys distdata survey.area
#' @docType data
#' @keywords datasets
NULL
