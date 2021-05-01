#' Density surface modelling
#'
#' `dsm` implements spatial models for distance sampling data. Models for detectability can be fitted using packages `mrds` or `Distance`. `dsm` fits generalized additive models to spatially-referenced data. See Miller et al (2013) for an introduction.
#'
#' Further information on distance sampling methods and example code is available at <http://distancesampling.org/R/>.
#'
#' For help with distance sampling and this package, there is a Google Group
#' <https://groups.google.com/forum/#!forum/distance-sampling>.
#'
#' A example analyses are available at <http://examples.distancesampling.org>.
#'
#' @name dsm-package
#' @import mgcv nlme mrds
#' @references Hedley, S. and S. T. Buckland. 2004. Spatial models for line
#' transect sampling. JABES 9:181-199.
#'
#' Miller, D. L., Burt, M. L., Rexstad, E. A., Thomas, L. (2013), Spatial
#' models for distance sampling data: recent developments and future
#' directions. Methods in Ecology and Evolution, 4: 1001-1010. doi:
#' 10.1111/2041-210X.12105 (Open Access, available at
#' <http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12105/abstract>)
#'
#' Wood, S.N. 2006. Generalized Additive Models: An Introduction with R.
#' CRC/Chapman & Hall.
#'
NULL


#' Data format for DSM
#'
#' Two `data.frame`s must be provided to [`dsm`][dsm]. They are referred to as
#' `observation.data` and `segment.data`.
#'
#' The `segment.data` table has the sample identifiers which define the
#' segments, the corresponding effort (line length) expended and the
#' environmental covariates that will be used to model abundance/density.
#' `observation.data` provides a link table between the observations used in
#' the detection function and the samples (segments), so that we can aggregate
#' the observations to the segments (i.e., `observation.data` is a "look-up
#' table" between the observations and the segments).
#'
#' `observation.data` - the observation `data.frame` must have (at least) the following columns:
#'   * `object` unique object identifier
#'   * `Sample.Label` the identifier for the segment where observation occurred
#'   * `size` the size of each observed group (e.g., 1 if all animals occurred
#'  individually)
#'   * `distance` distance to observation
#'
#' One can often also use `observation.data` to fit a detection function (so
#' additional columns for detection function covariates are allowed in this
#' table).
#'
#' `segment.data`: the segment `data.frame` must have (at least) the following
#' columns:
#'   * `Effort` the effort (in terms of length of the segment)
#'   * `Sample.Label` identifier for the segment (unique!)
#'   * ??? environmental covariates, for example location (projected latitude
#'   and longitude), and other relevant covariates (sea surface temperature,
#'   foliage type, altitude, bathymetry etc).
#'
#' @section Multiple detection functions:
#'
#' If multiple detection functions are to be used, then a column named `ddfobj`
#' must be included in `observation.data` and `segment.data`. This lets the
#' model know which detection function each observation is from. These are
#' numeric and ordered as the `ddf.obj` argument to [`dsm`][dsm], e.g.,
#' `ddf.obj=list(ship_ddf, aerial_ddf)` means ship detections have `ddfobj=1`
#' and aerial detections have `ddfobj=2` in the observation data.
#'
#' @section Mark-recapture distance sampling models:
#'
#' When using `mrds` models that include mark-recapture components (currently independent observer and trial modes are supported) then the format of the observation data needs to be checked to ensure that observations are not duplicated. The `observer` column is also required in the `observation.data`.
#' * *Independent observer mode* only unique observations (unique object IDs) are required.
#'  * *Trial mode* only observations made by observer 1 are required.
#'
#' @name dsm-data
NULL

#' Pan-tropical spotted dolphins in the Gulf of Mexico
#'
#' Data from a combination of several NOAA shipboard surveys conducted on
#' pan-tropical spotted dolphins in the Gulf of Mexico. The data consist of 47
#' observations of groups of dolphins. The group size was recorded, as well as
#' the Beaufort sea state at the time of the observation. Coordinates for each
#' observation and bathymetry data were also available as covariates for the
#' analysis. A complete example analysis (and description of the data) is
#' provided at <http://distancesampling.org/R/vignettes/mexico-analysis.html>.
#'
#' @references Halpin, P.N., A.J. Read, E. Fujioka, B.D. Best, B. Donnelly,
#' L.J. Hazen, C. Kot, K. Urian, E. LaBrecque, A. Dimatteo, J. Cleary, C. Good,
#' L.B. Crowder, and K.D. Hyrenbach. 2009. OBIS-SEAMAP: The world data center
#' for marine mammal, sea bird, and sea turtle distributions. Oceanography
#' 22(2):104-115
#'
#' NOAA Southeast Fisheries Science Center. 1996. Report of a Cetacean Survey
#' of Oceanic and Selected Continental Shelf Waters of the Northern Gulf of
#' Mexico aboard NOAA Ship Oregon II (Cruise 220)
#'
#' @name mexdolphins
#' @aliases obsdata segdata preddata pred.polys distdata survey.area
#' @docType data
#' @keywords datasets
NULL
