#' Fit a density surface model to segment-specific estimates of abundance
#' or density.
#'
#' Fits a density surface model (DSM) to detection adjusted counts from a spatially-referenced distance sampling analysis. \code{\link{dsm}} takes observations of animals, allocates them to segments of line (or strip transects) and optionally adjusts the counts based on detectability using a supplied detection function model. A generalized additive model, generalized mixed model or generalized linear model is then used to model these adjusted counts based on a formula involving environmental covariates.
#'
#' The response (LHS of `formula`) can be one of the following (with restrictions outlined below):
#' \tabular{ll}{
#'   \code{count} \tab count in each segment\cr
#'   \code{abundance.est} \tab estimated abundance per segment, estimation is via a Horvitz-Thompson estimator.\cr
#'   \code{density.est} \tab density per segment\cr
#'  }
#'
#' The offset used in the model is dependent on the response:
#' \tabular{ll}{
#'   count \tab area of segment multiplied by average probability of detection in the segment\cr
#'   estimated count \tab area of the segment\cr
#'   density \tab zero\cr
#'  }
#'
#' The \code{count} response can only be used when detection function covariates only vary between segments/points (not within). For example, weather conditions (like visibility or sea state) or foliage cover are usually acceptable as they do not change within the segment, but animal sex or behaviour will not work. The \code{abundance.est} response can be used with any covariates in the detection function.
#'
#' In the density case, observations can be weighted by segment areas via the \code{weights=} argument. By default (\code{weights=NULL}), when density is estimated the weights are set to the segment areas (using \code{segment.area} or by calculated from detection function object metadata and \code{Effort} data). Alternatively \code{weights=1} will set the weights to all be equal. A third alternative is to pass in a vector of length equal to the number of segments, containing appropriate weights.
#'
#' A example analyses are available at \url{http://examples.distancesampling.org}.
#'
#' @section Units:
#'
#' It is often the case that distances are collected in metres and segment lengths are recorded in kilometres. \code{dsm} allows you to provide a conversion factor (\code{convert.units}) to multiply the areas by. For example: if distances are in metres and segment lengths are in kilometres setting \code{convert.units=1000} will lead to the analysis being in metres. Setting \code{convert.units=1/1000} will lead to the analysis being in kilometres. The conversion factor will be applied to `segment.area` if that is specified.
#'
#' @section Large models:
#'
#' For large models, \code{engine="bam"} with \code{method="fREML"} may be useful. Models specified for \code{bam} should be as \code{gam}. READ \code{\link{bam}} before using this option; this option is considered EXPERIMENTAL at the moment. In particular note that the default basis choice (thin plate regression splines) will be slow and that in general fitting is less stable than when using \code{gam}. For negative binomial response, theta must be specified when using \code{bam}.
#'
#' @param formula formula for the surface. This should be a valid \code{\link{glm}}/\code{\link{gam}}/\code{\link{gamm}} formula. See "Details", below, for how to define the response.
#' @param ddf.obj result from call to \code{\link{ddf}} or \code{\link[Distance]{ds}}. If multiple detection functions are required a \code{list} can be provided. For strip/circle transects where it is assumed all objects are observed, see \code{\link{dummy_ddf}}. Mark-recapture distance sampling (\code{mrds}) models of type \code{io} (independent observers) are allowed.
#' @param segment.data segment data, see \code{\link{dsm-data}}.
#' @param observation.data observation data, see \code{\link{dsm-data}}.
#' @param engine which fitting engine should be used for the DSM (\code{\link{glm}}/\code{\link{gam}}/\code{\link{gamm}}/\code{\link{bam}}).
#' @param convert.units conversion factor to multiply the area of the segments by. See 'Units' below.
#' @param family response distribution (popular choices include \code{\link{quasipoisson}}, \code{\link{Tweedie}}/\code{\link{tw}} and \code{\link{negbin}}/\code{\link{nb}}). Defaults to \code{quasipossion}.
#' @param group if \code{TRUE} the abundance of groups will be calculated rather than the abundance of individuals. Setting this option to \code{TRUE} is equivalent to setting the size of each group to be 1.
#' @param control the usual \code{control} argument for a \code{gam}; \code{keepData} must be \code{TRUE} for variance estimation to work (though this option cannot be set for GLMs or GAMMs.
#' @param availability an estimate of availability bias. For count models used to multiply the effective strip width (must be a vector of length 1 or length the number of rows in \code{segment.data}); for estimated abundance/estimated density models used to scale the response (must be a vector of length 1 or length the number of rows in \code{observation.data}). Uncertainty in the availability is not handled at present.
#' @param segment.area if `NULL` (default) segment areas will be calculated by multiplying the `Effort` column in `segment.data` by the (right minus left) truncation distance for the `ddf.obj` or by `strip.width`. Alternatively a vector of segment areas can be provided (which must be the same length as the number of rows in `segment.data`) or a character string giving the name of a column in `segment.data` which contains the areas. If \code{segment.area} is specified it takes precident.
#' @param weights weights for each observation used in model fitting. The default, \code{weights=NULL}, weights each observation by its area (see Details). Setting a scalar value (e.g. \code{weights=1}) all observations are equally weighted.
#' @param method The smoothing parameter estimation method. Default is \code{"REML"}, using Restricted Maximum Likelihood. See \code{\link{gam}} for other options. Ignored for \code{engine="glm"}.
#' @param \dots anything else to be passed straight to \code{\link{glm}}/\code{\link{gam}}/\code{\link{gamm}}/\code{\link{bam}}.
#' @return a \code{\link{glm}}/\code{\link{gam}}/\code{\link{gamm}} object, with an additional element, \code{ddf} which holds the detection function object.
#' @author David L. Miller
# @seealso
#' @references Hedley, S. and S. T. Buckland. 2004. Spatial models for line transect sampling. JABES 9:181-199.
#'
#' Miller, D. L., Burt, M. L., Rexstad, E. A., Thomas, L. (2013), Spatial models for distance sampling data: recent developments and future directions. Methods in Ecology and Evolution, 4: 1001-1010. doi: 10.1111/2041-210X.12105 (Open Access, available at http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12105/abstract)
#'
#' Wood, S.N. 2006. Generalized Additive Models: An Introduction with R. CRC/Chapman & Hall.
#' @export
#' @importFrom stats quasipoisson
#' @importFrom utils packageVersion
#'
#' @examples
#' \dontrun{
#' library(Distance)
#' library(dsm)
#'
#' # load the Gulf of Mexico dolphin data (see ?mexdolphins)
#' data(mexdolphins)
#'
#' # fit a detection function and look at the summary
#' hr.model <- ds(distdata, max(distdata$distance),
#'                key = "hr", adjustment = NULL)
#' summary(hr.model)
#'
#' # fit a simple smooth of x and y to counts
#' mod1 <- dsm(count~s(x,y), hr.model, segdata, obsdata)
#' summary(mod1)
#'
#' # predict over a grid
#' mod1.pred <- predict(mod1, preddata, preddata$area)
#'
#' # calculate the predicted abundance over the grid
#' sum(mod1.pred)
#'
#' # plot the smooth
#' plot(mod1)
#'}
dsm <- function(formula, ddf.obj, segment.data, observation.data,
                engine="gam", convert.units=1,
                family=quasipoisson(link="log"), group=FALSE,
                control=list(keepData=TRUE), availability=1,
                segment.area=NULL, weights=NULL, method="REML",
                ...){

  stopifnot(engine %in% c("gam","bam","glm","gamm"))

  # if we have a model fitted using Distance, then just pull out the
  # ddf component
  if(is.null(ddf.obj)){
    stop("NULL detection functions no longer supported, see ?dummy_ddf")
  }
  # if we don't have one detection function, but the model was
  # fitted using Distance, then just pull out the ddf component
  if(all(class(ddf.obj)!="list")){
    if(all(class(ddf.obj)=="dsmodel")){
      ddf.obj <- ddf.obj$ddf
    }
  }else{
    if(length(ddf.obj) == 1){
      ddf.obj <- ddf.obj[[1]]
    }
    for(i in seq_along(ddf.obj)){
      if(all(class(ddf.obj[[i]])=="dsmodel")){
        ddf.obj[[i]] <- ddf.obj[[i]]$ddf
      }
    }
  }


  ## check the formula
  response <- as.character(formula)[2]
  # throw an error if we have one of the deprecated responses
  if(response %in% c("presence", "D", "density", "Dhat", "N", "Nhat", "n")){
    stop(paste("Response", response, "is deprecated, see ?dsm for details."))
  }
  possible.responses <- c("density.est",
                          "count",
                          "abundance.est")
  if(!(response %in% possible.responses)){
    stop(paste("Model must be one of:",
               paste(possible.responses, collapse=", ")))
  }

  ## check that the necessary columns exist in the data
  # NB this doesn't return anything just throws an error if something
  #    bad happens
  check.cols(ddf.obj, segment.data, observation.data, segment.area)


  ## build the data
  dat <- make.data(response, ddf.obj, segment.data, observation.data,
                   group, convert.units, availability, segment.area, family)

  ## if we are not modelling density, then add in the offset
  ##  to the formula
  if(!(response %in% c("density.est"))){
    formula <- as.formula(paste(c(as.character(formula)[c(2, 1, 3)],
                                "+ offset(off.set)"), collapse=""))
  }else{
    # set the weights if we are doing density estimation
    if(is.null(weights)){
      weights <- dat$segment.area
    }else if(length(weights)==1){
      weights <- rep(weights, nrow(dat))
    }
  }

  # if we're using a gamm engine, warn if using an old version of mgcv
  if(engine == "gamm" && dsm_env$old_mgcv){
      message("You are using mgcv version < 1.7-24, please update to at least 1.7-24 to avoid fitting problems.")
  }

  # GLMs and GAMMs don't support keeping the data
  if(engine %in% c("glm", "gamm")){
    # unsupported
    control$keepData <- NULL
  }

  ## run the engine
  args <- list(formula = formula,
               family  = family,
               data    = dat,
               weights = weights,
               control = control,
               method  = method,
               ...)

  # glm doesn't understand method= in the way we think of it
  if(engine == "glm"){
    args$method <- NULL
  }

  fit <- withCallingHandlers(do.call(engine, args),
                             warning=matrixnotposdef.handler)

  ## save knots
  if("knots" %in% names(match.call())){
    fit$knots <- get(as.character(match.call()$knots))
  }

  ## add the detection function object into the gam/gamm/glm object
  if(engine == "gamm"){
    fit$gam$ddf <- ddf.obj
    fit$gam$data <- dat
    # yucky way to get dsm.var/dsm.var.prop to work because gamm()
    #  doesn't store the call()
    fit$gam$gamm.call.list <- list(formula = formula,
                                   family  = family,
                                   data    = dat,
                                   control = control)
  }else{
    fit$ddf <- ddf.obj
  }

  class(fit) <- c("dsm",class(fit))

  ## return the model
  return(fit)

}
