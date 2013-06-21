#' Fit a density surface model to segment-specific estimates of abundance
#' or density.
#'
#' Given a detection function analysis, construct a density surface model (DSM)
#' based on environmental covariates.
#'
#' The response can be one of the following:
#' \tabular{ll}{
#'              \code{N}, \code{abundance} \tab count in each segment\cr
#'              \code{Nhat}, \code{abundance.est} \tab estimated abundance per segment, estimation is via a Horvitz-Thompson estimator. This should be used when there are covariates in the detection function.\cr
#'              \code{presence} \tab interpret the data as presence/absence (reember to change the \code{family} argument to \code{binomial()}\cr
#'              \code{D}, \code{density} \tab density per segment\cr
#'  }
#'
#' @param formula formula for the surface. This should be a
#'   valid \code{\link{glm}}/\code{\link{gam}}/\code{\link{gamm}} formula.
#' @param ddf.obj result from call to \code{\link{ddf}} or \code{\link{ds}}.
#'   If \code{ddf.obj} is \code{NULL} then strip transects are assumed.
#' @param segment.data segment data, see \code{\link{dsm-data}}.
#' @param observation.data observation data, see \code{\link{dsm-data}}.
#' @param engine which model should be used for the DSM (\code{\link{glm}}/
#'   \code{\link{gam}}/code{\link{gamm}}).
#' @param convert.units value to alter length to width for calculation
#'   of the offset.
#' @param family response distribution (popular choices include \code{\link{quasipoisson}}, \code{\link{Tweedie}} and \code{\link{negbin}}. Defaults to \code{quasipossion}.
#' @param \dots anything else to be passed straight to \code{\link{gam}}.
#' @param group should group abundance/density be modelled rather than individual abundance/density? This effectively sets the \code{size} column in \code{observation.data} to be 1.
#' @param control the usual \code{control} argument for a \code{gam}, \code{keepData} must be \code{TRUE} or variance estimation will not work.
#' @param availability an availability bias used to scale the counts/estimated  counts by. If we have \code{N} animals in a segment, then \code{N/availability} will be entered into the model. Uncertainty in the availability is not handled at present.
#' @param gamma parameter to \code{gam()} set to a value of 1.4 (from advice in Wood (2006)) such that the \code{gam()} is inclined to not 'overfit.'.
#' @param strip.width if \code{ddf.obj}, above, is \code{NULL}, then this is where the strip width is specified. Note that this is the total width, i.e. right truncation minus left truncation.
#' @return a \code{\link{glm}}/\code{\link{gam}}/\code{\link{gamm}} object, with an additional element, \code{ddf} which holds the detection function object.
#' @author David L. Miller
# @seealso
#' @references Hedley, S. and S. T. Buckland. 2004. Spatial models for line transect sampling. JABES 9:181-199.
#'
#' Wood, S.N. 2006. Generalized Additive Models: An Introduction with R. CRC/Chapman & Hall.
#' @export
#'
#' @examples
#'
#' library(Distance)
#' library(dsm)
#'
#' # load the Gulf of Mexico dolphin data (see ?mexdolphins)
#' data(mexdolphins)
#'
#' # fit a detection function and look at the summary
#' hr.model <- ds(mexdolphins$distdata, max(mexdolphins$distdata$distance), key = "hr", adjustment = NULL)
#' summary(hr.model)
#'
#' # fit a simple smooth of x and y
#' mod1<-dsm(N~s(x,y), hr.model, mexdolphins$segdata, mexdolphins$obsdata)
#' summary(mod1)
#'
#' # create an offset (in metres)
#' # each prediction cell is 444km2
#' off.set <- 444*1000*1000
#'
#' # predict over a grid
#' mod1.pred <- predict(mod1, mexdolphins$preddata, off.set)
#'
#' # calculate the predicted abundance over the grid
#' sum(mod1.pred)
#'
#' # plot the smooth
#' plot(mod1)
dsm <- function(formula, ddf.obj, segment.data, observation.data,
                engine="gam", convert.units=1,
                family=quasipoisson(link="log"), group=FALSE, gamma=1.4,
                control=list(keepData=TRUE), availability=1, strip.width=NULL,
                ...){

  # if we have a model fitted using Distance, then just pull out the
  # ddf component
  if(!is.null(ddf.obj)){
    if(all(class(ddf.obj)=="dsmodel")){
      ddf.obj <- ddf.obj$ddf
    }
  }

  ## check the formula
  response <- as.character(formula)[2]
  possible.responses <- c("D","density",
                          "N","Nhat","abundance","abundance.est",
                          "presence")
  if(!(response %in% possible.responses)){
    stop(paste("Model must be one of:",
               paste(possible.responses,collapse=", ")))
  }


  # if we are not modelling density, then add in the offset
  if(!(response %in% c("D","density","presence"))){
    formula <- as.formula(paste(c(as.character(formula)[c(2,1,3)],
                                "+ offset(off.set)"),collapse=""))
  }

  ## build the data
  dat <- make.data(response, ddf.obj, segment.data, observation.data,
                   group, convert.units, availability, strip.width)

  ## run the engine
  if(engine == "gam"){
    fit <- withCallingHandlers(gam(formula,family=family, data=dat, gamma=gamma,
               control=control, ...), warning = matrixnotposdef.handler)
  }else if(engine == "gamm"){
    # warn if using an old version of mgcv
    mgcv.version <- as.numeric(strsplit(as.character(packageVersion("mgcv")),
                                        "\\.")[[1]])
    if(mgcv.version[1]<1 | (mgcv.version[2]<7 |
                            (mgcv.version[2]==7 & mgcv.version[3]<24))){
      message("You are using mgcv version < 1.7-24, please update to at least 1.7-24 to avoid fitting problems.")
    }

    # unsupported
    control$keepData <- NULL
    fit <- withCallingHandlers(gamm(formula,family=family, data=dat,
                                    gamma=gamma,control=control, ...),
                               warning = matrixnotposdef.handler)
  }else if(engine == "glm"){
    fit <- glm(formula,family=family, data=dat, ...)
  }else{
    stop("engine must be one of 'gam', 'gamm' or 'glm'")
  }

  ## save knots
  if("knots" %in% names(match.call())){
    fit$knots <- get(as.character(match.call()$knots))
  }

  ## add the detection function object into the gam/gamm/glm object
  if(engine == "gamm"){
    fit$gam$ddf <- ddf.obj
    fit$gam$data <- dat
    fit$gam$gamma <- gamma
    # yucky way to get dsm.var/dsm.var.prop to work because gamm()
    #  doesn't store the call()
    fit$gam$gamm.call.list <- substitute(list(...))
    fit$gam$gamm.call.list$formula <- formula
    fit$gam$gamm.call.list$family <- family
    fit$gam$gamm.call.list$data <- dat
    fit$gam$gamm.call.list$gamma <- gamma
    fit$gam$gamm.call.list$control <- control
  }else{
    fit$ddf <- ddf.obj
    fit$gamma <- gamma
  }

  class(fit) <- c("dsm",class(fit))

  ## return the model
  return(fit)

}
