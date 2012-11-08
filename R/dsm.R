#' Fit a density surface model to segment-specific estimates of abundance
#' or density.
#'
#' Given a detection function analysis, construct a density surface model (DSM)
#' based on environmental covariates.
#'
#' @param formula formula for the surface. This should be a
#'   valid \code{\link{glm}}/\code{\link{gam}}/code{\link{gamm}} formula.
#' @param ddf.obj result from call to \code{\link{ddf}} or \code{\link{ds}}.
#'   If \code{ddf.obj} is \code{NULL} then strip transects are assumed.
#' @param segment.data segment data, see \code{\link{dsm-data}}.
#' @param observation.data observation data, see \code{\link{dsm-data}}.
#' @param engine which model should be used for the DSM (\code{\link{glm}}/
#'   \code{\link{gam}}/code{\link{gamm}}).
#' @param convert.units value to alter length to width for calculation
#'   of the offset.
#' @param family response distribution (popular choices include
#'   \code{\link{quasipoisson}}, \code{\link{Tweedie}} and \code{\link{negbin}}.
#' @param \dots anything else to be passed straight to \code{\link{gam}}.
#' @param group should group abundance/density be modelled rather than
#'  individual abundance/density? This effectively sets the \code{size} column
#'  in \code{observation.data} to be 1.
#' @param control the usual \code{control} argument for a \code{gam},
#'  \code{keepData} must be \code{TRUE} or variance estimation will not work.
#' @return a \code{\link{glm}}/\code{\link{gam}}/code{\link{gamm}} object, with
#'  an additional element, \code{ddf} which holds the detection function object.
#' @param gamma parameter to \code{gam()} set to a value of 1.4 (from advice in
#   Wood (2006)) such that the \code{gam()} is inclined to not 'overfit.'.
#'
#' @author David L. Miller
# @seealso
#' @references Hedley, S. and S. T. Buckland. 2004. Spatial models for line transect sampling. JABES 9:181-199.
#'
#' Wood, S.N. 2006. Generalized Additive Models: An Introduction with R. CRC/Chapman & Hall.
#' @export
# @keywords
# @examples
dsm <- function(formula, ddf.obj, segment.data, observation.data,
                engine="gam", convert.units=1,
                family=quasipoisson(link="log"), group=FALSE, gamma=1.4, 
                control=list(keepData=TRUE), ...){


  ## check the formula
  response <- as.character(formula)[2]
  possible.responses <- c("D","density","N","Nhat","abundance","abundance.est")
  if(!(response %in% possible.responses)){
    stop(paste("Model must be one of:",
               paste(possible.responses,collapse=", ")))
  }


  # if we are not modelling density, then add in the offset
  if(!(response %in% c("D","density"))){
    formula <- as.formula(paste(c(as.character(formula)[c(2,1,3)],
                                "+ offset(off.set)"),collapse=""))
  }

  ## build the data
  dat <- make.data(response, ddf.obj, segment.data, observation.data,
                   group, convert.units)

  ## run the engine
  if(engine == "gam"){
    fit <- gam(formula,family=family, data=dat, gamma=gamma, 
               control=control, ...)
  }else if(engine == "gamm"){
    fit <- gamm(formula,family=family, data=dat, ...)
  }else if(engine == "glm"){
    fit <- glm(formula,family=family, data=dat, ...)
  }else{
    stop("engine must be one of 'gam', 'gamm' or 'glm'")
  }

  ## add the detection function object into the gam/gamm/glm object
  fit$ddf <- ddf.obj

  class(fit) <- c("dsm",class(fit))

  ## return the model
  return(fit)

}
