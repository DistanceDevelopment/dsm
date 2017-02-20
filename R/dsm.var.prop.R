#' Prediction variance propogation for DSMs
#'
#' To ensure that uncertainty from the detection function is correctly propagated to the final variance estimate of abundace, this function uses a method first detailed in Williams et al (2011).
#'
#' The idea is to refit the spatial model but including an extra random effect. This random effect has zero mean and hence to effect on point estimates. Its variance is the Hessian of the detection function. Variance estimates then incorporate detection function uncertainty. Further mathematical details are given in the paper in the references below.
#'
#' Many prediction grids can be supplied by supplying a list of \code{data.frame}s to the function.
#'
#' Note that this routine is only useful if a detection function has been used in the DSM. It cannot be used when the \code{Nhat}, \code{abundance.est} responses are used.
#'
#' Based on (much more general) code from Mark Bravington and Sharon Hedley.
#'
#' Note that this routine simply calls \code{\link{dsm_varprop}}. If you don't require multiple prediction grids, the other routine will probably be faster.
#'
#' @inheritParams dsm.var.gam
#' @return a list with elements
#'         \tabular{ll}{\code{model} \tab the fitted model object\cr
#'                      \code{pred.var} \tab variance of each region given
#'                      in \code{pred.data}.\cr
#'                      \code{bootstrap} \tab logical, always \code{FALSE}\cr
#'                      \code{pred.data} \tab as above\cr
#'                      \code{off.set} \tab as above\cr
#'                      \code{model}\tab the fitted model with the extra term\cr
#'                      \code{dsm.object} \tab the original model, as above\cr
#'                      \code{model.check} \tab simple check of subtracting the coefficients of the two models to see if there is a large difference\cr
#'                      \code{deriv} \tab numerically calculated Hessian of the offset\cr
#'                      }
#' @author Mark V. Bravington, Sharon L. Hedley. Bugs added by David L. Miller.
#' @references
#' Williams, R., Hedley, S.L., Branch, T.A., Bravington, M.V., Zerbini, A.N. and Findlay, K.P. (2011). Chilean Blue Whales as a Case Study to Illustrate Methods to Estimate Abundance and Evaluate Conservation Status of Rare Species. Conservation Biology 25(3), 526-535.
#' @export
#' @importFrom stats as.formula update
#' @importFrom numDeriv grad
#' @examples
#' \dontrun{
#'  library(Distance)
#'  library(dsm)
#'
#'  # load the Gulf of Mexico dolphin data (see ?mexdolphins)
#'  data(mexdolphins)
#'
#'  # fit a detection function and look at the summary
#'  hr.model <- ds(distdata, max(distdata$distance),
#'                 key = "hr", adjustment = NULL)
#'  summary(hr.model)
#'
#'  # fit a simple smooth of x and y
#'  mod1 <- dsm(count~s(x, y), hr.model, segdata, obsdata)
#'
#'  # Calculate the variance
#'  # this will give a summary over the whole area in mexdolphins$preddata
#'  mod1.var <- dsm.var.prop(mod1, preddata, off.set=preddata$area)
#' }
dsm.var.prop <- function(dsm.obj, pred.data, off.set,
                         seglen.varname='Effort', type.pred="response") {

  ## pre-checking...
  # die if we have a gamm
  if(any(class(dsm.obj)=="gamm")){
    stop("GAMMs are not supported.")
  }

  # check that there are no covariates in the df model
  if(length(unique(dsm.obj$ddf$fitted)) > 1){
    stop("Covariate detection functions are not currently supported within dsm.var.prop.")
  }

  # break if we use the wrong response
  if(!(as.character(dsm.obj$formula)[2] %in% c("N", "count"))){
    stop("Variance propagation can only be used with count as the response.")
  }

  # if there is no ddf object, then we should stop!
  # thanks to Adrian Schiavini for spotting this
  if(is.null(dsm.obj$ddf)){
    stop("No detection function in this analysis, use dsm.var.gam")
  }

  ## end of checks

  ## data setup
  # if all the offsets are the same then we can just supply 1 and rep it
  if(length(off.set)==1){
    if(is.null(nrow(pred.data))){
      off.set <- rep(list(off.set), length(pred.data))
    }else{
      off.set <- rep(off.set, nrow(pred.data))
    }
  }

  # make sure if one of pred.data and off.set is not a list we break
  # if we didn't have a list, then put them in a list so everything works
  if(is.data.frame(pred.data) & is.vector(off.set)){
    pred.data <- list(pred.data)
    off.set <- list(off.set)
  }else if(is.list(off.set)){
    if(length(pred.data)!=length(off.set)){
      stop("pred.data and off.set don't have the same number of elements")
    }
  }

  # push the offsets into the data...
  for(i in seq_along(pred.data)){
    pred.data[[i]]$off.set <- off.set[[i]]
  }
  ## end data setup

  ## run varprop and calculate stuff!
  #varp <- dsm_varprop(dsm.obj, pred.data[[1]])
  #fit.with.pen <- varp$refit

# storage
vpred <- length(pred.data)
preddo <- list(length(pred.data))
varp <- list()

  # loop over the prediction grids
  for(ipg in seq_along(pred.data)){
    varp <- dsm_varprop(dsm.obj, pred.data[[1]])
    vpred[ipg] <- varp$var
    preddo[[ipg]] <- sum(varp$pred)
  }

  # Diagnostic from MVB
  # check that the fitted model isn't too different, used in summary()
  model.check <- summary(fitted(varp$refit) - fitted(dsm.obj))

  result <- list(pred.var = vpred,
                 bootstrap = FALSE,
                 var.prop = TRUE,
                 pred.data = pred.data,
                 pred = preddo,
                 off.set = off.set,
                 model = varp$refit,
                 dsm.object = dsm.obj,
                 model.check = model.check,
                 #deriv = firstD,
                 seglen.varname = seglen.varname,
                 type.pred=type.pred
                )

  class(result) <- "dsm.var"

  return(result)
}
