#' Prediction variance estimation assuming independence
#'
#' If one is willing to assume the the detection function and spatial model are
#' independent, this function will produce estimates of variance of predictions
#' of abundance, using the result that squared coefficients of variation will
#' add.
#'
#' @param dsm.obj a model object fitted by [`dsm`][dsm].
#' @param pred.data either: a single prediction grid or list of prediction
#' grids. Each grid should be a `data.frame` with the same columns as the
#' original data.
#' @param off.set a a vector or list of vectors with as many elements as there
#' are in `pred.data`. Each vector is as long as the number of rows in the
#' corresponding element of `pred.data`. These give the area associated with
#' each prediction cell. If a single number is supplied it will be replicated
#' for the length of `pred.data`.
#' @param seglen.varname name for the column which holds the segment length
#' (default value `"Effort"`).
#' @param type.pred should the predictions be on the "response" or "link"
#' scale? (default `"response"`).
#' @return a `list` with elements
#'   * `model` the fitted model object
#'   * `pred.var` variance of the regions given in `pred.data`.
#'   * `bootstrap` logical, always `FALSE`
#'   * `model` the fitted model with the extra term
#'   * `dsm.object` the original model (`dsm.obj` above)
#' @author David L. Miller
#' @importFrom stats coef vcov
#' @export
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
#'  mod1.var <- dsm.var.gam(mod1, preddata, off.set=preddata$area)
#' }
dsm.var.gam <- function(dsm.obj, pred.data, off.set,
                        seglen.varname='Effort', type.pred="response"){

  # strip dsm class so we can use gam methods
  class(dsm.obj) <- class(dsm.obj)[class(dsm.obj) != "dsm"]

  # if we have a gamm, then just pull out the gam object
  if(any(class(dsm.obj) == "gamm")){
    dsm.obj <- dsm.obj$gam
    is.gamm <- TRUE
  }

  # only support log links
  if(dsm.obj$family$link != "log"){
    stop("Only models with a log link are supported")
  }

  # if all the offsets are the same then we can just supply 1 and rep it
  if(length(off.set) == 1){
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
#    pred.data[[1]] <- pred.data
#    off.set[[1]] <- off.set
  }else if(is.list(off.set)){
    if(length(pred.data)!=length(off.set)){
      stop("pred.data and off.set don't have the same number of elements")
    }
  }

  # depending on whether we have response or link scale predictions...
  if(type.pred=="response"){
    tmfn <- dsm.obj$family$linkinv
  }else if(type.pred=="link"){
    tmfn <- identity
  }

  # grab the coefficients
  cft <- coef(dsm.obj)
  preddo <- list(length(pred.data))
  dpred.db <- matrix(0, length(pred.data), length(cft))

  # loop over the prediction grids
  for(ipg in seq_along(pred.data)){
    ### fancy lp matrix stuff
    # set the offset to be zero here so we can use lp
    pred.data[[ipg]]$off.set<-rep(0, nrow(pred.data[[ipg]]))

    lpmat <- predict(dsm.obj, newdata=pred.data[[ipg]], type='lpmatrix')
    lppred <- lpmat %**% cft

    # if the offset is just one number then repeat it enough times
    if(length(off.set[[ipg]]) == 1){
      this.off.set <- rep(off.set[[ipg]], nrow(pred.data[[ipg]]))
    }else{
      this.off.set <- off.set[[ipg]]
    }

    preddo[[ipg]] <-  this.off.set %**% tmfn(lppred)


    # NB previously this was done numerically but things didn't work
    #    when se ~=0 and actual 0s were produced. Use analytical expression
    #    for the log case here
    if(type.pred=="link"){
      dpred.db[ipg, ] <- this.off.set %**% lppred*lpmat
    }else{
      dpred.db[ipg, ] <- this.off.set %**% (tmfn(lppred)*lpmat)
    }
  }

  # "'vpred' is the covariance of all the summary-things." - MVB
  # so we want the diagonals if length(pred.data)>1
  # A B A^tr
  vpred <- dpred.db %**% tcrossprod(vcov(dsm.obj), dpred.db)

  if(is.matrix(vpred)){
    vpred <- diag(vpred)
  }

  result <- list(pred.var       = vpred,
                 bootstrap      = FALSE,
                 pred           = preddo,
                 var.prop       = FALSE,
                 pred.data      = pred.data,
                 off.set        = off.set,
                 dsm.object     = dsm.obj,
                 seglen.varname = seglen.varname,
                 type.pred      = type.pred
                )

  class(result) <- "dsm.var"

  return(result)
}
