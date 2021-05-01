#' Predict from a fitted density surface model
#'
#' Make predictions of density or abundance outside (or inside) the covered
#' area.
#'
#' If `newdata` is not supplied, predictions are made for the data that built
#' the model. Note that the order of the results will not necessarily be the
#' same as the `segdata` (segment data) `data.frame` that was supplied to
#' [`dsm`][dsm].
#'
#' The area `off.set` is used if that argument is supplied, otherwise it will
#' look for the areas in the column named `off.set` in the `newdata`. Either
#' way the link function (usually `log`) will be applied to the offsets, so
#' there is no need to log them before passing them to this function.
#'
#' @param object a fitted [`dsm`][dsm] object
#' @param newdata spatially referenced covariates e.g. altitude, depth,
#' distance to shore, etc. Covariates in the `data.frame` must have names
#' *identical* to variable names used in fitting the model
#' @param off.set area of each of the cells in the prediction grid. Should be
#' in the same units as the segments/distances given to `dsm`. Replaces the
#' column in `newdata` called `off.set` if it is supplied. Ignored if `newdata`
#' is not supplied
#' @param type what scale should the results be on. The default is
#' `"response"`, see [`predict.gam`][mgcv::predict.gam] for an explanation of
#' other options (usually not necessary)
#' @param \dots any other arguments passed to [`predict.gam`][mgcv::predict.gam]
#' @return predicted values on the response scale by default (unless `type` is
#' specified, in which case see [`predict.gam`][mgcv::predict.gam]).
#' @export
#' @seealso [`predict.gam`][mgcv::predict.gam], [`dsm.var.gam`][dsm.var.gam],
#' [`dsm.var.prop`][dsm.var.prop]
#' @author David L. Miller
#' @importFrom stats predict
predict.dsm <- function(object, newdata=NULL, off.set=NULL,
                        type="response", ...){

  if("gamm" %in% class(object)){
    object <- object$gam
  }

  # if there is no newdata, then use the data stored in the model object
  if(is.null(newdata)){
    newdata <- object$data

    if(!is.null(off.set)){
      warning("Ignoring supplied off.set as newdata was not supplied")
    }
    if(c(object$formula[[2]]) %in% c("density", "density.est")){
      newdata$off.set <- 1
    }
  }else{

    # if no offset was supplied
    if(!any("off.set" %in% colnames(newdata)) &
       is.null(off.set)){
      stop("You must supply off.set in data or as an argument.")
    }
    # if off.set was supplied, then replace what's in newdata
    if(!is.null(off.set)){
      newdata$off.set <- off.set
    }
    # applied the link
    newdata$off.set <- object$family$linkfun(newdata$off.set)
  }

  # remove the dsm class
  class(object) <- class(object)[class(object)!="dsm"]

  # actually do the predict call
  result <- predict(object, newdata, type=type, ...)

  if(c(object$formula[[2]]) %in% c("density", "density.est")){
    ### if we have density do the predictions on response scale right
    # grab standard error
    se.fit <- list(...)$se.fit
    # only need do this for type="response" but need to make sure
    # that the standard errors are okay too
    if(type=="response" & (is.null(se.fit) || !se.fit)){
      result <- result*newdata$off.set
    }else if(type=="response" & (!is.null(se.fit) && se.fit)){
      result$fit <- result$fit*newdata$off.set
      result$se.fit <- result$se.fit*newdata$off.set
    }
  }

  return(result)
}
