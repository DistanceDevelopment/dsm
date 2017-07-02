#' Predict from a fitted density surface model
#'
#' Make predictions outside (or inside) the covered area.
#'
#' If \code{newdata} is not supplied, predictions are made for the data that built the model. Note that the order of the results will not necessarily be the same as the \code{segdata} (segment data) \code{data.frame} that was supplied (it will be sorted by the \code{Segment.Label} field).
#'
#' @param object a fitted \code{\link{dsm}} object as produced by \code{dsm()}.
#' @param newdata spatially referenced covariates e.g. altitude, depth, distance to shore, etc. Covariates in the \code{data.frame} must have names *identical* to variable names used in fitting the DSM.
#' @param off.set area of each of the cells in the prediction grid. Should be in the same units as the segments/distances given to \code{dsm}. Ignored if there is already a column in \code{newdata} called \code{off.set}.
#' @param type what scale should the results be on. The default is
#'  \code{"response"}, see \code{\link{predict.gam}} for an explanation of other options (usually not necessary).
#' @param \dots any other arguments passed to \code{\link{predict.gam}}.
#' @return predicted values on the response scale (density/abundance).
#' @export
#'
#' @seealso \code{\link{predict.gam}} \code{\link{dsm.var.gam}} \code{\link{dsm.var.prop}} \code{\link{dsm.var.movblk}}
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
  }else{

    # if we don't have a density model, then set the offset
    # this is already set if we're using the data that was in the model
    # object, so don't re-log and ignore the off.set specified
    # thanks to Megan Furguson for pointing this out!
    if(!(c(object$formula[[2]]) %in% c("D", "presence", "density"))){
      if(is.null(newdata$off.set)){
        if(is.null(off.set)){
          stop("You must supply off.set in data or as an argument.")
        }
        newdata$off.set <- off.set
      }

      # apply the link function
      linkfn <- object$family$linkfun
      newdata$off.set <- linkfn(newdata$off.set)
    }
  }


  # remove the dsm class
  class(object) <- class(object)[class(object)!="dsm"]

  # actually do the predict call
  result <- predict(object, newdata, type=type, ...)


  ## if we have density do the predictions on response scale right
  # grab standard error
  se.fit <- list(...)$se.fit
  if((c(object$formula[[2]]) %in% c("D", "density"))){
    # only need do this for type="response" but need to make sure
    # that the standard errors are okay too
    if(type=="response" & (is.null(se.fit) || !se.fit)){
      result <- result*off.set
    }else if(type=="response" & se.fit){
      result$fit <- result$fit*off.set
      result$se.fit <- result$se.fit*off.set
    }
  }

  return(result)
}
