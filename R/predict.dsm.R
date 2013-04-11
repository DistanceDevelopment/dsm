#' Predict from a fitted density surface model
#'
#' Make predictions outside (or inside) the covered area.
#'
#' @param object a fitted \code{\link{dsm}} object as produced by \code{dsm()}.
#' @param newdata spatially referenced covariates e.g.
#'                sea temperature, depth, distance to shore, etc.
#'                Note covariates in this dataframe must have names *identical*
#'                to variable names used in fitting the GAM.
#' @param off.set area of each of the cells in the prediction grid. Ignored if 
#'                if their is already a column in \code{newdata} called
#'                \code{off.set}.
#'                The link function is applied to the offset inside this 
#'                function.
#' @param type what scale should the results be on. The default is 
#'             \code{"response"} and is almost always what you want.
#' @param \dots any other arguments passed to \code{\link{predict.gam}} or
#'              \code{\link{predict.glm}}.
#' @return predicted values on the response scale (density/abundance).
#'
#' @S3method predict dsm
#' @method predict dsm
#' @aliases predict.dsm
#'
#' @author David L. Miller
predict.dsm <- function(object, newdata=NULL, off.set=NULL,
                        type="response",...){

  if("gamm" %in% class(object)){
    object <- object$gam
  }

  if(is.null(newdata)){
    newdata <- object$data
  }

  # if we don't have a density model, then set the offset
  if(!(c(object$formula[[2]]) %in% c("D","presence","density"))){
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

  # remove the dsm class
  class(object) <- class(object)[class(object)!="dsm"]

  # actually do the predict call
  result<-predict(object, newdata, type=type,...)

  return(result)
}
