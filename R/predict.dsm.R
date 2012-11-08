#' Predict from a fitted density surface model
#'
#' Make predictions outside (or inside) the covered area.
#'
#' @param object: a fitted \code{\link{dsm}} object as produced by \code{gam()}.
#' @param newdata spatially referenced covariates e.g.
#'                sea temperature, depth, distance to shore, etc.
#'                Note covariates in this dataframe must have names *identical*·
#'                to variable names used in fitting the GAM.
#' @param off.set area of each of the cells in the prediction grid. Ignored if 
#'                if their is already a column in \code{newdata} called
#'                \code{off.set}.
#'                The link function is applied to the offset inside this 
#'                function.
#' @param \dots any other arguments passed to \code{\link{predict.gam}} or
#'              \code{\link{predict.glm}}.
#' @return predicted values on the response scale (density/abundance).
#'
#' @S3method predict dsm
#' @method predict dsm
#' @aliases predict.dsm
#'
#' @author David L. Miller
predict.dsm <- function(object, newdata, off.set=NULL, ...){


  if(is.null(newdata$off.set)){
    newdata$off.set <- off.set
  }

  # apply the link function
  linkfn <- object$family$linkfun
  newdata$off.set <- linkfn(newdata$off.set)

  # remove the dsm class
  class(object) <- class(object)[class(object)!="dsm"]

  # actually do the predict call
  result<-predict(object, newdata, type="response",...)

  return(result)
}
