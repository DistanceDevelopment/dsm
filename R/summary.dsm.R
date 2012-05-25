#' Summary for a fitted dsm model
#'
#' Gives a brief summary of a fitted \code{dsm} object. For more specific
#' information on the detection function and spatial model, call 
#' \code{summary()} on the \code{$ddf} and \code{$result} elements of the 
#' \code{dsm} object.
#'
#' @S3method summary dsm
#' @method summary dsm
#' @aliases summary.dsm
#' 
#' @param object a \code{dsm} object
#' @param \dots unused arguments for S3 compatibility
#' @return a summary object
#' 
#' @author David L. Miller
#'
summary.dsm<-function(object,...){

  # storage
  sinfo<-list()


  # what did we fit?
  sinfo$model.spec<-object$model.spec

  # summary for the detection function
  sinfo$ddf<-summary(object$ddf)
  sinfo$ddf$model.description<-ddf.model.description(object$ddf)

  # summary for the spatial model
  if(sinfo$model.spec$model == "GLM"){
    sinfo$glm<-summary(object$result)
  }else if(sinfo$model.spec$model == "GAM"){
    # for speed set freq=TRUE since we aren't going to use the Bayesian
    # confidence intervals
    sinfo$gam<-summary(object$result,freq=TRUE)
  }

  class(sinfo)<-"summary.dsm"

  return(sinfo)

}
