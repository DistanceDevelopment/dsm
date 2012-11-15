#' Summarize a fitted density surface model
#'
#' Gives a brief summary of a fitted \code{dsm} object.
#'
#' @S3method summary dsm
#' @method summary dsm
#' @aliases summary.dsm
#' 
#' @param object a \code{dsm} object
#' @param \dots other arguments passed to \code{\link{summary.gam}}.
#' @return a summary object
#'
#' @seealso dsm
#' @author David L. Miller
#'
summary.dsm<-function(object,...){
  class(object) <- class(object)[class(object)!="dsm"]
  return(summary(object))
}