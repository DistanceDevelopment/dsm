#' Summarize a fitted density surface model
#'
#' Gives a brief summary of a fitted [`dsm`][dsm] object.
#'
#' @method summary dsm
#' @aliases summary.dsm
#'
#' @param object a model fitted by [`dsm`][dsm]
#' @param \dots other arguments passed to [`summary.gam`][mgcv::summary.gam]
#' @return a summary object
#' @export
#'
#' @seealso [`dsm`][dsm]
#' @author David L. Miller
summary.dsm <- function(object, ...){
  class(object) <- class(object)[class(object)!="dsm"]
  NextMethod("summary", object)
}
