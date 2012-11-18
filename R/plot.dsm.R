#' Plot a density surface model.
#'
#' See \code{\link{plot.gam}}.
#'
#' @S3method plot dsm
#' @method plot dsm
#' @aliases plot.dsm
#'
#' @param x a \code{dsm} object
#' @param \dots other arguments passed to \code{\link{plot.gam}}.
#' @return a plot!
#'
#' @seealso dsm plot.gam
#' @author David L. Miller
#'
plot.dsm<-function(x,...){
  class(x) <- class(x)[class(x)!="dsm"]
  return(plot(x,...))
}
