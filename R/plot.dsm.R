#' Plot a density surface model.
#'
#' See [`plot.gam`][mgcv::plot.gam].
#'
#' @aliases plot.dsm
#'
#' @param x a model fitted by [`dsm`][dsm]
#' @param \dots other arguments passed to [`plot.gam`][mgcv::plot.gam].
#' @return a plot!
#' @export
#'
#' @seealso [`dsm`][dsm] [`plot.gam`][mgcv::plot.gam]
#' @author David L. Miller
#' @importFrom graphics plot
#'
plot.dsm <- function(x, ...){
  class(x) <- class(x)[class(x)!="dsm"]
  return(plot(x, ...))
}
