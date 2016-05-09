#' Print a description of a density surface model variance object
#'
#' This method only provides a short summary, see \code{\link{summary.dsm_varprop}}.
#'
#' @param x a \code{dsm} variance object
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return NULL
#' @export
#' @author David L. Miller
#' @seealso \code{\link{summary.dsm_varprop}}
#' @keywords utility
print.dsm.var<-function(x, ...){

  print(summary(x))

}
