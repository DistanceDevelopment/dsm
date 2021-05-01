#' Print a description of a density surface model variance object
#'
#' This method only provides a short summary, see
#' [`summary.dsm_varprop`][summary.dsm_varprop].
#'
#' @param x a `dsm` variance object
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return NULL
#' @export
#' @author David L. Miller
#' @seealso [`summary.dsm_varprop`][summary.dsm_varprop]
#' @keywords utility
print.dsm_varprop <- function(x, ...){
  print(summary(x))
}
