#' Print summary of density surface model variance object
#'
#' See \code{\link{summary.dsm_varprop}} for information.
#'
#' @param x a summary of \code{dsm} variance object
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return NULL
#' @export
#' @author David L. Miller
#' @seealso \code{\link{summary.dsm.var}}
#' @keywords utility
print.summary.dsm_varprop<-function(x, ...){

  cat("Summary of uncertainty in a density surface model calculated\n")
  cat(" by variance propagation.\n")
  cat("\nProbability of detection in fitted model and variance model\n")
  print(x$varprop_diagnostic)

  cat("\n")

  ## calculate the CI around the abundance estimate

    # this doesn't transform N, only se(N)
    # this probably should do the following:
    # lower =
    #  qlnorm(alpha/2, log(x$pred.est) - 0.5*log(x$cv^2+1),
    #         sqrt(log(x$cv^2+1)))
    # upper =
    #  qlnorm(1-alpha/2, log(x$pred.est) - 0.5*log(x$cv^2+1),
    #         sqrt(log(x$cv^2+1)))
  unconditional.cv.square <- x$cv^2
  asymp.ci.c.term <- exp(qnorm(1-x$alpha/2) *
                         sqrt(log(1+unconditional.cv.square)))
  asymp.tot <- c(x$pred.est / asymp.ci.c.term,
                 x$pred.est,
                 x$pred.est * asymp.ci.c.term)

  names(asymp.tot) <- c(paste0(x$alpha/2*100, "%"),
                        "Mean",
                        paste0((1-x$alpha/2)*100,"%"))

  cat("Approximate asymptotic confidence interval:\n")
  print(asymp.tot)
  cat("(Using log-Normal approximation)\n")

  cat("\n")
  cat("Detection function CV          :", paste(round(x$detfct.cv, 4),
                                                collapse=", "), "\n")
  cat("\n")
  cat("Point estimate                 :", x$pred.est,"\n")
  cat("Standard error                 :", sqrt(x$var),"\n")

  cat("Coefficient of variation       :", round(x$cv,4),"\n")
  cat("\n")

  invisible()
}
