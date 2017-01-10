#' Print summary of density surface model variance object
#'
#' See \code{\link{summary.dsm.var}} for information.
#'
#' @param x a summary of \code{dsm} variance object
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return NULL
#' @export
#' @author David L. Miller
#' @seealso \code{\link{summary.dsm.var}}
#' @keywords utility
print.summary.dsm.var<-function(x, ...){

  if(x$bootstrap){
    cat("Summary of bootstrap uncertainty in a density surface model\n")

    if(!x$ds.uncertainty){
      cat("Detection function uncertainty incorporated using the delta method.\n\n")

    }else{
      cat("Detection function uncertainty incorporated into boostrap.\n\n")
    }

    cat("Boxplot coeff     :",x$boxplot.coef,"\n")
    cat("Replicates        :",x$n.boot,"\n")
    cat("Outliers          :",x$boot.outliers,"\n")
    cat("Infinites         :",x$boot.infinite,"\n")
    cat("NAs               :",x$boot.NA,"\n")
    cat("NaNs              :",x$boot.NaN,"\n")
    cat("Usable replicates : ",x$boot.usable,
                              " (",100*(1-x$trim.prop),"%)\n",sep="")

    # Note that when we don't include detection function uncertainty in
    # our bootstrap then we need to do the log-Normal approx with the
    # CV calculated from the bootstrap, as the quantiles don't include
    # any of the uncertainty in the detection function.
    if(!x$ds.uncertainty){
      # asymptotic CI
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

      cat("Approximate asymptotic bootstrap confidence interval:\n")
      print(asymp.tot)
      cat("(Using log-Normal approximation)\n")

    }else{
      # DSU percentile CI
      cat("\nPercentile bootstrap confidence interval and median:\n")
      print(x$quantiles)
    }


  }else if(!x$bootstrap){
    cat("Summary of uncertainty in a density surface model calculated\n")
    if(x$varprop){
      cat(" by variance propagation.\n")
      cat("\nQuantiles of differences between fitted model and variance model\n")
      print(x$saved$model.check)
    }else{
      cat(" analytically for GAM, with delta method\n")

    }

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

  }

  cat("\n")
  cat("Point estimate                 :", x$pred.est,"\n")
  # print the individual CVs if we used the delta method
  if(x$bootstrap){
    if(!x$ds.uncertainty){
      if(!is.null(x$detfct.cv)) cat("CV of detection function       :",x$detfct.cv,"\n")
      cat("CV from bootstrap              :", round(x$bootstrap.cv,4),"\n")
      cat("Total standard error           :", x$se,"\n")
      cat("Total coefficient of variation :", round(x$cv,4),"\n")
    }else{
      cat("Standard error                 :", x$se,"\n")
      cat("Coefficient of variation       :", round(x$cv,4),"\n")
    }
  }else{
    if(x$varprop | is.null(x$saved$dsm.object$ddf)){
      cat("Standard error                 :", x$se,"\n")
      cat("Coefficient of variation       :", round(x$cv,4),"\n")
    }else{
      if(!is.null(x$detfct.cv)) cat("CV of detection function       :",x$detfct.cv,"\n")
      cat("CV from GAM                    :", round(x$gam.cv,4),"\n")
      cat("Total standard error           :", x$se,"\n")
      cat("Total coefficient of variation :", round(x$cv,4),"\n")
    }
  }
  cat("\n")


  if(!is.null(x$subregions)){
    i<-1
    for(subreg in x$subregions){
      cat("Subregion",i,"\n")
      print(subreg)
      cat("\n")
      i <- i+1
    }
  }

  invisible()
}
