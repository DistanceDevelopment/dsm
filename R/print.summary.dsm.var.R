#' Print summary of density surface model variance object
#' 
#' See \code{\link{summary.dsm.var}} for information.
#' 
#' @S3method print summary.dsm.var
#' @aliases print.summary.dsm.var
#' @method print summary.dsm.var
#' @param x a summary of \code{dsm} variance object
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return NULL
#' @author David L. Miller
#' @seealso \code{\link{summary.ds.var}}
#' @keywords utility
print.summary.dsm.var<-function(x,...){

  if(x$bootstrap){
    cat("Summary of bootstrap uncertainty in a density surface model\n")

    if(!x$ds.uncertainty){
      cat("Detection function uncertainty incorporated using the delta method.\n\n")

    }else{
      cat("Detection function uncertainty incorporated into boostrap.\n\n")
    }

    cat("Replicates        :",x$n.boot,"\n")
    cat("Outliers          :",x$boot.outliers,"\n") 
    cat("Infinites         :",x$boot.infinite,"\n") 
    cat("NAs               :",x$boot.NA,"\n")       
    cat("NaNs              :",x$boot.NaN,"\n")      
    cat("Usable replicates : ",x$boot.finite,
                              " (",100-(100*x$trim.prop),"%)\n",sep="")

    cat("\nPercentile method confidence intervals for density surface component only:\n")
    print(x$quantiles)

    cat("\n\n")
    cat("Point estimate :", x$pred.est,"\n",
        "Standard error :", x$se,"\n")

    cat("\n       Estimated CV for density surface model",
              round(x$boot.cv,4),"\n")

  }


  cat("\n")

  invisible()
}
