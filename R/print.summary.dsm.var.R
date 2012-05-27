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
    cat("Usable replicates : ",x$boot.usable,
                              " (",100-(100*x$trim.prop),"%)\n",sep="")

    cat("\nPercentrile bootstrap confidence interval and median:\n")
    if(!x$ds.uncertainty){
      cat("(Spatial model component only.)\n")
    }
    print(x$quantiles)


  }else if(!x$bootstrap){
    cat("Summary of uncertainty in a density surface model calculated\n")
    cat(" by variance propagation.\n")


  #model 
  #deriv
    cat("\nQuantiles of differences between fitted model and variance model\n")
    print(x$saved$model.check)

  }

  cat("\n\n")
  cat("Point estimate           :", x$pred.est,"\n")
  cat("Standard error           :", x$se,"\n")
  cat("Coefficient of variation :", round(x$cv,4),"\n")

  cat("\n")

  invisible()
}
