#' Print a description of a density surface model variance object
#' 
#' This method only provides a short summary, use the 
#' \code{\link{summary.dsm.var}} method for information.
#' 
#' @S3method print dsm.var
#' @aliases print.dsm.var
#' @method print dsm.var
#' @param x a \code{dsm} variance object
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return NULL
#' @author David L. Miller
#' @seealso \code{\link{summary.dsm.var}}
#' @keywords utility
print.dsm.var<-function(x,...){

  if(x$bootstrap){
    cat("Result of bootstrap uncertainty estimate\n")

    if(!x$ds.uncertainty){
      cat("Detection function uncertainty incorporated using the delta method.\n\n")

    }else{
      cat("Detection function uncertainty incorporated into boostrap.\n\n")
    }

    cat("Replicates        :",x$n.boot,"\n")
    #cat("Usable replicates : ",x$boot.usable,
    #                          " (",100-(100*x$trim.prop),"%)\n",sep="")

  }else if(!x$bootstrap){
    cat("Summary of uncertainty in a density surface model calculated\n")
    cat(" by variance propagation.\n")

  }

  cat("\n")

  invisible()
}
