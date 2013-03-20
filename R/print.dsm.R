#' Print a description of a density surface model object
#'
#' This method just gives a short description of the fitted model. Use the
#' \code{\link{summary.dsm}} method for more information.
#'
#' @S3method print dsm
#' @aliases print.dsm
#' @method print dsm
#' @param x a \code{dsm} object
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return NULL
#' @author David L. Miller
#' @seealso \code{\link{summary.ds}}
#' @keywords utility
print.dsm<-function(x,...){

  # the code here is chopped together from mgcv and mrds

  # if we fit a gamm then we should just grab teh gam bit for this
  if("gamm" %in% class(x)){
    x<-x$gam
  }

  ### General information

  cat("\nDensity surface model\n")
  cat("Response : ",as.character(x$formula)[2] , "\n")

  if(as.character(x$formula)[2]!="presence"){
    cat("\nDetection function : ",ddf.model.description(x$ddf),"\n")
  }

  cat("\nFormula: ",as.character(x$formula)[2],"~",
                    as.character(x$formula)[-c(1,2)],"\n")
  cat("\n")

  invisible()
}
