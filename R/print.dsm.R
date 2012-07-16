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

  ### General information

  cat("\nDensity surface model\n")
  cat("Response : ", x$model.spec$response, "\n")
  cat("Offset   : ", x$model.spec$offset, "\n")

  cat("\nDetection function : ",ddf.model.description(x$ddf),"\n")

  ### GAM things...
  if(x$model.spec$model=="GAM"){
    cat("\nSummary of GAM\n")
    cat("\nFormula: ",as.character(x$model.spec$formula)[2],"~",
                      as.character(x$model.spec$formula)[-c(1,2)],"\n")
  }else if(x$model.spec$model=="GLM"){
    cat("Print not implemented for GLMs at the moment\n")
  }
  
  cat("\n")

  invisible()
}
