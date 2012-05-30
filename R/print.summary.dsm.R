#' Print summary of density surface model object
#' 
#' See \code{\link{summary.dsm}} for information.
#' 
#' @S3method print summary.dsm
#' @aliases print.summary.dsm
#' @method print summary.dsm
#' @param x a summary of \code{dsm} object
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return NULL
#' @author David L. Miller
#' @seealso \code{\link{summary.ds}}
#' @keywords utility
print.summary.dsm<-function(x,...){

  # the code here is chopped together from mgcv and mrds

  ### General information

  cat("\nSummary of density surface model\n")
  cat("Response : ", x$model.spec$response, "\n")
  cat("Offset   : ", x$model.spec$offset, "\n")

  ### Detection function information
  cat("\nSummary of detection function\n")
  cat("\nForm: ",x$ddf$model.description,"\n")
  cat("Number of observations : ", x$ddf$n,"\n")
  cat("Distance range         : ", x$ddf$left, " - ",x$ddf$width,"\n")
  cat("Detection function AIC : ", x$ddf$aic, "\n")

  cat("\n")

  # information on average p
  parameters <- data.frame(Estimate=c(x$ddf$average.p))
  row.names(parameters) <- c("Average p")
  if(!is.null(x$ddf$average.p.se)){
    parameters$SE <- c(x$ddf$average.p.se)
    parameters$CV <- parameters$SE/parameters$Estimate
  }
  print(parameters)

  # was monotonicity enforced?
  if(x$ddf$mono & x$ddf$mono.strict){
    cat("\nStrict monotonicity constraints were enforced.\n")
  }else if(x$ddf$mono){
    cat("\nMonotonicity constraints were enforced.\n")
  }

  ### GAM things...
  if(x$model.spec$model=="GAM"){
    cat("\nSummary of GAM\n")
    cat("\nFormula: ",as.character(x$model.spec$formula),"\n")

    cat("\n")

    cat("Number of segments                   :",x$model.spec$n.segs,"\n")
    cat("Number of segments with observations :",x$model.spec$n.segs.withdata,
          "(",100*x$model.spec$n.segs.withdata/x$model.spec$n.segs,"%)\n")
    cat("R-sq.(adj)             :",formatC(x$gam$r.sq,digits=3,width=5),"\n")
    if(length(x$gam$dev.expl)>0){
      cat("Deviance explained     :",
          formatC(x$gam$dev.expl*100,digits=3,width=4),"%\n",sep="")
    }
    if(!is.null(x$gam$method)&&!(x$gam$method%in%c("PQL","lme.ML","lme.REML"))){
      cat(x$gam$method," score            : ",
          formatC(x$gam$sp.criterion,digits=5),"\n",sep="")
    }

  }else if(x$model.spec$model=="GLM"){
    cat("Model summary not implemented for GLMs at the moment\n")
  }
  
  cat("\n")

  invisible()
}
