#' Summarize the variance of a density surface model
#'
#' Gives a brief summary of a fitted [`dsm_varprop`][dsm_varprop] variance
#' object.
#'
#' @param object a `dsm.var` object
#' @param alpha alpha level for confidence intervals (default 0.05 to give a
#' 95% confidence internal)
#' @param ... unused arguments for S3 compatibility
#' @return a summary object
#' @export
#'
#' @seealso [`dsm_varprop`][dsm_varprop], [`summary.dsm.var`][summary.dsm.var]
#' @author David L. Miller
#' @importFrom stats qnorm update.formula
summary.dsm_varprop <- function(object, alpha=0.05, ...){

  if(all(is.na(object$pred)) & all(is.na(object$var)) & all(is.na(object$ses))){
    stop("Cannot get summary when newdata=NULL, rerun dsm_varprop with prediction data")
  }

  # storage
  sinfo <- list()
  # save the alpha value for cis
  sinfo$alpha <- alpha

  sinfo$varprop <- TRUE
  sinfo$saved<-object

  # abundance estimate
  sinfo$pred.est <- sum(object$pred, na.rm=TRUE)
  sinfo$var <- object$var

  # calculate the CV for the whole model
  sinfo$cv <- sqrt(sinfo$var)/sinfo$pred.est

  if(!all(class(object$old_model$ddf) == "list")){
    object$old_model$ddf <- list(object$old_model$ddf)
  }

  sinfo$detfct.cv <- c()

  for(i in seq_along(object$old_model$ddf)){

    # detection function CV too
    ddf.summary <- summary(object$old_model$ddf[[i]])

    cvp.sq <- (ddf.summary$average.p.se/
               ddf.summary$average.p)^2
    sinfo$detfct.cv <- c(sinfo$detfct.cv, sqrt(cvp.sq))
  }

  # save model check diagnostic
  sinfo$varprop_diagnostic <- varprop_check(object)

  class(sinfo) <- "summary.dsm_varprop"
  return(sinfo)
}
