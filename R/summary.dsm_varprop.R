#' Summarize the variance of a density surface model
#'
#' Gives a brief summary of a fitted \code{dsm_varprop} variance object.
#'
#' @aliases summary.dsm.var
#'
#' @param object a \code{dsm.var} object
#' @param alpha alpha level for confidence intervals (default 0.05 to give a 95\% confidence internal)
#' @param \dots unused arguments for S3 compatibility
#' @return a summary object
#' @export
#'
#' @seealso dsm_varprop
#' @author David L. Miller
#' @importFrom stats qnorm
summary.dsm_varprop<-function(object, alpha=0.05, ...){

  # storage
  sinfo<-list()
  # save the alpha value for cis
  sinfo$alpha <- alpha

  sinfo$varprop <- TRUE
  sinfo$saved<-object

  # abundance estimate
  sinfo$pred.est <- sum(object$pred, na.rm=TRUE)
  sinfo$var <- object$var

  # calculate the CV for the whole model
  sinfo$cv <- sqrt(sinfo$var)/sinfo$pred.est

  # detection function CV too
  ddf.summary <- summary(object$old_model$ddf)

  cvp.sq <- (ddf.summary$average.p.se/
             ddf.summary$average.p)^2
  sinfo$detfct.cv <- sqrt(cvp.sq)

  # save model check diagnostic
  sinfo$model.check <- summary(fitted(object$refit) - fitted(object$old_model))

  class(sinfo) <- "summary.dsm_varprop"
  return(sinfo)
}
