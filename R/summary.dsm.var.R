#' Summarize the varianc of a density surface model
#'
#' Gives a brief summary of a fitted \code{dsm} variance object. 
#'
#' @S3method summary dsm.var
#' @method summary dsm.var
#' @aliases summary.dsm.var
#' 
#' @param object a \code{dsm.var} object
#' @param alpha alpha level for confidence intervals
#' @param \dots unused arguments for S3 compatibility
#' @return a summary object
#' 
#' @author David L. Miller
#'
summary.dsm.var<-function(object, alpha=0.05, ...){

  sinfo<-list()

  if(object$bootstrap){

    #  short.var=short.var
    #  study.area.total=study.area.total
    pred.data <- object$pred.data

    sinfo$block.size <- object$block.size 
    sinfo$n.boot <- object$n.boot
    sinfo$bootstrap <- TRUE
    sinfo$ds.uncertainty <- object$ds.uncertainty

    # bootstrap abundances
    bootstrap.abund <- object$study.area.total

    # estimate from prediction
    mod1.pred <- dsm.predict(object$dsm.object,
                             newdata=pred.data,
                             off=object$off.set)
    pred.est <- sum(mod1.pred,na.rm=TRUE)

    # delta method, if necessary
    if(!object$ds.uncertainty){

      ddf.object<-object$dsm.object$ddf
      ddf.summary<-summary(ddf.object)

      # average p standard error
      sinfo$average.p.se <- ddf.summary$average.p.se

      ## calculate the variance via the delta method
      # find the cv squared of the p
      cvp.sq <- (ddf.summary$average.p.se/
                 ddf.summary$average.p)^2

      # cv squared of the Ns from the bootstrap
      cvNbs.sq <- (sqrt(trim.var(bootstrap.abund[is.finite(bootstrap.abund)]))/
                   pred.est)^2
#                   mean(bootstrap.abund[
#                       is.finite(bootstrap.abund)],na.rm=TRUE))^2

      # save the s.e. of N from bootstrap
      trimmed.variance <- trim.var(bootstrap.abund[
                                          is.finite(bootstrap.abund)])
      sinfo$N.bs.se <- sqrt(trimmed.variance)

      # cv of N
      cvN <- sqrt(cvp.sq+cvNbs.sq)
      sinfo$boot.cv <- cvN

      # variance (delta method)
      sinfo$var <- (cvN*pred.est)^2
      sinfo$se <- sqrt(sinfo$var)
    }else{
      # if we used detection function uncertainty

      # variance of the bootstrap abundances is the variance
      trimmed.variance <- trim.var(bootstrap.abund[
                                          is.finite(bootstrap.abund)])

      sinfo$var <- trimmed.variance
      sinfo$se <- sqrt(trimmed.variance)

      sinfo$boot.cv <- cvN
    }

    # general bootstrap stuff
    sinfo$trim.prop <- attr(trimmed.variance, "trim.prop")
    sinfo$boot.outliers <- attr(trimmed.variance, "outliers")
    sinfo$boot.infinite <- sum(is.infinite(bootstrap.abund))
    sinfo$boot.finite <- sum(!is.infinite(bootstrap.abund))
    sinfo$boot.NA <- sum(is.na(bootstrap.abund))
    sinfo$boot.NaN <- sum(is.nan(bootstrap.abund))

# don't use this at the moment
#    sinfo$boot.median <- median(bootstrap.abund[
#                                    is.finite(bootstrap.abund)],
#                                na.rm=TRUE)

    sinfo$quantiles <- quantile(bootstrap.abund[is.finite(bootstrap.abund)], 
                                c((1-alpha)/2, 1-((1-alpha)/2)),na.rm=TRUE)


    

  }else if(object$bootstrap==FALSE){
    # varprop stuff

  }

  class(sinfo) <- "summary.dsm.var"

  return(sinfo)

}
