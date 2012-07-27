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
#' @param boxplot.coef the value of \code{coef} used to calculate the outliers
#'        see \code{\link{boxplot}}.
#' @param \dots unused arguments for S3 compatibility
#' @return a summary object
#' 
#' @author David L. Miller
#'
summary.dsm.var<-function(object, alpha=0.05, boxplot.coef=1.5, ...){

  # storage
  sinfo<-list()

  if(object$bootstrap){
    # grab the predicted values
    mod1.pred <- dsm.predict(object$dsm.object,
                             newdata=object$pred.data,
                             off=object$off.set)
    sinfo$pred.est <- sum(mod1.pred,na.rm=TRUE)

    #  short.var=short.var
    #  study.area.total=study.area.total

    sinfo$block.size <- object$block.size 
    sinfo$n.boot <- object$n.boot
    sinfo$bootstrap <- TRUE
    sinfo$ds.uncertainty <- object$ds.uncertainty

    # bootstrap abundances
    bootstrap.abund <- object$study.area.total


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
      # save that
      sinfo$detfct.cv <- sqrt(cvp.sq)


      # save the s.e. of N from bootstrap
      trimmed.variance <- trim.var(bootstrap.abund[is.finite(bootstrap.abund)],
                                   boxplot.coef=boxplot.coef)
      sinfo$N.bs.se <- sqrt(trimmed.variance)

      # cv squared of the Ns from the bootstrap
      cvNbs.sq <- (sinfo$N.bs.se/sinfo$pred.est)^2
      # save that
      sinfo$bootstrap.cv <- sqrt(cvNbs.sq)

      # cv of N
      cvN <- sqrt(cvp.sq+cvNbs.sq)
      sinfo$cv <- cvN

      # variance (delta method)
      sinfo$var <- (cvN*sinfo$pred.est)^2
      sinfo$se <- sqrt(sinfo$var)
    }else{
      # if we used detection function uncertainty

      # variance of the bootstrap abundances is the variance
      trimmed.variance <- trim.var(bootstrap.abund[is.finite(bootstrap.abund)],
                                   boxplot.coef=boxplot.coef)

      sinfo$var <- trimmed.variance
      sinfo$se <- sqrt(trimmed.variance)

      sinfo$cv <- sinfo$se/sinfo$pred.est
    }

    ### general bootstrap stuff

    # how many duds did we have?
    sinfo$boxplot.coef <- boxplot.coef
    sinfo$trim.prop <- attr(trimmed.variance, "trim.prop")
    sinfo$trim.ind <- attr(trimmed.variance, "trim.ind")
    sinfo$boot.outliers <- attr(trimmed.variance, "outliers")
    sinfo$boot.infinite <- sum(is.infinite(bootstrap.abund))
    sinfo$boot.finite <- sum(!is.infinite(bootstrap.abund))
    sinfo$boot.NA <- sum(is.na(bootstrap.abund))
    sinfo$boot.NaN <- sum(is.nan(bootstrap.abund))
    sinfo$boot.usable <- sinfo$boot.finite - sinfo$boot.outliers 

    # grab the %ile c.i.s at alpha, 1-alpha and also median
    sinfo$quantiles <- quantile(bootstrap.abund[sinfo$trim.ind], 
                                c(alpha, 0.5, 1-alpha),na.rm=TRUE)
    attr(sinfo$quantiles,"names")[2] <- "Median"
    

  }else if(object$bootstrap==FALSE){
    # varprop stuff
    sinfo$saved<-object
    sinfo$bootstrap <- object$bootstrap

    # what if we had multiple areas (ie this is from a CV plot?)
    if(all(dim(object$pred.var)==1,1)){
      sinfo$se <- sqrt(object$pred.var)
    }else{
      # re run the variance calculation, putting everything together
      pd<-c()
      off<-c()
      for(i in 1:length(object$pred.data)){
        pd<-rbind(pd,object$pred.data[[i]])
        off<-rbind(off,object$off.set[[i]])
      }
      object$pred.data <- pd
      object$off.set <- as.vector(off)

      # grab the predicted values
      mod1.pred <- dsm.predict(object$dsm.object,
                               newdata=object$pred.data,
                               off=object$off.set)
      sinfo$pred.est <- sum(mod1.pred,na.rm=TRUE)

      var.prop <- dsm.var.prop(object$dsm.obj, object$pred.data, object$off.set,
                               object$seglen.varname, object$type.pred)

      
      sinfo$se <- sqrt(var.prop$pred.var)

    }

    sinfo$cv <- sinfo$se/sinfo$pred.est
  }

  class(sinfo) <- "summary.dsm.var"

  return(sinfo)

}
