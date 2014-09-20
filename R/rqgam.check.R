#' Randomised quantile residuals check plots for GAMs/DSMs
#'
#' Function operates as \code{\link{gam.check}} but using randomised quantile
#' residuals, a la Dunn and Smyth (1996). Checks of \code{k} are not computed,
#' these need to be done using \code{\link{gam.check}}.
#'
#' @param gam.obj a \code{gam}, \code{glm} or \code{dsm} object.
#' @param ... arguments passed on to all plotting functions
#' @return just plots!
#'
#' @author Based on code provided by Natalie Kelly, bugs added by Dave Miller
#' @export
#' @import statmod
#'
#' @examples
#'
#' library(Distance)
#' library(dsm)
#'
#' # load the Gulf of Mexico dolphin data (see ?mexdolphins)
#' data(mexdolphins)
#'
#' # fit a detection function and look at the summary
#' hr.model <- ds(mexdolphins$distdata, max(mexdolphins$distdata$distance),
#'                key = "hr", adjustment = NULL)
#'
#' # fit a simple smooth of x and y
#' mod1<-dsm(N~s(x,y), hr.model, mexdolphins$segdata, mexdolphins$obsdata)
#' rqgam.check(mod1)
rqgam.check<-function(gam.obj,...){

  # for negbin need to set $theta
  if(grepl("Negative Binomial",gam.obj$family$family)){
    gam.obj$theta <- gam.obj$family$getTheta()
  }


  # layout stuff
  opar <- par(mfrow=c(2,2))

  # grab the randomised quantile residuals
  # requires statmod package

  # need to do the right thing for mgcv's Tweedie
  if(grepl("^Tweedie",mod1$family$family)){
    qres <- qres.tweedie(gam.obj)
  }else{
    qres <- qresid(gam.obj)
  }


  # values of the linear predictor
  linpred <- napredict(gam.obj$na.action, gam.obj$linear.predictors)

  ## normal Q-Q plot
  qqnorm(qres,ylab="Randomised quantile residuals",...)

  ## resids vs. linear pred
  plot(linpred, qres,main="Resids vs. linear pred.",
         xlab="linear predictor",ylab="Randomised quantile residuals",...)

  ## histogram
  hist(qresid(gam.obj), main="Histogram of residuals",
       xlab="Randomised quantile residuals",...)


  ## Response vs. Fitted Values
  plot(fitted(gam.obj), gam.obj$y,
       main="Response vs. Fitted Values",
       xlab="Fitted Values", ylab="Response",...)
  #lines(lowess(gam.obj$fitted.values, gam.obj$model[,1]), col = 2)


  #residuals versus leverage
  #plot(gam.obj$hat, qresid(gam.obj), main="QRES versus leverage", xlab="Leverage", ylab="QRES")
  #lines(lowess(gam.obj$hat, qresid(gam.obj)), col = 2)
  #plot(gam.obj$fitted.values, sqrt(abs(qres)),
  #     main="Scale-Location (QRES)",
  #     xlab="Fitted", ylab="Sqrt(Abs(QRES))" )
  ##lines(lowess(gam.obj$fitted.values, sqrt(abs(qresid(gam.obj)))), col = 2)


  par(opar)
}
