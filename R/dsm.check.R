#' Diagnostic checks for \code{dsm} models.
#'
#' Produces four plots: quantile-quantile plot, scale-location plot
#' (optionally with LOESS line), fit to the residuals and empirical 
#' variogram.
#'
#' @param dsm.obj object resulting from a call to \code{\link{dsm.fit}}.
#' @param type the type of residuals to use for all of the plots.         

#' @param rep argument to be passed to \code{\link{qq.gam}} (default 0).
#' @param level argument to be passed to \code{\link{qq.gam}} (default 0.9).
#' @param rl.col argument to be passed to \code{\link{qq.gam}} (default 2).
#' @param rep.col argument to be passed to \code{\link{qq.gam}} (default 
#' "gray80").
#' @param loess should the LOESS smooth through the scale-location plot be
#' shown? (default \code{TRUE}).
#' @param gam.k complexity of the \code{\link{gam}} to be fitted to the 
#' residuals, see \code{\link{choose.k}} for more information (default 30).
#' @param vario.max maximum distance for the variogram; points further than 
#' this distance apart will be ignored.
#' @param ... other arguments to be passed to \code{\link{qq.gam}}.

#' @return a plot!
#' @author David L. Miller
#' @export
#
# analogous to gam.check for gam() objects, but dsm specific
# much of this is butchered from Simon Wood's gam.check()
# some other bits taken from Marra, Miller and Zanin (2011)
# what to produce?
# - qq
# - scale-location (with LOESS?)
# - fit to residuals 
# - empirical variogram
dsm.check<-function(dsm.obj,type=c("deviance","pearson","response"),
          ## arguments passed to qq.gam() {w/o warnings !}:
          rep=0, level=.9, rl.col=2, rep.col="gray80", loess=TRUE, 
          # What is k for the gam?
          gam.k=30, vario.max=100,...){

  # TODO
  # better way of selecting k? -- uniquecombs()?

  ### first, pull out the GAM part of the model
  if(all(class(dsm.obj)!="gam")){
    model<-dsm.obj$result
  }else{
    model<-dsm.obj
  }

  ### pre-amble, get things in order
  type <- match.arg(type)
  old.par<-par(mfrow=c(2,2))
  fitted.vals<-fitted(model)
  resids<-residuals(model,type=type)


  ### QQ-plot
  qq.gam(model, rep=rep, level=level, type=type, 
         rl.col=rl.col, rep.col=rep.col, ...)

  ### scale-location plot
  sl.dat<-data.frame(x=fitted.vals,y=abs(resids))
  plot(sl.dat,las=1,
       main="Scale-location plot",
       ylab="Abs. value of residuals",
       xlab="Predicted values",cex=0.3)

  if(loess){
    # loess fit..
    loe<-loess(y~x,data=sl.dat)
    nd<-seq(min(sl.dat$x,na.rm=T),max(sl.dat$x,na.rm=T))
    pred<-predict(loe,newdata=nd,by=0.01) 
    lines(nd,pred,col="grey")
  }

  ### fit to residuals -- check for residual spatial variation
  new.dat<-data.frame(x=model$data$x,
                      y=model$data$y,
                      z=resids)
  b<-gam(z~s(x,y,k=gam.k)-1,data=new.dat)
  vis.gam(b,plot.type="contour",main="Fit to residuals",
          asp=1,view=c("x","y"),type="response") 


  ### variogram
  # save a lot of things!
  all.dat<-c()

  for(tranid in unique(model$data$Transect.Label)){
    ind <- model$data$Transect.Label==tranid

    # let's assume that Sample.Label has the format
    #  Transect.Label-number
    this.dat <- data.frame(Effort = model$data$Effort[ind],
                           Sample.Label = model$data$Sample.Label[ind],
                           resids = resids[ind])

    # create a column of the just number of the segment within the transect
    this.dat$Sample.Number <- as.numeric(sub("\\d+-","",
                                         as.character(this.dat$Sample.Label)))
    # sort the data.frame by that
    this.dat <- this.dat[order(this.dat$Sample.Number),]

    # now pretend that that each segment lies along y=0, with distances
    # between the segments as the sum of the efforts.
    fake.dat <- data.frame(x=cumsum(this.dat$Effort),
                           resids=this.dat$resids)

    all.dat <- rbind(all.dat, fake.dat)
  }

  # calculate the acf...
  acf.fit <- acf(all.dat,plot=FALSE,type="covariance")

  plot(acf.fit$lag[,,2][,1], acf.fit$acf[,,2][,2],type="l", 
       xlab="lag",ylab="correlation", main="Autocorrelogram")

  #plot(all.dat)

  par(old.par)

  invisible()
#return(acf.fit)
#return(all.dat)
}
