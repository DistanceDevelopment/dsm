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
#' this distance apart will be ignored, see \code{\link{variog}} (default 100).
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
  qq.gam(model, rep=rep, level=level, type=type, rl.col=rl.col, rep.col=rep.col, ...)

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
  b<-gam(z~s(x,y,k=gam.k),data=new.dat)
  vis.gam(b,plot.type="contour",main="Fit to residuals",
          asp=1,view=c("x","y")) 

  ### variogram
  plot(1:10,1:10,type="n")
  text(5,5,label="coming soon")

  # need geoR for this!

  # overall variogram
#  coords<-matrix(0,length(model$model$x),2)
#  coords[,1]<-model$data$x
#  coords[,2]<-model$data$y
#  gb<-list(data=residuals(model,type="d"),coords=coords)
#  vg<-variog(gb,max.dist=vario.max,messages=FALSE)
#  #vg.env<-variog.mc.env(gb, obj.var = vg,messages=FALSE)
#  #plot(vg,envelope=vg.env,type="l",main="Emprical variogram",xlim=c(0,50),ylim=c(0,1))
#  plot(vg,type="l",main="Emprical variogram",xlim=c(0,50),ylim=c(0,1))


  #all.vg<-c()

  #for(tranid in unique(model$data$Transect.Label)){

  #  ind <- model$data$Transect.Label==tranid
  #  coords <- cbind(model$data$x[ind],model$data$y[ind])

  #  gb <- list(data=residuals(model,type="d")[ind],coords=coords)
  #  vg <- variog(gb,max.dist=vario.max,messages=FALSE)
  #  vg$x <- vg$u
  #  vg$y <- vg$v
  #  lines(vg,type="l",col=rgb(190,190,190,100,maxColorValue=255))

  #  all.vg<-rbind(all.vg,vg$y)

  #}

  #lines(x=vg$x,y=colMeans(all.vg),col="red")





  ##### OLD CODE
  ## taken from the Red Book
  #coords<-matrix(0,length(model$model$x),2)
  #coords[,1]<-model$data$x
  #coords[,2]<-model$data$y

  #gb<-list(data=residuals(model,type="d"),coords=coords)
  #vg<-variog(gb,max.dist=vario.max)

  #vg.env<-variog.mc.env(gb, obj.var = vg)

  ## plot the variogram
  #plot(vg,envelope=vg.env,type="l",main="Emprical variogram")
  ##### / OLD CODE

  par(old.par)
}
