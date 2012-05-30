#' Generate data from a fitted detection function
#'
#' When \code{ds.uncertainty} is \code{TRUE}, this procedure generates data
#' from the fitted detection function (assuming that it is correct).
#'
#' @param ds.object a fitted detection function object (as returned by a call
#'        to \code{ddf.ds()}.
#'
#'
#' @note This function changes the random number generator seed. To avoid any
#'       potential side-effects, use something like:
#'       \code{seed <- get(".Random.seed",envir=.GlobalEnv)}
#'       \code{### Run some code}
#'       \code{assign(".Random.seed",seed,envir=.GlobalEnv)}
#' @export
#' @author David L. Miller
generate.ds.uncertainty<-function(ds.object){

  n.ds.samples<-length(
                ds.object$data$distance[
                        ds.object$data$distance<=ds.object$ds$aux$width])

  # how many samples do we have so far?
  n.samps<-0
  dists<-c()
  
  # create an object to hold the parameters
  pars<-list()
  pars$scale<-ds.object$ds$aux$ddfobj$scale$parameters
  if(!is.null(ds.object$ds$aux$ddfobj$shape$parameters)){
    pars$shape<-ds.object$ds$aux$ddfobj$shape$parameters
  }
  if(!is.null(ds.object$ds$aux$ddfobj$adjustment$parameters)){
    pars$adjustment<-ds.object$ds$aux$ddfobj$adjustment$parameters
  }
  
  # make sure that a model gets fitted
  dud.df<-TRUE

  while(dud.df){
    # if we just have a half-normal key function then we can
    # directly simulate...
    if(ds.object$ds$aux$ddfobj$type=="hn" &
        is.null(ds.object$ds$aux$ddfobj$adjustment$parameters)){
    
      dists<-abs(rnorm(n.ds.samples,mean=0,sd=exp(pars$scale)))
    
    # otherwise we need to do some rejection sampling
    }else{
      # need to call out to mrds to get the data and model objects
      # into the correct format
      xmat <- mrds:::process.data(new.dists,ds.object$meta.data,
                                  check=FALSE)$xmat
      ddfobj <- mrds:::create.ddfobj(ds.object$call$dsmodel,xmat,
                          ds.object$meta.data,pars)

      while(n.samps < n.ds.samples){
    
        # generate some new distances
        new.dists<-data.frame(distance=runif(n.ds.samples-n.samps)*
                                          ds.object$ds$aux$width,
                              detected=rep(1,n.ds.samples-n.samps),
                              object=1:(n.ds.samples-n.samps))
    
        # generate acceptance probability
        U<-runif(n.ds.samples-n.samps)
    
        # do the rejection...
        # (evaluate the -log(L) then backtransform per-observation)
        # ONLY line transect at the moment!!
        inout <- exp(-mrds:::flt.lnl(ds.object$par,ddfobj,
                      misc.options=list(width=ds.object$ds$aux$width,
                                        int.range=ds.object$ds$aux$int.range,
                                        showit=FALSE, doeachint=TRUE,
                                        point=ds.object$ds$aux$point,
                                        integral.numeric=TRUE),TCI=FALSE))>U
        dists<-c(dists,new.dists$distance[inout])
    
        n.samps<-length(dists)
      }
    }
    # make sure that we got the right number
    dists<-dists[1:n.ds.samples]
    dists<-data.frame(distance=dists,
                      detected=rep(1,length(dists)),
                      object=1:length(dists))

    # fit the model to the new data
    ddf.call <- ds.object$call
    ddf.call$data <- dists
    ddf.call$meta.data <- ds.object$meta.data
    ddf.fitted <- try(eval(ddf.call))

    if(all(class(ddf.fitted)!="try-error")){
      dud.df <- FALSE
    }

  }


  # return the offset
  # in the future this could be the offset from a MCDS model too
  #return(rep(fitted(ddf.fitted,compute=TRUE,esw=TRUE)[1],length(N.round)))
  return(fitted(ddf.fitted,compute=TRUE,esw=TRUE)[1])

}