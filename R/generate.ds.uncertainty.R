#' Generate data from a fitted detection function
#'
#' When using [`dsm.var.movblk`][dsm.var.movblk] if `ds.uncertainty=TRUE`, this
#' procedure generates data from the fitted detection function (assuming that
#' it is correct).
#'
#' @param ds.object a fitted detection function object (as returned by a call
#' to [`ddf.ds`][mrds::ddf.ds]).
#'
#'
#' @note This function changes the random number generator seed. To avoid any
#' potential side-effects, use something like: `seed <-
#' get(".Random.seed",envir=.GlobalEnv)` before running code and
#' `assign(".Random.seed",seed,envir=.GlobalEnv)` after.
#' @export
#' @author David L. Miller
#' @importFrom stats runif rnorm fitted
generate.ds.uncertainty <- function(ds.object){

  n.ds.samples <- length(ds.object$data$distance[ds.object$data$distance <=
                                                 ds.object$ds$aux$width])

  # how many samples do we have so far?
  n.samps <- 0
  dists <- c()

  # create an object to hold the parameters
  pars <- list()
  pars$scale <- ds.object$ds$aux$ddfobj$scale$parameters
  if(!is.null(ds.object$ds$aux$ddfobj$shape$parameters)){
    pars$shape <- ds.object$ds$aux$ddfobj$shape$parameters
  }
  if(!is.null(ds.object$ds$aux$ddfobj$adjustment$parameters)){
    pars$adjustment <- ds.object$ds$aux$ddfobj$adjustment$parameters
  }

  # make sure that a model gets fitted
  dud.df <- TRUE

  while(dud.df){
    # if we just have a half-normal key function then we can
    # directly simulate...
    if(ds.object$ds$aux$ddfobj$type=="hn" &
       is.null(ds.object$ds$aux$ddfobj$adjustment$parameters)){

      dists <- abs(rnorm(n.ds.samples, mean=0, sd=exp(pars$scale)))

    # otherwise we need to do some rejection sampling
    }else{

      # since rejection sampling is time consuming, generate lots of
      # samples at once, we re-scale the number by the inverse of the 
      # ratio accepted. The first time over, let's make that 10x
      n.mult <- 10

      width <- ds.object$ds$aux$width

      while(n.samps < n.ds.samples){

        # how many samples should we take?
        this.n.samps <- n.mult*(n.ds.samples-n.samps)

        # generate some new distances
        xmat <- data.frame(distance = runif(this.n.samps, 0, width),
                           detected = rep(1, this.n.samps),
                           object   = 1:(this.n.samps),
                           binned   = rep(FALSE, this.n.samps))

        # create a ddf object
        ddfobj <- mrds::create.ddfobj(as.formula(ds.object$dsmodel), xmat,
                                      ds.object$meta.data,pars)

        # generate acceptance probability
        U <- runif(this.n.samps)

        ### do the rejection...
        ### (evaluate the -log(L) then backtransform per-observation)
        ### ONLY line transect at the moment!!
        ##inout <- U < exp(-mrds:::flpt.lnl(ds.object$par,ddfobj,
        ##              misc.options=list(width=width,
        ##                                int.range=ds.object$ds$aux$int.range,
        ##                                showit=FALSE, doeachint=TRUE,
        ##                                point=ds.object$ds$aux$point,
        ##                                integral.numeric=TRUE)))

        inout <- U <= mrds::detfct(xmat$distance, ddfobj,
                                   standardize=FALSE, width=width)

        dists <- c(dists, xmat$distance[inout])

        n.samps <- length(dists)

        # update the number of extra samples we make by inverting the ratio
        # of accepted to generated this round
        n.mult <- ceiling(1/(sum(inout)/this.n.samps))
      }
    }

    # make sure that we got the right number
    dists <- dists[1:n.ds.samples]
    dists <- data.frame(distance = dists,
                        detected = rep(1, length(dists)),
                        object   = seq_len(length(dists)))

    # fit the model to the new data
    ddf.call <- ds.object$call
    ddf.call$data <- dists
    ddf.call$meta.data <- ds.object$meta.data
    ddf.call$dsmodel <- as.formula(ds.object$dsmodel)
    ddf.fitted <- try(with(ds.object, eval(ddf.call)))

    # if it all went well, then set dud.df to FALSE and quit the loop
    if(all(class(ddf.fitted)!="try-error")){
      dud.df <- FALSE
    }else{
    # otherwise forget everything and start again
      n.samps <- 0
      dists <- c()
    }
  }

  # return the offset
  return(fitted(ddf.fitted, compute=TRUE, esw=TRUE)[1])
  # (in the future this could be the offset from a MCDS model too)
}
