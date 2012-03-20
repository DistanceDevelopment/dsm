#' Variance propogation for DSM models
#'
#' Based on code from Mark Bravington and Sharon Hedley
#' @export

varprop<-function(fitobj, inputobj, pred.grids, 
# first arg is dsm object, second is mrds
    pred.wts=lapply( pred.grids, function( x) rep( 1/nrow(x), nrow( x))),
    vmethod='delta', seglen.varname='Effort', 
    type.pred=c( 'response', 'link'), ...) {
###############

  tweakParams <- function(object, params) {
    if(missing(params)){
      return(object$par)
    }
    object$par <- params
    object$ds$aux$ddfobj <- mrds:::assign.par(object$ds$aux$ddfobj,params)
    return(object)
  }
  p0 <- tweakParams(inputobj) # returns the parameters to numderiv

  callo <- fitobj$call  
  fo2data <- fitobj$data
  if( is.null( fo2data)) {
    fo2data <- eval( callo$data, parent.frame()) # hit 'n' hope
  }


  funco <- function(p){
    ipo <- tweakParams(inputobj, p)
    ret <- log(2 *unique(predict(ipo, esw=TRUE, compute=TRUE)$fitted)* 
            fo2data[[seglen.varname]])
    return(ret)
  }

  firstD <- numderiv( funco, p0)
  
  if( all( firstD==0))
    warning( 'Doffset/Dpars==0... really??!!') # will permit this, for testing!
  
  formo <- fitobj$formula
  dmat.name <- '.D1'
  names.to.avoid <- unique( c( all.names( formo), names( fo2data)))
  while( dmat.name %in% names.to.avoid)
    dmat.name <- '.' %&% dmat.name
  fo2data[[ dmat.name]] <- firstD
  
  formo[[3]] <- call( '+', formo[[3]], as.symbol( dmat.name))
  
  # Tell it to penalize...
  #paraterm <- list( list( getHessian( inputobj)))
  paraterm<-list(list(inputobj$hess))
  names( paraterm) <- dmat.name
  callo$formula <- formo # should use pmatch
  callo$paraPen <- c( callo$paraPen, paraterm)
  callo$data <- fo2data # ideally not, as above
  
  fit.with.pen <- eval( callo, parent.frame())

  ####
  # check that it doesn't change the model fit much!
  #scatn( 'Comparison of fitted values when offset is flexible:')
  #print( summary( fitted( fit.with.pen) - fitted( fitobj)))
  #scatn( "Should check param ests or fitted vals in 'fit.with.pen' are similar to original...")

  cft <- coef( fit.with.pen)
  preddo <- numeric( length( pred.grids))
  real.preds <- list()

  names( preddo) <- names( pred.grids) # if any
  dpred.db <- matrix( 0, length( pred.grids), length( cft))
  
  type.pred <- match.arg( type.pred) 
  { tmfn; dtmfn} %<-% switch( type.pred,
    'response'=  
      list( fitobj$family$linkinv, function( eta) sapply( eta, numderiv, f=tmfn)),
    'link'= 
      list( identity, function( eta) 1)
  )
  

  vmethod <- match.arg( vmethod)
  for( ipg in seq_along( pred.grids)) {
    if(is.matrix(firstD)){
      pred.grids[[ ipg]][[ dmat.name]] <- matrix( 0, nrow( pred.grids[[ ipg]]), ncol( firstD))
    }else{
      pred.grids[[ ipg]][[ dmat.name]] <- rep(0, nrow( pred.grids[[ ipg]]))
    }


    # fancy lp matrix stuff
    lpmat <- predict( fit.with.pen, newdata=pred.grids[[ ipg]], type='lpmatrix')
    lppred <- lpmat %**% cft
    if( vmethod=='delta') {
      preddo[[ ipg]] <- pred.wts[[ ipg]] %**% tmfn( lppred)
      dpred.db[ ipg,] <- pred.wts[[ ipg]] %**% (dtmfn( lppred) * lpmat)
    }
  } 

  # "'vpred' is the covariance of all the summary-things." - MVB  
  vpred <- dpred.db %**% tcrossprod( vcov( fit.with.pen), dpred.db) # A B A^tr 
  # inverse link

  return(list(gam=fit.with.pen,vpred=vpred,preddo=preddo))
}

# from Mark Bravington's handy2
numderiv<-function (f, x0, eps = 1e-04, TWICE. = TRUE, param.name = NULL, 
    ..., SIMPLIFY = TRUE) 
{
    if (is.null(param.name)) 
        ff <- function(x, ...) f(x, ...)
    else ff <- function(x, ...) {
        ll <- c(list(x), list(...))
        names(ll)[1] <- param.name
        do.call("f", ll)
    }
    f0 <- ff(x0, ...)
    n <- length(x0)
    m <- matrix(0, length(f0), n)
    for (i in 1:n) {
        this.eps <- eps * if (x0[i] == 0) 
            1
        else x0[i]
        m[, i] <- (ff(x0 + this.eps * (1:n == i), ...) - f0)/this.eps
    }
    if (!is.null(dim(f0))) 
        dim(m) <- c(dim(f0), n)
    if (TWICE.) {
        mc <- match.call()
        mc$eps <- -eps
        mc$TWICE. <- FALSE
        m <- 0.5 * (m + eval(mc, sys.frame(sys.parent())))
    }
    if (any(dim(m) == length(m)) && SIMPLIFY) 
        m <- c(m)
    return(m)
}

# from handy2
"%<-%"<-function( a, b){
  # a must be of the form '{thing1;thing2;...}'
  a <- as.list( substitute( a))[-1]
  e <- sys.parent()
  stopifnot( length( b) == length( a))
  for( i in seq_along( a))
    eval( call( '<-', a[[ i]], b[[i]]), envir=e)
  NULL
}


# from mvbutils
"%**%"<-function (x, y)
{
    dimnames(x) <- NULL
    dimnames(y) <- NULL
    if (length(dim(x)) == 2 && length(dim(y)) == 2 && dim(x)[2] ==
        1 && dim(y)[1] == 1)
        return(c(x) %o% c(y))
    if ((!is.null(dim(x)) && any(dim(x) == 1)))
        dim(x) <- NULL
    if ((!is.null(dim(y)) && any(dim(y) == 1)))
        dim(y) <- NULL
    if (is.null(dim(x)) && is.null(dim(y))) {
        if (length(x) == length(y))
            x <- x %*% y
        else {
            if ((length(x) != 1) && (length(y) != 1))
                stop("lengths of x (" %&% length(x) %&% ") and y (" %&%
                  length(y) %&% ") are incompatible")
            else x <- x * y
        }
    }
    else x <- x %*% y
    if ((!is.null(dim(x)) && any(dim(x) == 1)))
        dim(x) <- NULL
    x
}
