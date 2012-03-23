#' Variance propogation for DSM models
#'
#' Rather than use a bootstrap to calculate the variance in a \code{dsm} model,
#' use the clever variance propogation trick from Williams et al. (2011).
#'
#' The idea is to refit the spatial model but including the Hessian of the 
#' offset as an extra term. This new model can then be used to calculate 
#' the variance of abundance estimates which incorporate detection function
#' uncertainty.
#' 
#' Based on (much more general) code from Mark Bravington and Sharon Hedley.
#'
#' @param dsm.obj an object returned from running \code{\link{dsm.fit}}.
#' @param pred.grids list of prediction grids. Each element should be a 
#'        \code{data.frame} with the same columns as the original data.
#' @param pred.area a list with as many elements as there are in 
#'        \code{pred.grids}. Each element is a vector as long as the number of
#'        rows in the corresponding element of \code{pred.grids}. These give
#'        the area associated with each prediction point. Default gives equal 
#'        weighting to each point such that the "area" sums to one.
#' @param seglen.varname name for the column which holds the segment length 
#'        (default value "Effort"). 
#' @param type.pred should the predictions be on the "response" or "link" scale?
#'        (default "response").
#' @return a list with elements
#'         \tabular{ll}{\code{model} \tab the fitted model object\cr
#'                      \code{var.pred} \tab variance of
#'                      }
#' @author Mark Bravington, Sharon Hedley. Bugs added by David L Miller.
#' @references 
#' Williams, R., Hedley, S.L., Branch, T.A., Bravington, M.V., Zerbini, A.N. and Findlay, K.P. (2011). Chilean Blue Whales as a Case Study to Illustrate Methods to Estimate Abundance and Evaluate Conservation Status of Rare Species. Conservation Biology 25(3), 526–535.
#' @export

varprop<-function(dsm.obj, pred.grids, 
    pred.area=lapply( pred.grids, function( x) rep( 1/nrow(x), nrow( x))),
    seglen.varname='Effort', type.pred="response", ...) {

  # pull out the ddf object
  ddf.obj <- dsm.obj$ddf
  # and the gam
  gam.obj <- dsm.obj$result

  # this function changes the parameters in the ddf object
  tweakParams <- function(object, params) {
    if(missing(params)){
      return(object$par)
    }
    object$par <- params
    object$ds$aux$ddfobj <- mrds:::assign.par(object$ds$aux$ddfobj,params)
    return(object)
  }

  # function to find the derivatives of -- the offset
  funco <- function(p){
    # set the parameters to be p
    ipo <- tweakParams(ddf.obj, p)
    # calculate the offset
    ret <- log(2 *unique(predict(ipo, esw=TRUE, compute=TRUE)$fitted)* 
            fo2data[[seglen.varname]])
    return(ret)
  }

  # pull out the data and the call
  callo <- gam.obj$call  
  fo2data <- dsm.obj$data

  # find the derivatives
  p0 <- tweakParams(ddf.obj) # returns the parameters to numderiv
  firstD <- numderiv( funco, p0)
  
  # if the derivatives were zero, throw an error
  if(all(firstD==0)){
    stop('Doffset/Dpars==0... really??!!') 
  }

  # now construct the extra term...
  formo <- gam.obj$formula
  dmat.name <- '.D1'
  names.to.avoid <- unique( c( all.names( formo), names( fo2data)))
  while( dmat.name %in% names.to.avoid){
    dmat.name <- '.' %&% dmat.name
  }
  fo2data[[ dmat.name]] <- firstD 
  formo[[3]] <- call( '+', formo[[3]], as.symbol(dmat.name))
  # put it all together
  paraterm<-list(list(ddf.obj$hess))
  names(paraterm) <- dmat.name
  callo$formula <- formo 
  callo$family<-gam.obj$family
  callo$paraPen <- c(callo$paraPen, paraterm)
  callo$data <- fo2data 
  # run the model
  fit.with.pen <- eval(callo, parent.frame())

  #### Diagnostics from Mark
  # check that it doesn't change the model fit much!
  #scatn( 'Comparison of fitted values when offset is flexible:')
  #print( summary( fitted( fit.with.pen) - fitted( gam.obj)))
  #scatn( "Should check param ests or fitted vals in 'fit.with.pen' are similar to original...")
  #### /Diagnostics from Mark

  cft <- coef( fit.with.pen)
  preddo <- numeric( length( pred.grids))
  real.preds <- list()

  names( preddo) <- names( pred.grids) # if any
  dpred.db <- matrix( 0, length( pred.grids), length( cft))
  
  # depending on whether we have response or link scale predictions...
  if(type.pred=="response"){
      tmfn<-gam.obj$family$linkinv
      dtmfn<-function( eta){sapply( eta, numderiv, f=tmfn)}
  }else if(type.pred=="link"){ 
      tmfn<-identity
      dtmfn<-function( eta){1}
  }

  # loop over the prediction grids
  for( ipg in seq_along( pred.grids)) {
    # if we have a single paramter model (e.g. half-normal) need to be careful
    if(is.matrix(firstD)){
      pred.grids[[ ipg]][[ dmat.name]] <- matrix( 0, nrow( pred.grids[[ ipg]]), ncol( firstD))
    }else{
      pred.grids[[ ipg]][[ dmat.name]] <- rep(0, nrow( pred.grids[[ ipg]]))
    }

    # fancy lp matrix stuff
    lpmat <- predict( fit.with.pen, newdata=pred.grids[[ ipg]], type='lpmatrix')
    lppred <- lpmat %**% cft
    preddo[[ ipg]] <- pred.area[[ ipg]] %**% tmfn( lppred)
    dpred.db[ ipg,] <- pred.area[[ ipg]] %**% (dtmfn( lppred) * lpmat)
  } 

  # "'vpred' is the covariance of all the summary-things." - MVB  
  # so we want the diagonals if length(pred.grids)>1
  vpred <- dpred.db %**% tcrossprod( vcov( fit.with.pen), dpred.db) # A B A^tr 

  return(list(model=fit.with.pen,pred.var=vpred))#,pred=preddo))
}


####### this is all utility stuff below here, taken from Mark's packages

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

# from mvbutils
"%**%"<-function(x, y){
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
