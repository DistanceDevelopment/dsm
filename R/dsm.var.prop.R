#' Variance propogation for DSM models
#'
#' Rather than use a bootstrap to calculate the variance in a \code{dsm} model,
#' use the clever variance propogation trick from Williams et al. (2011).
#'
#' The idea is to refit the spatial model but including the Hessian of the 
#' offset as an extra term. Variance estimates using this new model can then 
#' be used to calculate the variance of abundance estimates which incorporate 
#' detection function uncertainty. Further mathematical details are given in 
#' the paper in the references below.
#'
#' Many prediction grids can be supplied by supplying a list of 
#' \code{data.frame}s to the function.
#' 
#' Based on (much more general) code from Mark Bravington and Sharon Hedley.
#'
#' @param dsm.obj an object returned from running \code{\link{dsm.fit}}.
#' @param pred.data either: a single prediction grid or list of prediction 
#'        grids. Each grid should be a \code{data.frame} with the same 
#'        columns as the original data.
#' @param off.set a a vector or list of vectors with as many elements as there 
#'        are in \code{pred.data}. Each vector is as long as the number of
#'        rows in the corresponding element of \code{pred.data}. These give
#'        the area associated with each prediction point. 
#' @param seglen.varname name for the column which holds the segment length 
#'        (default value "Effort"). 
#' @param type.pred should the predictions be on the "response" or "link" scale?
#'        (default "response").
#' @return a list with elements
#'         \tabular{ll}{\code{model} \tab the fitted model object\cr
#'                      \code{pred.var} \tab covariances of the regions given
#'                      in \code{pred.data}. Diagonal elements are the 
#'                      variances in order\cr
#'                      \code{bootstrap} \tab logical, always \code{FALSE}\cr
#'                      \code{pred.data} \tab as above\cr
#'                      \code{off.set} \tab as above\cr
#'                      \code{model}\tab the fitted model with the extra term\cr
#'                      \code{dsm.object} \tab the original model, as above\cr
#'                      \code{model.check} \tab simple check of subtracting the
#'                        coefficients of the two models to see if there is a
#'                        large difference\cr
#'                      \code{deriv} \tab numerically calculated Hessian of the
#'                        offset\cr.
#'                      }
#' @author Mark V. Bravington, Sharon L. Hedley. Bugs added by David L. Miller.
#' @references
#' Williams, R., Hedley, S.L., Branch, T.A., Bravington, M.V., Zerbini, A.N. and Findlay, K.P. (2011). Chilean Blue Whales as a Case Study to Illustrate Methods to Estimate Abundance and Evaluate Conservation Status of Rare Species. Conservation Biology 25(3), 526-535.
#' @export
dsm.var.prop<-function(dsm.obj, pred.data,off.set,
    seglen.varname='Effort', type.pred="response") {

  pred.data.save<-pred.data
  off.set.save<-off.set

  # make sure if one of pred.data and off.set is not a list we break
  # if we didn't have a list, then put them in a list so everything works
  if(is.data.frame(pred.data) & is.vector(off.set)){
    pred.data <- list(pred.data)
    off.set <- list(off.set)
#    pred.data[[1]] <- pred.data
#    off.set[[1]] <- off.set
  }else if(is.list(off.set)){
    if(length(pred.data)!=length(off.set)){
      stop("pred.data and off.set don't have the same number of elements")
    }
  }

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

  # function to find the derivatives of the offset
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
    dmat.name <- paste('.',dmat.name,sep="")
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

  # Diagnostic from Mark
  # check that the fitted model isn't too different, used in summary()
  model.check<-summary(fitted(fit.with.pen) - fitted(gam.obj))

  cft <- coef(fit.with.pen)
  #preddo <- numeric(length(pred.data))
  preddo <- list(length(pred.data))

  #names(preddo) <- names(pred.data) # if any
  dpred.db <- matrix(0, length(pred.data), length(cft))

  # depending on whether we have response or link scale predictions...
  if(type.pred=="response"){
      tmfn <- gam.obj$family$linkinv
      dtmfn <- function(eta){sapply(eta, numderiv, f=tmfn)}
  }else if(type.pred=="link"){
      tmfn <- identity
      dtmfn <- function(eta){1}
  }

  # loop over the prediction grids
  for( ipg in seq_along(pred.data)) {
    # if we have a single paramter model (e.g. half-normal) need to be careful
    if(is.matrix(firstD)){
      pred.data[[ipg]][[dmat.name]] <- matrix(0,
                                              nrow(pred.data[[ipg]]),
                                              ncol(firstD))
    }else{
      pred.data[[ipg]][[dmat.name]] <- rep(0, nrow(pred.data[[ipg]]))
    }

    ### fancy lp matrix stuff
    # set the offset to be zero here so we can use lp
    pred.data[[ipg]]$off.set<-rep(0,nrow(pred.data[[ipg]]))

    lpmat <- predict( fit.with.pen, newdata=pred.data[[ ipg]], type='lpmatrix')
    lppred <- lpmat %**% cft

    # if the offset is just one number then repeat it enough times 
    if(length(off.set[[ipg]])==1){
      this.off.set <- rep(off.set[[ipg]],nrow(pred.data[[ipg]]))
    }else{
      this.off.set <- off.set[[ipg]]
    }

    preddo[[ipg]] <-  this.off.set %**% tmfn(lppred)
    dpred.db[ipg,] <- this.off.set %**% (dtmfn(lppred)*lpmat)
    # explanation of the above line and why we find this derivative
    # BTW in 'varpred', there is a decoy option 'vmethod' which at the
    # moment has to be '"delta"', for how to deal with nonlinearity in the
    # "link" of 'fitobj'. Could be done by simulation instead, and that would
    # be more accurate (if you did enough). However, in my limited experience:
    # once you've got a CV so big that the delta-method doesn't work, then
    # your estimate is officially Crap and there is not much point in
    # expending extra effort to work out exactly how Crap!
  }

  # "'vpred' is the covariance of all the summary-things." - MVB
  # so we want the diagonals if length(pred.data)>1
  # A B A^tr
  vpred <- dpred.db %**% tcrossprod(vcov(fit.with.pen), dpred.db)

  result <- list(pred.var = vpred,
                 bootstrap = FALSE,
                 pred.data = pred.data.save,
                 pred = preddo,
                 off.set = off.set.save,
                 model = fit.with.pen,
                 dsm.object = dsm.obj,
                 model.check = model.check,
                 deriv = firstD,
                 seglen.varname=seglen.varname,
                 type.pred=type.pred
                )

  class(result) <- "dsm.var"

  return(result)
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
                stop(paste("lengths of x (",length(x),") and y (",
                  length(y),") are incompatible",sep=""))
            else x <- x * y
        }
    }
    else x <- x %*% y
    if ((!is.null(dim(x)) && any(dim(x) == 1)))
        dim(x) <- NULL
    x
}
