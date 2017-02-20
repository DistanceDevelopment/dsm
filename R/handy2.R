# functions from Mark V Bravington's "handy2" library
# all code belongs to Mark!

"gam.fixed.priors" <-
function( ..., fixed.priors=character(0), debugging=FALSE, scale.trace=1){
  gamf <- gam
  formals( gamf)[ names( formals( sys.function())) %except% '...'] <- 
      formals( sys.function()) %without.name% '...'
  mc <- as.list( match.call( definition=gamf, expand.dots=TRUE))
  mc[[1]] <- quote( gam)
  fp <- fixed.priors
  mc$fixed.priors <- NULL
  
  check.ok <- function( x, must.be) {
      if( !is.null( mc[[x]]) && mc[[x]] != must.be)
        warning( x %&% "will be set to " %&% must.be)
      mc[[ x]] <<- must.be
    }
  check.ok( 'fit', TRUE)
  check.ok( 'G', NULL)
  check.ok( 'method', 'REML')
  check.ok( 'optimizer', c( 'outer', 'newton'))
  
  paraPen <- eval( mc$paraPen, parent.frame())
  if( !all( fp %in% names( paraPen)))
stop( "Fixed priors must correspond to things in paraPen (for now)")

  for( ifp in fp)
    paraPen[[ ifp]]$sp <- 1e8
  mc$paraPen <- paraPen
  
  #   was: gam( form.with.paraPen, paraPen=prior.dpar.fixo, knots=knots, method='REML', fit=FALSE)
  mc$fit <- FALSE
  mc <- as.call( mc)
  G.fix <- eval( mc, parent.frame()) 

  # For initial scale, fit with prior cfts set almost to zero
  b0 <- gam( G=G.fix, method='REML', scale=-1) 
  if( debugging) {
    mc$fit <- TRUE
    spunknown <- b0$sp*0-1
    pf <- parent.frame()
  }
  
  # Set up for 1D optimization that will start in the right place
  phi <- (sqrt(5)-1)/2
  untrans <- function( tx) {
    # Ensure start val is at initial fit
    # 0 -> 0; 1 -> Inf; 1-phi -> b0$sig2
    b0$sig2 * (phi/(1-phi)) * ( tx / (1-tx))
  }
  
  b.para1 <- NULL
  gamcrit <- function( trans.current.scale) {
    current.scale <- untrans( trans.current.scale)

    # Fix prior sp "internally"; this is not pretty
    G.fix$lsp0[ fp] <- log( current.scale) # NOT log( 1/current.scale) ... 
    b.para1 <<- gam( G=G.fix, method='REML', scale=current.scale)

    if( debugging) {
      for( ifp in fp)
        mc$paraPen[[ ifp]]$sp <- 1/current.scale
      mc$scale <- current.scale
      b.direct1 <- eval( mc, pf)
      # gam( form=form.with.paraPen, knots=knots, paraPen=prior.dpar.fixo, method='REML', 
      #        scale=current.scale)

      # Am suspicious of whether 'sp' is treated consistently when included in 'paraPen' vs
      # ... set directly
      # Have to leave for now, since it's hard to set correct elts of sp
      # for( ifp in fp)
      #   mc$paraPen[[ ifp]]$sp <- NULL
      # mc$sp <- replace( spunknown, fp, 1/current.scale)
      # b.direct2 <- eval( mc, pf)
      #gam( form=form.with.paraPen, knots=knots, paraPen=prior.dpar, method='REML', 
      #        sp=c( dlinkE.dpar=1/current.scale, spunknown), scale=current.scale)
    }
  
    marg.lglk <- -b.para1$gcv.ubre # ?? is this the REML marg lglk ??
    if( scale.trace>0)
      cat( 'Current scale: ',  current.scale, 'Marg lglk: ', marg.lglk, '\n')
  return( marg.lglk) 
  }

  if( debugging)
    mtrace( gamcrit)
  opto <- optimize( NEG( gamcrit), interval=0:1)  
  
  b.para1$call$fixed.priors <- fp
  b.para1$call[[1]] <- quote( gam.fixed.priors)
return( b.para1)
  
  # Also need a version where the scale param is fixed but the other smoopars aren't, ...
  # ... e.g. Poisson + soap + paraPen
}

"NEG" <-
function( f) { 
  if( is.null( f))
return( f) # useful for

  if( is.primitive( f)) {
    fargs <- formals( args( f)) # primitives don't have formals
    argo <- lapply( names( fargs), as.name)
    gbod <- list( as.name( '-'), as.call( c( list( substitute( f)), argo)))
    g <- function() 0
    body( g) <- as.call( gbod)
    formals( g) <- fargs
    environment( g) <- .GlobalEnv
  } else {
    # f is normal function
    g <- f
    body( g) <- substitute( { 
      mc <- match.call()
      mc[[1]] <- f
      -eval( mc, parent.frame())
    }, list( f=f))
    formals( g) <- formals( f)
    environment( g) <- environment( f) 
  }
return( g)
}

"%without.name%" <-
function( x, what) {
  new.names <- names( x) %except% what
  x[ new.names]
}


"%except%" <-
function (vector, condition)
vector[match(vector, condition, 0) == 0]


numderiv <- function(f, x0, eps = 1e-04, TWICE. = TRUE, param.name = NULL,
                     ..., SIMPLIFY = TRUE){
  if(is.null(param.name)){
    ff <- function(x, ...) f(x, ...)
  }else{
    ff <- function(x, ...){
      ll <- c(list(x), list(...))
      names(ll)[1] <- param.name
      do.call("f", ll)
    }
  }
  f0 <- ff(x0, ...)
  n <- length(x0)
  m <- matrix(0, length(f0), n)
  for(i in 1:n){
    if(x0[i] == 0){
      this.eps <- eps
    }else{
      this.eps <- eps * x0[i]
    }
    m[, i] <- (ff(x0 + this.eps * (1:n == i), ...) - f0)/this.eps
  }
  if(!is.null(dim(f0))){
    dim(m) <- c(dim(f0), n)
  }
  if(TWICE.){
    mc <- match.call()
    mc$eps <- -eps
    mc$TWICE. <- FALSE
    m <- 0.5 * (m + eval(mc, sys.frame(sys.parent())))
  }
  if(any(dim(m) == length(m)) && SIMPLIFY){
    m <- c(m)
  }
  return(m)
}

# from mvbutils
"%**%"<-function(x, y){
  dimnames(x) <- NULL
  dimnames(y) <- NULL

  if(length(dim(x)) == 2 && length(dim(y)) == 2 && dim(x)[2] ==
     1 && dim(y)[1] == 1){
    return(c(x) %o% c(y))
  }

  if((!is.null(dim(x)) && any(dim(x) == 1))){
    dim(x) <- NULL
  }

  if((!is.null(dim(y)) && any(dim(y) == 1))){
    dim(y) <- NULL
  }
  if(is.null(dim(x)) && is.null(dim(y))){
    if(length(x) == length(y)){
      x <- x %*% y
    }else{
      if ((length(x) != 1) && (length(y) != 1)){
          stop(paste("lengths of x (",length(x),") and y (",
            length(y),") are incompatible",sep=""))
      }else{
        x <- x * y
      }
     }
  }else{
    x <- x %*% y
  }

  if((!is.null(dim(x)) && any(dim(x) == 1))){
    dim(x) <- NULL
  }
  return(x)
}
