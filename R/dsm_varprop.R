#' Variance propagation for density surface models
#'
#' Calculate the uncertainty in predictions from a fitted DSM, including uncertainty from the detection function.
#'
#' When we make predictions from a spatial model, we also want to know the uncertainty about that abundance estimate. Since density surface models are 2 (or more) stage models, we need to incorporate the uncertainty from the earlier stages (i.e. the detection function) into our "final" uncertainty estimate.
#'
#' This function will refit the spatial model but include the Hessian of the offset as an extra term. Variance estimates using this new model can then be used to calculate the variance of predicted abundance estimates which incorporate detection function uncertainty. Importantly this requires that if the detection function has covariates, then these do not vary within a segment (so, for example covariates like sex cannot be used).
#'
#' For more information on how to construct the prediction grid \code{data.frame}, \code{newdata}, see \code{\link{predict.dsm}}.
#'
#' This routine is only useful if a detection function has been used in the DSM.
#'
#' Note that we use \code{Vc} here (see \code{\link{gamObject}}), which is the variance-covariance matrix for the spatial model, corrected for smoothing parameter uncertainty.
#'
#' @return a list with elements
#' \tabular{ll}{\code{old_model} \tab fitted model supplied to the function as \code{model}\cr
#'              \code{refit} \tab refitted model object, with extra term\cr
#'              \code{pred} \tab point estimates of predictions at \code{newdata}\cr
#'              \code{var} \tab total variance calculated over all of \code{newdata}\cr
#'              \code{ses} \tab standard error for each prediction cell in \code{newdata}\cr
#'  }
#' @author David L. Miller, based on code from Mark V. Bravington and Sharon L. Hedley.
#' @references
#' Williams, R., Hedley, S.L., Branch, T.A., Bravington, M.V., Zerbini, A.N. and Findlay, K.P. (2011). Chilean Blue Whales as a Case Study to Illustrate Methods to Estimate Abundance and Evaluate Conservation Status of Rare Species. Conservation Biology 25(3), 526-535.
#'
#'
#' @param model a fitted \code{\link{dsm}}
#' @param newdata the prediction grid
#' @param trace for debugging, see how the scale parameter estimation is going
#' @export
#'
#' @examples
#' \dontrun{
#'  library(Distance)
#'  library(dsm)
#'
#' # load the Gulf of Mexico dolphin data (see ?mexdolphins)
#' data(mexdolphins)
#' attach(mexdolphins)
#'
#' # fit a detection function and look at the summary
#' df <- ds(distdata, 7000, key = "hr", adjustment = NULL, formula=~beaufort)
#'
#' # fit a simple smooth of x and y
#' obsdata <- obsdata[obsdata$distance<=7000,]
#' mod1 <- dsm(N~s(x,y), df, segdata, obsdata, method="REML")
#'
#' # Calculate the variance
#' preddata$off.set <- preddata$area
#' mod1.varp <- dsm_varprop(mod1, preddata)
#'
#' # this will give a summary over the whole area in mexdolphins$preddata
#'
#' # detach the data
#' detach("mexdolphins")
#' }
dsm_varprop <- function(model, newdata, trace=FALSE){

  # die if the link isn't log
  if(model$family$link != "log"){
    stop("log link must be used!")
  }

  # extract the link & invlink
  linkfn <- model$family$linkfun
  linkinvfn <- model$family$linkinv

  # die if we're not using REML
  if(model$method != "REML"){
    stop("REML must be used for smoothing parameter selection")
  }

  # extract the call
  this_call <- as.list(model$call)
  # remvoe the function
  this_call[1] <- NULL

  # extract the detection function
  ddf <- model$ddf

  # function to differentiate
  mu_fn <- function(par, linkfn, ddf, data, ds_newdata){
    # set the detection function parameters to be par
    ddf$par <- par

    # calculate mu (effective strip width)
    mu <- predict(ddf, newdata=ds_newdata, esw=TRUE, compute=TRUE)$fitted

    # calculate log mu
    ret <- linkfn(mu)
    return(ret)
  }

  # extract the formula
  ds_formula <- ddf$ds$aux$ddfobj$scale$formula

  # if we don't have covariates then just setup the
  # data for predict.ds to be a distance of 0, which will be
  # ignored anyway
  if(ds_formula=="~1"){
    ds_newdata <- data.frame(distance=rep(0, nrow(model$data)))
  }else{
    # otherwise need the covars that are in the data (that we saved
    # with keepData=TRUE :))
    ds_newdata <- model$data[, all.vars(as.formula(ds_formula)), drop=FALSE]
    ds_newdata$distance <- 0
  }

  # probably a lot of duplicated stuff in the above, so let's just
  # pass the unique combinations
  u_ds_newdata <- mgcv::uniquecombs(ds_newdata)

  # find the derivatives of log(mu)
  firstD <- as.matrix(numderiv(mu_fn, ddf$par, linkfn=linkfn, ddf=ddf,
                               data=model$data, ds_newdata=u_ds_newdata))

  # repopulate with the duplicates back in
  firstD <- firstD[attr(u_ds_newdata, "index"), , drop=FALSE]
  # calculate the offset 2 * effective strip width * line length
  #  -- log(mu) + log(2*effort)
  firstD <- firstD + linkfn(2 * model$data$Effort)

  # put that in the data
  dat <- model$data
  dat[["XX"]] <- firstD

  # do a clever thing above to ensure we don't clobber
  # some other variable name

  ## build the model
  # update the formula to include the new term
  this_call$formula <- update.formula(this_call$formula,
                         paste0(paste(as.character(this_call$formula)[c(2,1,3)],
                                 collapse=" "), " + XX"))
  # update data
  this_call$data <- dat
  # add paraPen bit
  # hessian needs to be 2nd REAL Hessian not DS thing
  # that lives in $hessian
  opt_details <- attr(ddf$ds,"details")
  if(is.matrix(opt_details)){
    hess <- opt_details[nrow(opt_details),]$nhatend
  }else{
    hess <- opt_details$nhatend
  }
  if(any(is.na(hess))){
    # fall back to DS use if things are bad
    hess <- ddf$hessian
  }
  this_call$paraPen <- c(this_call$paraPen, list(XX=list(hess)))
  # tell gam.fixed.priors what to look for
  this_call$fixed.priors <- "XX"

  # set the trace on the scale parameter estimation
  this_call$scale.trace <- trace

  ## refit the model
  refit <- do.call("gam.fixed.priors", this_call)

  ## now do some predictions

  # make some zeros for the paraPen term so they have no mean effect
  newdata[["XX"]] <- matrix(0, nrow(newdata), ncol(firstD))

  # get everyone's favourite matrix
  Lp <- predict(refit, newdata=newdata, type="lpmatrix")
  # predictions on the link scale
  pred <- Lp %*% coef(refit)
  pred <- newdata$off.set * linkinvfn(pred)

  # get variance-covariance with smoothing parameter uncertainty
  Vc <- refit$Vc

  # this is why we can only use log link
  dNdbeta <- pred%**%Lp

  # make a sandwich
  var_p <- dNdbeta %*% Vc %*% dNdbeta

  # apply the link function to the offset
  newdata$off.set <- linkfn(newdata$off.set)
  # get the standard errors
  ses <- predict(refit, newdata=newdata, type="response", se.fit=TRUE)$se.fit

  # what should we return?
  ret <- list(old_model = model,
              refit     = refit,
              pred      = pred,
              var       = var_p,
              ses       = ses)

  class(ret) <- "dsm_varprop"

  return(ret)
}
