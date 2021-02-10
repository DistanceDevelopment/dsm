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
#' This routine is only useful if a detection function with covariates has been used in the DSM.
#'
#' Note that we can use \code{var_type="Vc"} here (see \code{\link{gamObject}}), which is the variance-covariance matrix for the spatial model, corrected for smoothing parameter uncertainty. See Wood, Pya & S{\"a}fken (2016) for more information.
#'
#' Models with fixed scale parameters (e.g., negative binomial) do not require an extra round of optimisation.
#'
#' @section Diagnostics:
#' The summary output from the function includes a simply diagnostic that shows the average probability of detection from the "original" fitted model (the model supplied to this function; column \code{Fitted.model}) and the probability of detection from the refitted model (used for variance propagation; column \code{Refitted.model}) along with the standard error of the probability of detection from the fitted model (\code{Fitted.model.se}), at the unique values of any factor covariates used in the detection function (for continuous covariates the 5%, 50% and 95% quantiles are shown). If there are large differences between the probabilities of detection then there are potentially problems with the fitted model, the variance propagation or both. This can be because the fitted model does not account for enough of the variability in the data and in refitting the variance model accounts for this in the random effect.
#'
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
#' Bravington, M.V., Miller, D.L. and Hedley, S.L. (2019) Reliable variance propagation for spatial density surface models. https://arxiv.org/abs/1807.07996
#'
#' Williams, R., Hedley, S.L., Branch, T.A., Bravington, M.V., Zerbini, A.N. and Findlay, K.P. (2011). Chilean Blue Whales as a Case Study to Illustrate Methods to Estimate Abundance and Evaluate Conservation Status of Rare Species. Conservation Biology 25(3), 526-535.
#'
#' Wood, S.N., Pya, N. and S{\"a}fken, B. (2016) Smoothing parameter and model selection for general smooth models. Journal of the American Statistical Association, 1-45.
#'
#' @param model a fitted \code{\link{dsm}}
#' @param newdata the prediction grid
#' @param trace for debugging, see how the scale parameter estimation is going
#' @param var_type which variance-covariance matrix should be used (\code{"Vp"} for variance-covariance conditional on smoothing parameter(s), \code{"Vc"} for unconditional). See \code{\link{gamObject}} for an details/explanation. If in doubt, stick with the default, \code{"Vp"}.
#' @export
#'
# @examples
# \dontrun{
# library(Distance)
# library(dsm)
#
# # load the Gulf of Mexico dolphin data (see ?mexdolphins)
# data(mexdolphins)
#
# # fit a detection function
# df <- ds(distdata, max(distdata$distance),
#          key = "hn", adjustment = NULL)
#
# # fit a simple smooth of x and y
# mod1 <- dsm(count~s(x, y), df, segdata, obsdata, family=tw())
#
# # Calculate the variance
# preddata$off.set <- preddata$area
# mod1.varp <- dsm_varprop(mod1, preddata)
# summary(mod1.varp)
# # this will give a summary over the whole area in mexdolphins$preddata
# }
dsm_varprop <- function(model, newdata, trace=FALSE, var_type="Vp"){

  # die if the link isn't log
  if(model$family$link != "log"){
    stop("log link must be used!")
  }

  # check for valid var_type
  if(!(var_type %in% c("Vp","Vc"))){
    stop("var_type must be \"Vp\" or \"Vc\"")
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
  # remove the function
  this_call[1] <- NULL

  # extract the detection function(s)
  ddf <- model$ddf
  if(all(class(ddf)!="list")){
    ddf <- list(ddf)
    # work around for dsms fitted in previous versions
    if(is.null(model$data$ddfobj)){
      model$data$ddfobj <- 1
    }
  }

  # get detection function info
  parskel <- list()
  for(i in seq_along(ddf)){
    parskel[[i]] <- ddf[[i]]$par
  }
  npars <- lapply(parskel, length)

  # new function to differentiate
  mu_fn <- function(par, linkfn, ddf, data, ds_newdata, skel){

    # put the parameter flesh back on the list bones?
    fleshy <- relist(par, skel)

    ret <- rep(0, nrow(data))

    # loop over the detection functions we have
    for(i in seq_along(ddf)){
      this_ddf <- ddf[[i]]


      # if we don't have a real detection function
      if("fake_ddf" %in% class(this_ddf)){
        next
      }
      # set the detection function parameters to be par
      this_ddf <- set_ddf_par(fleshy[[i]], this_ddf)

      # for io we need new data from both observers
      if(this_ddf$method == "io"){
        ds_newdata[[i]] <- rbind(ds_newdata[[i]], ds_newdata[[i]])
        ds_newdata[[i]]$observer <- c(rep(1, nrow(ds_newdata[[i]])/2),
                                      rep(2, nrow(ds_newdata[[i]])/2))
        # note that predict.io will return a vector of the correct
        # length
      }
      # calculate probability of detection
      this_p <- predict(this_ddf, newdata=ds_newdata[[i]],
                        compute=TRUE)$fitted

      # repopulate with the duplicates back in
      this_p <- this_p[attr(ds_newdata[[i]], "index"), drop=FALSE]

      # what is the width?
      if(this_ddf$method == "io"){
        this_width <- this_ddf$ds$meta.data$width
      }else{
        this_width <- this_ddf$meta.data$width
      }

      # calculate offset
      if(this_ddf$meta.data$point){
        # calculate log effective circle area
        # nb. predict() returns effective area of detection for points
        ret[data$ddfobj==i] <- linkfn(this_p * pi * this_width^2 *
                                      data$Effort[data$ddfobj==i])
      }else{
        # calculate log effective strip width
        ret[data$ddfobj==i] <- linkfn(2 * this_width * this_p *
                                      data$Effort[data$ddfobj==i])
      }
    }
    return(ret)
  }

  # now form the data structures we need to calculate the derivatives
  # using the function above...
  ds_newdata <- list()
  u_ds_newdata <- list()
  for(i in seq_along(ddf)){

    # get this detection function
    this_ddf <- ddf[[i]]
    # get all the covariates in this model
    df_vars <- all_df_vars(this_ddf)

    # if we don't have a real detection function
    if("fake_ddf" %in% class(this_ddf)){
      ds_newdata[[i]] <- NA
      u_ds_newdata[[i]] <- NA
      next
    }

    # if we don't have covariates
    # then just setup the data for predict.ds to be a distance of 0,
    # which will be ignored anyway
    if(length(df_vars)==0){
      ds_newdata[[i]] <- data.frame(distance=rep(0, sum(model$data$ddfobj==i)))
    }else{
      # otherwise need the covars that are in the data (that we saved
      # with keepData=TRUE :))
      ds_newdata[[i]] <- model$data[model$data$ddfobj==i, ]
      ds_newdata[[i]] <- ds_newdata[[i]][ , df_vars, drop=FALSE]
      ds_newdata[[i]]$distance <- 0
    }

    if(this_ddf$method == "io"){
      ds_newdata[[i]]$observer <- 1
    }


    # probably a lot of duplicated stuff in the above, so let's just
    # pass the unique combinations
    # inside mu_fn will return the right length
    u_ds_newdata[[i]] <- mgcv::uniquecombs(ds_newdata[[i]])

  }

  # find the derivatives of log(mu)
  pars <- unlist(parskel)
  firstD <- numderiv(mu_fn, pars, linkfn=linkfn, ddf=ddf,
                     data=model$data, ds_newdata=u_ds_newdata, skel=parskel)
  if(!is.matrix(firstD)){
    firstD <- matrix(firstD, ncol=length(pars))
  }

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

  # now form the paraPen bit
  hess <- matrix(0, length(pars), length(pars))
  ii <- 1

  # loop over our detection functions, extracting the hessian
  for(i in seq_along(ddf)){
    this_ddf <- ddf[[i]]

    if("fake_ddf" %in% class(this_ddf)){
      next
    }
#    opt_details <- attr(this_ddf$ds, "details")
#    if(is.matrix(opt_details)){
#      this_hess <- opt_details[nrow(opt_details), ]$nhatend
#    }else{
#      this_hess <- opt_details$nhatend
#    }
#    if(any(is.na(this_hess))){
#      # fall back to DS use if things are bad
#      this_hess <- this_ddf$hessian
#    }

    this_hess <- get_hessian(this_ddf)

    # drop that matrix into the big hessian
    hess[ii:(ii+nrow(this_hess)-1), ii:(ii+ncol(this_hess)-1)] <- this_hess

    ii <- ii+nrow(this_hess)
  }

  this_call$paraPen <- c(this_call$paraPen, list(XX=list(hess)))
  # tell gam.fixed.priors what to look for
  this_call$fixed.priors <- "XX"

  # set the trace on the scale parameter estimation
  this_call$scale.trace <- trace

  # is the scale estimated for this model?
  this_call$scale.estimated <- model$scale.estimated

  ## refit the model
  refit <- do.call("gam.fixed.priors", this_call)
  refit$data <- dat
  refit$ddf <- ddf
  class(refit) <- c("dsm", class(refit))

  ## now do some predictions

  # make some zeros for the paraPen term so they have no mean effect
  newdata[["XX"]] <- matrix(0, nrow(newdata), ncol(firstD))

  # get everyone's favourite matrix
  Lp <- predict(refit, newdata=newdata, type="lpmatrix")
  # predictions on the link scale
  pred <- Lp %*% coef(refit)
  pred <- newdata$off.set * linkinvfn(pred)

  # get variance-covariance matrix
  vc <- refit[[var_type]]

  # this is why we can only use log link
  dNdbeta <- t(pred)%*%Lp

  # make a sandwich
  var_p <- dNdbeta %*% vc %*% t(dNdbeta)

  # apply the link function to the offset
  # NB this is because refit is a gam not dsm object! If refit is dsm
  #    then this will get done in predict.dsm
  newdata$off.set <- linkfn(newdata$off.set)

  # if we are using Vc we need to set unconditional=TRUE
  if(var_type=="Vc"){
    uncond <- TRUE
  }else{
    uncond <- FALSE
  }
  # get the standard errors
  ses <- predict(refit, newdata=newdata, type="response", se.fit=TRUE,
                 unconditional=uncond)$se.fit

  # what should we return?
  ret <- list(old_model = model,
              refit     = refit,
              pred      = pred,
              var       = var_p,
              ses       = ses)

  class(ret) <- "dsm_varprop"

  return(ret)
}
