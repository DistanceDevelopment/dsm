#' @importFrom mgcv uniquecombs
#' @importFrom mrds DeltaMethod
varprop_check <- function(object){

  # calculate the difference between the prob. detection before we refitted
  # and after we did the refit -- big values == BAD

  # make a skeleton for the detection function parameters over all models
  parskel <- list()
  for(i in seq_along(object$old_model$ddf)){
    parskel[[i]] <- object$old_model$ddf[[i]]$par
  }

  # which parameters in the GAM are only in the refit?
  which.names <- !(names(coef(object$refit)) %in%
                   names(coef(object$old_model)))
  extra_gam_pars <- coef(object$refit)[which.names]
  # use these extra pars as the flesh for the skel
  ddf_corrections <- relist(extra_gam_pars, parskel)

  vp_diag <- list()

  for(ii in seq_along(object$old_model$ddf)){
    # get the data in order
    oddf <- object$old_model$ddf[[ii]]
    nd <- oddf$data

    # get all variables in this model
    df_vars <- all_df_vars(oddf)

    if(any(class(oddf)=="fake_ddf")){
      next
    }else if(length(df_vars) == 0){
      nd <- nd[1, "distance", drop=FALSE]
    }else{
      nd <- mgcv::uniquecombs(nd[, df_vars, drop=FALSE])

      # add distance column back in
      nd$distance <- 0
      nd$observer <- 1
      # faff to work out which summaries we need...
      numeric_ind <- lapply(nd, is.factor)
      numeric_ind[["distance"]] <- NULL
      numeric_ind[["observer"]] <- NULL
      numeric_ind <- which(!unlist(numeric_ind))
      for(i in seq_along(numeric_ind)){
        this_ind <- numeric_ind[i]

        quants <- quantile(nd[,this_ind], c(0.05, .5, .95))

        colname <- colnames(nd)[this_ind]
        nd[[colname]] <- NULL
        nd <- mgcv::uniquecombs(nd)

        quants <- quants[rep(1:3, nrow(nd))]

        nd <- nd[rep(1:nrow(nd), 3), , drop=FALSE]
        nd[[colname]] <- quants
      }

      nd <- mgcv::uniquecombs(nd)

      # for io models we need to have both observers there for the
      # prediction to work
      if(oddf$method == "io"){
        nd <- rbind(nd, nd)
        nd$observer <- c(rep(1, nrow(nd)/2), rep(2, nrow(nd)/2))
      }
    }
    rownames(nd) <- NULL

    # p's for old model
    old_p <- predict(oddf, compute=TRUE, newdata=nd)

    # old model standard errors
    predict_f <- function(par, model, newdata){
      model <- set_ddf_par(par, model)
      predict(model, compute=TRUE, newdata=nd)$fitted
    }
    old_p_se <- sqrt(diag(DeltaMethod(oddf$par, predict_f,
                                      solve(oddf$hessian),
                                      newdata=nd, model=oddf,
                                      delta=1e-8)$variance))

    # now need to construct the "new" refitted detection function
    # make a copy of the detection function
    fix_ddf <- oddf
    # correct the parameters for this model
    fix_ddf <- set_ddf_par(fix_ddf$par + ddf_corrections[[ii]], fix_ddf)

    # calculate detectabilities for new model
    new_p <- predict(fix_ddf, compute=TRUE, newdata=nd)

    # format data for return
    varprop_diagnostic <- data.frame("Fitted model"      = old_p$fitted,
                                     "Fitted model se"   = old_p_se,
                                     "Refitted model"    = new_p$fitted)
    nd$distance <- NULL
    varprop_diagnostic <- cbind.data.frame(nd, varprop_diagnostic)

    # sort by variable value
    if(nrow(varprop_diagnostic) > 1){
      varprop_diagnostic <- varprop_diagnostic[order(nd[, 1],
                                                     decreasing=FALSE), ]
    }
    # get rid of the observer column and shorten the table
    varprop_diagnostic$observer <- NULL
    varprop_diagnostic <- unique(varprop_diagnostic)

    # get a somewhat informative model description
    attr(varprop_diagnostic, "model_description") <- ddf.model.description(oddf)
    # store
    vp_diag[[ii]] <- varprop_diagnostic
  }

  return(vp_diag)
}
