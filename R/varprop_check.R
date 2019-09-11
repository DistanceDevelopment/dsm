#' @importFrom mgcv uniquecombs
#' @importFrom mrds DeltaMethod
varprop_check <- function(object){

  # calculate the difference between the prob. detection before we refitted
  # and after we did the refit -- big values == BAD

  # extract the "random effects" coefficients
  coefs <- coef(object$refit)
  coefs <- coefs[grepl("^XX\\d*", names(coefs))]

  if(!all(class(object$old_model$ddf)=="list")){
    object$old_model$ddf <- list(object$old_model$ddf)
  }

  # storage
  df <- c()
  varprop_diagnostic <- c()
  for(i in seq_along(object$old_model$ddf)){

    # get the data in order
    oddf <- object$old_model$ddf[[i]]
    nd <- oddf$data

    if("fake_ddf" %in% class(oddf)){
      next
    }else if(oddf$ds$aux$ddfobj$scale$formula == "~1"){
      nd <- nd[1, "distance", drop=FALSE]
    }else{
      vars <- all.vars(as.formula(oddf$ds$aux$ddfobj$scale$formula))
      nd <- mgcv::uniquecombs(nd[, vars, drop=FALSE])

      # add distance column back in
      nd$distance <- 1
      # faff to work out which summaries we need...
      numeric_ind <- lapply(nd, is.factor)
      numeric_ind[["distance"]] <- NULL
      numeric_ind <- which(!unlist(numeric_ind))
      for(i in seq_along(numeric_ind)){
        this_ind <- numeric_ind[i]

        quants <- quantile(nd[, this_ind], c(0.05, .5, .95))

        colname <- colnames(nd)[this_ind]
        nd[[colname]] <- NULL
        nd <- mgcv::uniquecombs(nd)

        quants <- quants[rep(1:3, nrow(nd))]

        nd <- nd[rep(1:nrow(nd), 3), , drop=FALSE]
        nd[[colname]] <- quants
      }

      nd <- mgcv::uniquecombs(nd)
    }
    rownames(nd) <- NULL

    # p's for old model
    old_p <- predict(oddf, compute=TRUE, newdata=nd)$fitted

    # old model standard errors
    predict_f <- function(par, model, newdata){
      model$par <- par
      predict(model, compute=TRUE, newdata=nd)$fitted
    }
    old_p_se <- sqrt(diag(DeltaMethod(oddf$par, predict_f,
                                      solve(oddf$hessian),
                                      newdata=nd, model=oddf,
                                      delta=1e-8)$variance))

    # make a new detection function with the corrected parameters
    fix_ddf <- object$old_model$ddf[[i]]
    fix_ddf$par <- fix_ddf$par + coefs[1:length(fix_ddf$par)]
    coefs <- coefs[-(1:length(fix_ddf$par))]

    # p's for new model
    new_p <- predict(fix_ddf, compute=TRUE, newdata=nd)$fitted

    df <- c(df, rep(i, length(new_p)))
    varprop_diagnostic <- rbind(varprop_diagnostic,
                                c(old_p, old_p_se, new_p))

  }

  # format the output
  varprop_diagnostic <- as.data.frame(varprop_diagnostic)
  names(varprop_diagnostic) <- c("Fitted model",
                                 "Fitted model se",
                                 "Refitted model")
  if(i>1){
    varprop_diagnostic <- cbind("Detection function" = df,
                                varprop_diagnostic)
  }
  return(varprop_diagnostic)
}
