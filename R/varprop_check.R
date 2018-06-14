#' @importFrom mgcv uniquecombs
#' @importFrom mrds DeltaMethod
varprop_check <- function(object){

  # calculate the difference between the prob. detection before we refitted
  # and after we did the refit -- big values == BAD

  # get the data in order
  oddf <- object$old_model$ddf
  nd <- oddf$data
  if(oddf$ds$aux$ddfobj$scale$formula == "~1"){
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

      quants <- quantile(nd[,this_ind], c(0.05, .5, .95))
      #quants <- median(nd[, this_ind])

      colname <- colnames(nd)[this_ind]
      nd[[colname]] <- NULL
      nd <- mgcv::uniquecombs(nd)

      quants <- quants[rep(1:3, nrow(nd))]

      nd <- nd[rep(1:nrow(nd), 3), , drop=FALSE]
      nd[[colname]] <- quants
      #nd <- cbind(nd, quants)
      #names(nd)[ncol(nd)] <- colname

    }

    nd <- mgcv::uniquecombs(nd)
  }
  rownames(nd) <- NULL

  # p's for old model
  old_p <- predict(object$old_model$ddf, compute=TRUE, newdata=nd)

  # old model standard errors
  predict_f <- function(par, model, newdata){
    model$par <- par
    predict(model, compute=TRUE, newdata=nd)$fitted
  }
  old_p_se <- sqrt(diag(DeltaMethod(object$old_model$ddf$par, predict_f,
                                    solve(object$old_model$ddf$hessian),
                                    newdata=nd, model=object$old_model$ddf,
                                    delta=1e-8)$variance))

  # get new detection function
  fix_ddf <- object$old_model$ddf
  which.names <- !(names(coef(object$refit)) %in% names(coef(object$old_model)))
  fix_ddf$par <- fix_ddf$par + coef(object$refit)[which.names]

  # p's for new model
  new_p <- predict(fix_ddf, compute=TRUE, newdata=nd)

  # format data for return
  varprop_diagnostic <- data.frame("Fitted model"      = old_p$fitted,
                                   "Fitted model se"   = old_p_se,
                                   "Refitted model"    = new_p$fitted)
  nd$distance <- NULL
  varprop_diagnostic <- cbind.data.frame(nd, varprop_diagnostic)
  # sort by variable value
  varprop_diagnostic <- varprop_diagnostic[order(nd[, 1], decreasing=FALSE), ]

  return(varprop_diagnostic)
}
