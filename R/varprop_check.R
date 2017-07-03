#' @importFrom mgcv uniquecombs
#' @importFrom mrds DeltaMethod
varprop_check <- function(object){

  # calculate the difference between the prob. detection before we refitted
  # and after we did the refit -- big values == BAD

  # get the data in order
  oddf <- object$old_model$ddf
  nd <- oddf$data
  if(oddf$ds$aux$ddfobj$scale$formula=="~1"){
    nd <- nd[1, "distance", drop=FALSE]
  }else{
    nd <- mgcv::uniquecombs(nd[, all.vars(as.formula(oddf$ds$aux$ddfobj$scale$formula)), drop=FALSE])
    nd$distance <- 1
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

  return(varprop_diagnostic)
}
