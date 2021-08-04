#' Using RMark models in density surface models
#'
#' Create a detection function or availability object from a fitted RMark model.
#'
#' Models fitted by RMark can be used in DSM.... CHAT HERE
#'
#' RMark models can't contain variables named `detected`, `distance` or `observer`.
#' Time is the observer id
#'
#' Obs should only contain *one* observation per detected group!
#'
#' @export
#' @param rmark_model model fitted by RMark
#' @param obs observation data, `data.frame` with columns `object` and `size`
#'        plus all covariates used in `rmark_model`.
#' @param width right truncation
#' @param left left truncation (default 0, no left truncation)
#' @param transect `"line"` or `"point"` transect
#' @author David L Miller
rmark_ddf <- function(rmark_model, obs, width, left=0, transect="line"){

  if(!is.numeric(obs$object)){
    stop("object should to be a numeric vector")
  }

  if(!(transect %in% c("line", "point"))){
    stop("transect should be \"line\" or \"point\"")
  }

  if(!("Huggins" %in% class(rmark_model))){
    stop("the RMark model must be a Huggins model")
  }

  # check obs data
  if(!all(c("object", "size") %in% names(obs))){
    stop("obs must have columns object and size")
  }

  df_obj <- list()

  # put object IDs in a data.frame...
  df_obj$data <- obs
  df_obj$data$detected <- rep(1, nrow(obs))
  df_obj$data$observer <- rep(1, nrow(obs))
  df_obj$data$distance <- rep(left, nrow(obs))

  # get predictions based on the RMark model and covariates in obs
  # replicate=TRUE will copy the data for each capture occasion (MR observer)
  preddata <- obs
  preds <- predict_real(model=rmark_model,
                        df=rmark_model$design.data$p,
                        parameter="p",
                        replicate=TRUE,
                        data=preddata)

  # set the fitted values
  pdot.triple <- function(p){
    1-prod(1-p)
  }
  # get the ps combined over the observers (times)
  df_obj$fitted <- aggregate(preds$real, list(preds$object), pdot.triple)$x
  names(df_obj$fitted) <- df_obj$data$object

  # truncation(s)
  df_obj$meta.data <- list()
  df_obj$meta.data$width <- width
  df_obj$meta.data$left <- left

  df_obj$meta.data$point <- FALSE
  if(transect == "point"){
    df_obj$meta.data$point <- TRUE
  }

  # save the rmark model
  df_obj$rmark_model <- rmark_model

  # make the method be "rmark"
  df_obj$method <- "rmark"

  class(df_obj) <- c("rmark_ddf", "ds", "ddf")
  return(df_obj)
}

#' Prediction for fake detection functions
#'
#' Prediction function for dummy detection functions. The function returns as
#' many 1s as there are rows in \code{newdata}. If \code{esw=TRUE} then the
#' strip width is returned.
#'
#' @export
#' @param object model object
#' @param newdata covariates needed to get a probabilities of detection
#' @param compute unused, compatability with [`mrds::predict`][mrds::predict]
#' @param int.range unused, compatability with [`mrds::predict`][mrds::predict]
#' @param esw should the strip width be returned?
#' @param \dots for S3 consistency
#' @author David L Miller
predict.rmark_ddf <- function(object, newdata=NULL, compute=FALSE,
                              int.range=NULL, esw=FALSE, ...){


  if(esw) stop("esw can't be predicted when using RMark models")

  #ret$fitted <- rep(1, nrow(newdata))
  newdata$detected <- 1
  newdata$observer <- 1
  newdata$distance <- object$meta.data$left

  # get predictions based on the RMark model and covariates in obs
  # replicate=TRUE will copy the data for each capture occasion (MR observer)
  preds <- predict_real(model=object$rmark_model,
                        df=object$rmark_model$design.data$p,
                        parameter="p",
                        replicate=TRUE,
                        data=newdata)

  # set the fitted values
  pdot.triple <- function(p){
    1-prod(1-p)
  }
  # get the ps combined over the observers (times)
  ret <- list()
  ret$fitted <- aggregate(preds$real, list(preds$object), pdot.triple)$x
  names(ret$fitted) <- newdata$object

  return(ret)
}

#' @export
print.fake_ddf <- function(x, ...){
  print(summary(x))
}

#' @export
summary.rmark_ddf <- function(object, ...){
  class(object) <- "summary.rmark_ddf"
  return(object)
}

#' @export
print.summary.rmark_ddf <- function(x, ...){
  cat("\nSummary for RMark dummy object \n")
  cat("Number of observations : ", length(x$fitted),"\n")
  cat("Distance range         : ", x$meta.data$left, " - ",
                                   x$meta.data$width,"\n")
  cat("\nModel : Huggins model\n\n")
  cat("AICc   : ", x$rmark_model$results$AICc, "\n")

  invisible()
}


