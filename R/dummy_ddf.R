#' Detection function objects when detection is certain
#'
#' Create a detection function object for strip/plot surveys for use in density
#' surface models.
#'
#' @export
#' @param object numeric vector of object identifiers, relating to the `object`
#' field in the observation data of the DSM.
#' @param size group size for each observation (default all groups size 1)
#' @param width right truncation
#' @param left left truncation (default 0, no left truncation)
#' @param transect `"line"` or `"point"` transect
#' @author David L Miller
dummy_ddf <- function(object, size=1, width, left=0, transect="line"){

  if(!is.numeric(object)){
    stop("object should to be a numeric vector")
  }
  if(!is.vector(object)){
    stop("object should be a numeric vector")
  }

  if(!(transect %in% c("line", "point"))){
    stop("transect should be \"line\" or \"point\"")
  }

  df_obj <- list()

  # put object IDs in a data.frame...
  df_obj$data <- data.frame(object   = object,
                            detected = rep(1, length(object)),
                            observer = rep(1, length(object)),
                            distance = rep(left, length(object)),
                            size     = size)
  # set the fitted values
  df_obj$fitted <- rep(1, length(object))
  names(df_obj$fitted) <- object

  # truncation(s)
  df_obj$meta.data <- list()
  df_obj$meta.data$width <- width
  df_obj$meta.data$left <- left

  df_obj$meta.data$point <- FALSE
  if(transect == "point"){
    df_obj$meta.data$point <- TRUE
  }

  # make the method be "dummy"
  df_obj$method <- "dummy"

  class(df_obj) <- c("fake_ddf", "ds", "ddf")
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
#' @param newdata how many 1s should we return?
#' @param compute unused, compatibility with [`mrds::predict`][mrds::predict]
#' @param int.range unused, compatibility with [`mrds::predict`][mrds::predict]
#' @param esw should the strip width be returned?
#' @param \dots for S3 consistency
#' @author David L Miller
predict.fake_ddf <- function(object, newdata=NULL, compute=FALSE,
                             int.range=NULL, esw=FALSE, ...){

  ret <- list()

  if(is.null(newdata)){
    newdata <- data.frame(dummy=object$fitted)
  }

  if(esw){
    ret$fitted <- rep(object$meta.data$width-object$meta.data$left,
                      nrow(newdata))
  }else{
    ret$fitted <- rep(1, nrow(newdata))
  }

  return(ret)
}

#' @export
print.fake_ddf <- function(x, ...){
  print(summary(x))
}

#' @export
summary.fake_ddf <- function(object, ...){
  class(object) <- "summary.fake_ddf"
  return(object)
}

#' @export
print.summary.fake_ddf <- function(x, ...){
  cat("\nSummary for dummy ds object \n")
  cat("Number of observations : ", nrow(x$data),"\n")
  cat("Distance range         : ", x$meta.data$left, " - ",
                                   x$meta.data$width,"\n")
  cat("\nModel : No detection function, strip transect\n\n")
  cat("AIC   : NA\n")

  invisible()
}


