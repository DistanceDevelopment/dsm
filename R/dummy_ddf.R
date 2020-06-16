#' Make a dummy detection function for strip transects
#'
#' Strip transects aren't really detection functions, but this allows you to use them like detection functions in \code{dsm}.
#'
#' @export
#' @param obs which observations are to be included? Numeric vector corresponding to object IDs
#' @param size group size of observed group (default all groups size 1)
#' @param width right truncation
#' @param left left truncation
#' @param transect "line" or "point" transect (character)
#' @author David L Miller
dummy_ddf <- function(obs, size=1, width, left=0, transect="line"){

  if(!is.numeric(obs)){
    stop("obs should to be a numeric vector")
  }
  if(!is.vector(obs)){
    stop("obs should be a numeric vector")
  }

  object <- list()

  # put object IDs in a data.frame...
  object$data <- data.frame(object   = obs,
                            detected = rep(1, length(obs)),
                            observer = rep(1, length(obs)),
                            distance = rep(left, length(obs)),
                            size     = size)
  # set the fitted values
  object$fitted <- rep(1, length(obs))
  names(object$fitted) <- obs

  # truncation(s)
  object$meta.data <- list()
  object$meta.data$width <- width
  object$meta.data$left <- left

  object$meta.data$point <- FALSE
  if(transect == "point"){
    object$meta.data$point <- TRUE
  }

  class(object) <- c("fake_ddf", "ds", "ddf")
  return(object)
}

#' Prediction for fake detection functions
#'
#' Dummy function to return the correct number of 1s.
#'
#' @export
#' @param object model object
#' @param newdata how many 1s should we return?
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

summary.fake_ddf <- function(object, ...){
  class(object) <- "summary.fake.ddf"
  return(object)
}

print.summary.fake.ddf <- function(x, ...){
  cat("\nSummary for dummy ds object \n")
  cat("Number of observations : ", nrow(x$data),"\n")
  cat("Distance range         : ", x$meta.data$left, " - ",
                                   x$meta.data$width,"\n")
  cat("\nModel : No detection function, strip transect\n\n")
  cat("AIC   : NA\n")

  invisible()
}


