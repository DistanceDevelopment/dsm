#' Handler to suppress the "matrix not positive definite" warning
#'
#' Internal function to suppress an annoying warnings from chol()
#'
#' See: https://stat.ethz.ch/pipermail/r-help/2012-February/302407.html
#' See: http://romainfrancois.blog.free.fr/index.php?post/2009/05/20/Disable-specific-warnings
#'
#' @param w a warning
#' @return not a warning if the warning was "matrix not positive definite"
#'
#' @author David L. Miller
matrixnotposdef.handler <- function(w){
  if(any(grepl("matrix not positive definite",w))){
    invokeRestart("muffleWarning")
  }
}
