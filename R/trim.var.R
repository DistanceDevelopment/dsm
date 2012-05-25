#' Trimmed variance
#'
#' Trim the variance estimates from the bootstrap. This is defined as the 
#' percentage defined as amount necessary to bring median and trimmed mean 
#' within 8% of each other these are defined as 'outliers'. 
#'
#' Code originally by Louise Burt.
#' 
#' @param untrimmed.bootstraps (usually the \code{$study.area.total} element
#'        of a returned \code{dsm} bootstrap object.
#'
#' @return trimmed variance
#'
#' @export
#'
#' @author Louise Burt
trim.var<-function(untrimmed.bootstraps){
  outliers <- boxplot.stats(untrimmed.bootstraps, coef=1.5)$out
  bootstrap.abund <-untrimmed.bootstraps[!(untrimmed.bootstraps %in% outliers)]

  ret <- var(bootstrap.abund)

  attr(ret,"trim.prop") <- length(outliers) / length(untrimmed.bootstraps)
  attr(ret,"untrimn") <- length(untrimmed.bootstraps)
  attr(ret,"outliers") <- length(outliers) 

  return(ret)
}
