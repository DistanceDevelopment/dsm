#' Observed versus expected diagnostics for fitted DSMs
#'
#' Given a covariate, calculate the observed and expected counts for each unique value of the covariate. This can be a useful goodness of fit check for DSMs.
#'
#' One strategy for model checking is to calculate observed and expected counts at different aggregations of the variable. If these match well then the model fit is good.
#'
#' @param model a fitted \code{dsm} model object
#' @param covar covariate to aggregate by (character)
#' @param cut vector of cut points to aggregate at. If not supplied, the unique values of \code{covar} are used.
#' @importFrom plyr ddply
#' @export
#' @author David L Miller, on the suggestion of Mark Bravington.
#' @export
#' @return \code{data.frame} with values of observed and expected counts.
#' @examples
#' library(Distance)
#' library(dsm)
#'
#' # example with the Gulf of Mexico dolphin data
#' data(mexdolphins)
#' hr.model <- ds(distdata, max(distdata$distance),
#'                key = "hr", adjustment = NULL)
#' mod1 <- dsm(count~s(x,y), hr.model, segdata, obsdata)
obs_exp <- function(model, covar, cut=NULL){

  # get data
  oe <- model$data
  # add in predictions
  oe$N <- predict(model)

  # do the aggregation if necessary
  if(!is.null(cut)){
    oe[[covar]] <- cut(oe[[covar]], breaks=cut)
  }

  # for each unique value of term, sum observed and expected
  oe <- plyr::ddply(oe, covar, function(x){
    data.frame(Observed = sum(x$count),
               Expected = sum(x$N))
  })

  # format the table
  cn <- oe[,1]
  oe <- t(oe[,2:3])
  colnames(oe) <- cn

  return(oe)
}
