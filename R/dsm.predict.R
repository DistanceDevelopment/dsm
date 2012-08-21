#' Produce estimated abundance at each cell in the prediction grid
#'
#' Function that extends the fitted density surface model from the covered 
#' region to the entire study region courtesy of predicting with covariates 
#' available in each of the cells of the prediction grid.
#'
#' @param model fitted model object (created by a call to \code{\link{gam}} or
#'                  \code{\link{dsm.fit}})
#' @param newdata spatially referenced covariates e.g.
#'                sea temperature, depth, distance to shore, etc.
#'                Note covariates in this dataframe must have names *identical* 
#'                to variable names used in fitting the GAM.
#' @param off.set offset for the cells in the prediction grid (\code{newdata}). 
#'                NB this should NOT have the link function applied already, 
#'                this happens below. (For example, it should _not_ be logged).
#' @return array of predicted values.

#' @note Presently, only density surface models fitted with \code{\link{mgcv}} 
#'       can be used for prediction with this incarnation of \code{dsm.predict}.
#' @author Eric Rexstad, David L. Miller 
#' @seealso predict.gam
#' @references Hedley, S. and S. T. Buckland. 2004. Spatial models for line transect sampling. JABES 9:181-199.
#'
#' Wood, S.N. 2006. Generalized Additive Models: An Introduction with R. CRC/Chapman & Hall.
#' @export
# Functions used:  predict.gam() from package mgcv
dsm.predict<-function(model, newdata=NULL, off.set=NULL){ 

  # do we have a dsm object or a gam?
  if("dsm" %in% class(model)){
    gam.model<-model$result
  }else{
    stop("Model must be the result of a call to dsm.fit!")
  }

  # is no new data are supplied, then we just predict at the data in 
  # the object.
  if(is.null(newdata)){
    prediction.grid<-model$data
  }else{
    # if we only got supplied 1 offset value, duplicate it enough times
    if(length(off.set)==1){
      off.set <- rep(off.set,nrow(newdata))
    }

    # apply the link function
    linkfn <- gam.model$family$linkfun
    off.set <- linkfn(off.set)

    # Append cell size of prediction grid to prediction grid  
    prediction.grid <- data.frame(newdata, off.set=off.set)
  }

  # actually do the predict call
  result<-predict(gam.model, newdata=prediction.grid, 
                  type="response", na.action=na.pass)

  return(result)
}
