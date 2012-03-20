#' Produce estimated abundance at each cell in the prediction grid
#'
#' Function that extends the fitted density surface model from the covered region
#' to the entire study region courtesy of predicting with covariates available in
#' each of the cells of the prediction grid.
#' @param gam.model fitted model object (created by previous call to \code{gam()}/
#'                  \code{dsm.fit()})
#' @param newdata spatially referenced (latitude and longitude) covariates e.g.,   
#'                sea temperature, depth, distance to shore, etc.
#'                Known as 'prediction.covariates.dat.r' in D6 DSM statement of 
#'                direction prediction grid data frame.
#'                
#' @param field logical indicating whether offset is provides as constant 
#'              (\code{FALSE}) or field in grid dataframe (\code{TRUE}).
#' @param off.set area of cells in prediction grid (\code{newdata}). 
#'                 Provided by Distance
#'                 Note covariates in this dataframe must have names *identical* 
#'                 to variable names used in fitting the \code{gam()}
#'                 If entered between double quotes, will be a field named in 
#'                 prediction grid, otherwise a numeric constant
#'			          Assume if a constant comes in, it will already have had the link 
#'                 transform applied, whereas if field, apply link fn.
#' @param silent  Set to true when used in conjunction with variance bootstraps; in 
#'                this case we don't wish to see the abundance estimate printed
#'                This identical nature of the names in the model object and dataframe
#'                is not currently enforced; however, some code along the line of
#'                \code{for (j in 1:length(names(newd))){
#'                   print(any(i<-grep(names(newd)[[j]], 
#'                     dimnames(attr(b$terms,"factors"))[[1]])))
#'                }}         
#'                would perform some checking of this
#'
#'                Failure to have all covariates in \code{newdata} results in message 
#'                of this sort:
#'                 Error in eval(expr, envir, enclos) : object "x2" not found
#'                   In addition: Warning message:
#'                   not all required variables have been supplied in newdata!
#'                   in: predict.gam(b, newd, type = "response") 
#'
#' @return predicted list consisting of:
#'          \tabular{ll}{\code{result} \tab one-dimensional array of predicted 
#'                       values (without se's)\cr
#'                       \code{call} \tab the call that generated this instance 
#'                       of the function}

#   Functions used:  predict.gam() from package mgcv
#' @note Presently, only density surface models fitted with \code{mgcv()} can be 
#' used for prediction with this incarnation of \code{dsm.predict()}.
#' Note the default arguments in \code{predict.gam()} \code{newdata.guaranteed=FALSE}
#' and \code{na.action=na.pass} (creation of NA as predicted value where a covariate 
#' is missing) seem believable; leave them at their defaults, *but* 
#' \code{type="response"} is not the default argument, and needs to be specified.
#' @author Eric Rexstad, David L. Miller 
# @seealso 
#' @references Hedley, S. and S. T. Buckland. 2004. Spatial models for line transect sampling. JABES 9:181-199.
#'
#' Wood, S.N. 2006. Generalized Additive Models: An Introduction with R. CRC/Chapman & Hall.
#' @export

dsm.predict<-function(gam.model, newdata=NULL, field=FALSE, off=NULL, 
                      silent=FALSE){
   #  Append cell size of prediction grid to prediction grid  if off.set 
   #  argument is a number, otherwise manufacture
   if (!field){ 
      prediction.grid <- data.frame(newdata, off.set=off)
   }else{
      # match.call()[] provides literal strings of the arguments passed to this 
      # function; I wish to build a string consisting of the string associated with 
      # the 3rd argument concatenated with the string of the 5th argument; in 
      # combination they give the name of the field in the dataframe containing 
      # cell size
      # if(gam.model$result$family$link=="log"){ # suggested by Louise 22 Sept 2008
      if(gam.model$family$link=="log"){
	      newdata$off.set<-eval(parse(text=paste("log(", match.call()[3], 
                               "$", match.call()[5], ")", sep=""))) 
      }else{
         newdata$off.set<-eval(parse(text=paste(match.call()[3], 
                               "$", match.call()[5], sep="")))
      }
      prediction.grid <- newdata 
   }
   # check whether surface was fitted using gam or glm, and call the 
   # appropriate 'predict'    
   if(gam.model$call[1]=="gam()"){
      result<-predict.gam(object=gam.model, newdata=prediction.grid, 
                          type="response", na.action=na.pass)
   }else{
      result<-predict.glm(object=gam.model, newdata=prediction.grid, 
                          type="response", na.action=na.pass)
   }

   if(!silent){
      cat("Total abundance in study area = ", sum(result, na.rm=TRUE), "\n")
   }

   result<-list(result=result, call=match.call())
}
