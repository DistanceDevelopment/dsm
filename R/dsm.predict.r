dsm.predict <- function(gam.model=gam.model, newdata=NULL, field=FALSE, off=NULL, silent=FALSE)
#
#    dsm.predict - purpose to predict abundance within study area at locations
#                  away from the transects surveyed.
#
#  Arguments:
#     gam.model - fitted model object (created by previous call to gam())
#     newdata   - know as 'prediction.covariates.dat.r' in D6 DSM statement of 
#                 direction prediction grid data frame, spatially referenced 
#                 (latitude and longitude) covariates e.g., sea temperature, depth, 
#                 distance to shore, etc.
#     field     - logical indicating whether offset is provides as constant (False) 
#                 or field in grid dataframe (True)
#     off.set   - area of cells in prediction grid (newdata), provided by Distance
#                 Note covariates in this dataframe must have names *identical* 
#                 to variable names used in fitting the gam
#                 If entered between double quotes, will be a field named in 
#                 prediction grid, otherwise a numeric constant
#			         Assume if a constant comes in, it will already have had the link 
#                 transform applied, whereas if field, apply link fn.
#     silent    - Set to true when used in conjunction with variance bootstraps; in 
#                 this case we don't wish to see the abundance estimate printed
#                 This identical nature of the names in the model object and dataframe
#                 is not currently enforced; however, some code along the line of
#                 for (j in 1:length(names(newd))){
#                    print(any(i<-grep(names(newd)[[j]], 
#                      dimnames(attr(b$terms,"factors"))[[1]])))
#                 }              
#                 would perform some checking of this

#               Failure to have all covariates in 'newdata' results in message 
#               of this sort:
#                 Error in eval(expr, envir, enclos) : object "x2" not found
#                   In addition: Warning message:
#                   not all required variables have been supplied in newdata!
#                   in: predict.gam(b, newd, type = "response") 
#
#   Value:
#         produces one-dimensional array of predicted values (without se's)
#         as well as the call that generated this instance of the function

#   Functions used:  predict.gam() from package mgcv
#                    Note the default arguments in predict.gam()
#                    newdata.guaranteed=FALSE and
#                    na.action=na.pass (creation of NA as predicted value 
#                                       where a covariate is missing)
#                    seem believable; leave them at their defaults, *but*
#                    type="response" is not the default argument, and needs to be 
#                    specified
{
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
