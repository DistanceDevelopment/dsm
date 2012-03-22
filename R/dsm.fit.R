#' Fits a GAM to segment-specific estimates of abundance
#' resulting from object-specific detection probabilities
#' 
#' Function constructs and then invokes a call to \code{gam()} and
#' returns the result of the fitting of the density surface model.
#' 
#' @param ddfobject result from call to \code{\link{ddf}}; might be usurpt by
#'   the \code{phat} argument below.
#'   If \code{ddfobject} is set to \code{NULL} when strip transects 
#'   are analyzed
#' @param phat if present, represents estimated detection probabilities 
#'   for each object present in the project database. This breaks the 
#'   obligation that detection functions for use in \code{dsm} need come
#'   from \code{mrds}. 
#' @param response response to be modelled; choices are:
#'   \tabular{ll}{\code{indiv} \tab - individual counts per segment\cr
#'                \code{group} \tab - group counts per segment\cr
#'                \code{indiv.est} \tab - estimated individual abundance per 
#'                  segment\cr
#'                \code{group.est} \tab - estimated group abundance per 
#'                  segment \cr
#'                \code{indiv.den} \tab individual density per segment\cr
#'                \code{group.den} \tab group density per segment\cr}
#' @param formula formula for the surface. This should be a
#'   valid \code{glm} or \code{gam} formula. In the GAM case, the \code{s}
#'   term should include basis definition (\code{bs} and \code{k} terms). If
#'   the soap film smoother is to be used (\code{bs="so"}), it must include
#'   \code{xt=list(bnd=bnd)} for the boundary, in this case the \code{k}
#'   parameter gives the complexity of the boundary smoother. See 
#'   \code{\link{smooth.construct.so.smooth.spec}} for further details.
#' @param model.defn a list comprised of
#'   \tabular{ll}{\code{function} \tab "gam" or "glm" \cr
#'                \code{family} \tab a \code{\link{family}} name (in 
#'                   quotes) for the response (link function definition 
#'                   follows) \cr
#'                \code{family.pars} \tab a list of named parameters 
#'                   for the family specified, can be \code{NULL}\cr
#'                \code{bnd} \tab Only required for soap film smoothers. 
#'                   A list with two elements, x and y which are the
#'                   vertices of the boundary polygon see
#'                   \code{\link{smooth.construct.so.smooth.spec}}\cr
#'                \code{knots} \tab Only required for soap film smoothers.
#'                   Internal knots for the soap film smoother.}
#' @param obsdata observation data, see \code{\link{dsm-data}}. 
#' @param segdata segment data, see \code{\link{dsm-data}}.
#' @param wghts weights, directly passed to \code{gam} or \code{glm}.
#' @param link link function, merged with \code{family} via 
#'   \code{eval(paste())}.
#' @param convert.units value to alter length to width for calculation 
#'   of the offset.
#' @return a list, consisting of:
#'   \tabular{ll}{\code{result} \tab object produced by the \code{gam} or
#'                  \code{glm} call.\cr
#'                \code{call.dsm} \tab the call to this function.
#'                \code{data} \tab the data object supplied in the call\cr
#'                \code{ddf} \tab the \code{\link{ddf}} object supplied}
#' @note Note that the gamma parameter to \code{gam()} is hardwired here; 
#'       set to a value of 1.4 (from advice in Wood (2006)) such that the 
#'       \code{gam()} is inclined to not 'overfit.'
#' @author Eric Rexstad, David L. Miller 
# @seealso 
#' @references Hedley, S. and S. T. Buckland. 2004. Spatial models for line transect sampling. JABES 9:181-199.
#'
#' Wood, S.N. 2006. Generalized Additive Models: An Introduction with R. CRC/Chapman & Hall.
#' @export 
# @keywords
# @examples
dsm.fit <- function(ddfobject, phat=NULL, response, formula,
                    model.defn=list(fn="gam",family="quasipoisson"), obsdata,
                    segdata, wghts=NULL, link='log',convert.units=1,...)
#
# History:
# This function has its orgins as perform.gam.fromddf (found in a txt 
# file dsm.R from Oct '05)
#
# This incarnation modifies that such that calling arguments conform 
# to that requested by the DISTANCE VB DSMNEngineInterface.MakeInputFile 
# by deleting region.table and sample.table and adding offset, link, and 
# weights as arguments
#
# Jan 2012, Dave Miller started updating and turning into a proper
#           R library.
#
{

#   #  This stolen from Laake
#   default.field.names <- function()
#   {
#   return(c("y.name<-'Latitude'","x.name<-'Longitude'","seglength.name<-'Effort'",
#            "segnum.name<-'Sample.Label'","distance.name<-'distance'","cluster.name<-'size'",
#            "esw.name<-'esw'","sightnum.name<-'object'","sampling.fraction.name<-'SF'"))
#   }
#   field.names <- default.field.names()
#   eval(parse(text=field.names))
   # above code just does this, but it is a bit obtuse
   # probably want to do something smart here...
  y.name<-'y'
  x.name<-'x'
  seglength.name<-'Effort'
  segnum.name<-'Sample.Label'
  distance.name<-'distance'
  cluster.name<-'size'
  esw.name<-'esw'
  sightnum.name<-'object'
  sampling.fraction.name<-'SF'

  # Truncate observations made at distances greater than the truncation width;
  # truncation value picked up from metadata contained within ddfobject
  # No truncation for strip transects
  if (!is.null(ddfobject)){
     obsdata<-obsdata[obsdata[,distance.name]<=ddfobject$meta.data$width,]
  }

  #  the following is borrowed (heavily) from dsm.count by Laake
  #  ER modification is to test for presence of phat argument and substitute 
  #     detection probabilities from phat if provided
  if(response=="indiv.est" | response=="group.est" | 
       response=="indiv.den" | response=="group.den"){
    if(!is.null(phat)){
      fitted.p<-phat
      object.data<-obsdata$sightnum.name
    }else{
      fitted.p<-ddfobject$fitted
      object.data<-names(ddfobject$fitted)
    }
    sig.prob<-data.frame(p=fitted.p, object=object.data)

    # Merge observations with sighting probabilities
    obsdata <- merge(obsdata, sig.prob, by=sightnum.name, all.x=T, sort=F)
    # Revision 10 Nov 2008; merge drops segments when detects are 
    # not made by primary see MLB and CODA
    obsdata <- obsdata[!is.na(obsdata$p), ]       

    # Check to see if any of the observations are missing a detection 
    # probability for designs of type 'trial' objects observed by 
    # trackers will not have computed detection probabilities, so, trap 
    # that type of design, interrogating the call to ddf (archived in 
    # ddfobject$call using an archane function 'languageEl'
    field.design<-substr(languageEl(ddfobject$call, which="method"),1,5)

    if(field.design!="trial" && any(is.na(obsdata$p))){
      cat("The following sighting numbers don't have a matching detection probability\n")
      print(obsdata[,sightnum.name][is.na(obsdata$p)])
      stop("Function terminated")
    }
  }

  # If response is group or group.est - then change all cluster values to 1
  # if density is the responnse, then the response variable needs to 
  # be # detected divided by area!!!
  if(response=="group" | response=="group.est"){
    obsdata[,cluster.name]<-rep(1,dim(obsdata)[1])
  }

  # Aggregate response values of the sightings over segments
  if(response=="indiv" | response=="group"){
    responsedata<-aggregate(obsdata[,cluster.name],
                            list(obsdata[,segnum.name]), sum)
    off.set <- "eff.area"

    # need to find the fitted ps to divide through by below
    # we don't have individual level covariates, so there is only
    # one p
  
    # did the user supply phat?
    if(!is.null(phat)){
      fitted.p<-unique(phat)
    }else if(any(names(obsdata)==esw.name)){
      fitted.p<-unique(obsdata[,esw.data]*ddfobject$width)
    }else{
      fitted.p<-unique(ddfobject$fitted)
    }    

    if(length(fitted.p)!=1){
      stop(paste("More than one unique probability of detection.\n",
                 "Are you sure this is an \"indiv\"/\"group\" analysis?"))
    }
  }else if (response=="indiv.est" | response=="group.est"){
    responsedata<-aggregate(obsdata[,cluster.name]/obsdata$p,
                            list(obsdata[,segnum.name]), sum)
    off.set<-"area"
  }else{
    responsedata <- aggregate(obsdata[,cluster.name]/obsdata$p,
                              list(obsdata[,segnum.name]), sum) 
    off.set <- "none"
  }
  
  # we'll just call the response N, whatever we're actually modelling
  names(responsedata)<-c(segnum.name,"N")

  # Next merge the response variable with the segment records and any 
  # response variable that is NA should be assigned 0 because these 
  # occur due to 0 sightings
  dat<-merge(segdata,responsedata,by=segnum.name,all.x=T)
  dat$N[is.na(dat$N)]<-0
  # With density, we need to transform response variable to a density 
  # by dividing by area    
  if (off.set=="none"){
    dat$N<-dat$N/2*dat[,seglength.name]*ddfobject$meta.data$width*convert.units
  }

  # when density is response, offset should be 1.
  dat$off.set<-switch(off.set,
                      eff.area=2*dat[,seglength.name]*fitted.p,
                      area=2*dat[,seglength.name]*ddfobject$meta.data$width,
                      none=1)
  # Altered 2 Feb 06 to use final argument passed into function 
  # from InputFileMaker
  if(!is.null(convert.units) & off.set!="none"){
    dat$off.set<-dat$off.set*convert.units
  }

  # Set offset as log of area or effective area
  if(off.set!="none"){
    dat$off.set<-log(dat$off.set)
  }

  # Create formula 
  if(response=="indiv.den" | response=="group.den"){
    formula<-as.formula(paste("N", deparse(formula,width.cutoff=500),
                              collapse=""))
  }else{
    formula<-as.formula(paste("N", deparse(formula,width.cutoff=500),
                                "+ offset(off.set)",collapse=""))
  }

  ###########################################
  ### Response distribution, link function etc
  # Paste link function argument together with family argument to 
  # present to gam/glm in the form:  family=quasipoisson(link="log")
  if(model.defn$family=="Tweedie"){
    # need to specify the Tweedie parameters
    if(is.null(model.defn$family.pars$p)){
       error("You must specify the p parameter to use the Tweedie family! See ?Tweedie.")
    }
    family.and.link<-eval(parse(text=paste(model.defn$family,
                         "(link='", link, "',p=",model.defn$family.pars$p,")",
                                           sep="")))
  }else if(model.defn$family=="quasi"){
    # specify the variance relationship for quasi
    if(is.null(model.defn$family.pars$variance)){
       error("You must specify the variance relation to use the quasi family! See ?family.")
    }
    family.and.link<-eval(parse(text=paste(model.defn$family,
                           "(link='", link, 
                        "',variance='",model.defn$family.pars$variance,"')",
                                           sep="")))
  }else{
    # if the family does not have extra parameters
    family.and.link<-eval(parse(text=paste(model.defn$family, 
                          "(link='", link, "')", sep="")))
  }

  if (!is.null(wghts)){
    wghts<-paste("dat$", wghts, sep="")
  }
  dat$area<-2*dat[,seglength.name]*ddfobject$meta.data$width*convert.units

  if(toupper(model.defn$fn)=="GAM"){
    # if we are doing soap film smoothing, we need make sure that we
    # don't mess up the knots
    if(grepl('bs = "so"',as.character(formula)[3])){

      # the boundary must be name bnd in the formula, 
      # otherwise this breaks
      bnd<-model.defn$bnd
      knots<-model.defn$knots
  
      if(is.null(wghts)){
        b<-try(gam(formula,family=family.and.link,data=dat,
               control=gam.control(keepData=TRUE),
               weights=NULL, gamma=1.4,...))
      }else{
        b<-try(gam(formula, family=family.and.link, data=dat,
               control=gam.control(keepData=TRUE),
               weights=eval(parse(text=wghts)),gamma=1.4,...))
      }

      # loop until the knots are okay but not more than 5 times
      max.wiggles<-1
      while((b[[1]]==
        "Error in check.knots(g) : Please (re)move problematic knots.\n") & 
        (max.wiggles<6)){

        # find the problem knots
        warning.mess<-names(last.warning)
        problem.knots<-as.numeric(
                  gsub("knot ([0-9]+) is in boundary cell of solution grid",
                               "\\1",warning.mess))

        # wiggle them
        for(i in problem.knots){
           this.knot<-knots[i,]
           this.knot<-this.knot+runif(2,-runif.scale,runif.scale)
           while(!inSide(bnd,this.knot[1],this.knot[2])){
              this.knot<-knots[i,]
              this.knot<-this.knot+runif(2,-runif.scale,runif.scale)
           }
           knots[i,]<-this.knot
        }

        # refit the model
        if(is.null(wghts)){
           b<-try(gam(formula,family=family.and.link,data=dat,
                  control=gam.control(keepData=TRUE),
                  weights=NULL, gamma=1.4,...))
        }else{
           b<-try(gam(formula, family=family.and.link, data=dat,
                  control=gam.control(keepData=TRUE),
                  weights=eval(parse(text=wghts)),gamma=1.4,...))
        }

         max.wiggles<-max.wiggles+1
      }

      # if we had too many wiggles above
      if(max.wiggles==5){
        stop("Too many soap knot failures! Try different knots!")
      }

    }else{
       # Non-soap model
       # Fit model  hardwiring gamma=1.4 per Wood (2006:254) who cites 
       # Kim and Gu(2004) for overfitting weights should be 'area' when 
       # density is response.
       if(is.null(wghts)){
         b<-gam(formula,family=family.and.link,data=dat,
                control=gam.control(keepData=TRUE),weights=NULL, gamma=1.4)
       }else{
         b<-gam(formula, family=family.and.link, data=dat,
                control=gam.control(keepData=TRUE),
                weights=eval(parse(text=wghts)),gamma=1.4)
       }
     }


  }else{
    # GLM case
    if(is.null(wghts)){
       b<-glm(formula,family=family.and.link,data=dat,control=glm.control(),
                                  weights=NULL,gamma=1.4)
    }else{
       b<-glm(formula,family=family.and.link,data=dat,control=glm.control(),
                                weights=eval(parse(text=wghts)),gamma=1.4)
    }
  }
  # Return model object
  ret<-list(result=b,call.dsm=match.call(),data=dat,ddf=ddfobject)
  class(ret)<-"dsm"

  return(ret)
}
