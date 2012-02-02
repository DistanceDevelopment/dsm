# fn.fit <- dsm.fit(cruise.segmented.detection.fn, response="indiv.est", formula=~s(seg.long)+s(seg.lat), model.defn=list(fn="gam",family="gaussian"), segdata=cruise.segmented.forddf$effort)
#===============================================================================================

dsm.fit <- function(ddfobject=ddfobject, phat=NULL, response=response, formula=formula,
                    model.defn=list(fn="gam",family="quasipoisson"), obsdata=obsdata,
                    segdata=segdata, wghts=NULL, link='log', off=NULL, convert.units=NULL)
#                    segdata=segdata, wghts=NULL, link='log', offset=offset)
#
# dsm.fit - fits a 'gam to segment-specific estimates of abundance
#           resulting from object-specific detection probabilities
#
#           This function has its orgins as perform.gam.fromddf (found in a txt file dsm.R from Oct '05)
#           This incarnation modifies that such that calling arguments conform to that requested
#           by the DISTANCE VB DSMNEngineInterface.MakeInputFile by deleting region.table and sample.table
#           and adding offset, link, and weights as arguments
#
#  Arguments:
#
#  ddfobject    - result from call to ddf; might be usurpt by phat argument below
#                   ddfobject is set to NULL when strip transects are analyzed
#  phat         - if present, represents estimated detection probabilities for each object present
#                 in the project database.  This breaks the obligation that detection functions for
#                 use in DSM need come from the MRDS engine.
#  response     - response to be modelled; choices are:
#                    group, indiv, group.est, indiv.est
#  formula      - formula for distribution surface
#  model.defn   - list comprised of
#                   function - gam or glm
#                   family - family of distribution for error term (link function follows)
#  obsdata      - object #, segment id, group size
#  segdata      - segment id and covariates data relevant to each segment (lat, long, depth, etc)
#  wghts        - direct pass-through argument to gam or glm
#  link         - direct pass-through argument subtly merged with family with eval(paste()) construct
#  offset       - direct pass-through argument **presently computed by this R-code**
#  convert.units - value to alter length to width for calculation of offset
#
#  Value:
#
#  dsm.fit      - object produced by mgcv
#
# Functions Used:  gam (from 'mgcv')
#
{

library(mgcv)
#
#  This stolen from Laake
default.field.names <- function()
{
return(c("y.name<-'Latitude'","x.name<-'Longitude'","seglength.name<-'Effort'",
         "segnum.name<-'Sample.Label'","distance.name<-'distance'","cluster.name<-'size'",
         "esw.name<-'esw'","sightnum.name<-'object'","sampling.fraction.name<-'SF'"))
}
field.names <- default.field.names()
eval(parse(text=field.names))

# variable.structure <- create.varstructure(model=ddfobject, region=region.table,
#                                          sample=sample.table, obs=obs.table)
# obsdata <- variable.structure$obs

#  Clumsy temporary fix because of field naming convention
#names(obsdata)[1] <- "Sig.no"
#names(obsdata)[5] <- "sample"
#
#   Truncate observations made at distances greater than the truncation width;
#       truncation value picked up from metadata contained within ddfobject
#       No truncation for strip transects
#
    if (!is.null(ddfobject)) obsdata<-obsdata[obsdata[,distance.name]<=ddfobject$meta.data$width,]

#  the following is borrowed (heavily) from dsm.count by Laake
#              ER modification is to test for presence of phat argument and substitute 
#                 detection probabilities from phat if provided
    if (response == "indiv.est" | response =="group.est" | response=="indiv.den" | response=="group.den")
    {
       if (!is.null(phat)) 
       {
          fitted.p <- phat
          object.data <- obsdata$sightnum.name
       }
       else
       {
          fitted.p <- ddfobject$fitted
          object.data <- names(ddfobject$fitted)
       }
       sig.prob <- data.frame(p = fitted.p, object = object.data)
#       old stmt: sig.prob <- data.frame(p = ddfobject$fitted, object = names(ddfobject$fitted))
#
#      Merge observations with sighting probabilities
#
       obsdata <- merge(obsdata, sig.prob, by = sightnum.name, all.x = T, sort = F)
#					Revision 10 Nov 2008; merge drops segments when detects are not made by primary see MLB and CODA
       obsdata <- obsdata[!is.na(obsdata$p), ]       
#
#      Check to see if any of the observations are missing a detection probability
#        for designs of type 'trial' objects observed by trackers will not have 
#        computed detection probabilities, so, trap that type of design, interrogating
#        the call to ddf (archived in ddfobject$call using an archane function 'languageEl'
#
       field.design <- substr(languageEl(ddfobject$call, which="method"),1,5)

       if(field.design != "trial" && any(is.na(obsdata$p)))
       {
          cat("The following sighting numbers don't have a matching detection probability\n")
          print(obsdata[,sightnum.name][is.na(obsdata$p)])
          stop("Function terminated")
       }
    }
#
#   If response is group or group.est - then change all cluster values to 1
#      if density is the responnse, then the response variable needs to be # detected divided by area!!!
    if (response == "group" | response =="group.est" )
    {
        obsdata[,cluster.name]<-rep(1,dim(obsdata)[1])
    }
#
#   Aggregate response values of the sightings over segments
#
    if (response == "indiv" | response =="group" )
    {
        responsedata <- aggregate(obsdata[,cluster.name],
                        list(obsdata[,segnum.name]), sum)
        off.set <- "eff.area"
    } else if (response == "indiv.est" | response == "group.est")
    {
        responsedata <- aggregate(obsdata[,cluster.name]/obsdata$p,
                        list(obsdata[,segnum.name]), sum)
        off.set <- "area"
    } else
    {
        responsedata <- aggregate(obsdata[,cluster.name]/obsdata$p,
                        list(obsdata[,segnum.name]), sum) 
        off.set <- "none"
    }
    
    names(responsedata)<-c(segnum.name,"N")
#
#   Next merge the response variable with the segment records and any response
#   variable that is NA should be assigned 0 because these occur due to 0 sightings
#
    dat<-merge(segdata,responsedata,by=segnum.name,all.x=T)
    dat$N[is.na(dat$N)]<-0
#          With density, we need to transform response variable to a density by dividing by area    
    if (off.set == "none") dat$N <- dat$N / 2 *dat[,seglength.name] * ddfobject$meta.data$width *convert.units
#          when density is response, offset should be 1.
    dat$off.set <- switch(off.set,
                   eff.area = 2 * dat[,seglength.name] * dat[,esw.name],
                   area = 2 * dat[,seglength.name] * ddfobject$meta.data$width,
                   none = 1)
#
#        Altered 2 Feb 06 to use final argument passed into function from InputFileMaker
    if(!is.null(convert.units) & off.set != "none")
       dat$off.set<-dat$off.set*convert.units
#
#
#   Set offset as log of area or effective area
#
    if (off.set != "none") dat$off.set<-log(dat$off.set)
#
#   Create formula 
#
    if (response=="indiv.den" | response=="group.den")
    {
           formula <- as.formula(paste("N", deparse(formula,width.cutoff=500),
                   collapse=""))
    } else
    {
            formula <- as.formula(paste("N", deparse(formula,width.cutoff=500),
                   "+ offset(off.set)",collapse=""))
    }
#
#   Paste link function argument together with family argument to present to gam/glm in the
#     form:  family=quasipoisson(link="log")
#
    family.and.link <- eval(parse(text=paste(model.defn$family, "(link='", link, "')", sep="")))
    if (!is.null(wghts)) wghts <- paste("dat$", wghts, sep="")
    dat$area = 2 * dat[,seglength.name] * ddfobject$meta.data$width * convert.units
#
#   Fit model  hardwiring gamma=1.4 per Wood (2006:254) who cites Kim and Gu(2004) for overfitting
#           weights should be 'area' when density is response.
    if (toupper(model.defn$fn) == "GAM") 
       {
        if (is.null(wghts)) {
          dsm.fit <- list(result=gam(formula, family=family.and.link, data = dat,
                                    control = gam.control(keepData=TRUE),
                                     weights = NULL, gamma=1.4), call.dsm=match.call())
        } else
        {
          dsm.fit <- list(result=gam(formula, family=family.and.link, data = dat,
                                     control = gam.control(keepData=TRUE),
                                     weights = eval(parse(text=wghts)), gamma=1.4), call.dsm=match.call())
        }
    } else
       {
        if (is.null(wghts)) {
          dsm.fit <- list(result=glm(formula, family=family.and.link, data = dat, control = glm.control(),
                                     weights = NULL, gamma=1.4), call.dsm=match.call())
        } else
        {
          dsm.fit <- list(result=glm(formula, family=family.and.link, data = dat, control = glm.control(),
                                     weights = eval(parse(text=wghts)), gamma=1.4), call.dsm=match.call())
        }
    }
#
#  Return model object
#
    dsm.fit
}

