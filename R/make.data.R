make.data <- function(response, ddfobject, segdata, obsdata, group,
                      convert.units){

  # probably want to do something smart here...
  seglength.name<-'Effort'
  segnum.name<-'Sample.Label'
  distance.name<-'distance'
  cluster.name<-'size'
  sightnum.name<-'object'

  # Truncate observations made at distances greater than the truncation width;
  # truncation value picked up from metadata contained within ddfobject
  # No truncation for strip transects
  if (!is.null(ddfobject)){
     obsdata<-obsdata[obsdata[,distance.name]<=ddfobject$meta.data$width,]
  }

  # Estimating group abundance/density
  if(group){
    obsdata[,cluster.name][obsdata[,cluster.name]>0] <- 1
  }


  ###### DETECTION FUNCTION COVARIATES?
  ### what kind of data are we working with?
  #  - are there covariates in the detection function?
  mcds<-TRUE
  fitted.p <- fitted(ddfobject)

  # what were the object identifiers from ddf?
  ddf.names <- names(fitted.p)

  # if we have a model with no covariates then we just pick the 1 value 
  #  we need
  if(ddfobject$ds$aux$ddfobj$scale$formula == "~1"){
    mcds <- FALSE
    fitted.p <- as.vector(unique(fitted.p))
  }
  ###################################################################

  # remove the segments with 0 observations
  obsdata <- obsdata[obsdata$object %in% ddf.names,]

  # Aggregate response values of the sightings over segments
  # for the density models
  if(response %in% c("D","density")){
### don't think this is correct
    if(group){
      responsedata <- aggregate(obsdata[,cluster.name]/fitted.p,
                                list(obsdata[,segnum.name]), sum)
    }else{
      responsedata <- aggregate(obsdata[,cluster.name],
                                list(obsdata[,segnum.name]), sum) 
    }
    off.set <- "none"
  }else if(response %in% c("N","abundance") & group & !mcds){
  # group abundance without covariates
    responsedata <- aggregate(obsdata[,cluster.name],
                              list(obsdata[,segnum.name]), sum)
    off.set <- "eff.area"

  }else if(response %in% c("N","abundance") & !mcds){
  # individual abundance without covariates
    off.set <- "eff.area"
    if(group){
      responsedata <- aggregate(obsdata[,cluster.name]/fitted.p,
                                list(obsdata[,segnum.name]), sum)
    }else{
      responsedata <- aggregate(obsdata[,cluster.name],
                                list(obsdata[,segnum.name]), sum)
    }
  }else if(mcds){
  # group abundance with covariates
  # individual abundance with covariates
    responsedata <- aggregate(obsdata[,cluster.name]/fitted.p,
                              list(obsdata[,segnum.name]), sum)
    off.set<-"area"
  }

  # name the response data columns
  names(responsedata)<-c(segnum.name,response)

  # Next merge the response variable with the segment records and any 
  # response variable that is NA should be assigned 0 because these 
  # occur due to 0 sightings
  dat <- merge(segdata,responsedata,by=segnum.name,all.x=T)
  dat[,response][is.na(dat[,response])] <- 0

  # calculate the "width" of the transect first, make sure we get it right
  # if we are doing left truncation
  width <- ddfobject$meta.data$width
  if(!is.null(ddfobject$meta.data$left)){
    width <- width - ddfobject$meta.data$left
  }

  # calculate the density (count/area)
  if(response %in% c("D","density")){
    dat[,response]<-dat[,response]/(2*dat[,seglength.name]*width)
  }else{
    # calculate the offset
    #   area we just calculate the area
    #   effective area multiply by p
    #   when density is response, offset should be 1.

    dat$off.set<-switch(off.set,
                        eff.area=2*dat[,seglength.name]*width*fitted.p,
                        area=2*dat[,seglength.name]*width,
                        none=1)

    # multiply up by conversion factor
    dat$off.set<-dat$off.set*convert.units

    # Set offset as log of area or effective area
    dat$off.set<-log(dat$off.set)
  }

  return(dat)
}
