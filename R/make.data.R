make.data <- function(response, ddfobject, segdata, obsdata, group,
                      convert.units, availability, strip.width, segment.area){

  # probably want to do something smart here...
  seglength.name<-'Effort'
  segnum.name<-'Sample.Label'
  distance.name<-'distance'
  cluster.name<-'size'

  # Truncate observations made at distances greater than the truncation width;
  # truncation value picked up from metadata contained within ddfobject
  # No truncation for strip transects
  if (!is.null(ddfobject)){
     obsdata <- obsdata[obsdata[,distance.name]<=ddfobject$meta.data$width,]
  }

  # Estimating group abundance/density
  if(group){
    obsdata[,cluster.name][obsdata[,cluster.name]>0] <- 1
  }

  # if we fitted a detection function
  if (!is.null(ddfobject)){
    # grab the probabilities of detection
    fitted.p <- fitted(ddfobject)

    # remove observations which were not in the detection function
    obsdata <- obsdata[obsdata$object %in% names(fitted.p),]
  }else{
    # strip transects or presence/absence data
    fitted.p <- 1
    if(is.null(strip.width) & (response != "presence")){
      stop("You must specify strip.width for strip transects!")
    }
  }

  # if the ps are all the same (count model) then just grab the 1 unique
  # value
  if(response %in% c("N","abundance")){
    fitted.p <- unique(fitted.p)
  }

  ## Aggregate response values of the sightings over segments
  if(response %in% c("D","density")){
    responsedata <- aggregate(obsdata[,cluster.name]/(fitted.p*availability),
                                list(obsdata[,segnum.name]), sum)
    off.set <- "none"
  }else if(response %in% c("N","abundance")){
    responsedata <- aggregate(obsdata[,cluster.name]/availability,
                              list(obsdata[,segnum.name]), sum)
    off.set <- "eff.area"
  }else if(response %in% c("Nhat","abundance.est")){
    responsedata <- aggregate(obsdata[,cluster.name]/(fitted.p*availability),
                              list(obsdata[,segnum.name]), sum)
    off.set<-"area"
  }else if(response == "presence"){
    responsedata <- aggregate(obsdata[,cluster.name],
                                list(obsdata[,segnum.name]), sum)
    responsedata$x[responsedata$x>0] <- 1
    responsedata$x[responsedata$x<1] <- 0
    off.set <- "none"
  }

  # name the response data columns
  names(responsedata)<-c(segnum.name,response)

  # Next merge the response variable with the segment records and any
  # response variable that is NA should be assigned 0 because these
  # occur due to 0 sightings
  dat <- merge(segdata,responsedata,by=segnum.name,all.x=T)
  dat[,response][is.na(dat[,response])] <- 0

  if(!is.null(segment.area)){

    # pull the column if segment.area is character
    if(is.character(segment.area)){
      segment.area <- dat[,segment.area]
    }

    dat$off.set <- switch(off.set,
                          eff.area=segment.area*fitted.p,
                          area=segment.area,
                          none=1)

    # if we have density then use that as the response
    if(response %in% c("D","density","Dhat","density.est")){
      dat[,response] <- dat[,response]/(segment.area*convert.units)
    }

  }else{
    # pull this from the detection function
    if (!is.null(ddfobject)){
      # calculate the "width" of the transect first, make sure we get it right
      # if we are doing left truncation
      width <- ddfobject$meta.data$width
      if(!is.null(ddfobject$meta.data$left)){
        width <- width - ddfobject$meta.data$left
      }
    }else{
    # or use strip.width if we have strip transects
      width <- strip.width
      # note we have to reset off.set
      if(response == "presence"){
        off.set <- "none"
      }else{
        off.set <- "area"
      }
    }

    # calculate the offset
    #   area we just calculate the area
    #   effective area multiply by p
    #   when density is response, offset should be 1 (and is ignored anyway)

    dat$off.set <- switch(off.set,
                          eff.area=2*dat[,seglength.name]*width*fitted.p,
                          area=2*dat[,seglength.name]*width,
                          none=1)

    # calculate the density (count/area)
    if(response %in% c("D","density","Dhat","density.est")){
      dat[,response] <- dat[,response]/(2*dat[,seglength.name]*
                                        width*convert.units)
    }
  }



  # multiply up by conversion factor
  dat$off.set <- dat$off.set*convert.units

  # Set offset as log of area or effective area
  dat$off.set <- log(dat$off.set)

  return(dat)
}
