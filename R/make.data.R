#' @importFrom stats aggregate
make.data <- function(response, ddfobject, segdata, obsdata, group,
                      convert.units, availability, segment.area,
                      family){

  # probably want to do something smart here...
  seglength.name<-'Effort'
  segnum.name<-'Sample.Label'
  distance.name<-'distance'
  cluster.name<-'size'

  # Estimating group abundance/density
  if(group){
    obsdata[, cluster.name][obsdata[, cluster.name] > 0] <- 1
  }

  # for single ddfs, make a list of 1
  if(any(c("ddf", "dsmodel")  %in% class(ddfobject))){
    ddfobject <- list(ddfobject)
    segdata$ddfobj <- 1
  }else{
    if(!any("ddfobj" %in% names(segdata))){
      stop("If there are multiple detection functions there must be a column named \"ddfobj\" in the segment frame, see ?\"dsm-data\"")
    }
  }

  # iterate over the list of ddfs
  full_obsdata <- c()
  segdata$segment.area <- NA
  segdata$width <- NA
  transect <- rep(NA, length(ddfobject))

  for(i in seq_along(ddfobject)){

    this_ddf <- ddfobject[[i]]

    # make this a mrds ddf object if we had a Distance one
    if(all(class(this_ddf)=="dsmodel")){
      this_ddf <- this_ddf$ddf
      ddfobject[[i]] <- this_ddf
    }
    # store the transect type for later
    transect[i] <- c("line", "point")[this_ddf$meta.data$point+1]

    # grab the probabilities of detection
    fitted.p <- fitted(this_ddf)

    # remove observations which were not in the detection function
    this_obsdata <- obsdata[obsdata$object %in% names(fitted.p), ]

    # Check that observations are between left and right truncation
    # warning only -- observations are excluded below
    # No truncation check for strip transects
    if(any(this_obsdata[, distance.name] > this_ddf$meta.data$width)){
      warning(paste("Some observations are outside of detection function", i,
              "truncation!"))
    }

    # reorder the fitted ps, making sure that
    # they match the ordering in obsdata
    fitted.p <- fitted.p[match(this_obsdata$object, names(fitted.p))]

    # set the segment area
    # calculate the "width" of the transect first, make sure we get it right
    # if we are doing left truncation
    width <- this_ddf$meta.data$width
    if(!is.null(this_ddf$meta.data$left)){
      width <- width - this_ddf$meta.data$left
    }

    # what if there are no matches? Perhaps this is due to the object
    # numbers being wrong? (HINT: yes.)
    if(nrow(this_obsdata) == 0){
      stop(paste("No observations in detection function", i,
                 "matched those in observation table. Check the \"object\" column."))
    }

    # put that all together

    # make a column for the detectabilities
    this_obsdata$p <- fitted.p

    # bind this to the full data
    full_obsdata <- rbind(full_obsdata, this_obsdata)

    # set the segment area for this detection function in the segdata
    segdata$width[segdata$ddfobj==i] <- width
    if(transect[i]=="point"){
      # here "Effort" is number of visits
      segdata$segment.area[segdata$ddfobj==i] <- pi *
        segdata$width[segdata$ddfobj==i]^2 *
        segdata[, seglength.name][segdata$ddfobj==i]
    }else{
    # line transects
      segdata$segment.area[segdata$ddfobj==i] <- 2 *
        segdata[, seglength.name][segdata$ddfobj==i] *
        segdata$width[segdata$ddfobj==i]
    }

  }
  # set the full observation data
  obsdata <- full_obsdata

  # set the segment area in the data
  if(!is.null(segment.area)){
    segdata$segment.area <- segment.area
  }


  ## Aggregate response values of the sightings over segments
  if(response %in% c("D","density","Dhat","density.est")){
    responsedata <- aggregate(obsdata[,cluster.name]/(fitted.p*availability),
                              list(obsdata[,segnum.name]), sum)
    off.set <- "none"
  }else if(response %in% c("N","count","n")){
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


  ## warn if any observations were not allocated
  responsecheck <- aggregate(obsdata[, cluster.name],
                             list(obsdata[, segnum.name]), sum)
  if(sum(obsdata[, cluster.name]) != sum(responsecheck[, 2])){
    message(paste0("Some observations were not allocated to segments!\n",
                   "Check that Sample.Labels match"))
  }

  # name the response data columns
  names(responsedata) <- c(segnum.name, response)

  # if the Sample.Labels don't match at all then we need to stop, nothing
  # can work as all the response values will be zero
  if(!any(segdata[, segnum.name] %in% responsedata[, segnum.name])){
    stop("No matches between segment and observation data.frame Sample.Labels!")
  }

  # Next merge the response variable with the segment records and any
  # response variable that is NA should be assigned 0 because these
  # occur due to 0 sightings
  dat <- merge(segdata, responsedata, by=segnum.name, all.x=TRUE)
  dat[, response][is.na(dat[, response])] <- 0

  # for the offsets with effective area, need to make sure that
  # the ps match the segments
  if(off.set == "eff.area"){
    dat$p <- NA
    for(i in seq_along(ddfobject)){
      this_ddf <- ddfobject[[i]]

      # if there are no covariates, and all the fitted ps are the same
      # then just duplicate that value enough times for the segments
      if(length(unique(fitted(this_ddf))) == 1){
        dat[dat$ddfobj == i, ]$p <- rep(fitted(this_ddf)[1],
                                        nrow(dat[dat$ddfobj == i, ]))
      }else if("fake_ddf" %in% class(this_ddf)){
        # predict for strip transect
        dat[dat$ddfobj == i, ]$p <- 1
      }else{

        if(this_ddf$method != "ds"){
          stop("Only \"ds\" models are supported!")
        }

        # extract formula
        df_formula <- as.formula(this_ddf$ds$aux$ddfobj$scale$formula)
        if(!is.null(this_ddf$ds$aux$ddfobj$shape$formula) &&
           this_ddf$ds$aux$ddfobj$shape$formula != "~1"){
          stop("Shape parameter formulae are not supported!")
        }

        # extract detection function variables
        df_vars <- all.vars(df_formula)

        # check these vars are in the segment table
        if(!all(df_vars %in% colnames(dat))){
          stop(paste0("Detection function covariates are not in the segment data",
                      "\n  Missing: ", df_vars[!(df_vars %in% colnames(dat))]))
        }

        # make a data.frame to predict for
        nd <- dat[dat$ddfobj == i, ][, df_vars, drop=FALSE]
        nd$distance <- 0

        dat[dat$ddfobj == i, ]$p <- predict(this_ddf, newdata=nd)$fitted
      }
    }
  }

  # check that none of the Effort values are zero
  if(any(dat[, seglength.name]==0)){
    stop(paste0("Effort values for segments: ",
                paste(which(dat[, seglength.name]==0), collapse=", "),
                " are 0."))
  }

  # correct segment area units
  dat$segment.area <- dat$segment.area*convert.units

  # calculate the offset
  #   area we just calculate the area
  #   effective area multiply by p
  #   when density is response, offset should be 1 (and is ignored anyway)
  dat$off.set <- switch(off.set,
                        eff.area = dat$segment.area*dat$p,
                        area     = dat$segment.area,
                        none     = 1)

  # calculate the density (count/area)
  if(response %in% c("D", "density", "Dhat", "density.est")){
    dat[, response] <- dat[, response]/(dat$segment.area)
  }


  # Set offset as log (or whatever link is) of area or effective area
  dat$off.set <- family$linkfun(dat$off.set)

  return(dat)
}
