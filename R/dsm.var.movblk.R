#' Variance estimation via parametric moving block bootstrap
#'
#' Estimate the variance in abundance over an area using a moving block
#' bootstrap. Two procedures are implemented, one incorporating detection
#' function uncertainty, one not.
#'
#'
#'
#'
#' @param n.boot number of bootstrap resamples.
#' @param dsm.object object returned from \code{\link{dsm.fit()}}.
#' @param pred.grid a \code{data.frame} that holds prediction points, must have
#'        the correct columns for other environmental covariates. 
#' @param ds.uncertainty incorporate uncertainty in the detection function?
#' @param samp.unit.name name sampling unit to resample (default 'Transect').
#' @param block.size number of segments in each block.
#' @param cell.size.name name of the column in pred.data corresponding to the
#'        size of the prediction cells (default NULL, see below).
#' @param cell.size scalar (or vector of length \code{nrow(pred.data)}) of the
#'        sizes of the prediction cells (default \code{NULL}, either 
#'        \code{cell.size} or \code{cell.size.name} MUST be specified).
#' @param bpfile path to a file to be used (usually by Distance) to get 
#'        generate a progress bar (default \code{NULL} -- no file written).
#' @export

## TODO
# documentation
# make up the missed replicates due to model errors?
# new sampler doesn't do point transects at the moment

# this used to be called param.movblk.variance
dsm.var.movblk <- function(n.boot, dsm.object, pred.data, 
                                  ds.uncertainty=FALSE,
                                  samp.unit.name='Transect',block.size, 
                                  cell.size.name=NULL, cell.size=NULL, 
                                  bpfile=NULL){

  # make sure either the size of the prediction cell or the name
  # of the corresponding column in pred.data is supplied
  if((is.null(cell.size) & is.null(cell.size.name))|
     (!is.null(cell.size) & !is.null(cell.size.name))){
    stop("You must supply either cell.size or cell.size.name!\n")
  }

  # construct cell size vector from field of prediction data
  if(is.null(cell.size)){
    cell.size<-pred.data[[cell.size.name]]
  }  

  # Initialize some data structures for use in variance computations
  study.area.total <- numeric(length=0)
  short.var <- data.frame(sumx=rep(0,nrow(pred.data)), 
                          sumx.sq=rep(0,nrow(pred.data)))

  # Sort out sampling unit for dsm object
  su.string <- paste("dsm.object$result$data$",samp.unit.name,sep="")
  dsm.object$result$data$sampling.unit <- eval(parse(text=su.string))

  # Following line removes any transect into which missing data was 
  # detected by the call to gam and recorded in '$na.action'
  # Consequence of this step should be that no sampling unit (transect)
  # can become part of the bootstrap when any segment has missing data.
  name.sampling.unit <- unique(dsm.object$result$data$sampling.unit)
  num.sampling.unit <- length(name.sampling.unit)

  # Reconstruct dsm model fitting command -- this is a call to gam()
  gam.call<-dsm.object$result$call
  gam.call$formula<-dsm.object$result$formula
  gam.call$family<-dsm.object$result$family

  # Get residuals 
  dsm.object$result$data$log.resids <- 
               resids.when.missing(model=dsm.object$result)

  # Sort out blocks for each sampling unit
  block.info <- block.info.per.su(block.size=block.size,
                                  data=dsm.object$result$data,
                                  name.su=name.sampling.unit)
  tot.num.blocks <- sum(block.info$num.block)
  num.blocks.required <- sum(block.info$num.req)
  block.vector <- 1:tot.num.blocks

  # Start bootstrapping
  for(i in 1:n.boot){
    # Compute proportion of bootstrapping completed, write it to file 
    # passed as argument 'bpfile'
    # This file will be read by Distance to present a progress bar.
    if(!is.null(bpfile)){
      progress <- round(i/n.boot, 2)* 100
      write(progress, file=bpfile, append=FALSE)
    }

    bs.blocks <- sample(block.vector, num.blocks.required, replace=TRUE)
    bs.resids <- generate.mb.sample(num.blocks.required, block.size, 
                                    bs.blocks, dsm.object$result$data, 
                                    block.info, num.sampling.unit)

    # Back transform to get bootstrap observations
    bs.samp <- dsm.object$result$data
    fit <- fitted(dsm.object$result)

    if(sum(is.na(fit)) > 0 ){
      stop(paste("Missing values detected in survey covariates,",
                 " cannot be used with moving block"))
    }
    bs.samp$N <- fit*exp(bs.resids)  

    # if we are incirporating detection function uncertainty, then 
    # need to resample the distances
    if(ds.uncertainty){

      #### This only deals with count data at the moment
      ####  => no individual level covariates
      #### Much of this can be put up top, doesn't need to be calculated
      #### each time

      # how many distances to generate? -- need to round
      n.ds.samples<-round(sum(bs.samp$N),0)
      bs.samp$N<-round(bs.samp$N,0)

      # how many samples do we have so far?
      n.samps<-0
      dists<-c()

      # create an object to hold the parameters
      pars<-list()
      ds.object<-dsm.object$ddf
      pars$scale<-ds.object$ds$aux$ddfobj$scale$parameters
      if(!is.null(ds.object$ds$aux$ddfobj$shape$parameters)){
        pars$shape<-ds.object$ds$aux$ddfobj$shape$parameters
      }
      if(!is.null(ds.object$ds$aux$ddfobj$adjustment$parameters)){
        pars$adjustment<-ds.object$ds$aux$ddfobj$adjustment$parameters
      }


      # do some rejection sampling
      while(n.samps < n.ds.samples){

        # generate some new distances
        new.dists<-data.frame(distance=runif(n.ds.samples-n.samps)*
                                          ds.object$ds$aux$width,
                              detected=rep(1,n.ds.samples-n.samps),
                              object=1:(n.ds.samples-n.samps))

        U<-runif(n.ds.samples-n.samps)

        # need to call out to mrds to get the data and model objects
        # into the correct format
        xmat <- mrds:::process.data(new.dists,ds.object$meta.data,
                                    check=FALSE)$xmat
        ddfobj <- mrds:::create.ddfobj(ds.object$call$dsmodel,xmat,
                            ds.object$meta.data,pars)

        # do the rejection...
        # (evaluate the -log(L) then backtransform per-observation
        # ONLY line transect at the moment!!
        inout <- exp(-mrds:::flt.lnl(ds.object$par,ddfobj,
                      misc.options=list(width=ds.object$ds$aux$width,
                                        int.range=ds.object$ds$aux$int.range,
                                        showit=FALSE, doeachint=TRUE,
                                        point=ds.object$ds$aux$point,
                                        integral.numeric=TRUE),TCI=FALSE))>U
        dists<-c(dists,new.dists$distance[inout])

        n.samps<-length(dists)
      }
      # make sure that we got the right number
      dists<-dists[1:n.ds.samples]
      dists<-data.frame(distance=dists,
                        detected=rep(1,length(dists)),
                        object=1:length(dists))

      # fit the new model
      ddf.call <- dsm.object$ddf$call
      ddf.call$data <- dists
      ddf.fitted<-eval(ddf.call)

      # find the offset(s)
      bs.samp$off.set<-rep(fitted(ddf.fitted,compute=TRUE,esw=TRUE)[1],nrow(bs.samp))

    }

    # Fit model to dsm bootstrap sample
    dud.replicate <- FALSE
    # Bit of a cheat because dud.replicate needs to be returned out 
    #  of the function(err) to be acted on 	
    # Handle chaos in gam fitting caused by pathological bootstrap resample
    #tryCatch(dsm.bs.model <- eval(parse(text=dsm.fit.command)), 
    #                              error=function(err){dud.replicate <<-TRUE})

    # put the boostrap data into the gam call
    gam.call$data<-bs.samp

    tryCatch(dsm.bs.model <- eval(gam.call, parent.frame()), 
                                  error=function(err){dud.replicate <<-TRUE})

    if(!dud.replicate){

      # is the offset supplied as a field in the predication data?
      if(!is.null(pred.data$off.set)){
        offset.field<-TRUE
        offset.field.off<-NULL
      }else{
        offset.field<-FALSE
        offset.field.off<-0
      }

      # Do prediction using newly fitted dsm model created from bootstrap sample
      dsm.predict.bs <- tryCatch(dsm.predict(dsm.bs.model,newdata=pred.data,
                                             field=offset.field,
                                             off=offset.field.off), 
                                 error= function(err) {rep(NA,length(fit))} )

      # Don't save all cell values for all reps, rather, 
      #  populate dataframe with machine formula components for each cell
      vector.cell.abundances <- dsm.predict.bs * cell.size
      short.var$sumx <- short.var$sumx + vector.cell.abundances
      short.var$sumx.sq <- short.var$sumx.sq + 
                           (vector.cell.abundances * vector.cell.abundances)
      study.area.total <- append(study.area.total, 
                                 sum(vector.cell.abundances, na.rm=TRUE))
    }else{
      study.area.total <- append(study.area.total, NA)
    }
  }


  result <- list(short.var=short.var, study.area.total=study.area.total)
}


# sort out the block information 
block.info.per.su <- function(block.size,data,name.su){
# block.size=block.size
#data=dsm.data
#name.su=unique.sampling.units
  unit <- NULL   
  unit$name <- name.su
  num.su <- length(name.su)
  for (i in 1:num.su) {
    unit$num[i] <- length(
                     data$sampling.unit[data$sampling.unit == unit$name[i]])
    unit$num.block[i] <- unit$num[i] - block.size + 1
    if(i == 1){
      unit$start.block[i] <- 1
    }else{
      unit$start.block[i] <- unit$start.block[i - 1] + unit$num.block[i - 1]
    }
  }
  unit$end.block <- cumsum(unit$num.block)
  # Get number of blocks required for each sampling unit  -- 
  #  samp units different sizes so need to get more observations than required 
  unit$num.req <- ceiling(unit$num/block.size)
  unit$num.req[unit$num.req==0] <- 1
  unit <- data.frame(unit)

  # returns a data frame, one row per sampling unit (e.g. transect)
  #  each row consists of:
  #  - name of the sampling unit  
  #  - number of segments in that unit
  #  - number of blocks available (= # segments - # blocks +1)
  #  - block # for start of segment 
  #  - block # for end of segment 
  #  - # of blocks needed for the sampling unit (> # needed) 
  
  return(unit)
}


# generates the vector of residuals which can be mapped back onto the data
generate.mb.sample <- function(num.blocks.required, block.size, which.blocks, 
                                dsm.data, unit.info, n.units){
# num.blocks.required
# block.len = len.block
# which.blocks = bs.blocks
# dsm.data
# unit.info = unit.info
# n.units = num.sampling.units

  bs <- NULL
  bs$block <- which.blocks
  bs <- data.frame(bs)

  for(i in 1:num.blocks.required){
    ## find the sampling unit that the block is in 

    # is the block a start or end point for sampling unit?
    j<-which(bs$block[i] == unit.info$start.block)
    if(length(j)==0){
      j<-which(bs$block[i] == unit.info$end.block)
    }
    # if not a start or end, then where is it?
    if(length(j)==0){
      find.block<-c(unit.info$start.block,bs$block[i])
      j<-which(sort(find.block)==bs$block[i])-1
    }

    bs$unit.name[i] <- as.character(unit.info$name[j])
    if(j == 1){
      bs$unit.block[i] <- bs$block[i]
    }else{
      bs$unit.block[i] <- bs$block[i] - unit.info$end.block[j - 1]
    }

    # pull out the data for this sampling unit
    x.unit <- dsm.data[dsm.data$sampling.unit == bs$unit.name[i], ]            
    x.unit <- data.frame(x.unit)

    # pull out the rows corresponding to this block
    # start.row is the block number and the end row is
    #  (block length) segments after that
    start.row <- bs$unit.block[i]
    end.row <- start.row + block.size - 1
    x.block <- x.unit[start.row:end.row, ]
    x.block <- data.frame(x.block)

    # create a new frame if necessary and store the data
    if (i == 1){
      bs.data <- x.block
    }else{
      bs.data <- rbind(bs.data, x.block)
    }
  }
    
  # Now need to map this onto data vector so same length 
  # (ie chopping off unwanted bits of blocks)

  temp <- bs.data$log.resids

  # storage
  bs.response <- c()

  # loop over the sample units
  for (j in 1:n.units) {
    # Get number of segments in the blocks
    tb <- unit.info$num.req[j] * block.size
    # grab enough residuals for this sample unit
    tran.response <- temp[1:unit.info$num[j]]
    # remove all of the ones we sampled (i.e. if we over sampled
    # make sure that we get rid of them too
    temp <- temp[-(1:tb)] 

    # store the result
    bs.response <- c(bs.response,tran.response) 
  } 
  return(bs.response)
}

# This function is substitute for "get.log.residuals.f"
# Modifications were needed to handle the possibility of missing values in 
# the fitted values. These missing cases need to correspond by placing missing
# values into the response variable (obs <-model$data$N)
resids.when.missing <- function(model=model.object){
  obs <- model$data$N

  # set to missing data for which there are missing values of covariates
  obs[model$na.action] <- NA  

  temp.fit <- fitted(model)
  fit <- numeric(length(obs))
  used.fit <- 1
  for (i in 1:length(obs)){
    if (!is.na(obs[i])){
      fit[i] <- temp.fit[used.fit]
      used.fit <- used.fit + 1
    }else{
      fit[i] <- NA
    }
  }

  # plus 0.001 (arbitrary) to avoid logging zero   
  res <- log(obs+0.001) - log(fit+0.001)  

  return(res)
}

