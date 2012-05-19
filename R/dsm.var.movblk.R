#' Variance estimation via parametric moving block bootstrap
#'
#' Estimate the variance in abundance over an area using a moving block
#' bootstrap. Two procedures are implemented, one incorporating detection
#' function uncertainty, one not.
#'
#' @param dsm.object object returned from \code{\link{dsm.fit()}}.
#' @param pred.data a \code{data.frame} that holds prediction points, must have
#'        the correct columns for other environmental covariates. It also MUST
#'        have a column called \code{cell.area} which gives the area for each
#'        prediction cell 
#' @param n.boot number of bootstrap resamples.
#' @param block.size number of segments in each block.
#' @param off.set offset to be used in the model, see \code{\link{offsets}}.
#'        Note that this is NOT logged.
#' @param ds.uncertainty incorporate uncertainty in the detection function?
#'        Note that this doesn't work for covariate models for the detection
#'        function at the moment.
#' @param samp.unit.name name sampling unit to resample (default 
#'        'Transect.Label').
#' @param bpfile path to a file to be used (usually by Distance) to 
#'        generate a progress bar (default \code{NULL} -- no file written).
#' @export

## TODO
# documentation
# make up the missed replicates due to model errors?
# detection function uncertainty:
#  * new sampler doesn't do point transects at the moment
#  * new sampler doesn't do covariates at the moment
#  * is the sample size correct?
#  * should we be calculating the offset in here?
#  * non-log link functions

# this used to be called param.movblk.variance
dsm.var.movblk <- function(dsm.object, pred.data, n.boot, block.size, 
                           off.set, ds.uncertainty=FALSE,
                           samp.unit.name='Transect.Label',
                           bpfile=NULL){

  # check the user didn't ask for individual level covars and detection
  # function uncertainty
  if(ds.uncertainty & (
          dsm.object$model.spec$response=="indiv.est"|
          dsm.object$model.spec$response=="group.est")){
    stop("Detection function uncertainty with covariates is not supported")
  }

  # Initialize storage
  study.area.total <- numeric(n.boot)
  short.var <- data.frame(sumx=rep(0,nrow(pred.data)), 
                          sumx.sq=rep(0,nrow(pred.data)))

  # Sort out sampling unit for dsm object
  dsm.object$result$data$sampling.unit <- 
                        dsm.object$result$data[[samp.unit.name]]

  # Following line removes any transect into which missing data was 
  # detected by the call to gam and recorded in '$na.action'
  # Consequence of this step should be that no sampling unit (transect)
  # can become part of the bootstrap when any segment has missing data.
  name.sampling.unit <- unique(dsm.object$result$data$sampling.unit)
  num.sampling.unit <- length(name.sampling.unit)

  # Get residuals 
  # this replaces call to -----> resids.when.missing(dsm.object$result)
  # need to do something about non-log links
  ######################################
  obs <- dsm.object$result$data$N
  obs[dsm.object$result$na.action]<-NA
  fit.vals <- rep(NA,length(obs))
  fit.vals[!is.na(obs)]<-fitted(dsm.object$result)
  # plus 0.001 (arbitrary) to avoid logging zero   
  dsm.object$result$data$log.resids <- log(obs+0.001) - log(fit.vals+0.001)
  ######################################

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

    # if we are incirporating detection function uncertainty, then 
    # need to resample the distances
    if(ds.uncertainty){

      # save the old probability of detection
      old.p<-fitted(dsm.object$ddf)[1]

      #### This only deals with count data at the moment
      ####  => no individual level covariates
      #### Much of this can be put up top, doesn't need to be calculated
      #### each time

      # call out to generate some ds data, and fit a model to that data,
      # asumming that the detection function model is correct
      seed <- get(".Random.seed",envir=.GlobalEnv) # messes with the seed...
      new.p<-generate.ds.uncertainty(dsm.object$ddf)
      assign(".Random.seed",seed,envir=.GlobalEnv) # recover the seed

      this.offset<-new.p*(exp(dsm.object$result$offset)/old.p)

      fit <- fitted(dsm.object$result)/exp(dsm.object$result$offset)*this.offset

      # replace the offset in the model 
      bs.samp$off.set<-log(this.offset)

      # replace the offset in the prediction grid
      off.set<-new.p*(off.set/old.p)
    }else{
    # if we're not doing detection function uncertainty
      fit <- fitted(dsm.object$result)
    }

    if(sum(is.na(fit)) > 0 ){
      stop(paste("Missing values detected in survey covariates,",
                 " cannot be used with moving block"))
    }
    bs.samp$N <- fit*exp(bs.resids)  

    ## Fit model to dsm bootstrap sample

    # Reconstruct dsm model fitting command -- this is a call to gam()
    gam.call<-dsm.object$result$call
    gam.call$formula<-dsm.object$result$formula
    gam.call$family<-dsm.object$result$family
    # put the boostrap data into the gam call
    gam.call$data<-bs.samp
    
    # Handle chaos in gam fitting caused by pathological bootstrap resample
    dsm.bs.model <- try(eval(gam.call)) 

    if(all(class(dsm.bs.model)!="try-error")){

      dsm.bs.model<-list(result=dsm.bs.model)
      class(dsm.bs.model)<-"dsm"

      # Do prediction using newly fitted dsm model created from bootstrap sample
      dsm.predict.bs <- try(dsm.predict(dsm.bs.model,
                                        newdata=pred.data,
                                        off=off.set))
      if(all(class(dsm.predict.bs)=="try-error")){
        dsm.predict.bs <- rep(NA,length(fit))
      }

      # Don't save all cell values for all reps, rather, 
      #  populate dataframe with machine formula components for each cell
      vector.cell.abundances <- dsm.predict.bs 
      short.var$sumx <- short.var$sumx + vector.cell.abundances
      short.var$sumx.sq <- short.var$sumx.sq + 
                           (vector.cell.abundances * vector.cell.abundances)
      study.area.total[i] <- sum(vector.cell.abundances, na.rm=TRUE)
    }else{
      study.area.total[i] <- NA
    }
  }

  result <- list(short.var=short.var, study.area.total=study.area.total)

  return(result)
}
