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
#'        Note that this should NOT be \code{log()}'d.
#' @param ds.uncertainty incorporate uncertainty in the detection function? See
#'        Details, below. Note that this feature is EXPERIMENTAL at the moment.
#' @param samp.unit.name name sampling unit to resample (default 
#'        'Transect.Label').
#' @param progress.file path to a file to be used (usually by Distance) to 
#'        generate a progress bar (default \code{NULL} -- no file written).
#' @param bs.file path to a file to store each boostrap round. This stores all
#'        of the bootstrap results rather than just the summaries, enabling
#'        outliers to be detected and removed. (Default \code{NULL}).
#' @param bar should a progress bar be printed to screen? (Default \code{TRUE}).
#'
#' @section Details: 
#'  Setting \code{ds.uncertainty=TRUE} will incorporate detection function 
#'  uncertainty directly into the bootstrap. This is done by generating 
#'  observations from the fitted detection function and then re-fitting a new 
#'  detection function (of the same form), then calculating a new effective 
#'  strip width. Rejection sampling is used to generate the observations (except
#'  in the half-normal case) so the procedure can be rather slow. Note that this
#'  is currently not supported with covariates in the detection function.
#'
#'  Setting \code{ds.uncertainty=FALSE} will incorporate detection function 
#'  uncertainty using the delta method. This assumes that the detection function
#'  and the spatial model are INDEPENDENT. This is probably not reasonable.
#'
#' @export

## TODO
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
                           progress.file=NULL, bs.file=NULL,bar=TRUE){

  # check the user didn't ask for individual level covars and detection
  # function uncertainty
  if(ds.uncertainty &
     dsm.object$ddf$ds$aux$ddfobj$scale$formula != "~1"){
#          dsm.object$model.spec$response=="indiv.est"|
#          dsm.object$model.spec$response=="group.est")){
    stop("Detection function uncertainty with covariates is not supported")
  }

  # Initialize storage
  study.area.total <- numeric(n.boot)
  short.var <- data.frame(sumx=rep(0,nrow(pred.data)), 
                          sumx.sq=rep(0,nrow(pred.data)))

  # save the original off.set that was supplied
  original.offset<-off.set

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
  if(dsm.object$model.spec$response %in% c("indiv.den","group.den")){
    obs <- dsm.object$result$data$D
  }else{
    obs <- dsm.object$result$data$N
  }
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

  # do we want to print a progress bar?
  if(bar){
    pb <- txtProgressBar(min=0,max=n.boot,style=3)
  }

  # first resample is the model predicitons themselves

  # if we supplied bs.file, then write detailed, per replicate
  # if it's the first time and the file already exists, modify the name 
  # use that
  if(!is.null(bs.file)){
    if(file.exists(bs.file)){
      a.number <- 1
      bs.file2 <- strsplit(bs.file,"\\.")[[1]]
      bs.file2.end <- paste('.',bs.file2[length(bs.file2)],sep="")
      bs.file2.start <- paste(bs.file2[1:(length(bs.file2)-1)],collapse=".")
      bs.file2 <- paste(bs.file2.start,"-",a.number,bs.file2.end,
                        collapse="",sep="")
      while(file.exists(bs.file2)){
        a.number <- a.number+1
        bs.file2 <- paste(bs.file2.start,"-",a.number,bs.file2.end,
                          collapse="",sep="")
      }
      warning(paste("Filename",bs.file,"was taken, writing to file",
                     bs.file2))
      bs.file<-bs.file2
    }
  }

  dsm.predict.bs <- try(dsm.predict(dsm.object,
                                    newdata=pred.data,
                                    off=off.set))
  if(all(class(dsm.predict.bs)=="try-error")){
    dsm.predict.bs <- rep(NA,length(short.var$sumx))
  }

  # Don't save all cell values for all reps, rather, 
  #  populate dataframe with machine formula components for each cell
  vector.cell.abundances <- dsm.predict.bs 
  short.var$sumx <- short.var$sumx + vector.cell.abundances
  short.var$sumx.sq <- short.var$sumx.sq + vector.cell.abundances^2
  study.area.total[1] <- sum(vector.cell.abundances, na.rm=TRUE)

  if(!is.null(bs.file)){
    # append to the file bs.file
    write.table(t(dsm.predict.bs),bs.file,append=TRUE,
                sep=",",col.names=FALSE)
  }

  off.set.save<-off.set

  # Start bootstrapping
  for(i in 2:n.boot){
    # Compute proportion of bootstrapping completed, write it to file 
    # passed as argument 'progress.file'
    # This file will be read by Distance to present a progress bar.
    if(!is.null(progress.file)){
      progress <- round(i/n.boot, 2)* 100
      write(progress, file=progress.file, append=FALSE)
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
      #seed <- get(".Random.seed",envir=.GlobalEnv) # messes with the seed...
      new.p<-generate.ds.uncertainty(dsm.object$ddf)
      #assign(".Random.seed",seed,envir=.GlobalEnv) # recover the seed

#cat("old=",old.p," new=",new.p,"\n")

      this.offset<-new.p*(exp(dsm.object$result$offset)/old.p)

      fit <- (fitted(dsm.object$result)/exp(dsm.object$result$offset))*
                                    this.offset

      # replace the offset in the model 
      bs.samp$off.set<-log(this.offset)

      # replace the offset in the prediction grid
      off.set<-new.p*(off.set.save/old.p)
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
    # if bnd or knots were used... 
    if(!is.null(gam.call$knots)){
      gam.call$knots <- dsm.object$model.spec$knots
    }
    if(!is.null(gam.call$bnd)){
      gam.call$bnd <- dsm.object$model.spec$bnd
    }

    # put the bootstrap data into the gam call
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
      short.var$sumx.sq <- short.var$sumx.sq + vector.cell.abundances^2
      study.area.total[i] <- sum(vector.cell.abundances, na.rm=TRUE)

      if(!is.null(bs.file)){
        # append to the file bs.file
        write.table(t(dsm.predict.bs),bs.file,append=TRUE,
                    sep=",",col.names=FALSE)
      }
    }else{
      study.area.total[i] <- NA
    }

    if(bar){
      setTxtProgressBar(pb, i)
    }

  }

  result <- list(short.var=short.var, 
                 study.area.total=study.area.total,
                 ds.uncertainty=ds.uncertainty,
                 bootstrap=TRUE,
                 pred.data=pred.data,
                 n.boot=n.boot, 
                 off.set=original.offset,
                 block.size=block.size,
                 bs.file=bs.file,
                 dsm.object=dsm.object
                )

  # give it a bit of class
  class(result)<-c("dsm.var")

  if(bar){
    cat("\n")
  }

  return(result)
}
