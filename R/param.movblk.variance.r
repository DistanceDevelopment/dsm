## TODO
# pretty-ify
# get rid of dead code
# split files
# re-write?
# add Tweedie stuff into get.model.details.f
# documentation
# add in detection function uncertainty stuff

#' Variance estimation via parametric moving block bootstrap
#'

param.movblk.variance <- function(n.boot, dsm.object, pred.object, 
                                  samp.unit.name='Transect',block.size, 
                                  field=FALSE, cell.name, bpfile=NULL){

#' @param n.boot number of bootstrap resamples.
#' @param dsm.object object returned from \code{\link{dsm.fit()}}.
#' @param pred.grid a \code{data.frame} that holds prediction points, must have
#'        the correct columns. 
#' @param samp.unit.name name sampling unit to resample (default 'Transect').
#' @param block.size number of segments in each block.
#' @param field=FALSE [[this and cell.name are linked, see below]]
#' @param cell.name=cell [[size of the cell to convert density to abundance]]
#' @param bpfile path to a file to be used (usually by Distance) to get 
#'        generate a progress bar (default \code{NULL} -- not used).

  #prediction.call <- pred.object$call
  # previously, prediction.call[3] held the name for the prediction data
  # let's just fix that now...
  pred.data<-pred.object
  prediction.call<-c("","","pred.data")

  # Tricky construction of cell size vector from field of prediction 
  # cell R object
  if(field){
    cell.string <- paste(prediction.call[3],"$",as.character(cell.name), sep="")
    cell.size <- eval(parse(text=cell.string))
  }else{
    cell.size <- cell.name
  }

  # Initialize some data structures for use in variance computations
  study.area.total <- numeric(length=0)
  short.var <- data.frame(sumx=rep(0,length(pred.object[[1]])), 
                          sumx.sq=rep(0,length(pred.object[[1]])))

  # Sort out sampling unit for dsm object
  su.string <- paste("dsm.object$result$data$",samp.unit.name,sep="")
  dsm.object$result$data$sampling.unit <- eval(parse(text=su.string))

  # Following line removes any transect into which missing data was 
  # detected by the call to gam and recorded in '$na.action'
  # Consequence of this step should be that no sampling unit (transect)
  # can become part of the bootstrap when any segment has missing data.
  name.sampling.unit <- unique(dsm.object$result$data$sampling.unit)
  num.sampling.unit <- length(name.sampling.unit)

  # Reconstruct dsm model fitting command - Note not using dsm.fit() 
  model.details <- get.model.details.f(dsm.object$result)
  dsm.fit.command <- get.dsm.fit.command.f(model.details)

  # Reconstruct prediction command
  pred.command <- paste("dsm.predict(gam.model=dsm.bs.model,newdata=",
                        prediction.call[3],",field=FALSE,off=0, silent=TRUE)",
                        sep="")

  # Get residuals 
  dsm.object$result$data$log.resids <- 
               resids.when.missing(model=dsm.object$result)

  # Sort out blocks for each sampling unit
  block.info <- block.info.per.su(block.size=block.size,
                                  data=dsm.object$result$data,
                                  name.su=name.sampling.unit)
  tot.num.blocks <- sum(block.info$num.block)
  k <- sum(block.info$num.req)
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

    bs.blocks <- sample(block.vector, k, replace=TRUE)
    bs.resids <- generate.mb.sample(k = k, block.size = block.size, 
                       which.blocks = bs.blocks, x = dsm.object$result$data, 
                       unit.info = block.info, n.units = num.sampling.unit)

    # Back transform to get bootstrap observations
    bs.samp <- dsm.object$result$data
    fit <- fitted(dsm.object$result)
    if(sum(is.na(fit)) > 0 ){
      stop("Missing values detected in survey covariates, cannot be used with moving block")
    }
    bs.samp$N <- fit*exp(bs.resids)  

    # Fit model to dsm bootstrap sample
    dud.replicate <- FALSE
    # Bit of a cheat because dud.replicate needs to be returned out 
    #  of the function(err) to be acted on 	
    # Handle chaos in gam fitting caused by pathological bootstrap resample
    tryCatch(dsm.bs.model <- eval(parse(text=dsm.fit.command)), 
                                  error=function(err){dud.replicate <<-TRUE})

    if(!dud.replicate){
      # Do prediction using newly fitted dsm model created from bootstrap sample
      dsm.predict.bs <- tryCatch(eval(parse(text=pred.command)), 
                                 error= function(err) {rep(NA,length(fit))} )

      # Don't save all cell values for all reps, rather, 
      #  populate dataframe with machine formula components for each cell
      vector.cell.abundances <- dsm.predict.bs$result * cell.size
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
# block.size=block.size,data=dsm.data,name.su=unique.sampling.units
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

  return(unit)
}


# generates the vector of residuals which can be mapped back onto the data
generate.mb.sample <- function(k, block.size, which.blocks, x, 
                               unit.info, n.units){
# k = required.num.blocks
# block.len = len.block
# which.blocks = bs.blocks
# x = data
# unit.info = unit.info
# n.units = num.sampling.units

  bs <- NULL
  bs$block <- which.blocks
  bs <- data.frame(bs)
  for(i in 1:k){
    for(j in 1:n.units){
      if(bs$block[i] >= unit.info$start.block[j] & 
         bs$block[i] <= unit.info$end.block[j]){
        bs$unit.name[i] <- as.character(unit.info$name[j])
        if(j == 1){
          bs$unit.block[i] <- bs$block[i]
        }else{
          bs$unit.block[i] <- bs$block[i] - unit.info$end.block[j - 1]
        }
      }
    }
    x.unit <- x[x$sampling.unit == bs$unit.name[i], ]            
    x.unit <- data.frame(x.unit)
    start.row <- bs$unit.block[i]
    end.row <- start.row + block.size - 1
    x.block <- x.unit[start.row:end.row, ]
    x.block <- data.frame(x.block)
    if (i == 1){
      bs.data <- x.block
    }else{
      bs.data <- rbind(bs.data, x.block)
    }
  }
    
  # Now need to map this onto data vector so same length 
  # (ie chopping off unwanted bits of blocks)
  temp <- bs.data$log.resids
  for (j in 1:n.units) {
    # Get number of observations in blocks
    tb <- unit.info$num.req[j] * block.size
    tran.resids <- temp[1:unit.info$num[j]]
    temp <- temp[-(1:tb)] 
    if (j==1){
      all.resids <- tran.resids
    }else{
      all.resids <- c(all.resids,tran.resids) 
    }
  } 
  return(all.resids)
}


# Extracts relevant details from model
get.model.details.f <- function(modelname=modelname){
  type <- class(modelname)[1]
  lnk <- modelname$family$link
  family <- modelname$family$family
  depvar <- all.vars(modelname$formula)[1]
  # Need to get specific model depending on model family 
  famstr <- substr(family,1,5)

  # if negative binomial need to get theta
  if (famstr=="Negat"){
    theta <- get.theta.f(family)
  }else{
    theta <- NULL
  }

  # if quasi need to get variance function
  if(famstr=="quasi" & nchar(family)==5){
    var <- modelname$family$varfun
  }else{ 
    var <- NULL
  }

  # Get formula (can't just use modelname$formula as composed of 3 bits)    
  form <- paste(modelname$formula[2],"~",modelname$formula[3],sep="")
  model <- list(type=type,lnk=lnk,family=family,famstr=famstr,
                theta=theta,var=var,form=form,depvar=depvar)

  return(model)
}


# Pastes together command for fitting glm/gam model
get.dsm.fit.command.f <- function(model=model.details){
  if (model$famstr == "Negat"){
    fam <- paste("negative.binomial(link=",model$lnk,",theta=",
                 model$theta,")",sep="")
  }else if(model$famstr=="quasi"&nchar(model$family)==5){ 
    fam <- paste(model$family,"(link=",model$lnk,",var=",model$var,")",sep="")
  }else{
    fam <- paste(model$family,"(link=",model$lnk,")",sep="")
  }
  command <- paste(model$type,"(",model$form,",data=bs.samp,family=",
                   fam,")",sep="")

  return(command)
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

  #return(res) ##?
}




# this was commented out anyway... replaced by redids.when.missing(), above?
# --------------------------------------------------------------------
#get.log.resids.f <- function(model=model.object) {
## Get residuals on log scale
#  obs <- model$data$N
#  fit <- fitted(model)
#  if ((min(obs)==0) | (min(fit)==0) | ((min(obs)==0)&(min(fit)==0))) {
#    obs <- obs + 1
#    fit <- fit + 1
#  }
#  res <- log(obs) - log(fit)
#res
#}
