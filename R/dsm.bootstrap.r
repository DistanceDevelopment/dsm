dsm.bootstrap.pre.post <- function(ddf.ID=ddf.ID, dsm.ID=dsm.ID, prediction.ID=prediction.ID,
                          method="movblock", replicates=replicates, blocksize=blocksize,
                          detectionrefit=FALSE, alpha=0.95, fname=fname,
                          field.flag=FALSE, cell.area=NULL, path=path, 
				  bpfile=bpfile, boxplot.coef=1.5)
{
#
#     Function:  pre-and post process information enroute to and returning from Louise's
#                bootstrapping routines
#
#   Arguments:
#     ddf.ID - ID number of ddf analysis to use, set to NULL if strip transect
#     dsm.ID - ID number of dsm model to use
#     prediction.ID - ID number of predict analysis to use for prediction grid
#     method - bootstrapping procedure (movingblock on nonparametric) to use
#     replicates - number of bootstrap replicates
#     blocksize - block size for moving block
#     detectionrefit - flag for including uncertainty in detection function i
#                      n estimates of precision [not currently employed]
#     alpha - magnitude of confidence interval
#     fname - file where cell-specific se's are written
#     field.flag - flag indicating whether prediction grid cell arera is fixed or field name (passed thru to bootstrap)
#     cell.area - size of each cell in prediction grid (used to convert density to abundance)
#     path - path to temporary directory to find file 'var.dat.r'
#     bpfile - filename where counter of progress is written by variance_par_mb
#	boxplot.coef - coefficient to determine how far out 'whiskers' extend from box; 
#				larger values make it more difficult for an observation to be an outlier, 
#				big numbers (e.g. 3, will reject few outliers, making large bootstrap variance)
#
#   Value:
#     vector of cell-by-cell standard errors over the prediction grid   written to 'fname'
#       list of
#               lower.abund.ci - lower bound for study area-wide CI, not adjusted for detection uncertainty
#               upper.abund.ci - upper bound for study area-wide CI, not adjusted for detection uncertainty
#               raw.bootstrap.values - all boostrap estimates of abundance not adjusted for detection uncertainty
#       and perhaps confidence interval across entire study area, and at sub-study area regions
#
#       Created 8/9 February, 3, 13 April, 27 June 2006 by Rexstad

#   Step 1, data preparation
#     Attach 'transect' to data used to build the detection function (ddf object)
  vardata.in <- read.table(file=paste(path,"var.dat.r",sep=""),header=TRUE, sep="\t")
#     rename the second column in var.dat.r to "Transect" as expected by Louise
  names(vardata.in)[2] <- "Transect"
  if (!is.null(ddf.ID)) ddfdata.name <- paste("ddf.dat.", ddf.ID, sep="")
#       Following statement does not work, because ddf data names segment labels "label", not
#           "Sample.label" as in the var.dat.r file; hence cannot merge upon the specified by=
#  ddf.data <- merge(eval(parse(text=ddfdata.name)), vardata.in, by="Sample.Label")

#     Attach 'transect' to data used for the dsm object

  dsmdata.name <- paste("dsm.", dsm.ID, "$result$data", sep="")
  ready.to.overwrite <- merge(eval(parse(text=dsmdata.name)), vardata.in, by="Sample.Label")
  new.dsm.object <- eval(parse(text=paste("dsm.",dsm.ID,sep=""))) 
  new.dsm.object$result$data <- ready.to.overwrite
#     Glue together names of ddf, dsm, and prediction R objects
  if (!is.null(ddf.ID)) ddf.object.name <- eval(parse(text=paste("ddf.", ddf.ID,sep="")))
  dsm.object.name <- eval(parse(text=paste("dsm.",dsm.ID,sep="")))  
  predict.object.name <- eval(parse(text=paste("dsm.predict.", prediction.ID,sep="")))
  
#     Call appropriate bootstrapping function, based on the value of 'method'
#       I'm thinking about not allowing writing to a file within the functions, or at least not
#         writing to the file reserved for containing the cell-by-cell SEs.

#      Set working directory to the temp directory defined by Distance
#         This will enable writing the bootstrap counter to simply write a filename, without
#         need to know the path, or pass the pathname to Louise's functions
  current.working.dir <- getwd()
  setwd(path)  
  if (method=="movblock") {
      bootstraps <- param.movblk.variance(dsm.object=new.dsm.object,pred.object=predict.object.name,
                       samp.unit.name='Transect',B=replicates, block=blocksize, field=field.flag, 
                       cell.name=cell.area, bpfile=bpfile)
#           values returned are: 'short.var' holding cell-specific sufficient stats, and
#                                'study.area.total' holding replicate-specific densities
      } else {
      print(cat("Sorry, we are unable to compute nonparametric bootstrap variance\n"))
#      bootstraps <- dsm.variance.nonpar.f(ddf.object=ddf.object.name, dsm.object=dsm.object.name,
#                            grid=grid.object.name, samp.unit.name='Transect', B=replicates)
      }
#     Send standard errors for each cell to a file for subsequent read into GIS layer  
#       having applied machine formula with sufficient statistics contained in the return
#       values from call to 'variance.par.mb.f'
  machine.var <- (bootstraps$short.var$sumx.sq - (bootstraps$short.var$sumx^2/replicates)) / (replicates-1)
  vector.cell.se <- sqrt(machine.var)
  
#   Done with the bootstrapping, restore correct working directory for remainder of function
  setwd(current.working.dir)
  write.table(vector.cell.se, file=fname, quote=FALSE, col.names=FALSE, sep="\t")
  
#   produce confidence interval for overall study area abundance
#
#         first step will be to trim largest x% of bootstraps; 
#               Percentage defined as amount necessary to bring median and trimmed mean within 8% of each other
#         these are defined as 'outliers,' and will be printed to note their existence, but
#         removed from all subsequent computations
 untrimmed.bootstraps <- bootstraps$study.area.total
 outliers <- boxplot.stats(untrimmed.bootstraps, coef=boxplot.coef)$out
 bootstrap.abund <- untrimmed.bootstraps[!(untrimmed.bootstraps %in% outliers)]
 trim <- length(outliers) / length(untrimmed.bootstraps)
  cat('\tBootstrap variance computations:  outlier removal\t\n')
  cat('         Using the method associated with boxplots (',
        boxplot.coef,' * interquartile range) above uppper quartile\n')
  cat('The following bootstrap abundance estimates were eliminated as upper ',
       round(trim*100,0), '% of bootstrap distribution\n')
  print(outliers, digits=2)
  
#  vector.abund.studyarea <- apply(bootstrap.abund, MARGIN=2, sum, na.rm=TRUE)  # sum over cells
  print(cat('\tResponse surface/Variance Plot:  Bootstrap distribution\t\n'))
  hist(bootstrap.abund[is.finite(bootstrap.abund)], xlab="Estimated values", 
       main=paste("Distribution of bootstrap estimates (after trimming largest", round(trim*100,0), "%)"))
  boot.infinite <- sum(is.infinite(bootstrap.abund))
  boot.NA <- sum(is.na(bootstrap.abund))
  boot.NaN <- sum(is.nan(bootstrap.abund))
  boot.outlier <- length(outliers)
  boot.median <- median(bootstrap.abund[is.finite(bootstrap.abund)], na.rm=TRUE)
  estimate.from.predict <- sum(predict.object.name$result, na.rm=TRUE)
  boot.cv <- sqrt(var(bootstrap.abund[is.finite(bootstrap.abund)], na.rm=TRUE)) / estimate.from.predict
  cat('\tBootstrap confidence interval for abundance within study area\t\n')
  cat('After trimming, there were ', sum(is.finite(bootstrap.abund)), " legitimate values, and ", 
             boot.infinite+boot.NA+boot.NaN+boot.outlier, " non-legitimate bootstrap replicates \n")
  cat("           ", boot.outlier, " Outliers, ", boot.infinite, " Infinites, ", 
             boot.NA, " NAs, and ", boot.NaN, " NaNs \n")
  cat("  Percentile method computed confidence intervals for density surface component only:\n")
  print(quantile(bootstrap.abund[is.finite(bootstrap.abund)], c((1-alpha)/2, 1-((1-alpha)/2)),na.rm=TRUE))
  
  cat("\n\n Point estimate ", estimate.from.predict, "   SE of bootstraps ", 
        sqrt(var(bootstrap.abund[is.finite(bootstrap.abund)], na.rm=TRUE)), 
        "\n       Estimated CV for density surface model ",   round(boot.cv,4), "\n")
#
#       approximate CV and CI incorporating uncertainty in detection probabilities
#
#         fixed bug ER 06/02/2008 15:31:33--"summary.ds" replaced by "summary" next 2 lines
   average.p <- summary(eval(parse(text=paste("ddf.", ddf.ID, sep=""))))$average.p
   average.p.se <- summary(eval(parse(text=paste("ddf.", ddf.ID, sep=""))))$average.p.se
   cv.detection <- average.p.se / average.p
   
   unconditional.cv.square <- cv.detection^2 + boot.cv^2
   asymp.ci.c.term <- exp(1.96*sqrt(log(1+unconditional.cv.square)))
   lci.asymp <- estimate.from.predict / asymp.ci.c.term
   uci.asymp <- estimate.from.predict * asymp.ci.c.term
   
   cat("\n\n  Approximate measure of precision incorporating detection function uncertainty\n")
   cat("      CV in detection probability = ", round(cv.detection,4), "\n")
   cat("      CV in overall estimate including density surface model AND detection probability = ", 
              round(sqrt(unconditional.cv.square),4), "\n")
   cat("      Confidence interval incorporating detection function uncertainty = ", lci.asymp, ", ",
       uci.asymp, "\n")

#   Other calculations to perform if there are subareas within studyarea of interest  
#       Problematic because although dsm.layers.X still resides within .Rdata, there no
#       longer remains the labels associated with them.  They were in a temporary file that evap.

#   Construct name of layer objects, from prediction.id
    layer.object.name <- eval(parse(text=paste("dsm.layers.", prediction.ID,sep="")))
  
#    There may be multiple subarea layers over which abundance is being estimated, 
#     and hence we need to perform  'aggregate' for each of those layers
    num.of.sublayers <- length(colnames(layer.object.name))
    if (num.of.sublayers != 2) {
     bootstrap.abund <- data.frame(bootstrap.abund, layer.object.name)  # I suggest substituting 'machine.var' for 'bootstrap.abund'
     for (sublayer in 2:num.of.sublayers)
    {
#     Test for layer in which all ids are zero (meaning this is the study area layer)
      if (sum(layer.object.name[,sublayer]) != 0)
        {
        subunit.sums <- aggregate(bootstrap.abund, 
                        list(ID.name=layer.object.name[,sublayer]), sum, na.rm=TRUE)
#    Eliminate the sum for cells with ID.name=0, meaning they were not in a subunit of interest
        subunit.sums <- subset(subunit.sums, ID.name!=0)
#    Eliminate 3 extraneous columns, i.e,  ID, ID.name and final column (with user-specified name) 
#           Trickery is needed to determine name of final column, and pass this to subset to remove
        subunit.columns <- length(colnames(subunit.sums))
        subunit.sums <- subset(subunit.sums, select=-c(ID.name, (subunit.columns-1):subunit.columns))  
#     Produce sub-unit specific quantiles
       cat('\tBootstrap confidence intervals for the sub-units specified\t\n')
       print(apply(subunit.sums[is.finite(subunit.sums)], MARGIN=1, quantile, 
              c((1-alpha)/2, 1-((1-alpha)/2)),na.rm=TRUE))
      }  # end of check regarding whether this is the study area layer  
    } # end of sublayer loop  
    }  #  if around sublayer loop if there is only 1
#     Results sent back from this function will be upper and lower bounds
#         of total study area abundance; this will be the DsmVarObject in VB
#         and these will be passed to the stats file
  sink()
  dsm.bootstrap.pre.post <- list(
       lower.abund.ci=quantile(bootstrap.abund[is.finite(bootstrap.abund)], (1-alpha)/2,na.rm=TRUE), 
       upper.abund.ci=quantile(bootstrap.abund[is.finite(bootstrap.abund)], 1-(1-alpha)/2,na.rm=TRUE),
       untrimmed=untrimmed.bootstraps,
       raw.bootstrap.values=bootstrap.abund)
 }                      
