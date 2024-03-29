#' Check for autocorrelation in residuals
#'
#' Once a DSM has been fitted to data, this function can be used to check for
#' autocorrelation in the residuals.
#'
#' @param dsm.obj a fitted dsm object.
#' @param Transect.Label label for the transect (default: `Transect.Label`).
#' Using different labels can be useful when transects are split over
#' geographical features or when transects are surveyed multiple times.
#' @param Segment.Label label for the segments (default: `Segment.Label`).The
#' result of calling [`order`][base::order] must make sense.
#' @param resid.type the type of residuals used, see
#' [`residuals.gam`][residuals.gam]. Defaults to `"scaled.pearson"` in the GAM
#' case and `"normalized"` in the GAMM case (which are equivalent).
#' @param fun the function to use, by default [`cor`][stats::cor], must take two
#' column vectors as arguments.
#' @param max.lag maximum lag to calculate at.
#' @param ylim user defined limits in y direction.
#' @param subset which subset of the data should the correlation function be
#' calculated on?
#' @param ... other options to pass to [`plot`][base::plot].
#'
#' @return a plot or a vector of `fun` applied at the lags.
#'
#' @section Details: Within each `Transect.Label`, segments will be sorted
#' according to their `Segment.Label`s. This may require some time to get right
#' for your particular data. If one has multiple surveys where transects are
#' revisited, for example, one may want to make `Transect.Label` a unique
#' transect-survey identifier. Neither label need to be included in the model,
#' they must just be present in the `$data` field in the model. This usually
#' means that they have to be in the segment data passed to [`dsm`][dsm].
#'
#'The current iteration of this function will only plot correlations nicely,
#'other things are up to you but you can get the function to return the data
#'(by assigning the result to an object).
#'
#' If there are NA values in the residuals then the correlogram will not be
#' calculated. This usually occurs due to `NA` values in the covariates (so the
#' smoother will not have fitted values there). Code like
#' `any(is.na(dsm.obj$data))` might be helpful.
#'
#' @importFrom graphics plot axis box legend lines abline segments
#' @importFrom stats cor residuals
#' @examples
#' \donttest{
#'  library(Distance)
#'  library(dsm)
#'
#'  # load the data, see ?mexdolphins
#'  data(mexdolphins)
#'
#'  # fit a model
#'  hr.model <- ds(distdata, truncation=6000,
#'                 key = "hr", adjustment = NULL)
#'  mod1 <- dsm(count~s(x,y), hr.model, segdata, obsdata)
#'
#'  # look at lag 1 differences up to a maximum of lag 9, using deviance
#'  # residuals
#'  dsm_cor(mod1, resid.type="deviance", max.lag=9,
#'          Segment.Label="Sample.Label")
#'}
#' @author David L. Miller
#' @export
dsm_cor <- function(dsm.obj, Transect.Label="Transect.Label",
                    Segment.Label="Segment.Label", max.lag=10,
                    resid.type="scaled.pearson",
                    fun=cor, ylim=c(0, 1), subset="all", ...){


  # only deal with the gam object
  if("gamm" %in% class(dsm.obj)){
    # pull the data out
    dat <- dsm.obj$gam$data
    #dsm.obj <- dsm.obj$gam
    dsm.obj <- dsm.obj$lme

    # residual type correction for gamm
    # "scaled.pearson" in GAM is the same as "normalized" in lme
    if(resid.type=="scaled.pearson"){
      resid.type <- "normalized"
    }
    # pull residuals
    resids <- residuals(dsm.obj, type=resid.type, level=1)
  }else{
    # pull the data out
    dat <- dsm.obj$data
    # pull residuals
    resids <- residuals(dsm.obj,type=resid.type)
  }


  # do some subsetting
  if(any(subset!="all")){
    dat <- dat[subset,]
    resids <- resids[subset]
  }

  # grab the labels
  if(is.character(Transect.Label)){
    tr.labs <- dat[[Transect.Label]]
  }else{
    tr.labs <- Transect.Label
  }
  if(is.null(tr.labs)){
    stop(paste0("No column called ",Transect.Label," in data"))
  }

  seg.labs <- dat[[Segment.Label]]
  if(is.null(seg.labs)){
    stop(paste0("No column called ",Segment.Label," in data"))
  }
  ### done with checking

  # storage
  lag.list <- list()

  tr.labs.u <- unique(tr.labs)

  # built the lags for all transects
  for(this.tr.lab in tr.labs.u){

    ind <- tr.labs == this.tr.lab
    these.resids <- resids[ind]

    # sort!
    dat2 <- dat[ind,]
    these.resids <- these.resids[order(dat2[[Segment.Label]])]

    # if there are NAs then something has gone wrong
    if(any(is.na(these.resids))){
      stop(paste0("There are NAs in residuals. Correlogram cannot be computed.",
                  "\nNB this is usually due to NA covariate values."))
    }

    # over all the lags
    for(lag in 1:max.lag){
      if(lag<length(these.resids)){
        # build an indicator
        lag.ind <- seq(lag,length(these.resids),1)[-1]
        lag.len <- length(lag.ind)

        lag.list[[as.character(lag)]] <- rbind(lag.list[[as.character(lag)]],
# testing
#                               cbind((1:length(these.resids))[1:lag.len],
#                                      lag.ind))
                               cbind(these.resids[1:lag.len],
                                     these.resids[lag.ind]))
      }
    }
  }

  # now work through the list we built, applying fun (correlations!)
  cors <- lapply(lag.list,function(x){fun(x[,1],x[,2])})

  # sample size for each line in the plot
  cors.samp.size <- c(length(resids),unlist(lapply(lag.list,length)))

  # assume cor for now...
  plot(x=c(0,as.numeric(names(lag.list))),
       y=c(1,unlist(cors)),ylim=ylim,
       xlab="Lag",ylab="Correlation",type="n",axes=FALSE,...)
  axis(1,at=c(0,as.numeric(names(lag.list))))
  axis(2)
  box()
  segments(x0=c(0,as.numeric(names(lag.list))),
           y0=rep(0,length(unlist(cors))+1),
           y1=c(1,unlist(cors)))
  # plot the zero line
  abline(h=0)

  # plot a 95% CI +/-2/sqrt(N)
  # where N is the sample size for that lag
  lines(x=c(-1,as.numeric(names(lag.list))),
        y=2/sqrt(cors.samp.size),lty=2)
  lines(x=c(-1,as.numeric(names(lag.list))),
        y=-2/sqrt(cors.samp.size),lty=2)


  invisible(cors)
}
