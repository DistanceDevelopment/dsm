#' Check for autocorrelation in residuals
#'
#' Once a DSM has been fitted to data, this function can be used to
#' check for autocorrelation in the residuals.
#'
#' @param dsm.obj a fitted dsm object.
#' @param Transect.Label label for the transect (default: \code{Transect.Label}). Using different labels can be useful when transects are split over geographical features.
#' @param resid.type the type of residuals used, see \code{\link{residuals.gam}}.
#' @param fun the function to use, by default \code{\link{cor}}, must take two column vectors as arguments.
#' @param max.lag maximum lag to calulate at.
#'
#' @return a plot or a vector of \code{fun} applied at the lags.
#'
#' @section Details: Note we assume that \code{Sample.Label} is in order within each transect, so sort data beforehand. The current iteration of this function will only plot correlations nicely, other things are up to you but you can get the function to return the data (by assigning the result to an object).
#'
#' @examples
#'
#'  library(Distance)
#'  library(dsm)
#'
#'  data(mexdolphins)
#'
#'  hr.model <- ds(mexdolphins$distdata, max(mexdolphins$distdata$distance), key = "hr", adjustment = NULL)
#'  mod1<-dsm(N~s(x,y), hr.model, mexdolphins$segdata, mexdolphins$obsdata)
#'
#'  dsm.cor(mod1,resid.type="d",max.lag=9)
#'
#' @author David L. Miller
#' @export
dsm.cor <- function(dsm.obj,Transect.Label="Transect.Label",max.lag=10, resid.type = c("deviance", "pearson","scaled.pearson","working", "response"),fun=cor){

  # pull the data out
  dat <- dsm.obj$data

  # pull residuals
  resids <- residuals(dsm.obj,type=resid.type)

  # grab the labels
  if(is.character(Transect.Label)){
    tr.labs <- dat[[Transect.Label]]
  }else{
    tr.labs <- Transect.Label
  }
  if(is.null(tr.labs)){
    stop(paste0("No column called ",Transect.Label," in data"))
  }

  ### done with checking

  # storage
  lag.list <- list()

  tr.labs.u <- unique(tr.labs)

  # built the lags for all transects
  for(this.tr.lab in tr.labs.u){

    ind <- tr.labs == this.tr.lab
    these.resids <- resids[ind]

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

  # now work through the list we built, applying fun

  cors <- lapply(lag.list,function(x){fun(x[,1],x[,2])})

  # assume cor for now...
  plot(x=c(0,as.numeric(names(lag.list))),
       y=c(1,unlist(cors)),
       xlab="Lag",ylab="Correlation",type="n",axes=FALSE)
  axis(1,at=c(0,as.numeric(names(lag.list))))
  axis(2)
  box()
  segments(x0=c(0,as.numeric(names(lag.list))),
           y0=rep(0,length(unlist(cors))+1),
           y1=c(1,unlist(cors)))
  # plot the zero line
  abline(h=0)
  # plot a 95% CI +/-2/sqrt(N)
  abline(h=2/sqrt(length(tr.labs)),lty=2)
  abline(h=-2/sqrt(length(tr.labs)),lty=2)

  invisible(cors)
}
