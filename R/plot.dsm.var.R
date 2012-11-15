#' Create plots of abundance uncertainty
#'
#' Note that the prediction data set must have \code{x} and \code{y} columns
#' even if these were not used in the model.
#'
#' @S3method plot dsm.var
#' @method plot dsm.var
#' @aliases plot.dsm.var
#' 
#' @param x a \code{dsm.var} object
#' @param poly a \code{list} or \code{data.frame} with columns \code{x} and 
#'        \code{y}, which gives the coordinates of a polygon to draw.
#' @param limits limits for the fill colours
#' @param breaks breaks for the colour fill
#' @param legend.breaks breaks as they should be displayed
#' @param xlab label for the \code{x} axis
#' @param ylab label for the \code{y} axis
#' @param observations should observations be plotted?
#' @param plot actually plot the map, or just return a \code{ggplot2} object?
#' @param boxplot.coef control trimming (as in \code{summary.dsm.var}), only
#'        has an effect if the bootstrap file was saved.
#' @param x.name name of the variable to plot as the x axis.
#' @param y.name name of the variable to plot as the y axis.
#' @param gg.grad optional \code{\link{ggplot}} gradient object.
#' @param \dots any other arguments
#' @return a plot 
#'
#' @section Details:
#'
#'  In order to get plotting to work with \code{\link{dsm.var.prop}} and
#'  \code{\link{dsm.var.gam}}, one must first format the data correctly since
#'  these functions are designed to compute very general summaries. One summary
#'  is calculated for each element of the list \code{pred} supplied to
#'  \code{dsm.var.prop} and \code{dsm.var.gam}.
#'
#'  For a plot of uncertainty over a prediction grid, \code{pred} (a
#'  \code{data.frame}), say, we can create the correct format by simply using
#'  \code{pred.new <- splot(pred,1:nrow(pred))}.
#' 
#' @author David L. Miller
#'
### TODO
# covariates in the detection function

plot.dsm.var<-function(x, poly=NULL, limits=NULL, breaks=NULL,
                       legend.breaks=NULL, xlab="x", ylab="y", 
                       observations=TRUE, plot=TRUE, boxplot.coef=1.5, 
                       x.name="x", y.name="y", gg.grad=NULL, ...){

  # I am exactly this lazy, sorry
  object <- x
  rm(x)

  # if the data isn't formatted correctly...
  if(length(object$pred)==1){
    stop("Looks like you're calling plot on a whole area summary, see ?plot.dsm.var for how to format your data correctly for plotting")
  }

  # if we did used the random effects trick, collapse everything down
  if(!object$bootstrap){
    pd<-c()
    off<-c()
    for(i in 1:length(object$pred.data)){
      pd<-rbind(pd,object$pred.data[[i]])
      off<-rbind(off,object$off.set[[i]])
    }
    object$pred.data <- pd
    object$off.set <- as.vector(off)
  }

  if(!all(c("width","height") %in% names(object$pred.data))){
      stop("No spatial data to create plot, need columns 'width' and 'height' in prediction data")
  }

  # predictions on the grid, from dsm.var.*
  mod.pred <- unlist(object$pred)

  if(object$bootstrap){

    #sinfo$block.size <- object$block.size 
    #sinfo$n.boot <- object$n.boot
    #sinfo$bootstrap <- TRUE
    #sinfo$ds.uncertainty <- object$ds.uncertainty

    # if we didn't save each bootstrap replicate
    if(is.null(object$bs.file)){
      # bootstrap cell abundances
      short.var <- object$short.var

      ### think about trimming here -- save everything and trim at this
      ### stage, will probably need to check that the var won't be too big

      n <- length(short.var$sumx.sq)

      cell.se<-sqrt((short.var$sumx.sq/(n-1))-((short.var$sumx/n)^2*(n/(n-1))))
      cell.cv <- cell.se/mod.pred

    }else{
      # if we did save each replicate...

      # load the data
      bs.save <- read.csv(object$bs.file,header=FALSE)

      # first col is just the ids
      bs.save <- bs.save[,-1]

      n <- ncol(bs.save)

      # don't do this because it's not consistant with what's in the summary()
      #cell.se <- sqrt(apply(trim.var,1,bs.save))
      #cell.cv <- cell.se/mod.pred

      # calculate the overall trimmed variance
      tv <- trim.var(object$study.area.total,boxplot.coef=boxplot.coef) 
      # there is an attribute that is an indicator of which to keep
      trim.ind <- attr(tv,"trim.ind")
      # keep those
      bs.save <- bs.save[trim.ind,]

      # calculate the variance
      cell.se <- apply(bs.save,2,sd)
      cell.cv <- cell.se/mod.pred

      rm(bs.save)
      gc()
    }

    # delta method, if necessary
    if(!object$ds.uncertainty){

      ddf.object<-object$dsm.object$ddf
      ddf.summary<-summary(ddf.object)

      ## calculate the variance via the delta method
      # find the cv squared of the p
      cvp.sq <- unique((ddf.summary$average.p.se/
                 ddf.summary$average.p))[1]^2

      # cv squared of the Ns from the bootstrap
      cell.cv.sq <- (cell.cv)^2

      # delta method
      cell.cv <- sqrt(cvp.sq+cell.cv.sq)
    }

  }else if(object$bootstrap==FALSE){
    # varprop stuff
    # pull out the standard errors
    if(is.null(dim(object$pred.var))){
      cell.se <- sqrt(object$pred.var)
    }else{
      cell.se <- sqrt(diag(object$pred.var))
    }

    cell.cv <- cell.se/mod.pred
  }

  # if we have some limits plot them
  if(is.null(limits)){
    limits <- c(min(cell.cv),max(cell.cv))
  }

  if(is.null(breaks)){
    breaks <- quantile(cell.cv)
    names(breaks) <- NULL
    #breaks <- seq(min(cell.cv),max(cell.cv),len=5)
    breaks <- round(breaks,2)
  }
  if(is.null(legend.breaks)){
    legend.breaks <- breaks
  }

  # put all the data together
  plotdata<-cbind(object$pred.data,cell.cv)

  # what are the names of the variables to be plotted on the x and y
  # axis? Rename!
  plotdata$x <- plotdata[[x.name]]
  plotdata$y <- plotdata[[y.name]]

  # build the plot
  gg.opts <- theme(panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  panel.background=element_blank(),
                  legend.key=element_blank())
  p <- ggplot(plotdata) + gg.opts
  p <- p + geom_tile(aes(x=x, y=y, fill=cell.cv, width=width, height=height))
  p <- p + coord_equal()

  if(is.null(gg.grad)){
    p <- p + scale_fill_gradient(low="white", high="black",
                                 limits=limits,
                                 rescaler = function(x, ...) x, oob = identity,
                                 breaks=legend.breaks)
  }else{
    p <- p + gg.grad
  }

  if(!is.null(poly)){
    poly$x <- poly[[x.name]]
    poly$y <- poly[[y.name]]
    p <- p+geom_path(aes(x=x, y=y),data=poly)
  }

  if(observations){
    object$dsm.object$data$x <- object$dsm.object$data[[x.name]]
    object$dsm.object$data$y <- object$dsm.object$data[[y.name]]
    object$dsm.object$ddf$data$x <- object$dsm.object$ddf$data[[x.name]]
    object$dsm.object$ddf$data$y <- object$dsm.object$ddf$data[[y.name]]

    p <- p + geom_line(aes(x=x, y=y,group=Transect.Label),
                        data=object$dsm.object$data)
    p <- p + geom_point(aes(x, y, size=size), data=object$dsm.object$ddf$data,
                        colour="blue",alpha=I(0.7))
    p <- p + labs(fill="CV",x=xlab,y=ylab, size="Counts")
  }else{
    p <- p + labs(fill="CV",x=xlab,y=ylab)
  }

  if(plot){
    # plot!
    print(p)
    invisible()
  }else{
    return(p)
  }
}
