#' Create plots of abundance uncertainty
#'
#' Note that the prediction data set must have \code{x} and \code{y} columns
#' even if these were not used in the model.
#'
#' @S3method plot dsm.var
#' @method plot dsm.var
#' @aliases plot.dsm.var
#' 
#' @param object a \code{dsm.var} object
#' @param poly a \code{list} or \code{data.frame} with columns \code{x} and 
#'        \code{y}, which gives the coordinates of a polygon to draw.
#' @param \dots any arguments that can usually be passed to a 
#' @param limits limits for the fill colours
#' \code{\link{plot}} method 
#' @return a plot 
#' 
#' @author David L. Miller
#'
### TODO
# plot observations
# plot transect lines
# covariates in the detection function


plot.dsm.var<-function(object, poly=NULL, ..., limits=NULL){

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

  if(!all(c("x","y","width","height") %in% names(object$pred.data))){
      stop("No spatial data to create plot, need columns 'x' and 'y' in prediction data")
  }

  # estimate from prediction
  mod.pred <- dsm.predict(object$dsm.object,
                          newdata=object$pred.data,
                          off=object$off.set)

  if(object$bootstrap){

    #sinfo$block.size <- object$block.size 
    #sinfo$n.boot <- object$n.boot
    #sinfo$bootstrap <- TRUE
    #sinfo$ds.uncertainty <- object$ds.uncertainty

    # bootstrap cell abundances
    short.var <- object$short.var

    ### think about trimming here -- save everything and trim at this
    ### stage, will probably need to check that the var won't be too big

    n <- length(short.var$sumx.sq)

    cell.se<-sqrt((short.var$sumx.sq/(n-1))-((short.var$sumx/n)^2*(n/(n-1))))
    cell.cv <- cell.se/mod.pred

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

  # put all the data together
  plotdata<-cbind(object$pred.data,cell.cv)

  # build the plot
  gg.opts <- opts(panel.grid.major=theme_blank(),
                  panel.grid.minor=theme_blank(),
                  panel.background=theme_rect())
  p <- ggplot(plotdata) + gg.opts
  p <- p + labs(fill="CV")
  p <- p+geom_tile(aes(x=x, y=y, fill=cell.cv, width=width, height=height))
  p <- p + coord_equal()
  p <- p + scale_fill_gradientn(colours=heat_hcl(200),
                                limits=limits)

  if(!is.null(poly)){
    p <- p+geom_path(aes(x=x, y=y),data=poly)
  }

  print(p)

  invisible()

}
