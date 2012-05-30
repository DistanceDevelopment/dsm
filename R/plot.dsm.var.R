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

    if(is.null(limits)){
      limits <- c(min(cell.cv),max(cell.cv))
    }

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

    #### general bootstrap stuff

    ## how many duds did we have?
    #sinfo$trim.prop <- attr(trimmed.variance, "trim.prop")
    #sinfo$trim.ind <- attr(trimmed.variance, "trim.ind")
    #sinfo$boot.outliers <- attr(trimmed.variance, "outliers")
    #sinfo$boot.infinite <- sum(is.infinite(bootstrap.abund))
    #sinfo$boot.finite <- sum(!is.infinite(bootstrap.abund))
    #sinfo$boot.NA <- sum(is.na(bootstrap.abund))
    #sinfo$boot.NaN <- sum(is.nan(bootstrap.abund))
    #sinfo$boot.usable <- sinfo$boot.finite - sinfo$boot.outliers 

    ## grab the %ile c.i.s at alpha, 1-alpha and also median
    #sinfo$quantiles <- quantile(bootstrap.abund[sinfo$trim.ind], 
    #                            c(alpha, 0.5, 1-alpha),na.rm=TRUE)
    #attr(sinfo$quantiles,"names")[2] <- "Median"


    

  }else if(object$bootstrap==FALSE){
    # varprop stuff
    #sinfo$saved<-object
    #sinfo$bootstrap <- object$bootstrap
    #sinfo$se <- sqrt(object$pred.var)

    #sinfo$cv <- sinfo$se/sinfo$pred.est

    cat("Plotting for variance propagation not implemented.\n\n")
  }

  invisible()

}
