# error messages when dead functions are used

#' Prediction variance propagation for DSMs
#'
#' This function is deprecated, use [dsm_var_prop].
#' @export
#' @inheritParams dsm_var_prop
dsm.var.prop <- function(dsm.obj, pred.data, off.set,
                         seglen.varname='Effort', type.pred="response") {
  stop("dsm.var.prop is deprecated, use dsm_var_prop")
}

#' Prediction variance estimation assuming independence
#'
#' This function is deprecated, use [dsm_var_gam].
#' @export
#' @inheritParams dsm_var_prop
dsm.var.gam <- function(dsm.obj, pred.data, off.set,
                        seglen.varname='Effort', type.pred="response"){
  stop("dsm.var.gam is deprecated, use dsm_var_gam")
}

#' Variance estimation via parametric moving block bootstrap
#'
#' This function is deprecated, use [dsm_var_movblk].
#' @export
#' @inheritParams dsm_var_movblk
dsm.var.movblk  <- function(dsm.object, pred.data, n.boot, block.size,
                           off.set, ds.uncertainty=FALSE,
                           samp.unit.name='Transect.Label',
                           progress.file=NULL, bs.file=NULL,bar=TRUE){
  stop("dsm.var.movblk is deprecated, use dsm_var_movblk")
}

#' Check for autocorrelation in residuals
#'
#' This function is deprecated, use [dsm_cor].
#' @export
#' @inheritParams dsm_cor
dsm.cor <- function(dsm.obj, Transect.Label="Transect.Label",
                    Segment.Label="Segment.Label", max.lag=10,
                    resid.type="scaled.pearson",
                    fun=cor, ylim=c(0, 1), subset="all", ...){
  stop("dsm.cor is deprecated, use dsm_cor")
}

#' Randomised quantile residuals check plot for GAMs/DSMs
#'
#' This function is deprecated, use [rqgam_check].
#' @export
#' @inheritParams rqgam_check
rqgam.check <- function(gam.obj, ...){
  stop("rqgam.check is deprecated, use rqgam_check")
}

#' Visualise concurvity between terms in a GAM
#'
#' This function is deprecated, use [vis_concurvity].
#' @export
#' @inheritParams vis_concurvity
vis.concurvity <- function(model, type="estimate"){
  stop("vis.concurvity is deprecated, use vis_concurvity")
}
