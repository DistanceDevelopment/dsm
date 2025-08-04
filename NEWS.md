# dsm 2.3.4

- fix bug where if predictions were NA (e.g., predictions outside the soap boundary), all variances were NA due to propagation (#4)
- fix CRAN NOTES regarding documentation links and special characters.

# dsm 2.3.3

* To improve consistency in functions in the package, some functions have changed from . separation to _. For function arguments some separators have been switched from _ to .. An error is now thrown when the "old" arguments/functions using . are used. This error will be removed in dsm 2.3.4.
  - dsm.var.prop -> dsm_var_prop
  - dsm.var.gam -> dsm_var_gam
  - dsm.var.movblk -> dsm_var_movblk
  - dsm.cor -> dsm_cor
  - rqgam.check -> rqgam_check
  - vis.concurvity -> vis_concurvity
  - dsm_varprop
    - var_type -> var.type

# dsm 2.3.2

* fixed bug in offset calculation in dsm_varprop se estimation where offsets were logged twice. Thanks to Megan Ferguson for spotting this.
* documentation now in rmarkdown format

- dsm 2.3.1

* transect type (line or point) is now determined from the detection function rather than using the transect= argument
* multiple detection functions are now possible, use a list() of fitted detection functions supplied to ddf.obj= argument to dsm(). dsm.var.gam, dsm.var.prop and dsm_varprop are compatible. Vignette coming soon. Thanks to Dave Fifield and Ewan Wakefield for extensive testing and suggestions. See ?"dsm-data" for more information on how to set up your data.
* strip transects now handled using the dummy_ddf() function.
* presence~ response type removed
* Response options "D", "density", "Dhat", "N", "Nhat" and "n" deprecated. Now only "count", "abundance.est" and "density.est" are allowed.
* Make predict.dsm work with tibbles when off.set is not a column in newdata (error is thrown in previous code's check of this). Thanks to Dave Fifield for reporting this issue.
* Fix bug with inconsistent predictions when density.est is used as the response. Code in predict.dsm has been simplified and documentation updated. Thanks to Dave Fifield for reporting this issue.
* mrds independent observer (io) and trial models are allowed for the detection function.Variance can be estimated via dsm_varprop (with the usual limitation that covariates may only vary at the segment level) or dsm.var.gam. This feature is currently experimental. Thanks to Doug Sigourney and Megan Ferguson for their help and suggestions. See ?"dsm-data" for more information on how to set up your data.
* An example of multiddf and mrds analysis in dsm is at https://github.com/densitymodelling/nefsc_fin_mrds_dsm
* availability= argument now works differently for count models. availability must now be the same length as number of segments (nrow(segment.data)) and will multiply the offset (effective area) rather than the number of animals observed to be more coherent with how detectability is handled in this case.
* dsm_varprop no longer requies newdata to be supplied. If newdata isn't present just the refitted and original models are returned, with other entries NA.

# dsm 2.3.0

* fixed dsm_varprop (gam.fixed.priors) for response distributions with fixed scale parameters. Negative binomial results before now may have been incorrect!

# dsm 2.2.18

* Clarification of the predict.dsm return documentation. Thanks to Dave Fifield for the suggestion.
* off.set=NULL lead to predictions of 0 when density was used. Now just returns density. Thanks to Dave Fifield for the catch.


# dsm 2.2.17

* Now Suggests package "sp" (needed for some parts of the mexdolphin example data)

# dsm 2.2.16

* check that any Sample.Labels match up, if none do then throw an error. Thanks to Madhura Davate for bringing this issue to my attention.
* gamma option removed, since it now is not ignored by gam(). If using method="GCV.Cp" then gamma=1.4 should still be used (but generally we advise using REML rather than GCV for smoothing parameter selection)
* dsm now requires mgcv version 1.8-23 or higher.
* dsm.var.prop (and dsm_varprop) now throw an error when there are no covariates in the detection function as it can't do anything useful
* new function obs_exp is a diagnostic function to compare observed and expected counts at different covariate aggregations

# dsm 2.2.15

* when mexdolphins is loaded using data(), the data no longer need to be attach()ed
* dsm.var.prop throws an error when count or N are not used as the response.
* dsm.var.prop now calls dsm_varprop, which keeps the scale parameter constant when refitting models for variance propagation. dsm_varprop provides a lower-level interface to variance propagation and is probably faster and more flexible.
* add vis.concurvity to visually assess concurvity (colinearity extended to smooth terms) in the model
* new model diagnostic function: plot_pred_by_term -- plots the effect of each smooth in the model spatially.
* remove "abundance" option for the response, use "count"
* dsm.var.prop and dsm_varprop have better diagnostics for checking the variance model hasn't done something weird
* dsm.var.movblk discards replicates where the model didn't converge
* fix potential issue with negative binomial response is used with variance propagation (due to scale parameterisations). dsm_varprop now refits nb() models as negbin() but we advise users to look into this more deeply. Thanks to Eric Rexstad for spotting this.
* fix bug in predict.dsm() where predictions for "density" response were not including the offset. Thanks to Len Thomas for spotting this.

# dsm 2.2.14

* nasty bug in predict.dsm fixed: was logging the offset twice when newdata was not provided, leading to very odd checking plots. This also effected the results of predict() when newdata= was not specified. Thanks to Megan Ferguson for alerting us to this.
* add support for point transects(!)
* rqgam.check() now produced ONLY the linear predictor vs. randomised quantile residuals plot to avoid confusion in interpreting the other plots
* REML is now default method for smoothing parameter selection

# dsm 2.2.13

* dsm.var.gam now has an example
* dsm.var.prop now checks that the detection function model does not have covariates

# dsm 2.2.12

* The alpha argument now works for analytical variance reported using dsm.var.prop and dsm.var.gam print.summary and summary methods. Thanks to Jason Roberts and Eric Rexstad for spotting this issue.
* dsm now uses numDeriv::grad to calculate derivatives for variance estimation -- this should be much faster. Thanks to Jason Roberts for this fix.
* dsm now requires mrds version 2.1.15 or higher, due to internal changes
* fixed issue where Effort column was required even when segment.area was specified. Thanks to Greg Lollback for spotting this issue.
* ensured that ordering of segments matched the ordering of the observations so that aggregation of observations are correct. Thanks to Dave Fifield for reporting this bug and suggesting a fix.
* check.cols no longer checks to see if there is a column called "distance" when strip transects are being used

# dsm 2.2.11

* dsm.var.prop uses update() rather than messing around with eval()
* print.summary.dsm.var didn't work for GAM variance when no detection function was specified, now fixed (thanks to Natalia Dellabianca for spotting this)
* Removed non-working plot() code in variance examples --  see the Mexico pantropical spotted dolphins for a proper treatment of plotting using projections (http://distancesampling.org/R/vignettes/mexico-analysis.html)
* plot.var.gam no longer crashes by default when distance data doesn't have columns "x" and "y"

# dsm 2.2.10

* corrected offsets in examples in dsm, dsm.var.prop, dsm.var.movblk
* updated mexdolphins data with lat/long reversal in prediction data and non-list storage
* stopped use of rqgam.check from being used with non-negbin/Tweedie responses (thanks to Phil Bouchet for highlighting this issue).

# dsm 2.2.9

* Moving block bootstrap now works when there is no detection function (suggested by Laura Williamson)
* Mexico dolphin data has had the areas of the prediction grids and segments re-projected. Thanks to Phil Bouchet for pointing out this error. Context: https://groups.google.com/d/msg/distance-sampling/oj1j8mSmnTc/3cOAAJgBw9MJ


# dsm 2.2.8

* Speed up and reduced memory usage in variance calculations, thanks to Filipe Dias for asking about this.
* Better error when no observations are matched to the detection function. Thanks to Ricardo Lima for highlighting this.
* Check for missing columns had incorrect logic, now corrected and with tests. Thanks to Julia Migné for spotting this.

# dsm 2.2.7

* Plotting check plots using random quantile residuals (rqgam.check) didn't work for Tweedie or tw family models, now fixed. Thanks to Laura Mannocci and Jason Roberts for pointing this out.
* Plotting CV used matrix() to coerce data into a data.frame -- thus converting all elements to the simplest type. Errors ensued. Thanks to Adrian Schiavini for finding this bug!

# dsm 2.2.6

* dsm.var.prop and dsm.var.gam now allow you to supply a single number as an offset, then copy it as many times as you need. Thanks to Filipe Dias for suggesting this.

# dsm 2.2.5

* dsm.var.prop will now throw an error if there is no detection function in the dsm that it's trying to calculate the variance for (thanks to Adrián Schiavini for finding this)
* variance/CV plotting threw an error when no detection function was supplied but observation=TRUE. Now error is sidestepped (again, thanks to Adrián Schiavini for finding this).

# dsm 2.2.4

* changed stop() to warning() when data not in the detection function are included in the observation data frame, as they are automatically excluded later, bug found by Eric Rexstad and Anna Cucknell

# dsm 2.2.3

* moving block boostrap was using observed abundances as the residuals and therefore gave wrong answers (Eric Rexstad found this bug)

# dsm 2.2.2

* version number reporting now works (Len Thomas found this bug)
* more documentation
* by default density and presence models are now weighted by area of each segment

# dsm 2.2.1

* switch over to new testthat spec
* added new response variables (Dhat, density.est, count, n), suggestion from Len Thomas

# dsm 2.1.4

* check that all observed distances are within width will now throw an error

# dsm 2.1.3

* CRAN-compatibility fixes

# dsm 2.1.2

* updates to dsm.cor
* more unified plotting
* variance estimation tweaks for non-dsm models

# dsm 2.0.6

* removed warning about if() comparison length in make.soapgrid need to pass bnd as a list()
* dsm.cor now handles gamm objects
* dsm.cor allows one to check residual autocorrelation in your models
* presence/absence data now no longer needs a dummy detection function object

# dsm 2.0.5

* started NEWS file
* fixed bug where dsm.var.gam() didn't return predictions so plotting the result would fail
* silly typo in docs
* check if we're using mgcv version >= 1.7-24 for GAMM fitting and warn if not

