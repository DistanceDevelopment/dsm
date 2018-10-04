library(dsm)
library(Distance)
library(testthat)

par.tol<-1e-5

context("test inputs")
# load the Gulf of Mexico dolphin data
data(mexdolphins)

# fit a detection function
suppressMessages(hn.model <- ds(distdata, max(distdata$distance),
                                adjustment = NULL))

test_that("formula specs",{

  ## models for count
  count.reml <- 936.0362722
  count.count <- dsm(count~s(x,y), hn.model, segdata, obsdata)
  expect_equal(unname(count.count$gcv.ubre), count.reml, tolerance=par.tol, check.attributes=FALSE)

  ## models for abund.est
  abund.est.gcv <- 992.0182117
  abund.est.abund.est<-dsm(abundance.est~s(x,y), hn.model, segdata,
                           obsdata)
  expect_equal(unname(abund.est.abund.est$gcv.ubre),abund.est.gcv,
                                                  tolerance=par.tol, check.attributes=FALSE)

  ## models for density
  D.reml <- -2812.33881
  density.density<-dsm(density~s(x,y), hn.model, segdata,
                       obsdata,
                       weights=rep(1,nrow(segdata)))
  expect_equal(unname(density.density$gcv.ubre), D.reml,
               tolerance=par.tol, check.attributes=FALSE)

  density.density.est<-dsm(density.est~s(x,y), hn.model, segdata,
                           obsdata,
                           weights=rep(1,nrow(segdata)))
  expect_equal(unname(density.density.est$gcv.ubre), D.reml,
               tolerance=par.tol, check.attributes=FALSE)

  # check that Effort is not zero
  mex_zero_effort <- segdata
  mex_zero_effort$Effort[c(1,5,10)] <- 0
  expect_error(dsm(abundance.est~s(x,y), hn.model, mex_zero_effort,
                           obsdata))

})

test_that("Missing columns cause errors",{

  seg <- segdata
  obs <- obsdata

  for(mcov in c("object","Sample.Label","size","distance")){
    obs_missing <- obsdata
    obs_missing[[mcov]] <- NULL
    expect_error(dsm(count~s(x,y), hn.model, seg, obs_missing),
                 paste0("Column(s) \"",mcov,
                        "\" not found in observation.data."),
                 fixed=TRUE)
  }

  for(mcov in c("Effort","Sample.Label")){
    seg_missing <- seg
    seg_missing[[mcov]] <- NULL
    expect_error(dsm(count~s(x,y), hn.model, seg_missing, obs),
                 paste0("Column(s) \"",mcov,
                        "\" not found in segment.data."),
                 fixed=TRUE)
  }

  # with segment area specified we only have a problem with Sample.Label
  seg_missing <- seg
  seg_missing[["Sample.Label"]] <- NULL
  expect_error(dsm(count~s(x,y), hn.model, seg_missing, obs),
               paste0("Column(s) \"Sample.Label\" not found in segment.data."),
               fixed=TRUE)

})

test_that("Error thrown when Sample.Labels don't match up",{

  seg <- segdata
  obs <- obsdata

  # now nothing will match
  seg$Sample.Label <- paste0(seg$Sample.Label, "-wat")

  expect_error(dsm(count~s(x,y), hn.model, seg, obs),
               "No matches between segment and observation data.frame Sample.Labels!",
               fixed=TRUE)

})


