library(dsm)
library(Distance)
library(testthat)

par.tol<-1e-5

context("test inputs")
# load the Gulf of Mexico dolphin data
data(mexdolphins)
attach(mexdolphins)

# fit a detection function
suppressMessages(hn.model <- ds(distdata, max(distdata$distance),
                                adjustment = NULL))

test_that("formula specs",{

  ## models for count
  count.gcv <- 42.9169051
  count.N<-dsm(N~s(x,y), hn.model, segdata, obsdata)
  expect_equal(unname(count.N$gcv.ubre), count.gcv,tolerance=par.tol)

  count.n<-dsm(n~s(x,y), hn.model, segdata, obsdata)
  expect_equal(unname(count.n$gcv.ubre), count.gcv,tolerance=par.tol)

  count.count<-dsm(count~s(x,y), hn.model, segdata,
                   obsdata)
  expect_equal(unname(count.count$gcv.ubre), count.gcv,tolerance=par.tol)

  count.abundance<-dsm(abundance~s(x,y), hn.model, segdata,
                       obsdata)
  expect_equal(unname(count.abundance$gcv.ubre), count.gcv,tolerance=par.tol)


  ## models for abund.est
  abund.est.gcv <- 57.3159048
  abund.est.Nhat<-dsm(Nhat~s(x,y), hn.model, segdata,
                      obsdata)
  expect_equal(unname(abund.est.Nhat$gcv.ubre), abund.est.gcv,tolerance=par.tol)
  abund.est.abund.est<-dsm(abundance.est~s(x,y), hn.model, segdata,
                           obsdata)
  expect_equal(unname(abund.est.abund.est$gcv.ubre),abund.est.gcv,
                                                  tolerance=par.tol)
  abund.est.abundance<-dsm(abundance.est~s(x,y), hn.model, segdata,
                       obsdata)
  expect_equal(unname(abund.est.abundance$gcv.ubre), abund.est.gcv,
                                                   tolerance=par.tol)


  ## models for density
  D.gcv <- 1.660703e-07
  density.D<-dsm(D~s(x,y), hn.model, segdata, obsdata,
                 weights=rep(1,nrow(segdata)))
  expect_equal(unname(density.D$gcv.ubre), D.gcv,tolerance=par.tol)

  density.density<-dsm(density~s(x,y), hn.model, segdata,
                       obsdata,
                       weights=rep(1,nrow(segdata)))
  expect_equal(unname(density.density$gcv.ubre), D.gcv,tolerance=par.tol)

  density.Dhat<-dsm(Dhat~s(x,y), hn.model, segdata,
                    obsdata,
                    weights=rep(1,nrow(segdata)))
  expect_equal(unname(density.Dhat$gcv.ubre), D.gcv,tolerance=par.tol)

  density.density.est<-dsm(density.est~s(x,y), hn.model, segdata,
                           obsdata,
                           weights=rep(1,nrow(segdata)))
  expect_equal(unname(density.density.est$gcv.ubre), D.gcv,tolerance=par.tol)

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
    expect_error(dsm(N~s(x,y), hn.model, seg, obs_missing),
                 paste0("Column(s) \"",mcov,
                        "\" not found in observation.data."),
                 fixed=TRUE)
  }

  for(mcov in c("Effort","Sample.Label")){
    seg_missing <- seg
    seg_missing[[mcov]] <- NULL
    expect_error(dsm(N~s(x,y), hn.model, seg_missing, obs),
                 paste0("Column(s) \"",mcov,
                        "\" not found in segment.data."),
                 fixed=TRUE)
  }

  # with segment area specified we only have a problem with Sample.Label
  seg_missing <- seg
  seg_missing[["Sample.Label"]] <- NULL
  expect_error(dsm(N~s(x,y), hn.model, seg_missing, obs),
               paste0("Column(s) \"Sample.Label\" not found in segment.data."),
               fixed=TRUE)

})



detach("mexdolphins")
