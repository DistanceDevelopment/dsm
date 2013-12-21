library(dsm)
library(Distance)
library(testthat)

par.tol<-1e-6

context("test inputs")

test_that("formula specs",{

  # load the Gulf of Mexico dolphin data
  data(mexdolphins)

  # fit a detection function and look at the summary
  hn.model <- ds(mexdolphins$distdata, max(mexdolphins$distdata$distance),
                 adjustment = NULL)

  ## models for count
  count.gcv <- 42.98519
  count.N<-dsm(N~s(x,y), hn.model, mexdolphins$segdata, mexdolphins$obsdata)
  expect_that(count.N$gcv.ubre, equals(count.gcv,tolerance=par.tol))
  count.n<-dsm(n~s(x,y), hn.model, mexdolphins$segdata, mexdolphins$obsdata)
  expect_that(count.n$gcv.ubre, equals(count.gcv,tolerance=par.tol))
  count.count<-dsm(count~s(x,y), hn.model, mexdolphins$segdata,
                   mexdolphins$obsdata)
  expect_that(count.count$gcv.ubre, equals(count.gcv,tolerance=par.tol))
  count.abundance<-dsm(abundance~s(x,y), hn.model, mexdolphins$segdata,
                       mexdolphins$obsdata)
  expect_that(count.abundance$gcv.ubre, equals(count.gcv,tolerance=par.tol))


  ## models for abund.est
  abund.est.gcv <- 57.40737
  abund.est.Nhat<-dsm(Nhat~s(x,y), hn.model, mexdolphins$segdata,
                      mexdolphins$obsdata)
  expect_that(abund.est.Nhat$gcv.ubre, equals(abund.est.gcv,tolerance=par.tol))
  abund.est.abund.est<-dsm(abundance.est~s(x,y), hn.model, mexdolphins$segdata,
                           mexdolphins$obsdata)
  expect_that(abund.est.abund.est$gcv.ubre,equals(abund.est.gcv,
                                                  tolerance=par.tol))
  abund.est.abundance<-dsm(abundance.est~s(x,y), hn.model, mexdolphins$segdata,
                       mexdolphins$obsdata)
  expect_that(abund.est.abundance$gcv.ubre, equals(abund.est.gcv,
                                                   tolerance=par.tol))


  ## models for density
  D.gcv <- 1.660703e-07
  density.D<-dsm(D~s(x,y), hn.model, mexdolphins$segdata, mexdolphins$obsdata)
  expect_that(density.D$gcv.ubre, equals(D.gcv,tolerance=par.tol))
  density.density<-dsm(density~s(x,y), hn.model, mexdolphins$segdata, mexdolphins$obsdata)
  expect_that(density.density$gcv.ubre, equals(D.gcv,tolerance=par.tol))
  density.Dhat<-dsm(Dhat~s(x,y), hn.model, mexdolphins$segdata,
                  mexdolphins$obsdata)
  expect_that(density.Dhat$gcv.ubre, equals(D.gcv,tolerance=par.tol))
  density.density.est<-dsm(density.est~s(x,y), hn.model, mexdolphins$segdata, mexdolphins$obsdata)
  expect_that(density.density.est$gcv.ubre, equals(D.gcv,tolerance=par.tol))


})
