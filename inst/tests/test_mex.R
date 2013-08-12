library(dsm)
library(testthat)

lnl.tol<-1e-4
par.tol<-1e-6

context("Mexico pantropical dolphin data")

test_that("Do we get the same results?",{

  # load the Gulf of Mexico dolphin data
  data(mexdolphins)

  # fit a detection function and look at the summary
  hn.model <- ds(mexdolphins$distdata, max(mexdolphins$distdata$distance),
                 adjustment = NULL)

  ddf.par <- 8.626542
  names(hn.model$ddf$par) <- NULL
  expect_that(hn.model$ddf$par, equals(ddf.par,tolerance=par.tol))

  # fit a simple smooth of x and y
  mod1<-dsm(N~s(x,y), hn.model, mexdolphins$segdata, mexdolphins$obsdata)
  summary(mod1)


  expect_that(mod1$gcv.ubre, equals(42.98519,tolerance=par.tol))


  expect_that(dsm.cor(mod1,resid.type="d",max.lag=9),
              throws_error("No column called Segment.Label in data"))
})
