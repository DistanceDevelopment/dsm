# test that we can use glms as models
library(dsm)
library(Distance)
library(testthat)

lnl.tol<-1e-4
par.tol<-1e-6

context("Do GLMs work")

# load the Gulf of Mexico dolphin data
data(mexdolphins)

# fit a detection function and look at the summary
hn.model <- suppressMessages(ds(distdata,
                                max(distdata$distance),
                                adjustment = NULL))

test_that("Do we get the same results?",{

  # fit a simple smooth of x and y
  mod1 <- dsm(count~x+y+depth, hn.model, segdata, obsdata, engine="glm")
  #summary(mod1)

  res_coef <- c(-1.445646175e+01, -6.034068174e-07,
                2.279085802e-06, 6.007334200e-04)
  names(res_coef) <- c("(Intercept)", "x", "y", "depth")
  expect_equal(coef(mod1), res_coef, tolerance=par.tol)


})

test_that("Density weighting",{

  # compare when we set the weights
  mod1 <- dsm(count~x+y, hn.model, segdata, obsdata, engine="glm")
  mod1.w <- dsm(density ~ x + y + depth, hn.model, segdata, obsdata,
                weights=mod1$data$segment.area, engine="glm")

  res_coef <- c(-1.445646196e+01, -6.034067589e-07,
                2.279085648e-06, 6.007334023e-04)
  names(res_coef) <- c("(Intercept)", "x", "y", "depth")
  expect_equal(coef(mod1.w), res_coef, tolerance=par.tol)

})
