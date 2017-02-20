library(dsm)
library(Distance)
library(testthat)

cv.tol<-1e-5
N.tol<-1e-4

## NB the 444km^2 for the prediction grid is INCORRECT but
## serves us fine for the purpose of these tests.
context("GAM variance")

# load the Gulf of Mexico dolphin data
data(mexdolphins)

# fit a detection function and look at the summary
hn.model <- suppressMessages(ds(distdata,
                                max(distdata$distance),
                                adjustment = NULL))

# fit a simple smooth of x and y
mod1 <- dsm(N~s(x,y), hn.model, segdata, obsdata)


# run the moving block bootstrap for 2 rounds
set.seed(1123)
mod1.var <- dsm.var.prop(mod1, preddata, off.set=preddata$area)

test_that("mexdolphins - results for s(x,y)",{
  # CV
  expect_equal(summary(mod1.var)$cv,
               0.1742568601, tol=cv.tol)
  # var
  expect_equal(mod1.var$pred.var,
               15336119.57, tol=N.tol)
  # test that the CIs are right
  expect_output(print(summary(mod1.var)),
                "2.5%     Mean    97.5% \\n16012.09 22473.35 31541.89")
})

test_that("different CIs work",{
  expect_output(print(summary(mod1.var, alpha=0.1)),
                "5%     Mean      95% \\n16908.97 22473.35 29868.86")
  expect_output(print(summary(mod1.var, alpha=0.02)),
                "1%     Mean      99% \\n15028.91 22473.35 33605.33")

})

## With no detection function
test_that("mexdolphins - works for NULL detection function",{
  mod1_nodf <-dsm(N~s(x,y), NULL, segdata, obsdata,
                  strip.width=8000)
  set.seed(1123)
  mod1.var <- dsm.var.gam(mod1_nodf, preddata, off.set=preddata$area)

  expect_equal(summary(mod1.var)$cv,
               0.1639757, tol=cv.tol)

  # throw an error if you want detection function uncertainty with no
  # detection function
  expect_error(dsm.var.prop(mod1_nodf, preddata, off.set=preddata$area),
               "No detection function in this analysis, use dsm.var.gam")

})

test_that("varprop doesn't work for estimated abundance", {

  mod1_Nhat <- dsm(abundance.est~s(x,y), hn.model, segdata, obsdata)
  expect_error(dsm.var.prop(mod1_Nhat, preddata, off.set=preddata$area),
               "Variance propagation can only be used with count as the response.")
})

#set.seed(1123)
#preddata$off.set <- preddata$area
#mod1.var <- dsm_varprop(mod1, preddata)
#
#test_that("mexdolphins - results for s(x,y)",{
#  # CV
#  expect_equal(summary(mod1.var)$cv,
#               0.2168384, tol=cv.tol)
#  # var
#  expect_equal(mod1.var$pred.var,
#               23747034.5355, tol=N.tol)
#  # test that the CIs are right
#  expect_output(print(summary(mod1.var)),
#                "2.5%     Mean    97.5% \\n14764.25 22473.39 34207.84")
#})
