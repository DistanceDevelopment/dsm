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
attach(mexdolphins)

# fit a detection function and look at the summary
hn.model <- suppressMessages(ds(distdata,
                                max(distdata$distance),
                                adjustment = NULL))

# fit a simple smooth of x and y
mod1<-dsm(N~s(x,y), hn.model, segdata, obsdata)

# run the moving block bootstrap for 2 rounds
set.seed(1123)
mod1.var <- dsm.var.prop(mod1, preddata, off.set=preddata$area)

test_that("mexdolphins - results for s(x,y)",{
  expect_that(summary(mod1.var)$cv,
              equals(0.2131066,tol=cv.tol))
})

## With no detection function
test_that("mexdolphins - works for NULL detection function",{
  mod1_nodf <-dsm(N~s(x,y), NULL, segdata, obsdata,
                  strip.width=8000)
  set.seed(1123)
  mod1.var <- dsm.var.gam(mod1_nodf, preddata, off.set=preddata$area)

  expect_that(summary(mod1.var)$cv,
              equals(0.1624528, tol=cv.tol))

  # throw an error if you want detection function uncertainty with no
  # detection function
  expect_error(dsm.var.prop(mod1_nodf, preddata, off.set=preddata$area),
               "No detection function in this analysis, use dsm.var.gam")

})
detach("mexdolphins")
