library(dsm)
library(Distance)
library(testthat)

cv.tol<-1e-5
N.tol<-1e-4


## NB the 444km^2 for the prediction grid is INCORRECT but
## serves us fine for the purpose of these tests.

context("Moving block bootstrap")

# load the Gulf of Mexico dolphin data
data(mexdolphins)

# fit a detection function and look at the summary
hn.model <- suppressMessages(ds(distdata,
                                max(distdata$distance),
                                adjustment = NULL))

# fit a simple smooth of x and y
mod1<-dsm(count~s(x,y), hn.model, segdata, obsdata)

# run the moving block bootstrap for 2 rounds
set.seed(1123)
mod1.movblk <- dsm.var.movblk(mod1, preddata, n.boot = 2,
                              block.size = 3, off.set = 444*1000*1000,bar=FALSE)

test_that("mexdolphins - bootstrap results for s(x,y)",{

  expect_equal(mod1.movblk$study.area.total,
              c(36787.71, 19407.48), tol=N.tol)

  expect_equal(summary(mod1.movblk)$cv[1],
               0.3574841, tol=cv.tol)

  expect_equal(summary(mod1.movblk)$bootstrap.cv[1],
               0.3340702, tol=cv.tol)
})

## With no detection function
test_that("mexdolphins - bootstrap works for strip transects",{
  fdf <- dummy_ddf(obsdata$object, obsdata$size, 8000)
  mod1_nodf <-dsm(count~s(x,y), fdf, segdata, obsdata)
  set.seed(1123)
  mod1.movblk_nodf <- dsm.var.movblk(mod1_nodf, preddata, n.boot=2,
                              block.size=3, off.set=444*1000*1000, bar=FALSE)

  expect_equal(summary(mod1.movblk_nodf)$cv[1],
               0.3340702, tol=cv.tol)

  # throw an error if you want detection function uncertainty with no
  # detection function
  expect_error(dsm.var.movblk(mod1_nodf, preddata, n.boot=2,
                              block.size=3, off.set=444*1000*1000, bar=FALSE,
                              ds.uncertainty=TRUE),
               "Cannot incorporate detection function uncertainty with no detection function!")

})
