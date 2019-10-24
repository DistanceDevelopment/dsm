library(dsm)
library(Distance)
library(testthat)

lnl.tol<-1e-4
par.tol<-1e-6

context("prediction")

# load the Gulf of Mexico dolphin data
data(mexdolphins)

# fit a detection function and look at the summary
hn.model <- suppressMessages(ds(distdata,
                                max(distdata$distance),
                                adjustment = NULL))

# check predictions when density is the response
test_that("predictions from density",{

  # fit a simple smooth of x and y
  mod1 <- dsm(density~s(depth), hn.model, segdata, obsdata)

  fake_dat <- mod1$data
  fake_dat$off.set <- NULL
  expect_equal(predict(mod1), predict(mod1, fake_dat, off.set=1))

  # check for error if we don't supply newdata but do supply off.set
  expect_warning(predict(mod1, off.set=1),
               "Ignoring supplied off.set as newdata was not supplied")


  # check you get different answers from different offsets
  expect_equal(2*predict(mod1, fake_dat, off.set=2),
               predict(mod1, fake_dat, off.set=4))
})


# and for count
test_that("predictions for count",{

  # fit a simple smooth of x and y
  mod1 <- dsm(count~s(depth), hn.model, segdata, obsdata)

  fake_dat <- mod1$data
  fake_dat$off.set <- NULL
  expect_equal(predict(mod1), predict(mod1, fake_dat, off.set=fake_dat$segment.area))

})


