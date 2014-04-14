library(dsm)
library(Distance)
library(testthat)

lnl.tol<-1e-4
par.tol<-1e-6


set.seed(1123)

context("Moving block bootstrap")

# load the Gulf of Mexico dolphin data
data(mexdolphins)

# fit a detection function and look at the summary
hn.model <- suppressMessages(ds(mexdolphins$distdata,
                                max(mexdolphins$distdata$distance),
                                adjustment = NULL))

# fit a simple smooth of x and y
mod1<-dsm(N~s(x,y), hn.model, mexdolphins$segdata, mexdolphins$obsdata)

# run the moving block bootstrap for 2 rounds
mod1.movblk <- dsm.var.movblk(mod1, mexdolphins$preddata, n.boot = 2,
                              block.size = 3, off.set = 444*1000*1000,bar=FALSE)

test_that("mexdolphins - bootstrap results for s(x,y)",{

  expect_that(mod1.movblk$study.area.total,
              equals(c(37141.16, 20659.74),tol=par.tol))

  expect_that(summary(mod1.movblk)$cv[1],
              equals(0.338609,tol=par.tol))

  expect_that(summary(mod1.movblk)$bootstrap.cv[1],
              equals(0.3137792,tol=par.tol))
})

