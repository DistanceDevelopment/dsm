lnl.tol<-1e-4
par.tol<-1e-6

context("Mexico pantropical dolphin data")

library(Distance)

# load the Gulf of Mexico dolphin data
data(mexdolphins)

# fit a detection function and look at the summary
hn.model <- suppressMessages(ds(distdata,
                                max(distdata$distance),
                                adjustment = NULL))

test_that("Do we get the same results?",{


  ddf.par <- 8.626542
  names(hn.model$ddf$par) <- NULL
  expect_equal(hn.model$ddf$par, ddf.par, tolerance=par.tol)

  # fit a simple smooth of x and y
  mod1 <- dsm(count~s(x, y), hn.model, segdata, obsdata)
  #summary(mod1)


  expect_equal(unname(mod1$gcv.ubre), 936.0362722, tolerance=par.tol, check.attributes=FALSE)


  expect_error(dsm_cor(mod1, resid.type="d", max.lag=9),
               "No column called Segment.Label in data")

  # predict(model) shoudld be the same as fitted(model)
  # DLM: is this true though??
  #expect_equal(as.vector(predict(mod1)), as.vector(fitted(mod1)), tolerance=par.tol)
})

test_that("Density weighting",{

  ## compare when we set the weights
  #mod1.w <- dsm(D~s(x,y), hn.model, segdata, obsdata,
  #              weights=mod1$data$segment.area/sum(mod1$data$segment.area))
  #expect_equal(fitted(mod1),fitted(mod1.w),tolerance=par.tol)

  # setting weights to 1 or another constant
  # compare when we set the weights
  mod1.w1 <- dsm(density.est~s(x,y), hn.model, segdata, obsdata,
                 weights=rep(1,nrow(segdata)))
  # compare when we set the weights
  mod1.w2 <- dsm(density.est~s(x,y), hn.model, segdata, obsdata,
                 weights=rep(100,nrow(segdata)))

  expect_equal(fitted(mod1.w1),fitted(mod1.w2),tolerance=par.tol)

  # scalar input of weights (same as weighting all as 1, or 10)
  mod1.ws1 <- dsm(density.est~s(x,y), hn.model, segdata, obsdata,
                weights=1)

  expect_equal(fitted(mod1.ws1),fitted(mod1.w2),tolerance=par.tol)

})



test_that("Density predictions",{
  # test predictions
  mod1 <- dsm(density.est~s(x,y), hn.model, segdata, obsdata, family=gaussian())

  # off.set=1 should be the same as newdata$off.set=1 and not supplying offset
  preddata$off.set <- 1
  expect_equal(predict(mod1, preddata), predict(mod1, preddata, offset=1))

})
