\name{rqgam.check}
\alias{rqgam.check}
\title{Randomised quantile residuals check plots for GAMs/DSMs}
\usage{
rqgam.check(gam.obj, ...)
}
\arguments{
  \item{gam.obj}{a \code{gam}, \code{glm} or \code{dsm}
  object.}

  \item{...}{arguments passed on to all plotting function}
}
\value{
just plots!
}
\description{
Function operates as \code{\link{gam.check}} but using
randomised quantile residuals, a la Dunn and Smyth (1996).
Checks of \code{k} are not computed, these need to be done
using \code{\link{gam.check}}.
}
\examples{
library(Distance)
library(dsm)

# load the Gulf of Mexico dolphin data (see ?mexdolphins)
data(mexdolphins)

# fit a detection function and look at the summary
hr.model <- ds(mexdolphins$distdata, max(mexdolphins$distdata$distance),
               key = "hr", adjustment = NULL)

# fit a simple smooth of x and y
mod1<-dsm(N~s(x,y), hr.model, mexdolphins$segdata, mexdolphins$obsdata)
rqgam.check(mod1)
}
\author{
Based on code provided by Natalie Kelly, bugs added by Dave
Miller
}

