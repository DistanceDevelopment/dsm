\name{dsm}
\alias{dsm}
\title{Fit a density surface model to segment-specific estimates of abundance
or density.}
\usage{
  dsm(formula, ddf.obj, segment.data, observation.data,
    engine = "gam", convert.units = 1,
    family = quasipoisson(link = "log"), group = FALSE,
    gamma = 1.4, control = list(keepData = TRUE),
    availability = 1, ...)
}
\arguments{
  \item{formula}{formula for the surface. This should be a
  valid
  \code{\link{glm}}/\code{\link{gam}}/\code{\link{gamm}}
  formula.}

  \item{ddf.obj}{result from call to \code{\link{ddf}} or
  \code{\link{ds}}.}

  \item{segment.data}{segment data, see
  \code{\link{dsm-data}}.}

  \item{observation.data}{observation data, see
  \code{\link{dsm-data}}.}

  \item{engine}{which model should be used for the DSM
  (\code{\link{glm}}/
  \code{\link{gam}}/code{\link{gamm}}).}

  \item{convert.units}{value to alter length to width for
  calculation of the offset.}

  \item{family}{response distribution (popular choices
  include \code{\link{quasipoisson}}, \code{\link{Tweedie}}
  and \code{\link{negbin}}. Defaults to
  \code{quasipossion}.}

  \item{\dots}{anything else to be passed straight to
  \code{\link{gam}}.}

  \item{group}{should group abundance/density be modelled
  rather than individual abundance/density? This
  effectively sets the \code{size} column in
  \code{observation.data} to be 1.}

  \item{control}{the usual \code{control} argument for a
  \code{gam}, \code{keepData} must be \code{TRUE} or
  variance estimation will not work.}

  \item{availability}{an availability bias used to scale
  the counts/estimated counts by. If we have \code{N}
  animals in a segment, then \code{N/availability} will be
  entered into the model. Uncertainty in the availability
  is not handled at present.}

  \item{gamma}{parameter to \code{gam()} set to a value of
  1.4 (from advice in Wood (2006)) such that the
  \code{gam()} is inclined to not 'overfit.'.}
}
\value{
  a \code{\link{glm}}/\code{\link{gam}}/\code{\link{gamm}}
  object, with an additional element, \code{ddf} which
  holds the detection function object.
}
\description{
  Given a detection function analysis, construct a density
  surface model (DSM) based on environmental covariates.
}
\details{
  The response can be one of the following: \tabular{ll}{
  \code{N}, \code{abundance} \tab count in each segment\cr
  \code{Nhat}, \code{abundance.est} \tab estimated
  abundance per segment, estimation is via a
  Horvitz-Thompson estimator. This should be used when
  there are covariates in the detection function.\cr
  \code{presence} \tab interpret the data as
  presence/absence (reember to change the \code{family}
  argument to \code{binomial()}\cr \code{D}, \code{density}
  \tab density per segment\cr }
}
\examples{
library(Distance)
library(dsm)

# load the Gulf of Mexico dolphin data (see ?mexdolphins)
data(mexdolphins)

# fit a detection function and look at the summary
hr.model <- ds(mexdolphins$distdata, max(mexdolphins$distdata$distance), key = "hr", adjustment = NULL)
summary(hr.model)

# fit a simple smooth of x and y
mod1<-dsm(N~s(x,y), hr.model, mexdolphins$segdata, mexdolphins$obsdata)
summary(mod1)

# create an offset (in metres)
# each prediction cell is 444km2
off.set <- 444*1000*1000

# predict over a grid
mod1.pred <- predict(mod1, mexdolphins$preddata, off.set)

# calculate the predicted abundance over the grid
sum(mod1.pred)

# plot the smooth
plot(mod1)
}
\author{
  David L. Miller
}
\references{
  Hedley, S. and S. T. Buckland. 2004. Spatial models for
  line transect sampling. JABES 9:181-199.

  Wood, S.N. 2006. Generalized Additive Models: An
  Introduction with R. CRC/Chapman & Hall.
}

