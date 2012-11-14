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
  \code{\link{glm}}/\code{\link{gam}}/code{\link{gamm}}
  formula.}

  \item{ddf.obj}{result from call to \code{\link{ddf}} or
  \code{\link{ds}}.  If \code{ddf.obj} is \code{NULL} then
  strip transects are assumed.}

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
  and \code{\link{negbin}}.}

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
  1.4 (from advice in}
}
\value{
  a \code{\link{glm}}/\code{\link{gam}}/code{\link{gamm}}
  object, with an additional element, \code{ddf} which
  holds the detection function object.
}
\description{
  Given a detection function analysis, construct a density
  surface model (DSM) based on environmental covariates.
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

