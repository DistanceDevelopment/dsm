\name{summary.dsm}
\alias{summary.dsm}
\title{Summary for a fitted dsm model}
\usage{
  \method{summary}{dsm} (object, ...)
}
\arguments{
  \item{object}{a \code{dsm} object}

  \item{\dots}{unused arguments for S3 compatibility}
}
\value{
  a summary object
}
\description{
  Gives a brief summary of a fitted \code{dsm} object. For
  more specific information on the detection function and
  spatial model, call \code{summary()} on the \code{$ddf}
  and \code{$result} elements of the \code{dsm} object.
}
\author{
  David L. Miller
}

