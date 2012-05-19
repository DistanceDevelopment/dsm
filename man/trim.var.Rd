\name{trim.var}
\alias{trim.var}
\title{Trimmed variance}
\usage{
  trim.var(untrimmed.bootstraps)
}
\arguments{
  \item{untrimmed.bootstraps}{(usually the
  \code{$study.area.total} element of a returned \code{dsm}
  bootstrap object.}
}
\value{
  trimmed variance
}
\description{
  Trim the variance estimates from the bootstrap. This is
  defined as the percentage defined as amount necessary to
  bring median and trimmed mean within 8% of each other
  these are defined as 'outliers'.
}
\details{
  Code originally by Louise Burt.
}
\author{
  Louise Burt
}

