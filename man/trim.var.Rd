\name{trim.var}
\alias{trim.var}
\title{Trimmed variance}
\usage{
trim.var(untrimmed.bootstraps, boxplot.coef = 1.5)
}
\arguments{
  \item{untrimmed.bootstraps}{(usually the
  \code{$study.area.total} element of a returned \code{dsm}
  bootstrap object.}

  \item{boxplot.coef}{the value of \code{coef} used to
  calculate the outliers see \code{\link{boxplot}}.}
}
\value{
trimmed variance
}
\description{
Trim the variance estimates from the bootstrap. This is
defined as the percentage defined as amount necessary to
bring median and trimmed mean within 8% of each other these
are defined as 'outliers'.
}
\details{
Code originally by Louise Burt.
}
\author{
Louise Burt
}

