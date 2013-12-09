\name{generate.mb.sample}
\alias{generate.mb.sample}
\title{Generate a vector of residuals to be mapped back onto the data}
\usage{
generate.mb.sample(num.blocks.required, block.size, which.blocks, dsm.data,
  unit.info, n.units)
}
\arguments{
  \item{num.blocks.required}{number of blocks that we
  need.}

  \item{block.size}{number of segments per block.}

  \item{which.blocks}{which blocks should be sampled.}

  \item{dsm.data}{the \code{$data} element of the result of
  a call to \code{dsm.fit}.}

  \item{unit.info}{result of calling
  \code{\link{block.info.per.su}}.}

  \item{n.units}{number of sampling units.}
}
\value{
vector of log-residuals
}
\description{
Not usually used on its own, called from within
\code{\link{dsm.var.movblk}}.
}

