\name{block.info.per.su}
\alias{block.info.per.su}
\title{Find the block information}
\usage{
  block.info.per.su(block.size, data, name.su)
}
\arguments{
  \item{block.size}{number of segments per block}

  \item{data}{data used to build the model}

  \item{name.su}{names of the (unique) sampling units (ie.
  transects)}
}
\value{
  a \code{data.frame} with the following columns
  \tabular{ll}{name \tab the sample unit name (e.g.
  transect label) \cr num.seg \tab number of segments in
  that transect \cr num.block \tab number of blocks
  available\cr start.block \tab block # for first block\cr
  end.block \tab block # for last block\cr num.req \tab
  number of blocks needed for the unit\cr }
}
\description{
  Takes the transect data and works out how many blocks of
  a given size (in segment terms) fit into each.
}

