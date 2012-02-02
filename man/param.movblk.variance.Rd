\name{param.movblk.variance}
\alias{param.movblk.variance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Parametric moving-block variance estimator }
\description{
  Function reconstructs the density surface model used to fit the original response variable.
  A sample is drawn using a moving block approach to allocation residuals, hence retaining
  some spatial correlation structure.  Prediction is performed after the fitting of the original
  model to data that have been permuted via the reallocation of residuals.  The resulting
  estimated response in each cell is stored in a matrix, and vector of replicate total abundances
  is stored for summary when all replicates are processed.
}
\usage{
param.movblk.variance(B = num.boot, dsm.object = dsm.object, 
                      pred.object = predict.object, samp.unit.name = "Transect",
			    block = num.segs.in.block, field = FALSE, 
			    cell.name = cell, bpfile = bpfile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{B}{ Number of bootstrap replicates }
  \item{dsm.object}{ Object containing the density surface model fitted }
  \item{pred.object}{ Object containing the name of the prediction grid }
  \item{samp.unit.name}{ Scale of sampling within which moving blocks will be rearranged }
  \item{block}{ Number of segments in each block }
  \item{field}{ Boolean indicating whether cell size is constant for all cells (FALSE) or
                if a field containing cell-specific size is provided (TRUE) }
  \item{cell.name}{ If field=TRUE, the name of the field containing cell size }
  \item{bpfile}{ Bootstrap progress file to update Distance as to progress }
}
\details{
  Note that memory is allocated (using \code{memory.limit}) to the full amount
  available on a Windows machine (or 2Gb) whichever is smaller.
}
\value{
  Matrix of cell-specific summary statistics is returned for producing measures of
  uncertainty at the cell level.  Also the total abundance for the entire study region
  is added to a vector of replicate estimates.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Louise Burt, RUWPA, \email{louise@mcs.st-and.ac.uk} \cr and Eric Rexstad \email{ericr@mcs.st-and.ac.uk} }
\note{  
To prevent storing the bootstrap replicate values (perhaps thousands) for each cell in the
prediction grid (perhaps many thousands), this code makes use of the 'machine formula' for
computing variance, by storing summary statistics (sum and sum of squares) for each cell of 
the grid.

The problem with this approach is that rogue values at the prediction cell level cannot be
filtered for the final bootstrap variance estimate as they are filtered for the overall study
region (\code{study.area.total}) estimate of variance.  The result is that cell-level 
measures of precision are quite large because they are influenced by occasional 'rogue' values.
}
\seealso{  \code{\link{dsm.bootstrap.pre.post}} }
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

