\name{dsm.bootstrap.pre.post}
\alias{dsm.bootstrap.pre.post}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Pre- and post-processing of moving blook boostrap replicates }
\description{
  Organize data to preform moving block bootstrap analysis.  When replicates are
  finished, remove outliers induced by rogue replicates.  Then print results
  summarizing number of ill-formed replicates, histogram of replicate values (after
  outlier removal), and combine uncertainty of density-surface model with 
  uncertainty of detection function fitting.
}
\usage{
dsm.bootstrap.pre.post(ddf.ID = ddf.ID, dsm.ID = dsm.ID, prediction.ID = prediction.ID,
                      method = "movblock", replicates = replicates, blocksize = blocksize, 
			    detectionrefit = FALSE, alpha = 0.95, fname = fname, 
			    field.flag = FALSE, cell.area = NULL, path = path, 
			    bpfile = bpfile, boxplot.coef = 1.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ddf.ID}{ Analysis number containing the fitted detection function.  If NULL, then detection probabilities were
  			assigned rather than computed, hence variability associated with fitting the detection function cannot
			  be incorporated into the overall measures of precision. }
  \item{dsm.ID}{ Analysis number containing the fitted density surface model. }
  \item{prediction.ID}{ Analysis number containing the predicted abundance for which precision estimate is required. }
  \item{method}{ Only legitimate variance computation method is parametric moving block }
  \item{replicates}{ Number of replicates for the bootstrap }
  \item{blocksize}{ Size (in segments) of each block for moving block bootstrap }
  \item{detectionrefit}{ Boolean indicating whether detection function should also be refitted for each bootstrap (disabled at present) }
  \item{alpha}{ 1-Type I error rate for confidence interval construction }
  \item{fname}{ Temporary file to contain matrix of cell-specific measures of precision }
  \item{field.flag}{Boolean indicating whether cell size is constant for all cells (FALSE) or
                if a field containing cell-specific size is provided (TRUE) }
  \item{cell.area}{ If field=TRUE, the name of the field containing cell size }
  \item{path}{ Path to temporary directory where  }
  \item{bpfile}{ Temporary filename where progress of replicates is written }
  \item{boxplot.coef}{ Outlier multiplier.  The larger the value, the fewer values are considered outliers. It measures the muliplier of the interquartile range added to the upper quartile (75\%) to establish the threshold for an outlier.  If you wish there to be no outliers excluded from bootstrap variance estimation, set this value to some large number. }
}
\details{
  Variances are also computed for each polygon of the small area estimation requested.
}
\value{
  List consisting of:
  \item{lower.abund.ci }{Lower confidence bound using percentile method, after trimming}
  \item{upper.abund.ci }{Upper confidence bound using percentile method, after trimming}
  \item{untrimmed }{Vector of study region total abundances, untrimmed}
  \item{raw.bootstrap.values }{Vector of study region total abundances, after trimming}
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Eric Rexstad, RUWPA, \email{ericr@mcs.st-and.ac.uk} }
%\note{ ~~further notes~~ 
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{  }

