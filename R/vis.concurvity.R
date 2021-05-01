#' Visualise concurvity between terms in a GAM
#'
#' Plot measures of how much one term in the model could be explained by
#' another. When values are high, one should consider re-running variable
#' selection with one of the offending variables removed to check for stability
#' in term selection.
#'
#' These methods are considered somewhat experimental at this time. Consult
#' [`concurvity`][mgcv::concurvity] for more information on how concurvity
#' measures are calculated.
#'
#' @param model fitted model
#' @param type concurvity measure to plot, see [`concurvity`][mgcv::concurvity]
#' @author David L Miller
#' @export
#' @importFrom graphics image layout
#' @examples
#' \dontrun{
#' library(Distance)
#' library(dsm)
#'
#' # load the Gulf of Mexico dolphin data (see ?mexdolphins)
#' data(mexdolphins)
#'
#' # fit a detection function and look at the summary
#' hr.model <- ds(distdata, max(distdata$distance),
#'                key = "hr", adjustment = NULL)
#'
#' # fit a simple smooth of x and y to counts
#' mod1 <- dsm(count~s(x,y)+s(depth), hr.model, segdata, obsdata)
#'
#' # visualise concurvity using the "estimate" metric
#' vis.concurvity(mod1)
#'}
vis.concurvity <- function(model, type="estimate"){

  # calculate concurvity for this model
  cc <- concurvity(model, full=FALSE)[[type]]

  # remove diagonal elements
  diag(cc) <- NA

  # setup plotting
  layout(matrix(1:2, ncol=2), widths=c(5,1))
  opar <- par(mar=c(5, 6, 5, 0) + 0.1)

  # main plot
  image(z=cc, x=1:ncol(cc), y=1:nrow(cc), ylab="", xlab="",
        axes=FALSE, asp=1, zlim=c(0,1))
  axis(1, at=1:ncol(cc), labels = colnames(cc), las=2, lwd=0)
  axis(2, at=1:nrow(cc), labels = rownames(cc), las=2, lwd=0)

  # legend
  opar <- par(mar=c(5, 0, 4, 3) + 0.1)
  image(t(matrix(rep(seq(0, 1, len=100), 2), ncol=2)),
        x=1:3, y=1:101, zlim=c(0,1), axes=FALSE, xlab="", ylab="")
  axis(4, at=seq(1,101,len=5), labels = round(seq(0,1,len=5),1), las=2)

  # reset graphical pars
  par(opar)
}
