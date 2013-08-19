#' Check column names exist
#'
#' Internal function to check that supplied `data.frames` have the correct columns.
#' @param ddf.obj a ddf object from `mrds`
#' @param segment.data segment data as defined in \code{\link{dsm}}
#' @param observation.data observation data as defined in \code{\link{dsm}}
#' @param strip.width strip width if strip transects are being used
#' @param segment.area area of segments
#' @return nothing, but throws an error if something went wrong
#' @author David Lawrence Miller

check.cols <- function(ddf.obj, segment.data, observation.data, strip.width,
                       segment.area){

  checks <-list(segment.data = c("Effort","Sample.Label"),
                observation.data = c("object","Sample.Label","size","distance"))

  for(i in 1:length(checks)){
    check.res <- checks[[i]] %in% names(get(names(checks)[[i]]))
    if(!any(check.res)){

      stop(paste0("Column(s) \"",
                  paste(checks[[i]][!check.res],collapse="\", \""),
                  "\" not found in ", names(checks)[[i]],
                  ".\n  Check ?\"dsm-data\"."))

    }
  }

  invisible()
}
