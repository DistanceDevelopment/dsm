#' Check column names exist
#'
#' Internal function to check that supplied `data.frames` have the correct
#' columns and checks that sample labels are all unique.
#' @param ddf.obj a `ddf` object from [`mrds`][mrds]
#' @param segment.data segment data as defined in [`dsm`][dsm]
#' @param observation.data observation data as defined in [`dsm`][dsm]
#' @param segment.area area of segments
#' @return nothing, but throws an error if something went wrong
#' @author David Lawrence Miller
check.cols <- function(ddf.obj, segment.data, observation.data, segment.area){

  ## check that the columns are there
  checks <- list(segment.data = c("Effort", "Sample.Label"),
                observation.data = c("object", "Sample.Label", "size",
                                    "distance"))

  # don't need to check distance if we don't have a detection function
  if(is.null(ddf.obj)){
    checks$observation.data <- checks$observation.data[checks$observation.data!=
                                                       "distance"]
  }

  if(!is.null(segment.area)){
    checks$segment.data <- checks$segment.data[checks$segment.data != "Effort"]
  }

  for(i in seq_len(length(checks))){
    check.res <- checks[[i]] %in% names(get(names(checks)[[i]]))
    if(any(!check.res)){

      stop(paste0("Column(s) \"",
                  paste(checks[[i]][!check.res],collapse="\", \""),
                  "\" not found in ", names(checks)[[i]],
                  ".\n  Check ?\"dsm-data\"."))
    }
  }

  ## check that Sample.Label is unique
  if(length(segment.data$Sample.Label)!=length(
                            unique(segment.data$Sample.Label))){
    warning("'Sample.Labels are non-unique in segment data!")
  }

  invisible()
}
