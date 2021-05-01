dsm_env <- new.env(parent = emptyenv())

.onAttach <- function(...){

  # uses packageStartupMessage which can then be
  # suppressed
  version <- utils::packageVersion("dsm")
  built <- utils::packageDescription("dsm",fields="Built")
  hello <- paste0("This is dsm ",version,"\nBuilt: ",built)
  packageStartupMessage(hello)

  # Check the MGCV version. If we do it here, there's no need to do it on every
  # dsm() call.
  mgcv.version <- as.numeric(strsplit(as.character(packageVersion("mgcv")),
                                      "\\.")[[1]])
  if(mgcv.version[1]<1 | (mgcv.version[2]<7 |
     (mgcv.version[2]==7 & mgcv.version[3]<24))){
    old_version <- TRUE
  } else {
    old_version <- FALSE
  }
  assign("old_mgcv", old_version, envir = dsm_env)
}
