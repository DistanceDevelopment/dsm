.onAttach<-function(library,pkgname)
{
  # this now conforms with new R conventions
  # uses packageStartupMessage which can then be
  # surpressed
  dsm.version()
}

# taken from mgcv
dsm.version <- function(){ 
  library(help=dsm)$info[[1]] -> linfo

  version <- linfo[pmatch("Version",linfo)]
  built <- linfo[pmatch("Built",linfo)]

  um <- strsplit(version," ")[[1]]
  version <- um[nchar(um)>0][2]

  um <- strsplit(built,"R")[[1]]
  built <- um[nchar(um)>0][2]

  hello <- paste("This is dsm ",version,"\nBuilt: R",built,sep="")
  packageStartupMessage(hello)
}

