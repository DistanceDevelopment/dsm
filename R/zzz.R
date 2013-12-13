.onAttach<-function(...){
  # uses packageStartupMessage which can then be
  # surpressed
  version <- utils::packageVersion("dsm")

  hello <- return(paste0("This is dsm ",version,"\n"))
  packageStartupMessage(hello)
}
