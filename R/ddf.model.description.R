# builds the description of the model to be printed
# only used internally
ddf.model.description <- function(model){

  if(any(class(model)=="fake_ddf")){
    return("No detection function")
  }

  if(model$method=="io"){
    mod.str <- "MR model: Independent observer, DS model:"
    model <- model$ds
  }else{
    mod.str <- ""
  }


  key <- switch(model$ds$aux$ddfobj$type,
                hn = "Half-normal",
                hr = "Hazard-rate")
  mod.str <- paste(mod.str, key, "key function")
  if(!is.null(model$ds$aux$ddfobj$adjustment)){
    adj.series <- switch(model$ds$aux$ddfobj$adjustment$series,
                        cos  = "cosine",
                        herm = "Hermite polynomial",
                        poly = "simple polynomial")
    mod.str <- paste(mod.str, "with", adj.series, "adjustment term")

    adj.order <- model$ds$aux$ddfobj$adjustment$order
    if(length(adj.order)>1){
      mod.str <- paste(mod.str, "s", sep="")
    }
    mod.str <- paste(mod.str, "of order", paste(adj.order, collapse=","))
  }

  return(mod.str)
}
