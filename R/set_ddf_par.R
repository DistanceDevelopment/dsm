# set the parameter values for model to be par
set_ddf_par <- function(par, model){
  if(model$method %in% c("io", "trial")){
    model$mr$mr$coefficients <- par[seq_along(model$mr$mr$coefficients)]
    model$ds$par <- par[(length(model$mr$mr$coefficients)+1):length(par)]
  }else{
    model$par <- par
  }

  return(model)
}

