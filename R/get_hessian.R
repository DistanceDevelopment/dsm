# extract the hessian from a ddf object
# hessian needs to be 2nd REAL Hessian from the optimisation, not
# "distance sampling" version (product of 1st derivs) that lives in $hessian
get_hessian <- function(this_ddf){

  hess <- matrix(0, length(this_ddf$par), length(this_ddf$par))

  first_ind <- 1

  # get the io hessian first if it exists
  if(this_ddf$method %in% c("io", "trial")){
    mr_hess <- this_ddf$mr$hessian
    hess[1:nrow(mr_hess), 1:ncol(mr_hess)] <- mr_hess
    first_ind <- 1 + nrow(mr_hess)
    ds_bit <- this_ddf$ds
  }else{
    ds_bit <- this_ddf
  }


  # extract the ds hessian
  opt_details <- attr(ds_bit$ds, "details")
  if(is.matrix(opt_details)){
    this_hess <- opt_details[nrow(opt_details), ]$nhatend
  }else{
    this_hess <- opt_details$nhatend
  }

  if(any(is.na(this_hess))){
    # fall back to DS use if things are bad
    return(this_ddf$hessian)
  }

  hess[first_ind:(first_ind+nrow(this_hess)-1),
       first_ind:(first_ind+nrow(this_hess)-1)] <- this_hess

  return(hess)
}
