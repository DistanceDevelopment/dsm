# extract all variables in a detection function
# works like all.vars but for ddf objects
all_df_vars <- function(this_ddf){

  if("fake_ddf" %in% class(this_ddf)){
    df_vars <- NULL
  }else if(this_ddf$method == "io"){
    df_vars <- c(all.vars(as.formula(this_ddf$ds$ds$aux$ddfobj$scale$formula)),
                 all.vars(as.formula(this_ddf$mr$model)))
    df_vars <- setdiff(df_vars, "distance")
  }else{
    # extract formulae
    df_formula <- as.formula(this_ddf$ds$aux$ddfobj$scale$formula)
    if(!is.null(this_ddf$ds$aux$ddfobj$shape$formula) &&
       this_ddf$ds$aux$ddfobj$shape$formula != "~1"){
      stop("Shape parameter formulae are not supported!")
    }

    # extract detection function variables
    df_vars <- all.vars(df_formula)
  }

  return(df_vars)
}
