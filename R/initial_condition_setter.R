#' Initial Condition Setter
#'
#' Sets the initial conditions for the nlm_zip() and em_zip() functions, based on inputs.
#'
#' @param data A data set, where the observed count vector and time vector (if variable) are the last two columns.
#' @param n_leaves An integer specifying the number of distinct variable combinations in the data set.
#' @param p_0 A numeric vector of initial values for the risk probabilities. If the vector is of length 1, this value will be repeated for each leaf. If NULL, then an automatic method (initial_method) will be used.
#' @param l_0 A numeric vector of initial values for the rates. If the vector is of length 1, this value will be repeated for each leaf. If NULL, then an automatic method (initial_method) will be used.
#' @param initial_method A character string indicating the method for automatically setting the initial conditions. The character string can be an element of c("mean",mle","mm"), with "mean" being the default when p_0, l_0 are not provided.
#' @param variable_time A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE).
#'
#' @return A list containing vectors of initial values for the risk probability and rat, as well as a logical value for whether mean is being used as the initial condition method (TRUE) or not (FALSE).
#' @examples
#' initial_condition_setter(knee_pain_obs,8)
#'
initial_condition_setter<-function(data,n_leaves,p_0=NULL,l_0=NULL,initial_method=NULL,variable_time=TRUE){
  if(length(p_0) != length(l_0)){
    if(min(length(p_0),length(l_0))==0){
      stop("if one of p_0, l_0 is provided, the other must too")
    }else{
      warning("Length of vectors of initial values for parameters don't match.")
    }
  }

  if(max(length(p_0),length(l_0))>0 & length(initial_method)>0 ){
    stop("Initial condition method can only be provided when neither p_0, l_0 are")
  }

  if(max(length(p_0),length(l_0),length(initial_method))==0){
    initial_method<-"mean"
    initial_mean<-TRUE
  }else{
    initial_mean<-FALSE
  }

  if(length(initial_method)>1){
    stop("Please provide only one initial method")
  }

  if(length(initial_method)>0){
    if(!(initial_method %in% c("mean","mme","mle") )){
      stop("Unknown initial condition method chosen - Please select either mean, mme or mle")
    }
    if(initial_method=="mme"){
      if(variable_time){
      mme<-mme_variable_time_zip(data)
      }else{
      mme<-uniform_time_zip(data,"mme")
      }
      l0_vec<-mme$lambda
      p0_vec<-mme$prob
    }else if(initial_method=="mle"){
      if(variable_time){
        mle<-mle_variable_time_zip(data)
      }else{
        mle<-uniform_time_zip(data,"mle")
      }
      l0_vec<-mle$lambda
      p0_vec<-mle$prob
    }else{
      p0_vec<-c()
      l0_vec<-c()
      initial_mean<-TRUE
    }
  }

  if(length(p_0)==0){
    p_0<-p0_vec
  }else if(length(p_0)==1){
    p_0<-rep(p_0,n_leaves)
  }else if(length(p_0)!=n_leaves){
    stop("Initial value vector for risk probability doesn't match number of leaves.")
  }

  if(length(l_0)==0){
    l_0<-l0_vec
  }else if(length(l_0)==1){
    l_0<-rep(l_0,n_leaves)
  }else if(length(l_0)!=n_leaves){
    stop("Initial value vector for rates doesn't match number of leaves.")
  }

  if(min(c(l_0,p_0))<0){
    stop("Negative parameters estimated as initial condition - consider an alternative method.")
  }

  return(list(l=l_0,p=p_0,initial_mean = initial_mean))

}
