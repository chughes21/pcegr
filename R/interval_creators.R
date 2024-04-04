#' The Gamma Highest Density Interval Creator
#'
#' This function calculates the highest density interval for the parameters of a Gamma posterior, as with a Poisson response variable.
#'
#' This function calculates the highest density interval for the estimated rates of the leaf stages in a PCEG or ZIPCEG.
#'
#' @param data A data set where the observed response vector and time vector (if applicable and variable) are the last two columns.
#' @param mod A StagedTree model fit to the data set, as produced by [pceg()].
#' @param ci A numeric value between 0 and 1 specifying the percentage confidence for the highest density interval.
#' @param N An integer specifying the number of iterations for the empirical functions used.
#' @param zip A logical value indicating whether the model specified is zero-inflated (TRUE) or not (FALSE).
#'
#' @return A series of highest density intervals for each parameter.
#' @export
#'
#' @examples
#' mod<-pceg(knee_pain_obs,2,TRUE,TRUE)
#' hdi_gamma_extractor(knee_pain_obs,mod,zip=FALSE)
hdi_gamma_extractor<-function(data,mod,ci=0.95,N=10000,zip=TRUE){

  variable_time <- mod$event.tree$variable.time

  n<-dim(data)[2]-1*variable_time+1*zip

  posterior<-mod$prior.distribution[[n]] + mod$data.summary[[n]]
  posterior<-posterior[which(!is.na(posterior[,1])),]
  shapes<-posterior[,1]
  scales<-1/(posterior[,2])
  m=length(shapes)
  output<-list()
  for(i in 1:m){
    marg_post<-distribution_gamma(N,shape=shapes[i],scale=scales[i])
    output[[i]]<-hdi(marg_post,ci=ci)
  }
  return(output)
}

#' The Beta Highest Density Interval Creator
#'
#' This function calculates the highest density interval for the parameters of a Dirichlet posterior, using a Beta distribution.
#' The Dirichlet posterior is present for the situations with a categorical response variable in any type of CEG.
#'
#' This function calculates the highest density interval for the estimated transition probabilities for the situation stages in a PCEG, ZIPCEG or CEG.
#'
#' @param data A data set, where the observed count vector and time vector (if variable) are the last two columns
#' @param mod A StagedTree model fit to the data set, as produced by [pceg()].
#' @param ci A numeric value between 0 and 1 specifying the perecentage confidence for the highest density interval.
#' @param N An integer specifying the number of iterations for the empirical functions used.
#' @param level_rel_final A non-positive integer indicating where the desired variable is relative to the response variable.
#' @param zip A logical value indicating whether the model specified is zero-inflated (TRUE) or not (FALSE).
#'
#' @return A series of highest density intervals for each parameter of interest.
#' @export
#'
#' @examples
#' mod<-pceg(knee_pain_obs,2,TRUE,TRUE)
#' hdi_beta_extractor(knee_pain_obs,mod,zip=FALSE)
hdi_beta_extractor<-function(data,mod,ci=0.95,N=10000,level_rel_final=-1,zip=TRUE){

  poisson_response<-mod$event.tree$poisson.response
  variable_time<-mod$event.tree$variable.time

  if(level_rel_final == 0 & poisson_response){
    stop("Beta Distribution not for Poisson Response - Use Gamma")
  }
  n=dim(data)[2]-1*variable_time+level_rel_final+1*zip
  p=nlevels(data[,n])
  if(level_rel_final==-1&zip){
    p<-2
  }
  posterior<-mod$prior.distribution[[n]] + mod$data.summary[[n]]

  ind<-which(!is.na(posterior[,1]))
  k<-length(ind)

  if(k>1){
    posterior<-posterior[ind,]
  }else{
    posterior<-matrix(posterior[ind,],nrow=1,byrow=TRUE)
  }

  output<-list()

  for(i in 1:k){
    v<-as.numeric(posterior[i,])
    m=sum(v)
    output[[i]]<-list()
    for(j in 1:p){
      marg_post<-distribution_beta(N,shape1=v[j],shape2=m-v[j])
      output[[i]][[j]]<-hdi(marg_post,ci=ci)
    }
  }
  return(output)
}
