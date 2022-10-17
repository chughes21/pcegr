#' The Gamma Highest Density Interval Creator
#'
#' This function calculates the highest density interval for the parameters of a Gamma posterior, as with a Poisson response variable.
#'
#' This function calculates the highest density interval for the estimated rates of the leaf stages in a PCEG or ZIPCEG. The expected value itself can be found using [value_extractor()].
#'
#' @param data A data set, where the observed count vector and time vector (if variable) are the last two columns
#' @param ceg A ceg model fit to the data set, as produced by CEG.AHC.POISS.
#' @param ci A numeric value between 0 and 1 specifying the perecentage confidence for the highest density interval.
#' @param N An integer specifying the number of iterations for the empirical functions used.
#' @param variable_time A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE).
#' @param zip A logical value indicating whether the model specified is zero-inflated (TRUE) or not (FALSE).
#'
#' @return A series of highest density intervals for each parameter.
#' @export
#'
#' @examples
#' mod<-pceg(knee_pain_obs,2,TRUE,TRUE)
#' hdi_gamma_extractor(knee_pain_obs,mod,zip=FALSE)
hdi_gamma_extractor<-function(data,ceg,ci=0.95,N=10000,variable_time=TRUE,zip=TRUE){
  summ<-value_extractor(data=data,ceg=ceg,variable_time = variable_time,zip=zip)
  shapes<-summ[,1]+summ[,3]
  scales<-1/(summ[,2]+summ[,4])
  n=length(shapes)
  output<-list()
  for(i in 1:n){
    posterior<-distribution_gamma(N,shape=shapes[i],scale=scales[i])
    output[[i]]<-hdi(posterior,ci=ci)
  }
  return(output)
}

#' The Beta Highest Density Interval Creator
#'
#' This function calculates the highest density interval for the parameters of a Dirichlet posterior, using a Beta distribution.
#' The Dirichlet posterior is present for the situations with a categorical response variable in any type of CEG.
#'
#' This function calculates the highest density interval for the estimated transition probabilities for the situation stages in a PCEG, ZIPCEG or CEG. The expected value itself can be found using [value_extractor()].
#'
#' @param data A data set, where the observed count vector and time vector (if variable) are the last two columns
#' @param ceg A ceg model fit to the data set, as produced by CEG.AHC.POISS.
#' @param ci A numeric value between 0 and 1 specifying the perecentage confidence for the highest density interval.
#' @param N An integer specifying the number of iterations for the empirical functions used.
#' @param level_rel_final A non-positive integer indicating where the desired variable is relative to the response variable.
#' @param poisson_response A logical value indicating whether the response variable is Poisson (TRUE) or categorical (FALSE).
#' @param variable_time A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE).
#' @param zip A logical value indicating whether the model specified is zero-inflated (TRUE) or not (FALSE).
#'
#' @return A series of highest density intervals for each parameter of interest.
#' @export
#'
#' @examples
#' mod<-pceg(knee_pain_obs,2,TRUE,TRUE)
#' hdi_beta_extractor(knee_pain_obs,mod,zip=FALSE)
hdi_beta_extractor<-function(data,ceg,ci=0.95,N=10000,level_rel_final=-1,poisson_response=TRUE,variable_time = TRUE,zip=TRUE){
  summ<-value_extractor(data=data,ceg=ceg,level_rel_final = level_rel_final,poisson_response = poisson_response,variable_time = variable_time,zip=zip)
  if(level_rel_final == 0 & poisson_response){
    stop("Beta Distribution not for Poisson Response - Use Gamma")
  }
  n=dim(data)[2]-1*variable_time+level_rel_final
  p=nlevels(data[,n])

  a<-summ[,1:p]
  prior<-summ[,(p+1):(2*p)]

  a_star<-matrix(a+prior,ncol=p)

  k<-length(a_star[,1])

  output<-list()

  for(i in 1:k){
    v<-a_star[i,]
    m=sum(v)
    output[[i]]<-list()
    for(j in 1:p){
      posterior<-distribution_beta(N,shape1=v[j],shape2=m-v[j])
      output[[i]][[j]]<-hdi(posterior,ci=ci)
    }
  }
  return(output)
}
