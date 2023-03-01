#' The variance inflation factor calculation for a level of a staged tree
#'
#' @param mod A StagedTree model fit to the data set, as produced by pceg() or zipceg().
#' @param var A positive integer value detailing which level of the tree to calculate the VIF for.
#'
#' @return A numeric value for the variance inflation factor.
#' @export
#'
#' @examples
#' mod1<-pceg(knee_pain_obs,2,TRUE,TRUE)
#' vif(mod1,3)
vif<-function(mod,var){

  if(var < 0){
   stop("The level of the tree must be positive")
  }

  if(var > length(mod$data.summary)){
    stop("The inputted level is not in the tree")
  }

  if((var == length(mod$data.summary)) &  mod$event.tree$poisson.response){
    stop("VIF is undefinded for Poisson response variables")
  }

  post<-mod$posterior.expectation[[var]]
  data<-mod$data.summary[[var]]
  prior<-mod$prior.distribution[[var]]

  ind<-which(!is.na(post[,1]))
  J<-length(ind)

  l1<-0

  for(j in ind){
    pv<-post[j,]
    dv<-data[j,]
    l1<-l1+sum(dv*log(pv))
  }

  if(J > 1){
  data_null<-colSums(data[ind,])
  prior_null<-colSums(prior[ind,])
  post_null<-(data_null+prior_null)/sum(data_null+prior_null)

  l0<-sum(data_null*log(post_null))
  }else{
    l0<-l1
  }
  R2<-1-(l1/l0)
  vif<-1/(1-R2)

  return(vif)
}
