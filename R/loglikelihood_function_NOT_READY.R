#Currently not ready due to missing the data dependent part of a Poisson
#Also could fold into the VIF

#' Log Likelihood Calculator
#'
#' Calculates the log likelihood and BIC for a StagedTree object.
#'
#' @param mod A StagedTree object.
#'
#' @return A list containing the log likelihood, number of stages and BIC
#' @export
#'
#' @examples
#' mod<-pceg(knee_pain_obs,2,TRUE,TRUE)
#' loglik(mod)
loglik<-function(mod){

  data.sum<-mod$data.summary
  posterior<-mod$posterior.expectation

  poisson.response<-mod$event.tree$poisson.response

  stages<-0
  result<-0
  n<-0

  for(i in 1:length(posterior)){

    if(poisson.response & i==length(posterior)){
      poisson<-TRUE
    }else{poisson<-FALSE}

    d<-data.sum[[i]]
    post<-posterior[[i]]

    ind<-which(!is.na(d[,1]))

    d<-d[ind,]
    if(poisson){
     post<-post[ind]
    }else{
    post<-post[ind,]
    }

    for(j in 1:length(ind)){
     if(poisson){

       stop("Function doesn't work for Poissons - yet")


       }else{
      if(length(ind)>1){

        temp<-sum(d[j,]*log(post[j,]))

       }else{

        temp<-sum(d*log(post))

       }
      }
    result<-result+temp
    if(i==1 || (length(ind)==1)){
      stages<-stages+length(post)-1
    }else{
    stages<-stages+dim(post)[2]-1
     }
    }

    n<-n+sum(d)

  }

  return(list(loglik=result,stages=stages,BIC=stages*log(n)-2*result))

}
