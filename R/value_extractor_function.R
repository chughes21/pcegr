#value extractor

#value extractor can calculate the posterior or sample estimate of the rates for a Poisson
#or probs for a multinomial

#' The Expected Value Extractor
#'
#' This function calculates the expected values of parameters based on the data set and a chosen pceg model.
#'
#' @param data A data set, where the observed response vector and time vector (if applicable and variable) are the last two columns
#' @param ceg A ceg model fit to the data set, as produced by pceg().
#' @param level_rel_final A non-positive integer indicating where the desired variable is relative to the response variable.
#' @param poisson_response A logical value indicating whether the response variable is Poisson (TRUE) or categorical (FALSE).
#' @param poisson_time_variable A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE), if applicable.
#' @param posterior A logical value indicating whether the estimates of the posterior (TRUE) or sample (FALSE) expected value should be calculated.
#' @param zip A logical value indicating whether the model specified is zero-inflated (TRUE) or not (FALSE).
#' @param true_value A numeric vector specifying the true values of the parameters for comparison, if known. If unknown, this will be NA.
#'
#' @return A matrix displaying the observed data, the prior values, and the expected values of the parameters.
#' @export
#'
#' @examples
value_extractor<-function(data,ceg,level_rel_final = 0,poisson_response=TRUE,poisson_time_variable=TRUE,posterior = TRUE, zip=TRUE, true_value = NA){

  if(!poisson_response & zip){
    stop("Zero Inflated Poisson Requires Poisson Response")
  }

  n<-dim(data)[2] - 1 - 1*poisson_time_variable#if there are variable times, they will be an extra column
  data_levels<-sapply(data[,1:n],nlevels)
  if(zip){
    data_levels<-c(data_levels,risk = 2)}
  n<-n+1*zip #the number of levels of the tree after including the probability of risk
  p<-prod(data_levels)

  n_level<-n + level_rel_final #the level we care about

  if(n_level>1){
    start<-sum(cumprod(data_levels)[1:(n_level-1)])+2 #1 for s0, 1 for next situation
    end<-sum(cumprod(data_levels)[1:n_level]) + 1 #+1 for s0
  }else if(n_level == 1){
    start<-2
    end<-1+data_levels[1]
  }else{
    start<-1
    end<-1
  }
  ind<-ceg$stages
  ind<-ind[ind>=start & ind<=end]

  true_value_input<-!all(is.na(true_value)) #could also do any(!is.na)

  p=length(ceg$data[[ind[1]]])

  k1<-3*p
  k2<- k1 +1*true_value_input

  out.mat<-matrix(nrow=length(ind),ncol=k2)

  j<-1


  for(i in seq_along(ind)){
    out.mat[i,1:p]<-ceg$data[[ind[i]]]
    if(posterior){
      out.mat[i,(p+1):(2*p)]<-ceg$prior[[ind[i]]]
    }else{
      out.mat[i,(p+1):(2*p)]<-0
    }
    if(level_rel_final == 0 & poisson_response){
      out.mat[i,k1]<-(out.mat[i,1]+out.mat[i,3])/(out.mat[i,2]+out.mat[i,4])
    }else{
      out.mat[i,(2*p+1):k1]<-(out.mat[i,1:p]+out.mat[i,(p+1):(2*p)])/sum(out.mat[i,1:(2*p)])
    }

    if(true_value_input){
      if(level_rel_final ==0 & i == 1 & zip){ #if it's the first leaf, so all the no risks
        out.mat[i,k2]<-NA
      }else{
        out.mat[i,k2]<-true_value[j]
        j=j+1
      }
    }

  }
  return(out.mat)
}
