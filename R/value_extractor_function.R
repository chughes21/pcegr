#' The Expected Value Extractor
#'
#' This function calculates the expected values of parameters based on the data set and a chosen PCEG model.
#'
#' This function is capable of calculating the expected rate for the leaf stages in a PCEG or ZIPCEG, or the expected transition probabilities of the situation stages in any type of CEG.
#'
#' @param data A data set, where the observed response vector and time vector (if applicable and variable) are the last two columns
#' @param ceg A ceg model fit to the data set, as produced by pceg().
#' @param level_rel_final A non-positive integer indicating where the desired variable is relative to the response variable.
#' @param poisson_response A logical value indicating whether the response variable is Poisson (TRUE) or categorical (FALSE).
#' @param variable_time A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE), if applicable.
#' @param posterior A logical value indicating whether the estimates of the posterior (TRUE) or sample (FALSE) expected value should be calculated.
#' @param zip A logical value indicating whether the model specified is zero-inflated (TRUE) or not (FALSE).
#' @param dec_place An integer value detailing how many decimal places the outputs should be rounded to. If NA, no rounding will occur.
#' @param true_value A numeric vector specifying the true values of the parameters for comparison, if known. If unknown, this will be NA.
#'
#' @return A matrix displaying the observed data, the prior values, and the expected values of the parameters.
#' @export
#'
#' @examples
#' mod<-pceg(knee_pain_obs,2,TRUE,TRUE)
#' value_extractor(knee_pain_obs,mod,zip=FALSE)
value_extractor<-function(data,ceg,level_rel_final = 0,poisson_response=TRUE,variable_time=TRUE,posterior = TRUE, zip=TRUE, dec_place = NA, true_value = NA){

  if(!poisson_response & zip){
    stop("Zero Inflated Poisson Requires Poisson Response")
  }

  if(!poisson_response & variable_time){
    stop("Variable Time Requires Poisson Response")
  }

  n<-dim(data)[2] - 1 - 1*variable_time#if there are variable times, they will be an extra column
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

  out.mat<-data.frame(matrix(nrow=length(ind),ncol=k2))

  j<-1

  num<-c(1:p)

  if(level_rel_final == 0 & poisson_response){
    colnames(out.mat)[1:p]<-c("y_bar","t_bar")
    colnames(out.mat)[(p+1):(2*p)]<-c("prior_a","prior_b")
    colnames(out.mat)[k1]<-"exp_value"
  }else{
    colnames(out.mat)[1:p]<-paste0("cat",num)
    colnames(out.mat)[(p+1):(2*p)]<-paste0("prior",num)
    colnames(out.mat)[(2*p+1):k1]<-paste0("exp_value",num)
  }

  if(true_value_input){
    colnames(out.mat)[k2]<-"true"
  }

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

  if(level_rel_final == 0 & poisson_response){
    out.mat<-out.mat[,-(2*p+1)]
  }

  if(!(is.na(dec_place))){
    out.mat<-round(out.mat,dec_place)
  }


  return(out.mat)
}

#' The Total Expected Value Extractor
#'
#' This function prints the expected values of parameters based on the data set and a chosen PCEG model for each variable.
#'
#' This function is the [value_extractor()] function repeated at each level.
#'
#' @param data A data set, where the observed response vector and time vector (if applicable and variable) are the last two columns
#' @param ceg A ceg model fit to the data set, as produced by pceg().
#' @param level_exclude A non-positive integer integer indicating which variables (if any) should be excluded from the output, is relative to the response variable.
#' @param poisson_response A logical value indicating whether the response variable is Poisson (TRUE) or categorical (FALSE).
#' @param variable_time A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE), if applicable.
#' @param posterior A logical value indicating whether the estimates of the posterior (TRUE) or sample (FALSE) expected value should be calculated.
#' @param zip A logical value indicating whether the model specified is zero-inflated (TRUE) or not (FALSE).
#' @param dec_place An integer value detailing how many decimal places the outputs should be rounded to. If NA, no rounding will occur.
#'
#' @return
#' @export
#'
#' @examples
total_value_extractor<-function(data,ceg,level_exclude = NA,poisson_response=TRUE,variable_time=TRUE,posterior = TRUE, zip=TRUE, dec_place = NA){
  n<-dim(data)[2]-1*variable_time

  ind<-c(-(n-1):0)

  if(!(is.na(level_exclude))){
    ind<-ind[-which(ind %in% level_exclude)]
  }

  for(i in ind){
    print(paste0("Level ",n+i,"- ",colnames(data)[i]))
    print(value_extractor(data,ceg,level_rel_final = i,poisson_response,variable_time,posterior,zip,dec_place))
  }
}


