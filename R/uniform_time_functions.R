#' The Uniform Time ZIP Estimator
#'
#' This function can estimate the risk proportions and rates for counts with uniform time, using either a maximum likelihood or methods of moments estimator.
#'
#' @param data A data set, where the observed count vector and time vector (if included) are the last two columns.
#' @param method A character string indicating the method for parameter estimation. The character string can be an element of c("mle","mm").
#' @param time_input A logical value indicating whether the time has been included in the data set (TRUE) or not (FALSE). As time is uniform, if it is included, then all times should be equal to 1.
#'
#' @return A list containing a matrix with the desired outputs for each evolution of the process, as well as vectors of the estimates for the rates and risk probabilities.
#' @export
#'
#' @examples
#' uniform_time_zip(knee_pain_obs[,1:4])
uniform_time_zip<-function(data,method = "mle",time_input = FALSE){
  n<-dim(data)[2] - 1 - 1*time_input #if there are times, they will be the last column.

  if(time_input){
    if(max(data[,n+2]>1 | min(data[,n+2])<1)){
      return("Error - Time should be uniform") #currently only works for uniform time
    }
  }

  if(!(method %in% c("mle","mm"))){
    stop("Unknown estimation method chosen - Please choose either mle or mm")
  }

  #below is also in other zip functions

  path_details<-refactored_tree_matrix(data,TRUE,time_input)
  data.use<-path_details$data_use
  n<-path_details$num_var
  p<-path_details$p
  tree_matrix<-path_details$tree_matrix

  output_matrix<-cbind(tree_matrix, p_hat=rep(0,p), l_hat=rep(0,p),n_zero = rep(0,p),n_pos=rep(0,p),y_bar=rep(0,p),t_bar = rep(0,p))

  l<-c()
  propor<-c()

  for(i in 1:p){
    v<-tree_matrix[i,]
    ind<-which(row.match(data.use[,1:n],v)==1 )
    data.temp<-data.use[ind,]
    y = data.temp[,n+1]
    m<-length(y)

    if(method == "mle"){
      y_plus<-y[y>0]
      n_plus<-length(y_plus)

      r_plus<-sum(y)/n_plus #sum(y) is same as sum(y_plus) cause difference is zeroes

      y_bar<-mean(y)

      lambda<-lambertW0(-r_plus*exp(-r_plus))+r_plus
      prob<-y_bar/lambda
    }

    if(method == "mm"){

    y_bar<-mean(y)
    s_sq<-var(y)

    lambda<-(y_bar^2+s_sq)/y_bar - 1
    prob<-(y_bar^2)/(y_bar^2+s_sq-y_bar)
    }

    output_matrix$l_hat[i]=lambda
    output_matrix$p_hat[i]=prob
    output_matrix$n_pos[i]=sum(y>0)
    output_matrix$n_zero[i]=m-output_matrix$n_pos[i]
    output_matrix$y_bar[i]=sum(y)
    output_matrix$t_bar[i]=m

    #maybe don't need these

    l[i]<-lambda
    propor[i]<-prob

  }

  return(list(summary=output_matrix,lambda = l, prob = propor))
}

#' The Methods of Moments Estimator under Variable Time
#'
#' Calculates the method of moments estimator (mme) for a ZIP distribution under variable time, by dividing the counts by observation time.
#'
#' @param data A data set, where the observed count vector and time vector (if included) are the last two columns.
#'
#' @return A list containing mmes for at-risk rate and risk proportion for each leaf
#' @export
#'
#' @examples
#' mme_variable_time_zip(knee_pain_obs)
#'
mme_variable_time_zip<-function(data){

  n<-dim(data)[2]
  data.temp<-data[,1:(n-2)]

  y<-data[,n-1]
  t<-data[,n]

  if(min(t)<=0){
    stop("Time must be positive")
  }

  z<-y/t

  data.temp<-cbind(data.temp,z=z)

  return(uniform_time_zip(data.temp,"mm",FALSE))
}

#' The Maximum Likelihood Estimator under Variable Time
#'
#' Calculates the maximum likelihood estimator (mle) for a ZIP distribution under variable time, by dividing the counts by observation time.
#'
#' @param data A data set, where the observed count vector and time vector (if included) are the last two columns.
#'
#' @return A list containing mles for at-risk rate and risk proportion for each leaf
#' @export
#'
#' @examples
#' mle_variable_time_zip(knee_pain_obs)
#'
mle_variable_time_zip<-function(data){

  n<-dim(data)[2]
  data.temp<-data[,1:(n-2)]

  y<-data[,n-1]
  t<-data[,n]

  if(min(t)<=0){
    stop("Time must be positive")
  }

  z<-y/t

  data.temp<-cbind(data.temp,z=z)

  return(uniform_time_zip(data.temp,"mle",FALSE))
}

