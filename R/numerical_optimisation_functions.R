#' The logit function
#'
#' This function maps a probability to a value on the real line, for the purposes of the [nlm_zip()] function.
#'
#' @param p A numeric value between 0 and 1
#'
#' @return A numeric value on the real line.
#'
glogit = function (p){ log(p/(1-sum(p)))}

#' The logistic function
#'
#' This function maps a value on the real line to a value between 0 and 1. It is the inverse of the [glogit()] function.
#'
#' @param p A numeric value on the real line.
#'
#' @return A numeric value between 0 and 1.
#'
glogitinv = function(p){ exp(p)/(1+sum(exp(p)))}

#' The zero-inflated Poisson likelihood for a single observation.
#'
#' This function calculates the likelihood for a single observation assuming a zero-inflated Poisson distribution
#'
#' @param p A numeric value between 0 and 1 specifying the risk probability.
#' @param lambda A numeric value specifying the rate of the Poisson process.
#' @param y An integer value for the observed event count.
#' @param t A positive numeric value for the observation time.
#'
#' @return A numeric value specifying the likelihood.
#'
f<-function(p, lambda, y, t){
  return(p*(dpois(y,lambda*t))+(1-p)*(y==0))
}

#' The zero-inflated Poisson log likelihood
#'
#' This function calculates the zero-inflated Poisson log likelihood for a count vector and time vector, given values of the parameters of a zero-inflated Poisson distribution.
#'
#' @param params A numeric vector specfiying the parameters for the distribution.
#' @param y An integer vector of observed event counts.
#' @param t A numeric vector of observation times.
#'
#' @return A numeric value specifying the zero-inflated Poisson log likelihood.
#' @export
#'
#' @examples
#' params<-c(0.5,1)
#' y<-c(1,2,0,0,0,1,5)
#' t<-c(0.1,0.5,0.2,0.25,0.1,0.3,2)
#' ziplike(params,y,t)

ziplike<-function(params, y, t){
  p <-glogitinv(params[1]) #for constraints, so pi is a prop
  lambda<-exp(params[-1]) #for constraints, so lambda is nonneg
  loglik<--sum(log(f(p,lambda,y,t)))
  loglik
}

#' The Numerical Optimisation Method
#'
#' @param data A data set, where the observed count vector and time vector (if variable) are the last two columns.
#' @param variable_time A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE).
#'
#' @return A list containing a matrix with the desired outputs for each evolution of the process, as well as vectors of the estimates for the rates and risk probabilities.
#' @export
#'
#' @examples
#' nlm_zip(knee_pain_obs)
nlm_zip<-function(data,variable_time = TRUE){
  #below is also in other zip functions

  path_details<-refactored_tree_matrix(data,TRUE,variable_time)
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

    if(variable_time){
      t<-data.temp[,n+2]
    }else{
      t<-rep(1,m)
    }

    if(sum(y)==0){
      p_0 <- 1
      l_0<-0
    }else{
      p_0<-glogit(sum(y>0)/m) #could maybe base this on l_0
      l_0<-log(sum(y)/(sum(t[y>0])))
    }

    est<-suppressWarnings(nlm(ziplike,c(p_0,l_0),y = y,t = t)$estimate) #leaf = i is just for when it's being printed really
    prob<-glogitinv(est[1])
    lambda<-exp(est[2])

    output_matrix$l_hat[i]=lambda
    output_matrix$p_hat[i]=prob
    output_matrix$n_pos[i]=sum(y>0)
    output_matrix$n_zero[i]=m-output_matrix$n_pos[i]
    output_matrix$y_bar[i]=sum(y)
    output_matrix$t_bar[i]=sum(t)

    #maybe don't need these

    l[i]<-lambda
    propor[i]<-prob

  }

  return(list(summary=output_matrix,lambda = l, prob = propor))
}
