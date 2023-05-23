#' The complete data log likelihood for a Poisson response
#'
#' This function calculates the complete data log likelihood for a count vector and time vector, given values of the parameters of a zero-inflated Poisson distribution.
#'
#' @param params A numeric vector of the parameters for the distribution.
#' @param y An integer vector of observed event counts.
#' @param t A numeric vector of observed times.
#'
#' @return A list containing the likelihood value and a vector of weights for the em algorithm
#'
cdll<-function(params,y,t){
  # p <-glogitinv(params[1]) #for constraints, so pi is a prop
  #  lambda<-exp(params[-1])
  p<-params[1]
  lambda<-params[-1]
  omega <- p*exp(-lambda*t)/(p*exp(-lambda*t)+(1-p)*(y==0))
  lik<-sum(omega*log(p) + (1-omega)*log(1-p)+ omega*log(dpois(y,lambda*t))) #if using nlm, use minus
  return(list(lik=lik,weights=omega))
}

#' The EM algorithm for a leaf
#'
#'This function carries out the Expectation-Maximisation algorithm for a leaf, given initial values of the parameters and the count and time vectors.
#'
#'For a specific leaf, which is a unique evolution of the process, the event counts and observed times are used, alongside the [cdll()] function and some initial parameter values,
#'to carry out the Expectation-Maximisation algorithm.
#'
#' @param p_0 A numeric initial value for the risk probability
#' @param l_0 A numeric initial value for the rate.
#' @param y An integer vector of observed event counts.
#' @param t A numeric vector of observed times.
#' @param max_iter An integer for the maximum number of iterations for the algorithm.
#' @param tol A numeric value which represents the minimum change in the complete data log likelihood needed to continue the algorithm.
#'
#' @return A list containing the estimates of the risk probability and the rate, and the number of iterations necessary
#' @export
#'
#' @examples
#' y<-c(1,0,2,0)
#' t<-c(0.5,0.7,0.8,1.2)
#' params<-c(0.5,1.5)
#' em(params[1],params[2],y,t)
em<-function(p_0,l_0,y,t,max_iter=10000,tol = 1e-10){

  n<-length(y)

  p<-p_0
  lambda<-l_0

  p_chain<-c(p)
  lambda_chain<-c(lambda)

  result<-cdll(c(p,lambda),y,t)
  lik<-result$lik
  omega<-result$weights

  diff<-1
  i<-0

  while(diff>tol & i <= max_iter){
    p<-sum(omega)/n
    lambda<-sum(omega*y)/sum(omega*t)
    temp<-cdll(c(p,lambda),y,t)
    omega<-temp$weights
    val<-temp$lik
    diff=val-lik
    lik<-val
    i<-i+1

    p_chain<-c(p_chain,p)
    lambda_chain<-c(lambda_chain,lambda)

  }

  return(list(p=p,lambda=lambda,p_chain=p_chain,lambda_chain=lambda_chain,iter=i))
}



#' The EM algorithm
#'
#' This function carries out the Expectation-Maximisation algorithm across all leaves.
#'
#' For each leaf, which represents a unique evolution of the process, the [em()] function is applied.
#'
#' @param data A data set, where the observed count vector and time vector (if variable) are the last two columns.
#' @param p_0 A numeric vector of initial values for the risk probabilities. If the vector is of length 1, this value will be repeated for each leaf.
#' @param l_0 A numeric vector of initial values for the rates. If the vector is of length 1, this value will be repeated for each leaf
#' @param variable_time A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE).
#' @param max_iter An integer for the maximum number of iterations for the algorithm.
#' @param tol A numeric which represents the minimum change in the complete data log likelihood needed to continue the algorithm.
#'
#' @return A list containing a matrix with the desired outputs for each evolution of the process, as well as vectors of the estimates for the rates and risk probabilities.
#' @export
#'
#' @examples
#' em_zip(knee_pain_obs)
em_zip<-function(data,p_0=0.5,l_0=1, variable_time = TRUE, max_iter = 10000,tol=1e-10){

  #below is also in other zip functions

  path_details<-refactored_tree_matrix(data,TRUE,variable_time)
  data.use<-path_details$data_use
  n<-path_details$num_var
  p<-path_details$p
  tree_matrix<-path_details$tree_matrix

  output_matrix<-cbind(tree_matrix, p_hat=rep(0,p), l_hat=rep(0,p),n_zero = rep(0,p),n_pos=rep(0,p),y_bar=rep(0,p),t_bar = rep(0,p))

  l<-list()
  propor<-list()

  if(length(p_0) != length(l_0)){
    warning("Length of vectors of initial values for parameters don't match.")
  }

  if(length(p_0)==1){
    p_0<-rep(p_0,p)
  }else if(length(p_0)!=p){
    stop("Initial value vector for risk probability doesn't match number of leaves.")
  }

  if(length(l_0)==1){
    l_0<-rep(l_0,p)
  }else if(length(l_0)!=p){
    stop("Initial value vector for rates doesn't match number of leaves.")
  }

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

    prob<-p_0[i]
    lambda<-l_0[i]

    result<-em(prob,lambda,y,t,max_iter=max_iter,tol=tol)
    lambda<-result$lambda
    prob<-result$p
    lambda_chain<-result$lambda_chain
    p_chain<-result$p_chain

    output_matrix$l_hat[i]=lambda
    output_matrix$p_hat[i]=prob
    output_matrix$n_pos[i]=sum(y>0)
    output_matrix$n_zero[i]=m-output_matrix$n_pos[i]
    output_matrix$y_bar[i]=sum(y)
    output_matrix$t_bar[i]=sum(t)

    propor[[i]]<-prob
    l[[i]]<-lambda

  }

  return(list(summary=output_matrix,lambda=l,prob=propor))
}
