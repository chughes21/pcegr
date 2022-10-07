#EM Algorithm Functions

#' The complete data log likelihood for a Poisson response
#'
#' This function calculates the complete data log likelihood for a count vector and time vector, given values of the parameters of a Poisson distribution.
#'
#' @param params The parameters for the distribution
#' @param y The observed count vector
#' @param t The observed time vector
#'
#' @return A list containing the likelihood value and a vector of weights for the em algorithm
#' @export
#'
#' @examples
#' y<-c(1,0,2,0)
#' t<-c(0.5,0.7,0.8,0.3)
#' params<-c(0.5,1.5)
#' cdll(params,y,t)
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
#' @param p_0 The initial value for the risk probability
#' @param l_0 The initial value for the rate
#' @param y The observed count vector
#' @param t The observed time vector
#' @param tol The tolerance to stop the algorithm
#' @param max_iter The maximum number of iterations for the algorithm
#'
#' @return A list containing the estimates of the risk probability and the rate, and the number of iterations necessary
#' @export
#'
#' @examples
em<-function(p_0,l_0,y,t,tol = 1e-10,max_iter=10000){

  n<-length(y)

  p<-p_0
  lambda<-l_0

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
  }

  return(list(p=p,lambda=lambda,iter=i))
}



#' The EM algorithm
#'
#' This function carries out the Expectation-Maximisation algorithm across all leaves.
#'
#' @param data The entire data set, where the observed count vector and time vector are the last two columns
#' @param p_0 The initial value for the risk value
#' @param l_0 The initial value for the rate. If NA
#' @param variable_time An indicator of whether the observed time is uniform or variable
#' @param max_iter The tolerance to stop the algorithm
#' @param tol The maximum number of iterations for the algorithm
#'
#' @return
#' @export
#'
#' @examples
em_zip<-function(data,p_0=NA,l_0=NA, variable_time = TRUE, max_iter = 10000,tol=1e-10){
  n<-dim(data)[2] - 1 - 1*variable_time #if there are variable times, they will be an extra column
  data_levels<-sapply(data[,1:n],nlevels)
  data.use<-path_refactor(data,n,data_levels)
  p<-prod(data_levels)
  Z<-lapply(data_levels,vec_from)#changed because if they don't have same number of levels for each, we get a list.
  #so we start with a list
  tree_matrix<-Z[[1]]
  for(i in 2:n){
    tree_matrix<-expand_grid(tree_matrix,Z[[i]]) #changed to list too
    colnames(tree_matrix)[1:i]<-colnames(data.use[,1:n])[1:i]
  }
  tree_matrix<-as.data.frame(tree_matrix)

  output_matrix<-cbind(tree_matrix, p_hat=rep(0,p), l_hat=rep(0,p),n_zero = rep(0,p),n_pos=rep(0,p),y_bar=rep(0,p),t_bar = rep(0,p))

  l<-c()
  propor<-c()

  if(is.na(p_0)){
    p_0<-rep(0.5,p)
  }else if(length(p_0)==1){
    p_0<-rep(p_0,p)
  }

  if(is.na(l_0)){
    l_0<-rep(1,p)
  }else if(length(l_0)==1){
    l_0<-rep(l_0,p)
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

    result<-em(prob,lambda,y,t,tol=tol,max_iter=max_iter)
    lambda<-result$lambda
    prob<-result$p

    output_matrix$l_hat[i]=lambda
    output_matrix$p_hat[i]=prob
    output_matrix$n_pos[i]=round(output_matrix$p_hat[i]*m)
    output_matrix$n_zero[i]=m-output_matrix$n_pos[i]
    output_matrix$y_bar[i]=sum(y)
    output_matrix$t_bar[i]=sum(t)

    propor[i]<-prob
    l[i]<-lambda

  }

  return(list(summary=output_matrix,lambda=l,prob=propor))
}
