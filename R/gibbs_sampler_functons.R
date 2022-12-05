#' The Path refactor function
#'
#' This function takes the categorical covariates from a data set and refactors each one from 0 to the number of levels - 1.
#'
#' @param data A data set where the covariates are factors.
#' @param n An integer for the number of columns to refactor, beginning from the start.
#' @param data_levels An integer vector of the number of levels for each covariate.
#'
#' @return A data set where the categorical covariates have been refactor to start at 0 and end at the number of levels - 1.
#'
#'
path_refactor<-function(data,n,data_levels){
  data_out<-data
  for(i in 1:n){
    n_levels=data_levels[i]
    l_max=n_levels-1
    levels_temp<-levels(data[,i])
    data_temp_check<-factor(data[,i],levels=unique(c(levels_temp,c(l_max:0))))
    data_temp_change<-data_temp_check
    for(j in 1:n_levels){
      data_temp_change[which(data_temp_check == levels_temp[j])]=l_max-j+1
    }
    data_out[,i]=factor(data_temp_change,levels=c(l_max:0))
  }
  return(data_out)
}

#' The Gibbs Sampler
#'
#' A function to carry out Gibbs sampling for each unique covariate combination in the data set.
#'
#' @param data A data set where the observed count vector and time vector (if variable) are the last two columns.
#' @param N An integer specifying the number of iterations for the Gibbs sampler.
#' @param variable_time A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE).
#' @param a A numeric vector specifying the shape hyperparameter of the Gamma priors for the Poisson rates. If the vector is of length 1, this value will be repeated for each leaf.
#' @param b A numeric vector specifying the rate hyperparameter of the Gamma priors for the Poisson rates. If the vector is of length 1, this value will be repeated for each leaf.
#' @param c A numeric vector specifying the alpha hyperparameter of the Beta priors for the risk probabilities. If the vector is of length 1, this value will be repeated for each leaf.
#' @param d A numeric vector specifying the beta hyperparameter of the Beta priors for the risk probabilities. If the vector is of length 1, this value will be repeated for each leaf.
#'
#'
#' @return A list containing a matrix with the desired outputs for each evolution of the process, as well as vectors of the estimates for the rates and risk probabilities.
#' @export
#'
#' @examples
#' gibbs_result<-gibbs_zip(knee_pain_obs)
#' gibbs_result$summary
gibbs_zip<-function(data,N = 1000, variable_time = TRUE, a = 1, b = 1, c = 1, d  =  1 ){

  #below is also in other zip functions

  path_details<-refactored_tree_matrix(data,variable_time)
  data.use<-path_details$data_use
  n<-path_details$num_var
  p<-path_details$p
  tree_matrix<-path_details$tree_matrix

  output_matrix<-cbind(tree_matrix, p_hat=rep(0,p), l_hat=rep(0,p),n_zero = rep(0,p),n_pos=rep(0,p),y_bar=rep(0,p),t_bar = rep(0,p))

  l<-list()
  propor<-list()

  if(length(a) != length(b)){
    warning("Length of vectors of prior hyperparameters for Gamma priors don't match.")
  }

  if(length(c) != length(d)){
    warning("Length of vectors of prior hyperparameters for Beta priors don't match.")
  }

  ind_even<-seq(2,2*p,2)

  if(length(a)==1){
    a<-rep(a,p)
  }else if(length(a)==2*p){
    a<-a[ind_even]
  }else{stop("Vector of shape hyperparameters doesn't match number of leaves.")}

  if(length(b)==1){
    b<-rep(b,p)
  }else if(length(b)==2*p){
    b<-b[ind_even]
  }else{stop("Vector of rate hyperparameters doesn't match number of leaves.")}

  if(length(c)==1){
    c<-rep(c,p)
  }else if(length(c)==2*p){
    c<-c[ind_even]
  }else{stop("Vector of alpha hyperparameters doesn't match number of leaves.")}

  if(length(d)==1){
    d<-rep(d,p)
  }else if(length(d)==2*p){
    d<-d[ind_even]
  }else{stop("Vector of beta hyperparameters doesn't match number of leaves.")}

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

    r<-rep(1,m)
    r[which(y==0)]<-0
    if(sum(y)>0){
    l_0<-sum(y)/(sum(t[y>0]))
    }else{
    l_0<-0
    }
    p_0<-sum(r)/m #could maybe base this on l_0
    lambda<-rep(l_0,N)
    prob<-rep(p_0,N)

    for(j in 2:N){
      r=(y==0)*(runif(m)<1/(1+(1-prob[j-1])/(prob[j-1]*exp(-lambda[j-1]*t))))+(y>0)
      lambda[j]=rgamma(1,a[i]+sum(y),b[i]+sum(r*t))
      prob[j]=rbeta(1,c[i]+sum(r),m-sum(r)+d[i])
    }
    output_matrix$l_hat[i]=mean(lambda)
    output_matrix$p_hat[i]=mean(prob)
    output_matrix$n_pos[i]=round(output_matrix$p_hat[i]*m)
    output_matrix$n_zero[i]=m-output_matrix$n_pos[i]
    output_matrix$y_bar[i]=sum(y)
    output_matrix$t_bar[i]=sum(t)

    l[[i]]<-lambda
    propor[[i]]<-prob

  }

  return(list(summary=output_matrix,lambda = l, prob = propor))
}
