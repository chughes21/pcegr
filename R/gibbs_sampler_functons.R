#' The vec_from function
#'
#' A function to create a descending vector of natural numbers
#'
#' @param x A positive integer
#'
#' @return A vector of integers from $x-1$ to 0.
#' @export
#'
#' @examples
#' vec_from(7)
vec_from<-function(x){
  return(c((x-1):0))
}



#' The Path refactor function
#'
#' This function takes the categorical covariates from a data set and refactors each one from 0 to the number of levels - 1.
#'
#' @param data A data set where the covariates are factors.
#' @param n An integer for the number of columns to refactor, beginning from the start.
#' @param data_levels An integer vector of the number of levels for each covariate.
#'
#' @return A data set where the categorical covariates have been refactor to start at 0 and end at the number of levels - 1.
#' @export
#'
#' @examples
#' N<-100
#' h<-sample(c("Short","Average","Tall"),size=N,replace=TRUE)
#' h<-factor(h,levels=c("Short","Average","Tall"))
#' w<-sample(c("Normal","Over"),size=N,replace=TRUE)
#' w<-factor(w,levels=c("Normal","Over"))
#' y<-rpois(N,lambda=2)
#' data<-data.frame(height=h,weight=w,y=y)
#' summary(data)
#' n<-2
#' data_levels<-sapply(data[,1:n],nlevels)
#' data.refactor<-path_refactor(data,n=2,data_levels)
#' summary(data.refactor)
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
#' @param a A numeric value for the shape hyperparameter of the Gamma prior for the Poisson rate.
#' @param b A numeric value for the rate hyperparameter of the Gamma prior for the Poisson rate.
#' @param c A numeric value for the alpha hyperparameter of the Beta prior for the risk probability.
#' @param d A numeric value for the beta hyperparameter of the Beta prior for the risk probability.
#'
#'
#' @return A list containing a matrix with the desired outputs for each evolution of the process, as well as vectors of the estimates for the rates and risk probabilities.
#' @export
#'
#' @examples
gibbs_zip<-function(data,N = 1000, variable_time = TRUE, a = 1, b = 1, c = 1, d  =  1 ){
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

  l<-list()
  propor<-list()

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
    l_0<-sum(y)/(sum(t[y>0]))
    p_0<-sum(r)/m #could maybe base this on l_0
    lambda<-rep(l_0,N)
    prob<-rep(p_0,N)

    for(j in 2:N){
      r=(y==0)*(runif(m)<1/(1+(1-prob[j-1])/(prob[j-1]*exp(-lambda[j-1]*t))))+(y>0)
      lambda[j]=rgamma(1,a+sum(y),b+sum(r*t))
      prob[j]=rbeta(1,c+sum(r),m-sum(r)+d)
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
