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
    output_matrix$n_pos[i]=round(output_matrix$p_hat[i]*m)
    output_matrix$n_zero[i]=m-output_matrix$n_pos[i]
    output_matrix$y_bar[i]=sum(y)
    output_matrix$t_bar[i]=m

    #maybe don't need these

    l[i]<-lambda
    propor[i]<-prob

  }

  return(list(summary=output_matrix,lambda = l, prob = propor))
}