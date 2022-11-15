quantile_band<-function(data,mod,signif = 0.05, limit=NA,shift = TRUE, poisson_response=TRUE,variable_time=TRUE,zip=TRUE){
  if(poisson_response == FALSE){
    stop("Quantile band plots only well-defined for count models.")
  }

  #a lot of the below is copied from chi_square - if that changes, so should this

  path_details<-refactored_tree_matrix(data,variable_time)
  data_use<-path_details$data_use
  n<-path_details$num_var
  p<-path_details$p
  tree<-path_details$tree_matrix
  data_levels<-path_details$data_levels

  output<-merge_separator(mod,n,p,tree,data_levels,zip)
  tree<-output$tree
  rates<-output$rates
  probs<-output$probs

  if(is.na(limit)){
    max_y<-TRUE
  }else if(length(limit)==1){
    limit<-rep(limit,p)
  }else if(length(limit) != p){
    stop("Limit vector of incorrect length")
  }

  leaves<-list()

  for(k in 1:p){

    v<-tree[k,c(1:n)]
    ind<-which(row.match(data_use[,1:n],v)==1 )
    lambda_stage<-tree$rate_stage[k]
    lambda<-rates[lambda_stage]

    if(zip){
      prop_stage<-tree$prob_stage[k]
      prop<-probs[prop_stage]
    }else{
      prop<-1
    }

    y<-data_use[ind,n+1]

    if(variable_time){
      t<-data_use[ind,n+2]
    }else{
      t<-rep(1,length(ind))
    }

    if(max_y){
      lim<-max(y)
    }else{
      lim<-limit[k]
    }

    count_vec<-numeric(length = lim+1)
    quant_vec<-matrix(nrow=lim+1,ncol=2)
    median_vec<-numeric(length = lim+1)

    for(j in 0:lim){
      count_vec[j+1]<-length(which(y == j))
      prob_vec<-f(prop,lambda,j,t)
      temp<-qpoibin(c(signif/2,1-signif/2,0.5),prob_vec)
      quant_vec[j+1,]<-temp[1:2]
      median_vec<-temp[3]
      if(shift){
        count_vec<-count_vec-median_vec
        quant_vec<-quant_vec-median_vec
      }
    }

    x<-c(0:lim)

    if(shift){
    data.temp<-data.frame(x,count = count_vec,left = quant_vec[,1],right = quant_vec[,2])
    }else{
    data.temp<-data.frame(x,count = count_vec,left = quant_vec[,1],right = quant_vec[,2],median = median_vec)
    }

    if(shift){
    leaves[[i]]<-ggplot(data=data.temp)+
    geom_line(mapping = aes(x = left, y = x),col="green")+geom_point(mapping = aes(x = left, y = x),col="green")+
    geom_line(mapping = aes(x = right, y = x),col="green")+geom_point(mapping = aes(x = right, y = x),col="green")+
    geom_line(mapping = aes(x = count, y = x),col="red")+geom_point(mapping = aes(x = count, y = x),col="red")+
    xlab("Shifted Observed Counts")+ylab("Event counts")+ggtitle("Shifted QUantile Band Plot")
    }else{
      leaves[[k]]<-ggplot(data=data.temp)+
        geom_line(mapping = aes(x = left, y = x),col="green")+geom_point(mapping = aes(x = left, y = x),col="green")+
        geom_line(mapping = aes(x = right, y = x),col="green")+geom_point(mapping = aes(x = right, y = x),col="green")+
        geom_line(mapping = aes(x = count, y = x),col="red")+geom_point(mapping = aes(x = count, y = x),col="red")+
        geom_line(mapping = aes(x = median, y = x),col="black")+
        xlab("Raw Observed Counts")+ylab("Event counts")+ggtitle("Raw Quantile Band Plot")
    }

  }
  print(do.call(ggarrange,list(plotlist=leaves,nrow=p/2,ncol=2)))

}
