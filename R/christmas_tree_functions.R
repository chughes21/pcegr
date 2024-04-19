#' The Counter Function
#'
#' Counts how many elements of a vector match a specific value
#'
#' @param x A value
#' @param v A vector
#'
#' @return An integer
counter<-function(x,v){
  return(length(which(v==x)))
}


#' The Quantile Band Function
#'
#' This function creates a quantile band plot, a diagnostic plot of the suitability of a model to count data.
#'
#' @param data A data set where the observed response vector and time vector (if applicable and variable) are the last two columns.
#' @param mod A StagedTree model fit to the data set, as produced by pceg() or zipceg()
#' @param signif A numeric value specifying the significance level for the quantiles.
#' @param limit An integer vector specifying the maximium number of possible counts to analyse. If it is a vector of length one, this value will be used for all leaves. If NA, this will be the maximum count recorded per leaf.
#' @param shift A logical value indicating whether the raw observed counts and quantiles should be used (FALSE), or whether they should be shifted by the median (TRUE).
#' @param max_per_plot An integer value specifying the maximum number of leaves that can be shown in a single plot.
#' @param plot.leaves A logical value indicating whether quantile band plots should be produced for each leaf (TRUE) or not (FALSE).
#' @param plot.overall A logical value indicating whether a quantile band plot should be produced for the model as a whole (TRUE) or not (FALSE).
#' @param zip A logical value indicating whether the model specified is zero-inflated (TRUE) or not (FALSE).
#'
#' @return A quantile-band plot for each leaf.
#' @export
#'
#' @examples
#' mod<-pceg(knee_pain_obs,2,TRUE,TRUE)
#' quantile_band(knee_pain_obs,mod,limit=10,plot.leaves=TRUE,plot.overall=TRUE,zip=FALSE)
#' zipmod<-zipceg(knee_pain_obs,equivsize=2,poisson_response=TRUE,variable_time=TRUE)
#' quantile_band(knee_pain_obs,zipmod,limit=10,plot.leaves=TRUE,plot.overall=TRUE,zip=TRUE)
quantile_band<-function(data,mod,signif = 0.05, limit=NA,shift = TRUE, max_per_plot = 8, plot.leaves = TRUE, plot.overall = TRUE, zip=FALSE){

  poisson_response<-mod$event.tree$poisson.response
  remove_risk_free<-mod$remove.risk.free.edges
  variable_time<-mod$event.tree$variable.time

  if(poisson_response == FALSE){
    stop("Quantile band plots only well-defined for count models.")
  }

  if(signif <= 0 || signif >= 1){
    stop("Significance level should be greater than 0 and less than 1.")
  }

  #a lot of the below is copied from chi_square - if that changes, so should this

  path_details<-refactored_tree_matrix(data,poisson_response,variable_time)
  data_use<-path_details$data_use
  n<-path_details$num_var
  p<-path_details$p
  tree<-path_details$tree_matrix

  posterior<-mod$posterior.expectation
  stage.struct<-mod$stage.structure

  n1<-n+1*poisson_response +1*zip

  rates<-parameter_extractor(stage.struct,posterior,n1,poisson_response,remove_risk_free)
  probs<-parameter_extractor(stage.struct,posterior,n1-1,poisson_response,remove_risk_free)

  if(zip & !(remove_risk_free)){
    len<-seq(2,2*p,by=2)
    rates<-rates[len]
  }

  probs<-probs[,2]

  max_y<-FALSE

  if(is.na(limit)){
    max_y<-TRUE
    limit_overall<-max(data[,n+1])
  }else if(length(limit)==1){
    limit_overall<-limit
    limit<-rep(limit,p)
  }else if(length(limit) != p){
    stop("Limit vector of incorrect length")
  }else if(length(limit) == p){
    limit_overall<-max(limit)
  }

  leaves<-list()

  count_vec_overall<-rep(0,limit_overall+1)

  if(p > max_per_plot){
    n_plot<-p%/%max_per_plot
    p_plot<-rep(max_per_plot, n_plot)
    rem<-p%%max_per_plot
    if(rem>0){
      p_plot<-c(p_plot,rem)
      n_plot<-n_plot+1
    }
    ind_plot_start<-c(1,cumsum(p_plot)+1)
    ind_plot_end<-ind_plot_start[-1]-1
    ind_plot_start<-ind_plot_start[-(n_plot+1)]
    print(paste0("Number of leaves greater than maximum number of plots per page - ",n_plot," pages needed, use back arrow to see previous plots"))
  }else{
    n_plot<-1
    p_plot<-p
    ind_plot_start<-1
    ind_plot_end<-p}

  for(k in 1:p){

    v<-tree[k,c(1:n)]
    ind<-which(row.match(data_use[,1:n],v)==1 )
    lambda<-rates[k]

    if(zip){
      prop<-probs[k]
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

    quant_vec<-matrix(nrow=lim+1,ncol=2)
    median_vec<-numeric(length = lim+1)

    x<-c(0:lim)

    lim_diff<-limit_overall-lim

    count_vec<-sapply(x,counter,v=y)
    prob_mat<-sapply(x,f,p=prop,lambda=lambda,t=t)
    temp<-apply(prob_mat,MARGIN = 2, FUN = qpoibin, qq = c(signif/2,1-signif/2,0.5),wts=NULL)

    quant_vec<-temp[1:2,]
    median_vec<-temp[3,]

    if(shift){
      count_vec<-count_vec-median_vec
      quant_vec<-quant_vec-matrix(rep(median_vec,2),nrow=2,byrow=TRUE)
    }

    if(shift){
    data.temp<-data.frame(x,count = count_vec,left = quant_vec[1,],right = quant_vec[2,])
    }else{
    data.temp<-data.frame(x,count = count_vec,left = quant_vec[1,],right = quant_vec[2,],median = median_vec)
    }

    if(shift){
    leaves[[k]]<-ggplot(data=data.temp)+
    geom_path(mapping = aes(x = left, y = x),col="green")+geom_point(mapping = aes(x = left, y = x),col="green")+
    geom_path(mapping = aes(x = right, y = x),col="green")+geom_point(mapping = aes(x = right, y = x),col="green")+
    geom_path(mapping = aes(x = count, y = x),col="red")+geom_point(mapping = aes(x = count, y = x),col="red")+
    xlab("Shifted Observed Counts")+ylab("Event counts")+ggtitle(paste0("Shifted Quantile Band Plot - Leaf ",k))
    }else{
      leaves[[k]]<-ggplot(data=data.temp)+
        geom_path(mapping = aes(x = left, y = x),col="green")+geom_point(mapping = aes(x = left, y = x),col="green")+
        geom_path(mapping = aes(x = right, y = x),col="green")+geom_point(mapping = aes(x = right, y = x),col="green")+
        geom_path(mapping = aes(x = count, y = x),col="red")+geom_point(mapping = aes(x = count, y = x),col="red")+
        geom_path(mapping = aes(x = median, y = x),col="black")+
        xlab("Raw Observed Counts")+ylab("Event counts")+ggtitle(paste0("Raw Quantile Band Plot - Leaf ",k))
    }

    count_vec_temp<-c(count_vec,rep(0,lim_diff))
    count_vec_overall<-count_vec_overall+count_vec_temp

    prob_mat_temp<-cbind(prob_mat,matrix(0,nrow=dim(prob_mat)[1],ncol=lim_diff))
    if(k==1){
      prob_mat_overall<-prob_mat_temp
    }else{
      prob_mat_overall<-rbind(prob_mat_overall,prob_mat_temp)
    }
  }
  if(plot.leaves){
  for(i in 1:n_plot){
  p.temp<-p_plot[i]
  start.temp<-ind_plot_start[i]
  end.temp<-ind_plot_end[i]
  ind.temp<-c(start.temp:end.temp)
  leaves.temp<-list()
  for(j in 1:p.temp){
    leaves.temp[[j]]<-leaves[[ind.temp[j]]]
  }
  print(do.call(ggarrange,list(plotlist=leaves.temp,nrow=ceiling(p.temp/2),ncol=2)))
  }}

  if(plot.overall){
    temp<-apply(prob_mat_overall,MARGIN = 2, FUN = qpoibin, qq = c(signif/2,1-signif/2,0.5),wts=NULL)

    quant_vec_overall<-temp[1:2,]
    median_vec_overall<-temp[3,]
    x<-c(0:limit_overall)

    if(shift){
     # count_vec_overall<-count_vec_overall-median_vec_overall
      quant_vec_overall<-quant_vec_overall-matrix(rep(median_vec_overall,2),nrow=2,byrow=TRUE)
    }

    if(shift){
      data.temp<-data.frame(x,count = count_vec_overall,left = quant_vec_overall[1,],right = quant_vec_overall[2,])
    }else{
      data.temp<-data.frame(x,count = count_vec_overall,left = quant_vec_overall[1,],right = quant_vec_overall[2,],median = median_vec_overall)
    }
    if(shift){
      ggplot(data=data.temp)+
        geom_path(mapping = aes(x = left, y = x),col="green")+geom_point(mapping = aes(x = left, y = x),col="green")+
        geom_path(mapping = aes(x = right, y = x),col="green")+geom_point(mapping = aes(x = right, y = x),col="green")+
        geom_path(mapping = aes(x = count, y = x),col="red")+geom_point(mapping = aes(x = count, y = x),col="red")+
        xlab("Shifted Observed Counts")+ylab("Event counts")+ggtitle("Shifted Quantile Band Plot - Overall")
    }else{
      ggplot(data=data.temp)+
        geom_path(mapping = aes(x = left, y = x),col="green")+geom_point(mapping = aes(x = left, y = x),col="green")+
        geom_path(mapping = aes(x = right, y = x),col="green")+geom_point(mapping = aes(x = right, y = x),col="green")+
        geom_path(mapping = aes(x = count, y = x),col="red")+geom_point(mapping = aes(x = count, y = x),col="red")+
        geom_path(mapping = aes(x = median, y = x),col="black")+
        xlab("Raw Observed Counts")+ylab("Event counts")+ggtitle("Raw Quantile Band Plot - Overall")
    }
  }

}

#mid_quantile_(signif, lambda, prop, t)

