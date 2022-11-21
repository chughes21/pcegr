#' The Quantile Band Function
#'
#' This function creates a quantile band plot, a diagnostic plot of the suitability of a model to count data.
#'
#' @param data A data set where the observed response vector and time vector (if applicable and variable) are the last two columns.
#' @param mod A StagedTree model fit to the data set, as produced by pceg() or zipceg()
#' @param signif A numeric value specifying the significance level for the quantiles.
#' @param limit An integer vector specifying the maximium number of possible counts to analyse. If it is a vector of length one, this value will be used for all leaves. If NA, this will be the maximum count recorded per leaf.
#' @param shift A logical value indicating whether the raw observed counts and quantiles should be used (FALSE), or whether they should be shifted by the median (TRUE).
#' @param poisson_response A logical value indicating whether the response variable is Poisson (TRUE) or categorical (FALSE).
#' @param variable_time A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE), if applicable.
#' @param zip A logical value indicating whether the model specified is zero-inflated (TRUE) or not (FALSE).
#'
#' @return A quantile-band plot for each leaf.
#' @export
#'
#' @examples
#' mod<-pceg(knee_pain_obs,2,TRUE,TRUE)
#' quantile_band(knee_pain_obs,mod,limit=10,zip=FALSE)
quantile_band<-function(data,mod,signif = 0.05, limit=NA,shift = TRUE, poisson_response=TRUE,variable_time=TRUE,zip=FALSE){
  if(poisson_response == FALSE){
    stop("Quantile band plots only well-defined for count models.")
  }

  if(signif <= 0 || signif >= 1){
    stop("Significance level should be greater than 0 and less than 1.")
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

  max_y<-FALSE

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
      median_vec[j+1]<-temp[3]

    }
    if(shift){
      count_vec<-count_vec-median_vec
      quant_vec<-quant_vec-median_vec
    }

    x<-c(0:lim)

    if(shift){
    data.temp<-data.frame(x,count = count_vec,left = quant_vec[,1],right = quant_vec[,2])
    }else{
    data.temp<-data.frame(x,count = count_vec,left = quant_vec[,1],right = quant_vec[,2],median = median_vec)
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

  }
  print(do.call(ggarrange,list(plotlist=leaves,nrow=p/2,ncol=2)))

}