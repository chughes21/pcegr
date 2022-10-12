#' The State Imputer
#'
#' This function imputes the latent risk states based on the data, the estimated parameters, and the choice of imputation method.
#'
#' @param data A data set, where the observed count vector and time vector (if variable) are the last two columns
#' @param lambda A numeric vector of Poisson rates used to impute the states.
#' @param prob A numeric vector of the risk probabilities used to impute the states.
#' @param variable_time A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE).
#' @param stoch A logical value indicating whether stochastic (TRUE) or deterministic (FALSE) imputation should be used.
#' @param all_risk A logical value indicating whether all observations should be considered at risk. This is necessary when converting a PCEG to a ZIPCEG for the purposes of calculating the Bayes Factor.
#'
#' @return A data set containing refactored covariates, imputed risk states, and observed counts and times.
#' @export
#'
#' @examples
state_imputer<-function(data,lambda=1,prob=1,variable_time=TRUE,stoch=TRUE,all_risk=FALSE){
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

  data.out<-cbind(data,state = factor(rep("Risk",length(data[,1])),levels=c("No Risk","Risk")))
  if(variable_time){
    data.out<-data.out[,c(1:n,n+3,n+1,n+2)]
  }else{
    data.out<-data.out[,c(1:n,n+2,n+1)] #n+2 here is n+1 above
  }

  for(i in 1:p){
    v<-tree_matrix[i,]
    ind<-which(row.match(data.use[,1:n],v)==1 )
    y<-data.use[ind,n+1]
    m<-length(y)
    t=data.use[ind,n+2]*variable_time+rep(1,m)*(1-variable_time)

    if(stoch){
      r=(y==0)*(runif(m)<1/(1+(1-prob[i])/(prob[i]*exp(-lambda[i]*t))))+(y>0)
      r[r==0]<-"No Risk"
      r[r==1]<-"Risk"
    }else{
      ind_zero<-which(y==0)
      t_zero<-t[ind_zero]
      ind_order<-order(t_zero,decreasing=TRUE)
      total_risk_free<-ceiling((1-prob[i])*m)
      ind_risk_free<-ind_order[1:total_risk_free]
      ind_risk_free<-ind_zero[ind_risk_free]
      r<-rep("Risk",m)
      r[ind_risk_free]<-"No Risk"
    }

    if(all_risk){
      r<-rep("Risk",length(ind))
    }

    data.out[ind,n+1]<-r
  }
  return(data.out)
}



#' The ZIPCEG function
#'
#' This function fits a ZIPCEG model to the chosen data set, based on chosen methods of parameter estimation and state imputation.
#'
#' @param data A data set where the observed count vector and time vector (if variable) are the last two columns.
#' @param method A character string indicating the method for parameter estimation. The character string can be an element of c("Gibbs","nlm","EM","mle","mm").
#' @param iter The number of iterations for the Gibbs sampler or Expectation-Maximisation algorithm.
#' @param equivsize A numeric value specifying the equivalent sample size for the prior, a measure of confidence in the prior.
#' @param poisson_response A logical value indicating whether the response variable is Poisson (TRUE) or categorical (FALSE).
#' @param variable_time A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE), if applicable.
#' @param stoch_imputation A logical value indicating whether stochastic (TRUE) or deterministic (FALSE) imputation should be used.
#' @param gamma_alpha A numeric value for the shape hyperparameter of the Gamma prior for the Poisson rate, if applicable.
#' @param gamma_beta A numeric value for the rate hyperparameter of the Gamma prior for the Poisson rate.
#' @param c A numeric value for the alpha hyperparameter of the Beta prior for the risk probability.
#' @param d A numeric value for the beta hyperparameter of the Beta prior for the risk probability.
#' @param p_0 A numeric initial value for the risk probability.
#' @param l_0 A numeric initial value for the rate.
#' @param tol A numeric which represents the minimum change in the complete data log likelihood needed to continue the Expectation-Maximisation algorithm.
#' @param var_disc An integer value specifying which variable to discretise. If 0, no discretisation is necessary.
#' @param disc_length An integer value specifying how many neighbours can be searched over for the purposes of variable discetisation. If 0, all other possible stages may be merged over.
#' @param restrict A logical value indicating whether variable discretisation should be restricted to stages with the same unfolding of the process (TRUE) or not (FALSE).
#' @param mirror A logical value indicating whether variable discretisation should be equivalent across each unfolding of the process (TRUE) or not (FALSE).
#' @param cat_limit An integer value specifying the minimum number of categories to the variable can be discretised to. If 0, there is no minimum number of categories.
#'
#' @return A list specifying a ZIPCEG model. The list contains: the prior for the final model, the data for the final model, the stage numbers for the final model,
#' the stage structure for the final model, the vector of likelihoods after each merging, the details of the stages merged at each step,
#' the comparison set of stages left to be merged, a list of the stages merged, and the log marginal likelihood of the final model.
#' @export
#'
#' @examples
zipceg<-function(data,method="Gibbs",iter = 10000, equivsize=2, poisson_response = TRUE,
                      variable_time = FALSE, stoch_imputation = TRUE, gamma_alpha =1, gamma_beta = 2, beta_c = 1, beta_d = 1,
                      p_0 = NA, l_0 = NA, tol=1e-10, var_disc = 0, disc_length = 0, restrict = FALSE, mirror = FALSE, cat_limit = 0){

  if(!(method %in% c("Gibbs","nlm","EM","mle","mm"))){
    stop("Unknown estimation method chosen - Please choose either Gibbs, nlm or EM")
  }

  if((method %in% c("mle","mm")) & variable_time){
    stop(paste0("The ",method," method is only available for uniform time"))
  }

  if((method %in% c("nlm","mle","mm")) & iter){
    warning(paste0("Nonzero amount of iterations chosen for ",method," method"))
  }

  if(stoch_imputation & !variable_time){
    warning("Stochastic Imputation unnecessary for uniform time - deterministic chosen instead")
    stoch_imputation <- FALSE }

  if(method=="Gibbs"){
    out<-gibbs_zip(data,N=iter,variable_time,gamma_alpha,gamma_beta,beta_c,beta_d)
    le<-unlist(lapply(out$lambda,mean))
    pe<-unlist(lapply(out$prob,mean))
  }

  if(method=="nlm"){
    out<-nlm_zip(data,variable_time)
    pe<-out$prob
    le<-out$lambda
  }

  if(method=="EM"){
    out<-em_zip(data,p_0=p_0,l_0=l_0,variable_time=variable_time,max_iter=iter,tol=tol)
    pe<-out$prob
    le<-out$lambda
  }

  if(method=="mle"){
    out<-mle_zip(data,time_input = variable_time)
    pe<-out$prob
    le<-out$lambda
  }

  if(method=="mm"){
    out<-mm_zip(data,time_input = variable_time)
    pe<-out$prob
    le<-out$lambda
  }

  data.adj<-state_imputer(data,lambda=le,prob=pe,variable_time = variable_time, stoch = stoch_imputation)

  return(pceg(data.adj,equivsize=equivsize, poisson_response=poisson_response, variable_time = variable_time, zip=TRUE, gamma_alpha=gamma_alpha, gamma_beta=gamma_beta, var_disc=var_disc, disc_length=disc_length, restrict=restrict, mirror=mirror, cat_limit=cat_limit))

}

#' The Iterative ZIPCEG function
#'
#' This function performs the ZIPCEG model selection a specified number of times, producing plots of the rates and stage ranks, and selecting the model with the largest log marginal likelihood.
#'
#' @param data A data set where the observed count vector and time vector (if variable) are the last two columns.
#' @param method A character string indicating the method for parameter estimation. The character string can be an element of c("Gibbs","nlm","EM","mle","mm").
#' @param iter_f An integer specifying the number of iterations for the parameter estimation method chosen, if necessary.
#' @param iter_total An integer specifying the number of iterations the model selection process should be performed for.
#' @param plot_ranks A logical value indicating whether the stage ranks, based on rates, should be plotted (TRUE) or not (FALSE).
#' @param plot_rates A logical value indicating whether the stage rates should be plotted (TRUE) or not (FALSE).
#' @param plot_probs A logical value indicating whether the risk probabilities should be plotted (TRUE) or not (FALSE).
#' @param hist A logical value indicating whether the plot should be in the form of a histogram (TRUE) or not (FALSE).
#' @param violin A logical value indicating whether the plot should be in the form of a violin plot (TRUE) or not (FALSE).
#' @param scatter A logical value indicating whether the plot should be in the form of a scatter plot (TRUE) or not (FALSE).
#' @param equivsize A numeric value specifying the equivalent sample size for the prior, a measure of confidence in the prior.
#' @param poisson_response A logical value indicating whether the response variable is Poisson (TRUE) or categorical (FALSE).
#' @param variable_time A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE), if applicable.
#' @param stoch_imputation A logical value indicating whether stochastic (TRUE) or deterministic (FALSE) imputation should be used.
#' @param gamma_alpha A numeric value for the shape hyperparameter of the Gamma prior for the Poisson rate, if applicable.
#' @param gamma_beta A numeric value for the rate hyperparameter of the Gamma prior for the Poisson rate.
#' @param c A numeric value for the alpha hyperparameter of the Beta prior for the risk probability.
#' @param d A numeric value for the beta hyperparameter of the Beta prior for the risk probability.
#' @param p_0 A numeric initial value for the risk probability.
#' @param l_0 A numeric initial value for the rate.
#' @param tol A numeric which represents the minimum change in the complete data log likelihood needed to continue the Expectation-Maximisation algorithm.
#' @param var_disc An integer value specifying which variable to discretise. If 0, no discretisation is necessary.
#' @param disc_length An integer value specifying how many neighbours can be searched over for the purposes of variable discetisation. If 0, all other possible stages may be merged over.
#' @param restrict A logical value indicating whether variable discretisation should be restricted to stages with the same unfolding of the process (TRUE) or not (FALSE).
#' @param mirror A logical value indicating whether variable discretisation should be equivalent across each unfolding of the process (TRUE) or not (FALSE).
#' @param cat_limit An integer value specifying the minimum number of categories to the variable can be discretised to. If 0, there is no minimum number of categories.
#'
#'
#' @return A list containing: a matrix of the estimated rates for each leaf, a matrix of the ranks for each leaf, a numeric value of the log marginal likelihood of the MAP model, and the MAP model itself.
#' @export
#'
#' @examples
zipceg.iter<-function(data, method = "Gibbs", iter_f = 10000, iter_total = 1, plot_ranks = TRUE, plot_rates = TRUE,
                           plot_probs = FALSE, hist = FALSE, violin = FALSE, scatter = FALSE, equivsize=2,
                           poisson_response = TRUE, variable_time = FALSE,stoch_imputation=TRUE, gamma_alpha = 1,
                           gamma_beta = 2, beta_c = 1, beta_d = 1,p0=NA,l0=NA,tol=1e-10,
                           var_disc = 0, disc_length = 0, restrict = FALSE, mirror = FALSE, cat_limit = 0){
  if(sum(hist,scatter,violin)>1 ){
    stop("Only 1 display option possible between histogram, violin and scatter") #default is lines
  }

  if((hist | violin) & plot_ranks){
    warning("Histogram and violin only well defined for rates. Ranks have been excluded")
  }

  if(plot_rates & plot_probs){
    stop("Can't Return Rates and State Probabilities - Choose one ")
  }

  n1<-dim(data)[2] - 1 - 1*variable_time#if there are variable times, they will be an extra column
  data_levels<-sapply(data[,1:n1],nlevels)
  data_levels<-c(data_levels,risk = 2)
  n2<-n1+1 #the number of levels after including the probability of risk
  p<-prod(data_levels)

  start_rates<-sum(cumprod(data_levels)[1:n1])+2 #+ 1 for s0, +1 for next situation
  end_rates<-sum(cumprod(data_levels))+1 #+1 for s0

  start_probs<- sum(cumprod(data_levels)[1:(n1-1)])+2 #the risk state prob situations
  end_probs<-start_rates-1

  num<-p/2

  stage_ranks<-matrix(nrow=iter_total,ncol=num)
  rates<-matrix(nrow=iter_total,ncol=num)
  probs<-matrix(nrow=iter_total,ncol=num)

  score<-numeric(iter_total)

  for(j in 1:iter_total){
    ceg.temp<-zipceg(data,method=method,iter = iter_f, equivsize=equivsize, poisson_response = poisson_response,
                          variable_time = variable_time , stoch_imputation = stoch_imputation,
                          gamma_alpha =gamma_alpha, gamma_beta = gamma_beta, beta_c = beta_c, beta_d = beta_d,
                          p_0 = p0, l_0 = l_0, tol=tol, var_disc = var_disc, disc_length = disc_length, restrict = restrict, mirror = mirror)
    risk_stages<-seq(from = start_rates+1, by =2,length.out = num)
    prob_stages<-seq(from = start_probs,to=end_probs)
    merged<-ceg.temp$merged
    merged_rates<-merged[,which(merged[3,]==n2)]
    merged_probs<-merged[,which(merged[3,]==n1)]
    data.mat<-matrix(nrow=num,ncol=8)
    data.mat[,1]<-1:num
    colnames<-c("leaf","y","t","prior_a","prior_b","rate","rate rank","prob")
    for(i in 1:num){
      data.mat[i,2:3]<-ceg.temp$data[[risk_stages[i]]]
      data.mat[i,4]<-data.mat[i,2]/data.mat[i,3]
    }
    data.mat[,5]<-rank(data.mat[,4])
    ind<-which(is.na(data.mat[,2]))
    for(k in ind){
      ref1<-merged[1,which(merged[2,]==risk_stages[k])]
      ref2<-which(risk_stages==ref1)
      data.mat[k,-1]=data.mat[ref2,-1]
    }
    rates[j,]<-data.mat[,4]
    stage_ranks[j,]<-data.mat[,5]

    score[j]<-ceg.temp$lik

    if(j==1){
      map_mod<-ceg.temp
      map_score<-map_mod$lik
    }
    if(j>1){
      if(score[j]>map_score){
        map_mod<-ceg.temp
        map_score<-score[j]
      }
    }

  }

  if(plot_ranks){

    data.temp<-data.frame(x=1:iter_total,leaf=stage_ranks)
    data.long<-melt(data.temp,id.vars="x")

    data.long<-mutate(data.long,label = if_else(x==iter_total,as.character(variable),NA_character_))

    if(scatter){
      rank_plot<-ggplot(data=data.long,aes(x=x,y=value,group=variable))+geom_point(aes(col=variable),position=position_jitter(width=0.1,height=0.01))+
        xlab("Iteration")+ylab("Rank")+
        ggtitle("Stage Rank by Iteration")+
        geom_label_repel(mapping=aes(label=label))
    }else{

      rank_plot<-ggplot(data=data.long,aes(x=x,y=value,group=variable))+geom_line(aes(col=variable),position=position_jitter(width=0.1,height=0.01))+
        xlab("Iteration")+ylab("Rank")+
        ggtitle("Stage Rank by Iteration")+
        geom_label_repel(mapping=aes(label=label))
    }
  }

  if(plot_rates){

    data.temp<-data.frame(x=1:iter_total,leaf=rates)
    data.long<-melt(data.temp,id.vars="x")

    data.long<-mutate(data.long,label = if_else(x==iter_total,as.character(variable),NA_character_))

    if(hist){
      xlim_l<-floor(min(data.temp[,-1]))
      xlim_u<-ceiling(max(data.temp[,-1]))

      leaves<-list()

      #in order to get leaves arranged vertically rather than horizontally

      index_mat<-data.frame(cbind(1:(num/2),((num/2+1):num)))
      index_hor<- c(t(index_mat))

      for(j in 1:num){
        i<-index_hor[j]
        y=data.temp[,i+1]
        m=length(y)
        data.temp2<-data.frame(x=c(1:m),y=y)
        leaves[[j]]<-  ggplot(data=data.temp2,mapping=aes(x=y))+geom_histogram(col=i,fill=i)+
          #geom_density()+
          #  geom_vline(col="red",xintercept=mean(y))+geom_vline(xintercept=true_lambda[i])+
          xlab("Rate") +ggtitle(paste("Leaf", i, sep = " "))+xlim(xlim_l,xlim_u)+ylim(0,iter_total)
      }

      rate_plot<-do.call(ggarrange,list(plotlist=leaves,nrow=num/2,ncol=2))
    }else if(violin){
      rate_plot<- ggplot(data=data.long,aes(x=variable,y=value))+geom_violin(aes(col=variable,fill=variable))+
        xlab("Leaf")+ylab("Rate")+
        ggtitle("Estimated rates for each leaf")+theme(legend.position="none")
    }else if(scatter){
      rate_plot<- ggplot(data=data.long,aes(x=x,y=value,group=variable))+geom_point(aes(col=variable),position=position_jitter(width=0.1,height=0.01))+
        xlab("Iteration")+ylab("Rate")+
        ggtitle("Rate by Iteration")+
        geom_label_repel(mapping=aes(label=label))
    }else{rate_plot<-ggplot(data=data.long,aes(x=x,y=value,group=variable))+geom_line(aes(col=variable),position=position_jitter(width=0.1,height=0.01))+
      xlab("Iteration")+ylab("Rate")+
      ggtitle("Estimated Rate by Iteration")+
      geom_label_repel(mapping=aes(label=label))
    }
  }

  if(hist){
    print(rate_plot)
  }else if(plot_ranks & plot_rates){
    print(ggarrange(rank_plot, rate_plot),nrow=2,ncol=1)
  }else if (plot_ranks){
    print(rank_plot)
  }else{
    print(rate_plot)
  }

  return(list(rates = rates, stage_ranks = stage_ranks,score=score,mod=map_mod))

}
