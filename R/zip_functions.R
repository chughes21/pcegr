#' The State Imputer
#'
#' This function imputes the latent risk states based on the data, the estimated parameters, and the choice of imputation method.
#'
#' @param data A data set, where the observed count vector and time vector (if variable) are the last two columns
#' @param lambda A numeric vector of the Poisson rates used to impute the states. If the vector is of length 1, this value will be repeated for each leaf.
#' @param prob A numeric vector of the risk probabilities used to impute the states.If the vector is of length 1, this value will be repeated for each leaf.
#' @param variable_time A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE).
#' @param stoch A logical value indicating whether stochastic (TRUE) or deterministic (FALSE) imputation should be used.
#' @param all_risk A logical value indicating whether all observations should be considered at risk. This is necessary when converting a PCEG to a ZIPCEG for the purposes of calculating the Bayes Factor.
#'
#' @return A data set containing refactored covariates, imputed risk states, and observed counts and times.
#' @export
#'
#' @examples state_imputer(knee_pain_obs)
state_imputer<-function(data,lambda=1,prob=0.5,variable_time=TRUE,stoch=TRUE,all_risk=FALSE){
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

  if(length(prob) != length(lambda)){
    warning("Vectors of parameters don't match.")
  }

  if(length(prob)==1){
    prob<-rep(prob,p)
  }else if(length(prob)!=p){
    stop("Vector of risk probabilities doesn't match number of leaves.")
  }

  if(length(lambda)==1){
    lambda<-rep(lambda,p)
  }else if(length(lambda)!=p){
    stop("Vector of rates doesn't match number of leaves.")
  }

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
#' This function takes the same inputs as the [pceg()], along with chosen methods of parameter estimation and state imputation, to fit a ZIPCEG model. As with [pceg()] it can perform variable discretisation methods, but is not to be used for vanilla CEGs.
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
#' the comparison set of stages left to be merged, and the log marginal likelihood of the final model.
#' @export
#'
#' @examples  zipceg(knee_pain_obs,"nlm",variable_time=TRUE)
zipceg<-function(data,method="Gibbs",iter = 10000, equivsize=2, poisson_response = TRUE,
                      variable_time = FALSE, stoch_imputation = TRUE, gamma_alpha =1, gamma_beta = 2, beta_c = 1, beta_d = 1,
                      p_0 = 0.5, l_0 = 1, tol=1e-10, var_disc = 0, disc_length = 0, restrict = FALSE, mirror = FALSE, cat_limit = 0){

  if(!(method %in% c("Gibbs","nlm","EM","mle","mm"))){
    stop("Unknown estimation method chosen - Please choose either Gibbs, nlm, EM, mle or mm")
  }

  if((method %in% c("mle","mm")) & variable_time){
    stop(paste0("The ",method," method is only available for uniform time"))
  }

  if((method %in% c("nlm","mle","mm")) & (iter!=10000)){
    warning(paste0("Nonzero, non-default amount of iterations chosen for ",method," method"))
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

  if(method %in% c("mle","mm")){
    out<-uniform_time_zip(data,method=method,time_input = variable_time)
    pe<-out$prob
    le<-out$lambda
  }

  data.adj<-state_imputer(data,lambda=le,prob=pe,variable_time = variable_time, stoch = stoch_imputation)

  return(pceg(data.adj,equivsize=equivsize, poisson_response=poisson_response, variable_time = variable_time, zip=TRUE, gamma_alpha=gamma_alpha, gamma_beta=gamma_beta, var_disc=var_disc, disc_length=disc_length, restrict=restrict, mirror=mirror, cat_limit=cat_limit))

}

#' The Iterative ZIPCEG function
#'
#' This function performs the ZIPCEG model selection a specified number of times and can produce plots and perform model selection.
#'
#' This function takes the same inputs as the [zipceg()] function, except carries out the model selection process a specified number of times. Using the results of each iteration, different types of plots can be used to display either the estimated risk probabilities or estimated rates. The default plot is a violin plot. The function will also select the _maximum_ _a_ _posteriori_ (MAP) model from these iterations.
#'
#' @param data A data set where the observed count vector and time vector (if variable) are the last two columns.
#' @param method A character string indicating the method for parameter estimation. The character string can be an element of c("Gibbs","nlm","EM","mle","mm").
#' @param iter_total An integer specifying the number of iterations the model selection process should be performed for.
#' @param iter_f An integer specifying the number of iterations for the parameter estimation method chosen, if necessary.
#' @param plot_rates A logical value indicating whether the stage rates should be plotted (TRUE) or not (FALSE).
#' @param plot_probs A logical value indicating whether the risk probabilities should be plotted (TRUE) or not (FALSE).
#' @param hist A logical value indicating whether the plot should be in the form of a histogram (TRUE) or not (FALSE).
#' @param violin A logical value indicating whether the plot should be in the form of a line plot (TRUE) or not (FALSE).
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
#' @return A list containing: a matrix of the estimated rates or risk probabilities for each leaf across iterations, a numeric value of the log marginal likelihood for the MAP model, and the MAP model itself.
#' @export
#'
#' @examples zipceg.iter(knee_pain_obs,"nlm",iter_total=100,plot_ranks=FALSE,violin=TRUE)
zipceg.iter<-function(data, method = "Gibbs", iter_total = 10, iter_f = 1000, plot_rates = TRUE,
                           plot_probs = FALSE, hist = FALSE, line = FALSE, scatter = FALSE, equivsize=2,
                           poisson_response = TRUE, variable_time = FALSE,stoch_imputation=TRUE, posterior = TRUE, gamma_alpha = 1,
                           gamma_beta = 2, beta_c = 1, beta_d = 1,p0=NA,l0=NA,tol=1e-10,
                           var_disc = 0, disc_length = 0, restrict = FALSE, mirror = FALSE, cat_limit = 0){
  if(sum(hist,scatter,line)>1 ){
    stop("Only 1 display option possible between histogram, line and scatter") #default is lines
  }

  if(plot_rates & plot_probs){
    stop("Only 1 plot option possible between rates and risk probabilities ")
  }

  path_details<-refactored_tree_matrix(data,variable_time)
  n<-path_details$num_var
  p<-path_details$p
  tree<-path_details$tree_matrix
  data_levels<-path_details$data_levels
  data_levels_zip<-c(data_levels,risk=2)

  n_zip<-n+1

  if(plot_rates){
  rates<-matrix(nrow=iter_total,ncol=p)
  }

  if(plot_probs){
  probs<-matrix(nrow=iter_total,ncol=p)
  }

  if(n>1){

    start_probs<-sum(cumprod(data_levels)[1:(n-1)])+2 #1 for s0, 1 for next situation
    end_probs<-sum(cumprod(data_levels)[1:n]) + 1 #+1 for s0

    start_rates<-sum(cumprod(data_levels_zip)[1:(n_zip-1)])+2 #1 for s0, 1 for next situation
    end_rates<-sum(cumprod(data_levels_zip)[1:n_zip]) + 1

  }else if(n == 1){
    start_probs<-2
    end_probs<-1+data_levels[1]

    start_rates<-end+1
    end_rates<- 1+ sum(cumprod(data_levels_zip)[1:n_zip])

  }else{
    start<-1
    end<-1

    start_zip<-2
    end_zip<-3 #this might not be right

  }

  score<-numeric(iter_total)

  for(i in 1:iter_total){
    ceg.temp<-zipceg(data,method=method,iter = iter_f, equivsize=equivsize, poisson_response = poisson_response,
                          variable_time = variable_time , stoch_imputation = stoch_imputation,
                          gamma_alpha =gamma_alpha, gamma_beta = gamma_beta, beta_c = beta_c, beta_d = beta_d,
                          p_0 = p0, l_0 = l_0, tol=tol, var_disc = var_disc, disc_length = disc_length, restrict = restrict, mirror = mirror)

    ind<-ceg.temp$stages
    ind_probs<-ind[ind>=start_probs & ind<=end_probs]
    ind_rates<-ind[ind>=start_rates & ind<=end_rates]

    rates.temp<-value_extractor(data,ceg.temp,level_rel_final=0,poisson_response=poisson_response,variable_time=variable_time,posterior=posterior,zip=TRUE)
    probs.temp<-value_extractor(data,ceg.temp,level_rel_final=-1,poisson_response,variable_time,posterior=posterior,zip=TRUE) #if zip=FALSE, we'll ignore this anyway

    rates.temp<-rates.temp[,dim(rates.temp)[2]] #if start inputting true values, this may be incorrect
    probs.temp<-probs.temp[,dim(probs.temp)[2]]

    seq_prob<-seq(start_probs,end_probs)
    seq_rate<-seq(start_rates+1,end_rates,by=2)
    if(i == 1){
    tree<-cbind(tree,prob_stage=seq_prob,rate_stage=seq_rate)
    }else{
      tree$prob_stage<-seq_prob
      tree$rate_stage<-seq_rate
    }

    merged<-ceg.temp$merged
    m<-max(merged[3,])
    m_prob<-m-1

    merged_rates<-merged[,which(merged[3,]==m)]
    merged_probs<-merged[,which(merged[3,]==m_prob)]

    merged_list_rates<-merged_list_extractor(merged_rates)
    merged_list_prob<-merged_list_extractor(merged_probs)

    if(plot_rates){

    k<-1

    stage_count<-0

    for(j in 1:length(ind_rates)){

      if(k <= length(merged_list_rates)){
        if(ind_rates[j]==merged_list_rates[[k]][1]){
          stage_comp<-merged_list_rates[[k]]
          replace<-which(tree$rate_stage %in% stage_comp)
          tree$rate_stage[replace]<-j
          k<-k+1
          stage_count<-stage_count+length(replace)
        }else{
          replace<-which(tree$rate_stage == ind_rates[j])
          tree$rate_stage[replace]<-j
          stage_count<-stage_count+length(replace)
        }
      }else{
        replace<-which(tree$rate_stage == ind_rates[j])
        tree$rate_stage[replace]<-j
        stage_count<-stage_count+length(replace)
      }
    }

    if(stage_count != p){
      stop("All stages not accounted for - fix")
    }

    rates[i,]<-rates.temp[tree$rate_stage]

    }

    if(plot_probs){
      k<-1

      stage_count<-0

      for(j in 1:length(ind_probs)){

        if(k <= length(merged_list_prob)){
          if(ind_probs[j]==merged_list_prob[[k]][1]){
            stage_comp<-merged_list_prob[[k]]
            replace<-which(tree$prob_stage %in% stage_comp)
            tree$prob_stage[replace]<-j
            k<-k+1
            stage_count<-stage_count+length(replace)
          }else{
            replace<-which(tree$prob_stage == ind_probs[j])
            tree$prob_stage[replace]<-j
            stage_count<-stage_count+length(replace)
          }
        }else{
          replace<-which(tree$prob_stage == ind_probs[j])
          tree$prob_stage[replace]<-j
          stage_count<-stage_count+length(replace)
        }
      }

      if(stage_count != p){
        stop("All stages not accounted for - fix")
      }

      probs[i,]<-probs.temp[tree$prob_stage]
    }

    score[i]<-ceg.temp$lik

    if(i==1){
      map_mod<-ceg.temp
      map_score<-map_mod$lik
    }
    if(i>1){
      if(score[i]>map_score){
        map_mod<-ceg.temp
        map_score<-score[i]
      }
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

      index_mat<-data.frame(cbind(1:(p/2),((p/2+1):p)))
      index_hor<- c(t(index_mat))

      for(j in 1:p){
        i<-index_hor[j]
        y=data.temp[,i+1]
        m=length(y)
        data.temp2<-data.frame(x=c(1:m),y=y)
        leaves[[j]]<-  ggplot(data=data.temp2,mapping=aes(x=y))+geom_histogram(col=i,fill=i)+
          #geom_density()+
          #  geom_vline(col="red",xintercept=mean(y))+geom_vline(xintercept=true_lambda[i])+
          xlab("Rate") +ggtitle(paste("Leaf", i, sep = " "))+xlim(xlim_l,xlim_u)+ylim(0,iter_total)
      }

      rate_plot<-do.call(ggarrange,list(plotlist=leaves,nrow=p/2,ncol=2))
    }else if(line){
      rate_plot<-ggplot(data=data.long,aes(x=x,y=value,group=variable))+geom_line(aes(col=variable),position=position_jitter(width=0.1,height=0.01))+
        xlab("Iteration")+ylab("Rate")+
        ggtitle("Estimated Rate by Iteration")+
        geom_label_repel(mapping=aes(label=label))
    }else if(scatter){
      rate_plot<- ggplot(data=data.long,aes(x=x,y=value,group=variable))+geom_point(aes(col=variable),position=position_jitter(width=0.1,height=0.01))+
        xlab("Iteration")+ylab("Rate")+
        ggtitle("Rate by Iteration")+
        geom_label_repel(mapping=aes(label=label))
    }else{
    rate_plot<- ggplot(data=data.long,aes(x=variable,y=value))+geom_violin(aes(col=variable,fill=variable))+
      xlab("Leaf")+ylab("Rate")+
      ggtitle("Estimated rates for each leaf")+theme(legend.position="none")
    }

  print(rate_plot)
  }

  if(plot_probs){

    data.temp<-data.frame(x=1:iter_total,leaf=probs)
    data.long<-melt(data.temp,id.vars="x")

    data.long<-mutate(data.long,label = if_else(x==iter_total,as.character(variable),NA_character_))

    if(hist){
      xlim_l<-floor(min(data.temp[,-1]))
      xlim_u<-ceiling(max(data.temp[,-1]))

      leaves<-list()

      #in order to get leaves arranged vertically rather than horizontally

      index_mat<-data.frame(cbind(1:(p/2),((p/2+1):p)))
      index_hor<- c(t(index_mat))

      for(j in 1:p){
        i<-index_hor[j]
        y=data.temp[,i+1]
        m=length(y)
        data.temp2<-data.frame(x=c(1:m),y=y)
        leaves[[j]]<-  ggplot(data=data.temp2,mapping=aes(x=y))+geom_histogram(col=i,fill=i)+
          #geom_density()+
          #  geom_vline(col="red",xintercept=mean(y))+geom_vline(xintercept=true_lambda[i])+
          xlab("Risk Prob") +ggtitle(paste("Leaf", i, sep = " "))+xlim(xlim_l,xlim_u)+ylim(0,iter_total)
      }

      prob_plot<-do.call(ggarrange,list(plotlist=leaves,nrow=p/2,ncol=2))
    }else if(line){
      prob_plot<-ggplot(data=data.long,aes(x=x,y=value,group=variable))+geom_line(aes(col=variable),position=position_jitter(width=0.1,height=0.01))+
        xlab("Iteration")+ylab("Risk Prob")+
        ggtitle("Estimated Risk Prob by Iteration for Each leaf")+
        geom_label_repel(mapping=aes(label=label))
    }else if(scatter){
      prob_plot<- ggplot(data=data.long,aes(x=x,y=value,group=variable))+geom_point(aes(col=variable),position=position_jitter(width=0.1,height=0.01))+
        xlab("Iteration")+ylab("Risk Prob")+
        ggtitle("Estimated Risk Prob by Iteration for Each leaf")+
        geom_label_repel(mapping=aes(label=label))
    }else{
    prob_plot<- ggplot(data=data.long,aes(x=variable,y=value))+geom_violin(aes(col=variable,fill=variable))+
      xlab("Leaf")+ylab("Risk Prob")+
      ggtitle("Estimated Risk Prob by Iteration for Each leaf")+theme(legend.position="none")
    }

    print(prob_plot)
  }

  if(plot_rates){
    output<-rates
  }else if(plot_probs){
    output<-probs
  }else{output<-NA}


  return(list(output=output, score=score,mod=map_mod))

}
