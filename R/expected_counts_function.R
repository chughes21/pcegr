#' The Merged Stage Finder
#'
#' This function finds which stages have merged with a given unique stage from the merging details of a CEG.
#'
#' @param unique_stage An integer value
#' @param M An integer matrix detailing the merging present in a CEG.
#'
#' @return An integer vector detailing the stages that have merged with the unique stage given
#'
merged_stage_finder<-function(unique_stage,M){
  M_unique<-as.matrix(M[,which(M[1,]==unique_stage)])
  stages_merged<-c(unique_stage,M_unique[2,])
  stages_merged<-as.vector(sort(stages_merged,decreasing=FALSE))
  return(stages_merged)
}

#' The Merged List Extractor
#'
#' This function finds a list of the stages which have merged given the merging details of a CEG.
#'
#' @param merged A matrix specifying the merging present in a CEG.
#'
#' @return A list of integer vectors detailing the stages which have merged.
#'
merged_list_extractor<-function(merged){
  level<-max(merged[3,])
  merged_stages<-as.matrix(merged[,which(merged[3,]==level)])
  stages_top<-unique(merged_stages[1,])
  stages_bottom<-unique(merged_stages[2,])

  stage_ref<-stages_bottom[which(sapply(stages_bottom,FUN=function(x) vec_in(x,stages_top))==TRUE)]

  num_ref<-length(stage_ref)

  merged_stages_adj<-merged_stages

  for(i in 1:num_ref){

    col_ref_top<-which(merged_stages_adj[1,]==stage_ref[i])
    col_ref_bottom<-which(merged_stages_adj[2,]==stage_ref[i])
    new_top_stage<-merged_stages[1,col_ref_bottom]
    merged_stages_adj[1,col_ref_top]<-new_top_stage

  }

  stages_top_adj<-as.vector(sort(unique(merged_stages_adj[1,]),decreasing=FALSE))

  output<-lapply(stages_top_adj,function(x) merged_stage_finder(x,merged_stages_adj))
  return(output)
}

#' The vec_in function
#'
#' This function checks whether a given value is in a vector.
#'
#' @param x A numeric value
#' @param v A numeric vector
#'
#' @return A logical value indicating whether the given value was in the vector (TRUE) or not (FALSE)
#'
vec_in<-function(x,v){
  check<-x %in% v
  return(check)
}

#' The merge_Separator function
#'
#' @param mod An object of the S3 class StagedTree
#' @param n An integer value
#' @param p An integer value
#' @param tree A matrix
#' @param data_levels A numeric vector
#' @param zip A logical value indicating whether the model specified is zero-inflated (TRUE) or not (FALSE).
#'
#' @return A list containing a matrix and two vectors.
merge_separator<-function(mod,n,p, tree,data_levels, zip=FALSE){

  if(zip){
    data_levels_zip<-c(data_levels,risk=2)
  }else{
    data_levels_zip<-data_levels
  }

  n_zip<-n+1*zip

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
  ind<-mod$stages
  ind_probs<-ind[ind>=start_probs & ind<=end_probs]
  ind_rates<-ind[ind>=start_rates & ind<=end_rates]

  rates<-mod$posterior.expectation[[n_zip+1]]
  probs<-mod$posterior.expectation[[n_zip]] #if zip=FALSE, we'll ignore this anyway

  probs<-probs[,dim(probs)[2]]

  seq_prob<-seq(start_probs,end_probs)
  if(zip){
    seq_rate<-seq(start_rates+1,end_rates,by=2)
    seq_temp<-seq_rate-(start_rates-1)
    rates<-rates[c(1,seq_temp)]
  }else{
    seq_rate<-seq(start_rates,end_rates)
  }

  if(zip){
    tree<-cbind(tree,prob_stage=seq_prob,rate_stage=seq_rate)
  }else{
    tree<-cbind(tree,rate_stage=seq_rate)
  }

  merged<-mod$merged
  m<-max(merged[3,])
  m_prob<-m-1

  merged_rates<-merged[,which(merged[3,]==m)]
  merged_probs<-merged[,which(merged[3,]==m_prob)]

  merged_list_rates<-merged_list_extractor(merged_rates)
  merged_list_prob<-merged_list_extractor(merged_probs)

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

  if(zip){

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

  }

  rates<-rates[!(is.na(rates))]
  probs<-probs[!(is.na(probs))]

  return(list(tree = tree, rates = rates, probs = probs))

}

#' The Parameter Extractor Function
#'
#' @param stage_struct A list detailing the stage structure from a StagedTree object.
#' @param posterior A list detailing the posterior expectations of parameters from a StagedTree object.
#' @param var An integer value indicating which variable should have its parameters extracted.
#' @param poisson_response A logical value indicating whether the response variable is Poisson (TRUE) or categorical (FALSE).
#' @param remove_risk_free A logical value indicating whether the risk free leaves and edges should be removed (TRUE) or not (FALSE).
#'
#' @return A numeric matrix or vector (for Poisson rates) detailing the individual posterior expectations for each situation for the chosen variable.
parameter_extractor<-function(stage_struct, posterior, var, poisson_response = TRUE, remove_risk_free = TRUE){

  if(remove_risk_free & !poisson_response){
    stop("Risk and Risk free requires Poisson Response ")
  }

  n<-length(posterior)

  post<-as.matrix(posterior[[var]])
  solution<-stage_struct[[var]]

  if(var>1){
  m<-dim(post)[1]

  rrf_ind<-remove_risk_free & (n == var)

  if(rrf_ind){
    len<-seq(from=2,to=m,by=2)
    m<-m/2
  }else{
    len<-c(1:m)
  }

  if(m != length(solution)){
    stop("Posterior and Stage Structure have different lengths - check inputs again")
  }

  output<-as.matrix(post[len,])

  for(i in len){
    if(rrf_ind){
      i<-i/2
    }
    ind<-solution[[i]]
    j<-min(ind)
    if(rrf_ind){
      j<-2*j
    }
    if(!is.na(j)){
    #output[ind,]<-matrix(rep(post[j,],times=length(ind)),byrow=TRUE,nrow=length(ind)) used to do it this way
    #I'm not sure why I used min - when doing model_combining it breaks it - maybe for rrf?
    output[ind,]<-matrix(rep(post[i+i*rrf_ind,],times=length(ind)),byrow=TRUE,nrow=length(ind))
    }
  }

  if(poisson_response & (n==var)){
    output<-as.vector(output)
  }
  }else{
    output<-post
  }
  return(output)
}



#' The Chi Square calculator function.
#'
#' This function tests the goodness of fit of a pceg model to a data set using a Chi squared calculation.
#'
#' For each leaf stage in a pceg model created by [pceg()], this function calculates the observed event counts up to some integer limit, and then calculates
#' the expected event counts for that leaf based on the model estimates of the parameters for the Poisson distribution. Then, for each leaf $i$ and each event count $j$,
#' a chi-square contribution is calculated using the formula \eqn{\frac{(O_{ij}-E_{ij})^2}{E_{ij}}} and output as a matrix.
#'
#' @param data A data set, where the observed count vector and time vector (if variable) are the last two columns
#' @param mod A StagedTree model fit to the data set, as produced by pceg() or zipceg().
#' @param stages A logical value indicating whether the calculations should be based on stages (TRUE) or leaves (FALSE).
#' @param limit An integer where the number of event counts greater than or equal to this integer are grouped together.
#' @param min_exp An integer specifying the minimum expected count necessary for its chi square contribution to be considered.
#' @param zip A logical value indicating whether the model specified is zero-inflated (TRUE) or not (FALSE).
#' @param dec_place An integer value detailing how many decimal places the outputs should be rounded to. If NA, no rounding will occur.
#'
#' @return A list of three matrices and three numeric values. The first matrix is the observed count matrix, the second is the expected count matrix, and the third is the chi-squared contribution matrix. All three matrices are in ascending order of rate and risk probability. The first numeric value is the sum of the chi square contributions. The second numeric value is the degrees of freedom, the product of the amount of stages for risk probabilities and rates and the upper limit on counts. The third numeric value is the p-value for a Chi-squared distribution with the observed Chi-squared statistics at the calculated degrees of freedom.
#' @export
#'
#' @examples
#' mod1<-pceg(knee_pain_obs,2,TRUE,TRUE)
#' chi_sq_calculator(knee_pain_obs,mod1,zip=FALSE)
#'
#' mod2<-zipceg(knee_pain_obs,"nlm",variable_time=TRUE)
#' chi_sq_calculator(knee_pain_obs,mod2)
chi_sq_calculator<-function(data,mod,stages = TRUE, limit=4,min_exp=5,zip=FALSE, dec_place = NA){

  poisson_response<-mod$event.tree$poisson.response
  remove_risk_free<-mod$remove.risk.free.edges
  variable_time<-mod$event.tree$variable.time

  if(limit<=0){
    stop("Please input a valid count limit greater than 0")
  }

  if(min_exp<0){
    stop("Please input a non-negative minimum expectation")
  }

  if(!poisson_response & variable_time){
    stop("Variable Time Requires Poisson Response")
  }

  if(!poisson_response & zip){
    stop("Zero Inflated Poisson Requires Poisson Response")
  }

  #below is copied into quantile_band - if this changes, so should that
  #make into own function

  path_details<-refactored_tree_matrix(data,poisson_response,variable_time)
  data_use<-path_details$data_use
  n<-path_details$num_var
  p<-path_details$p
  tree<-path_details$tree_matrix

  #a lot of the below is copied into quantile_band - if this changes, so should that

  posterior<-mod$posterior.expectation
  stage_struct<-mod$stage.structure

  n1<-n+1*poisson_response +1*zip

  rates<-parameter_extractor(stage_struct,posterior,n1,poisson_response,remove_risk_free)
  probs<-parameter_extractor(stage_struct,posterior,n1-1,poisson_response,remove_risk_free)

  if(zip & !(remove_risk_free)){
    len<-seq(2,2*p,by=2)
    rates<-rates[len]
  }

  probs<-probs[,2]

  obs.mat<-matrix(nrow=p,ncol=limit+1)
  exp.mat<-matrix(nrow=p,ncol=limit+1)

  min.exp.mat<-matrix(nrow=p,ncol=limit+1)

  x<-c(0:(limit-1))

  unique_rates<-sort(unique(rates),decreasing = FALSE)
  if(zip){
    unique_probs<-sort(unique(probs),decreasing = FALSE)
  }else{
    unique_probs<-1
  }

  n_stages_rates<-length(unique_rates)
  n_stages_probs<-length(unique_probs)

  n_stages_combo<-n_stages_rates*n_stages_probs

  stage_combo_index<-matrix(data=NA,nrow=p,ncol=2)

  for(k in 1:p){
    v<-tree[k,c(1:n)]
    ind<-which(row.match(data_use[,1:n],v)==1 )
    lambda<-rates[k]
    stage_combo_index[k,2]<-which(unique_rates==lambda)

    if(zip){
      prop<-probs[k]
      stage_combo_index[k,1]<-which(unique_probs==prop)
    }else{
      prop<-1
      stage_combo_index[k,1]<-1
    }

    y<-data_use[ind,n+1]

    if(variable_time){
      t<-data_use[ind,n+2]
    }else{
      t<-rep(1,length(ind))
    }

    obs.mat[k,1:limit]<-sapply(x,counter,v=y)
    exp.mat[k,1:limit]<-colSums(sapply(x,f,p=prop,lambda=lambda,t=t))

    obs.mat[k,limit+1]<-length(which(y >= limit))
    exp.mat[k,limit+1]<-length(ind)-sum(exp.mat[k,1:limit])

  }

  if(stages){

    exp.mat.new<-matrix(data=NA,nrow=n_stages_combo,ncol=limit+1)
    obs.mat.new<-matrix(data=NA,nrow=n_stages_combo,ncol=limit+1)

    count<-1
    for(i in 1:n_stages_probs){
      for(j in 1:n_stages_rates){
        combo<-c(i,j)
        ind_temp<-which((!colSums(t(stage_combo_index)!=combo))==TRUE)

        if(length(ind_temp)>1){
          exp.mat.new[count,]<-colSums(exp.mat[ind_temp,])
          obs.mat.new[count,]<-colSums(obs.mat[ind_temp,])
          count<-count+1
        }else if(length(ind_temp)==1){
          exp.mat.new[count,]<-exp.mat[ind_temp,]
          obs.mat.new[count,]<-obs.mat[ind_temp,]
          count<-count+1
        }
      }
    }
    exp.mat<-exp.mat.new
    obs.mat<-obs.mat.new
  }

  df<-n_stages_combo*limit

  ind.na<-which(is.na(obs.mat[,1]))
  if(length(ind.na)>0){
  obs.mat<-obs.mat[-ind.na,]
  exp.mat<-exp.mat[-ind.na,]
  n_stages_combo<-n_stages_combo-length(ind.na)
  }

  min.exp.mat<-exp.mat<min_exp

  if(sum(min.exp.mat)>0){
    warning("Some expected counts below minimum expected count allowed - these have been excluded from the calculation, consider decreasing count limit")
  }

  chi.mat<-(exp.mat-obs.mat)^2/exp.mat

  chi.mat[min.exp.mat]<-0

  if(!is.na(dec_place)){
    exp.mat<-round(exp.mat,dec_place)
    chi.mat<-round(chi.mat,dec_place)
  }

  v1<-paste0(limit,"+")
  v<-c(x,v1)

  colnames(obs.mat)<-v
  colnames(exp.mat)<-v
  colnames(chi.mat)<-v

  if(stages){
    u<-paste0("Stage",1:n_stages_combo)
  }else{
    u<-paste0("Leaf",1:p)
  }

  rownames(obs.mat)<-u
  rownames(exp.mat)<-u
  rownames(chi.mat)<-u

  chi_sq<-sum(chi.mat)
  p_value<-pchisq(chi_sq,df)

  return(list(obs=obs.mat,exp=exp.mat,chi.mat=chi.mat,chi_sq = chi_sq,df=df,p_value=p_value))

}

