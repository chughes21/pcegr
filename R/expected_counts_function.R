#' The Merged Stage Finder
#'
#' A function to find which stages have merged with a given unique stage from the merging details of a CEG.
#'
#' @param unique_stage An integer value
#' @param M An integer matrix detailing the merging present in a CEG.
#'
#' @return An integer vector detailing the stages that have merged with the unique stage given
#'
#' @examples
merged_stage_finder<-function(unique_stage,M){
  M_unique<-as.matrix(M[,which(M[1,]==unique_stage)])
  stages_merged<-c(unique_stage,M_unique[2,])
  stages_merged<-as.vector(sort(stages_merged,decreasing=FALSE))
}

#' The Merged List Extractor
#'
#' A function to find a list of the stages which have merged given the merging details of a CEG.
#'
#' @param merged A matrix specifying the merging present in a CEG.
#'
#' @return A list of intger vectors detailing the stages which have merged.
#'
#' @examples
merged_list_extractor<-function(merged){
  level<-max(merged[3,])
  merged_stages<-merged[,which(merged[3,]==level)]
  stages_top<-unique(merged_stages[1,])
  stages_bottom<-unique(merged_stages[2,])

  stage_ref<-stages_bottom[which(sapply(stages_bottom,FUN=function(x) vec_in(x,stages_top))==TRUE)]

  num_ref<-length(stage_ref)

  merged_stages_adj<-merged_stages

  for(i in 1:num_ref){

    col_ref_top<-which(merged_stages[1,]==stage_ref[i])
    col_ref_bottom<-which(merged_stages[2,]==stage_ref[i])
    new_top_stage<-merged_stages[1,col_ref_bottom]
    merged_stages_adj[1,col_ref_top]<-new_top_stage

  }

  stages_top_adj<-as.vector(sort(unique(merged_stages_adj[1,]),decreasing=FALSE))

  output<-lapply(stages_top_adj,function(x) merged_stage_finder(x,merged_stages_adj))
  return(output)
}

#' The vec_in function
#'
#' A function to check whether a given value is in a vector
#'
#' @param x A numeric value
#' @param v A numeric vector
#'
#' @return A logical value indicating whether the given value was in the vector (TRUE) or not (FALSE)
#'
#' @examples vec_in(7,c(1:9))
#' vec_in(10,c(1:9))
vec_in<-function(x,v){
  check<-x %in% v
  return(check)
}

#' The expected count function.
#'
#' A function to test the goodness of fit of a pceg model to a data set.
#'
#' For each leaf stage in a pceg model created by [pceg()], this function calculates the observed event counts up to some integer limit, and then calculates
#' the expected event counts for that leaf based on the model estimates of the parameters for the Poisson distribution. Then, for each leaf $i$ and each event count $j$,
#' a chi-square contribution is calculated using the formula $\frac{(O_{ij}-E_{ij})^2}{E_{ij}}$ and output as a matrix.
#'
#' @param data A data set, where the observed count vector and time vector (if variable) are the last two columns
#' @param ceg A ceg model fit to the data set, as produced by pceg().
#' @param limit An integer where the number of event counts greater than or equal to this integer are grouped together.
#' @param poisson_response A logical value indicating whether the response variable is Poisson (TRUE) or categorical (FALSE).
#' @param variable_time A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE).
#' @param posterior A logical value indicating whether the estimates of the parameters used should be the posterior estimate (TRUE) or sample estimate (FALSE).
#' @param zip A logical value indicating whether the model specified is zero-inflated (TRUE) or not (FALSE).
#'
#' @return A list of three matrices and a numeric value. The first matrix is the observed count matrix, the second is the expected count matrix, and the third is the chi-squared contribution matrix.The numeric value is the sum of the chi square contributions.
#' @export
#'
#' @examples
expected_count_calculator<-function(data,ceg,limit=4,poisson_response=TRUE,variable_time=TRUE,posterior = TRUE, zip=TRUE){

  if(!poisson_response & variable_time){
    stop("Variable Poisson Time Requires Poisson Response")
  }

  if(!poisson_response & zip){
    stop("Zero Inflated Poisson Requires Poisson Response")
  }

  path_details<-refactored_tree_matrix(data,variable_time)
  data_use<-path_details$data_use
  n<-path_details$num_var
  p<-path_details$p
  tree<-path_details$tree_matrix
  data_levels<-path_details$data_levels
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
  ind<-ceg$stages
  ind_probs<-ind[ind>=start_probs & ind<=end_probs]
  ind_rates<-ind[ind>=start_rates & ind<=end_rates]

  rates<-value_extractor(data,ceg,level_rel_final=0,poisson_response=poisson_response,variable_time=variable_time,posterior=posterior,zip=zip)
  probs<-value_extractor(data,ceg,level_rel_final=-1,poisson_response,variable_time,posterior=posterior,zip) #if zip=FALSE, we'll ignore this anyway

  rates<-rates[,dim(rates)[2]] #if start inputting true values, this may be incorrect
  probs<-probs[,dim(probs)[2]]

  seq_prob<-seq(start_probs,end_probs)
  if(zip){
    seq_rate<-seq(start_rates+1,end_rates,by=2)
  }else{
    seq_rate<-seq(start_rates,end_rates)
  }

  if(zip==TRUE){
    tree<-cbind(tree,prob_stage=seq_prob,rate_stage=seq_rate)
  }else{
    tree<-cbind(tree,rate_stage=seq_rate)
  }

  merged<-ceg$merged
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

  obs.mat<-matrix(nrow=p,ncol=limit+1)
  exp.mat<-matrix(nrow=p,ncol=limit+1)

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

    prob_sum<-0

    for(j in 0:(limit-1)){
      obs.mat[k,j+1]<-length(which(y == j))
      total_prob<-sum(f(prop,lambda,j,t))
      exp.mat[k,j+1]<-total_prob
      prob_sum<-prob_sum+total_prob

    }

    obs.mat[k,limit+1]<-length(which(y >= limit))
    exp.mat[k,limit+1]<-length(ind)-sum(exp.mat[k,1:limit])

    chi.mat<-(exp.mat-obs.mat)^2/exp.mat
  }

  return(list(obs=obs.mat,exp=exp.mat,chi.mat=chi.mat,chi_sq = sum(chi.mat)))

}
