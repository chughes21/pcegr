#' The Matrix to Vector Function
#'
#' Convert a matrix into a vector, row by row.
#'
#' @param mat A numeric matrix with m rows and n columns
#'
#' @return A numeric vector with mn elements, where each row of the matrix is concatenated to the next.
matrix_to_vector<-function(mat){
  return(c(t(mat)))
}


#' The Tree Repeater Function
#'
#' Repeat stage probabilities for each path they are a constituent of.
#'
#' @param vec A numeric vector of stage probabilities.
#' @param num_leaves An integer count of the number of leaves in the StagedTree object.
#'
#' @return A numeric vector the length of num_leaves where each element is the corresponding stage probability for that part of the path.
tree_repeater<-function(vec,num_leaves){
  i=num_leaves/length(vec)
  return(rep(vec,each=i))
}

#' The Path Probability Extractor
#'
#' Calculate the probability of traversing each path in a StagedTree object. For a StagedTree with a Poisson response, it calculates the probability of traversing the path until the leaf where the response is recorded.
#'
#' @param mod A StagedTree object.
#' @param precision An integer value specifying how many decimal places numerical checks should be done to.
#'
#' @return A list containing: a stages matrix, outlining the individual probability for each transition along a path, and a paths matrix, calculating the cumulative product of transition probabilities such that the final column is the probability of traversing each path.
#' @export
#'
#' @examples
#' mod<-pceg(knee_pain_obs,2,TRUE,TRUE)
#' path_prob_extractor(mod)
path_prob_extractor<-function(mod,precision = 5){

  post<-mod$posterior.expectation
  stage.struct<-mod$stage.structure
  poisson.response<-mod$event.tree$poisson.response
  remove.risk.free<-mod$remove.risk.free.edges
  num_var<-mod$event.tree$num.variable-1*poisson.response
  data_levels<-mod$event.tree$num.category[1:num_var]

  params<-c(list(post[[1]]),lapply(c(2:num_var),parameter_extractor,stage_struct = stage.struct, posterior = post, poisson_response = poisson.response,remove_risk_free = remove.risk.free))

  edges<-cumprod(data_levels)
  p<-max(edges)
  leaves<-p/edges

  params_vec<-lapply(params,matrix_to_vector)

  prob_mat<-sapply(params_vec,tree_repeater,num_leaves=p)
  path_prob_mat<-t(apply(prob_mat,1,cumprod))

  colnames(prob_mat)[1:num_var]<-paste0("Var",1:num_var)
  colnames(path_prob_mat)[1:num_var]<-paste0("Var",1:num_var)

  if(!poisson.response){
    colnames(prob_mat)[num_var]<-"Y"
    colnames(path_prob_mat)[num_var]<-"Y"
  }

  if(round(sum(path_prob_mat[,num_var]),precision)!=1){
    stop("Path probabilities don't sum to 1")
  }

  if(round(sum(colSums(path_prob_mat)-p/edges),precision)!=0){
    stop("Situation counts don't match")
  }

  return(list(stages = prob_mat, paths = path_prob_mat))
}


#' The Kullback-Leibler Divergence function
#'
#' Calculate the Kullback-Leibler (KL) divergence between a StagedTree object and a reference model, either another StagedTree object or a vector of path probabilities that sum to one.
#'
#' @param new_mod A StagedTree object.
#' @param old_mod A StagedTree object. Must be included if old_prob is not.
#' @param old_prob A numeric vector of path probabilities that sums to one. Must be included if old_mod is not.
#' @param precision An integer value specifying how many decimal places numerical checks should be done to.
#'
#' @return A numeric value for the KL divergence between the StagedTree object and reference model.
#' @export
#'
#' @examples
#' mod<-pceg(knee_pain_obs,2,TRUE,TRUE)
#' mod2<-pceg(knee_pain_obs,2,TRUE,TRUE,indep=c(2:3))
#' kl_divergence(mod,old_mod=mod2)
#' equiv_prob<-rep(1/8,8)
#' kl_divergence(mod,old_prob=equiv_prob)
kl_divergence<-function(new_mod,old_mod = NULL, old_prob = NULL,precision=5){

  old_mod_check<-length(old_mod)>0
  old_prob_check<-length(old_prob)>0

  poisson.response<-new_mod$event.tree$poisson.response
  num_var<-new_mod$event.tree$num.variable-1*poisson.response

  if(!old_mod_check & !old_prob_check){
    stop("At least one of a model or probability vector must be included")
  }

  if(old_mod_check & old_prob_check){
    stop("Only one of a model or probability vector must be included")
  }

  if(old_prob_check){
    if(round(sum(old_prob),precision)!=1){
      stop("Probability vector must sum to 1")
    }
  }

  if(old_mod_check){
    if(old_mod$event.tree$poisson.response != poisson.response){
      stop("Both models should have same response variable")
    }

    if(new_mod$event.tree$num.variable != old_mod$event.tree$num.variable){
      stop("Both models should have same number of variables")
    }
  }

  new_path_prob<-path_prob_extractor(new_mod,precision)$paths[,num_var]

  if(old_mod_check){
    old_path_prob<-path_prob_extractor(old_mod,precision)$paths[,num_var]
  }
  else{old_path_prob<-old_prob}

  p<-length(new_path_prob)
  if(p!=length(old_path_prob)){
    stop("Distributions do not act on same number of atoms (paths)")
  }
  ratio<-new_path_prob/old_path_prob
  lr<-log(ratio)
  plr<-new_path_prob*log(ratio)
  return(sum(plr))
}

