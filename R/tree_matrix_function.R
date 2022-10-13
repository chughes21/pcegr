#' The Tree Matrix Function
#'
#' A function to refactor the data and create a refactored Matrix that matches the covariates in the tree.
#'
#' @param data
#' @param variable_time
#'
#' @return A list containing the refactored data, the number of levels in the tree, the number of leaves in the tree, a vector of the levels for each variable, and the refactored tree matrix.
#' @export
#'
#' @examples refactored_tree_matrix(knee_pain_obs,TRUE)
refactored_tree_matrix<-function(data,variable_time){

  n<-dim(data)[2] - 1 - 1*variable_time#if there are variable times, they will be an extra column
  data_levels<-sapply(data[,1:n],nlevels)
  p<-prod(data_levels)
  names.mat<-expand.grid(l1)

  Z<-lapply(data_levels,vec_from)
  tree_matrix<-Z[[1]]
  for(i in 2:n){
    tree_matrix<-expand_grid(tree_matrix,Z[[i]])
    colnames(tree_matrix)[1:i]<-colnames(names.mat)[1:i]
  }
  tree_matrix<-as.data.frame(tree_matrix)

  data_refactor<-path_refactor(data,n,data_levels)

  return(list(data_use=data_use,num_var=n,p=p,data_levels=data_levels, tree_matrix=tree_matrix))

}