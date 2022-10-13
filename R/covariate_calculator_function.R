#' The Covariate Calculator Function
#'
#' @param data
#' @param variable_time
#'
#' @return
#' @export
#'
#' @examples
covariate_calculator<-function(data,variable_time){

  n<-dim(data)[2] - 1 - 1*variable_time#if there are variable times, they will be an extra column
  data_levels<-sapply(data[,1:n],nlevels)
  data_levels_list<-lapply(data[,1:n],levels)
  p<-prod(data_levels)
  l1<-lapply(data_levels,vec_from)
  l2<-lapply(data_levels,vec_to)
  X<-expand.grid(l1) #this is the different covariate combinations
  X<-cbind(X,nonzero=apply(X,1,nonzero_checker))
  X=X[order(X$nonzero,decreasing = TRUE),] #ignore the far left hand side labels
  # X_fact<-X #will define as factor in later for loop

  alpha_ind<-expand.grid(l2) #other way around because we will solve the alphas this way
  alpha_ind<-cbind(alpha_ind,nonzero=apply(alpha_ind,1,nonzero_checker)) #changed because what if the variable has value 2 (3 levels)
  alpha_ind<-cbind(alpha_ind,ind=rep(0,p)) #placeholder for the label, now need to calculate the labels
  #can maybe remove now
  possible_alpha<-c(0:(p-1))

  #also need to create a matrix, like X, but that matches the tree
  #this is to incorporate merging later
  #X and T are ordered for calculation, Z is ordered for display
  #do it using expand_grid because this varies the first variable slowest,
  #not fastest like expand.grid, so that it matches the tree

  Z<-lapply(data_levels,vec_from)#changed because if they don't have same number of levels for each, we get a list.
  #so we start with a list
  tree_matrix<-Z[[1]]
  for(i in 2:n){
    tree_matrix<-expand_grid(tree_matrix,Z[[i]]) #changed to list too
    colnames(tree_matrix)[1:i]<-colnames(X)[1:i]
  }
  tree_matrix<-as.data.frame(tree_matrix)

  #need to create new matrices to use in case variables not binary

  extra_levels<-data_levels-rep(2,n)
  n_diff<-sum(extra_levels)
  n1=n+n_diff

  data_refactor<-path_refactor(data,n,data_levels)

  for(i in 1:n){
    copy<-extra_levels[i]+1
    name_base<-names(tree_matrix)[i]
    name_levels<-data_levels_list[[i]]

    X_temp<-X[i]
    X_temp_keep<-as.data.frame(keeper(X_temp[,1],data_levels[i]-1))

    tree_temp<-tree_matrix[i]
    tree_temp_keep<-as.data.frame(keeper(tree_temp[,1],data_levels[i]-1))

    alpha_ind_temp<-alpha_ind[i]
    alpha_ind_temp_keep<-as.data.frame(keeper(alpha_ind_temp[,1],data_levels[i]-1))

    data_temp<-data_refactor[,i]
    data_temp_keep<-factor(keeper(data_temp,data_levels[i]-1),levels=c(1,0))
    data_temp_keep<-as.data.frame(data_temp_keep)

    if(copy>1){
      names(X_temp_keep)<-paste(name_base,name_levels[1],sep="")
      names(tree_temp_keep)<-paste(name_base,name_levels[1],sep="")
      names(alpha_ind_temp_keep)<-paste(name_base,name_levels[1],sep="")
      names(data_temp_keep)<-paste(name_base,name_levels[1],sep="")
    }else{
      names(X_temp_keep)<-name_base
      names(tree_temp_keep)<-name_base
      names(alpha_ind_temp_keep)<-name_base
      names(data_temp_keep)<-name_base
    }

    if(i == 1){
      new_X<-X_temp_keep
      new_tree_matrix<-tree_temp_keep
      new_alpha_ind<-alpha_ind_temp_keep
      new_data<-data_temp_keep
    }else{
      new_X<-cbind(new_X,X_temp_keep)
      new_tree_matrix<-cbind(new_tree_matrix,tree_temp_keep)
      new_alpha_ind<-cbind(new_alpha_ind,alpha_ind_temp_keep)
      new_data<-cbind(new_data,data_temp_keep)
    }

    if(copy > 1){
      for(j in 2:copy){
        X_temp_keep<-as.data.frame(keeper(X_temp[,1],data_levels[i]-j))
        tree_temp_keep<-as.data.frame(keeper(tree_temp[,1],data_levels[i]-j))
        alpha_ind_temp_keep<-as.data.frame(keeper(alpha_ind_temp[,1],data_levels[i]-j))
        data_temp_keep<-factor(keeper(data_temp,data_levels[i]-j),levels=c(1,0))
        data_temp_keep<-as.data.frame(data_temp_keep)

        names(X_temp_keep)<-paste(name_base,name_levels[j],sep="")
        names(tree_temp_keep)<-paste(name_base,name_levels[j],sep="")
        names(alpha_ind_temp_keep)<-paste(name_base,name_levels[j],sep="")
        names(data_temp_keep)<-paste(name_base,name_levels[j],sep="")

        new_X<-cbind(new_X,X_temp_keep)
        new_tree_matrix<-cbind(new_tree_matrix,tree_temp_keep)
        new_alpha_ind<-cbind(new_alpha_ind,alpha_ind_temp_keep)
        new_data<-cbind(new_data,data_temp_keep)
      }
    }
  }

  X<-cbind(new_X,nonzero=apply(new_X,1,nonzero_checker)) #keep around for later
  X_fact<-X

  alpha_ind<-cbind(new_alpha_ind,nonzero=apply(new_alpha_ind,1,nonzero_checker))
  alpha_ind=alpha_ind[order(alpha_ind$nonzero,decreasing=FALSE),] #ignore the far left hand side labels
  alpha_ind<-cbind(alpha_ind,ind=c(0:(p-1)))

  tree_matrix<-new_tree_matrix

  component_names<-c(p)

  possible_names<-names(alpha_ind)[1:n1]

  data_use<-cbind(new_data,data[,-c(1:n)])

  for(j in 1:p){
    ind<-which(alpha_ind[j,1:n1]==1)
    m=length(ind)
    if(m == 0){
      component_names[j]="Intercept"
    }else if (m ==  1){
      component_names[j]=possible_names[ind]
    }else{
      component_names[j]=do.call(paste, c(as.list(possible_names[ind]), sep = ":"))
    }
  }

  for(i in 1:n1){
    X_fact[,i]=factor(X_fact[,i],levels=c(1:0)) #need to be a factor for later
    tree_matrix[,i]=factor(tree_matrix[,i],levels=c(1:0))
  }


  # alpha_ind=alpha_ind[order(alpha_ind$ind),] #ignore the far left hand side labels

  #n

  alpha_cols<-list()
  alpha_cols[[1]]=NA
  for(i in 2:p){
    alpha_cols[[i]]=which(alpha_ind[i,1:n1]==1) #the nonzero covariates for each alpha component
  }

  T<-list() #T will be the paths to the lambdas, not in same order of the tree
  for(i in 1:p){
    T[[i]]<-c(0)
    k=2
    for(j in 2:p){
      v<-alpha_cols[[j]]
      l<-length(v)
      check<-rep(1,l)
      if(sum(X[i,v]-check)==0){ #This is why we keep X around
        T[[i]][k]=j-1
        k=k+1
      }
    }
  }


  return(list(data_use=data_use,num_var=n,num_binary_var=n1,num_alpha=p,data_levels=data_levels,
              covariate_matrix=X_fact,alpha_ind=alpha_ind,lambda_paths=T,tree_matrix=tree_matrix,component_names=component_names)) #return X_fact, not x

}
