#' The Prior Set Function
#'
#' This function sets the prior for the situations in the tree.
#'
#' This function sets the prior for the situations in the tree which have categorical response variables, based on the equivalent sample size and how zero edge counts should be treated.
#'
#' @param sample A list containing the edge counts of the situations in a data set.
#' @param struc_zero A logical value indicating whether zero edge counts in the sample should be considered as structural zeroes (TRUE) or sampling zeroes (FALSE) for the purposes of setting the prior.
#' @param numbvariables An integer value specifying how many covariates are in the data set, or levels in the event tree.
#' @param numbcat An integer vector specifying the number of categories for each covariate in the data set, or outgoing edges at each level of the event tree.
#' @param numb An integer vector specifying the number of situations at each level of the event tree.
#' @param equivsize A numeric value specifying the equivalent sample size for the prior, a measure of confidence in the prior.
#'
#' @return A list containing the prior hyperparameters for each situation and its outgoing edges in the event tree.
#'
prior_set<-function(sample,struc_zero,numbvariables,numbcat,numb,equivsize){
  check<-unlist(lapply(sample,zero_checker)) #for each situation, check how many edges have a zero count
  node_starts<-c(1,cumsum(numb[-numbvariables])+1) #start of each level
  node_ends<-cumsum(numb) #end of each level
  check_sum<-sum(check[1:(node_starts[numbvariables]-1)]) # only focus on non-leaves
  n_leaves<-as.integer(numbcat[numbvariables]*numb[numbvariables]) #number of leaves
  last_level<-c(node_starts[numbvariables]:node_ends[numbvariables]) #actually this is the last level before leaves
  check_last_level<-check[last_level]
  prior<-c()

  leaves<-seq(from=node_ends[numbvariables]+1,by=1,length=n_leaves)

  leaf_list<-c()

  for(i in 1:numbvariables){
    for(j in 1:numb[i]){
      leaf_list<-c(leaf_list,list(rbind(rep((n_leaves/(numbcat[i]*numb[i])),numbcat[i]))))
    }
  }

  node_list<-c()
  for(i in 1:numbvariables){
    node_list<-c(node_list,list(node_starts[i]:node_ends[i])) #list of all non-leaf, divided by level
  }

  paths<-matrix(nrow=n_leaves,ncol=numbvariables+1)
  paths[,1]<-1
  paths[,numbvariables+1]<-leaves
  if(numbvariables>1){
  for(k in 2:numbvariables){
    paths[,k]=rep(node_list[[k]],each=prod(numbcat[k:numbvariables]))
    }
  }

  paths_edit<-paths
  leaf_list_edit<-c()

  if(struc_zero==FALSE || check_sum == 0){
    ind_non_leaf<-NA
    struc_zero_leaf<-NA
    count<-1
    for(i in 1:numbvariables){
      for(j in 1:numb[i]){
        prior<-c(prior ,list(rbind(rep(equivsize*leaf_list[[count]][1]/n_leaves,numbcat[i])))) #could also use leaf_vec
        count<-count+1
      }
    }
  }else if(struc_zero==TRUE & check_sum > 0){
    ind<-which(check>0)
    ind_non_leaf<-ind[which(ind<=node_starts[numbvariables]-1)] #nodes that don't enter leaves
    ind_leaf<-ind[which(ind>node_starts[numbvariables]-1)] #nodes that enter leaf, not an actual leaf
    sampling_zero_leaf<-last_level[which(check_last_level==1)] #sampling zeroes into leaves have only one edge 0
    struc_zero_leaf<-last_level[which(check_last_level==2)] #structural zeroes into leaves have both edges 0
    for(i in struc_zero_leaf){
      ind_temp<-which(paths_edit[,numbvariables]==i)
      paths_edit<-paths_edit[-ind_temp,]
    }
    n_leaves_edit<-dim(paths_edit)[1]
    count<-1
    for(i in 1:numbvariables){
      for(j in 1:numb[i]){
        v<-numeric(numbcat[i])
        paths_zoom<-paths[which(paths[,i]==count),]
        paths_zoom_edit<-paths_edit[which(paths_edit[,i]==count),]
        edges<-unique(paths_zoom[,i+1]) #what should be the true edges
        for(k in 1:numbcat[i]){
          v[k]=length(which(paths_zoom_edit[,i+1]==edges[k]))
        }
        leaf_list_edit<-c(leaf_list_edit,list(v))
        prior<-c(prior,list(equivsize*leaf_list_edit[[count]]/n_leaves_edit))
        count<-count+1
      }
    }

  }
  return(list(prior=prior,no_merge=c(ind_non_leaf,struc_zero_leaf)))
}
