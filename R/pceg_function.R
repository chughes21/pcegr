#' The Poisson Score
#'
#' This function calculates the marginal likelihood for Poisson responses.
#'
#' @param data_sample A list where each element of the list contains the observed count vector and time vector for a unique unfolding of the process.
#' @param gamma_prior A list of numeric vectors where the numeric vectors are the prior hyperparameters for the Gamma distribution.
#'
#' @return A numeric value of the log marginal likelihood for the data.
#'
poisson_score<-function(data_sample, gamma_prior){
  p<-length(data_sample)
  post_a<-0
  post_b<-0
  y<-c()
  t<-0
  l<-0
  for(i in 1:p){
    y=data_sample[[i]][,1]
    t=data_sample[[i]][,2]
    prior_a<-gamma_prior[[i]][1]
    prior_b<-gamma_prior[[i]][2]

    if(prior_a > 0 & prior_b > 0){
      post_a<-prior_a + sum(y)
      post_b<-prior_b+sum(t)
      C<-0
      n<-length(y)
      if(n > 0){
        for(j in 1:n){
          C = C + y[j]*log(t[j]) - log(factorial(y[j]))
        }
      }
      l= l + prior_a*log(prior_b)-post_a*log(post_b) + lgamma(post_a)-lgamma(prior_a) + C
      #if prior is equal to posterior they should cancel, so could also exclude
    }

  }
  return(l)
}


#' The log Bayes factor Between nested models
#'
#' This function calculates the change in log marginal likelihood of a model when two stages are merged into one.
#'
#' @param sample_sum A list where each element is an integer sum of the observed event counts and observation time for each stage.
#' @param prior A list of numeric vectors where the numeric vectors are the hyperparameters of the Gamma prior distribution for each stage.
#' @param stage1 An integer value specifying the first stage to be merged.
#' @param stage2 An integer value specifying the second stage to be merged.
#'
#' @return A numeric value specifying the change in log marginal likelihood after the two stages are merged.
#'
bayes_factor<-function(sample_sum, prior, stage1, stage2){
  prior_a_1<-prior[[stage1]][1]
  prior_b_1<-prior[[stage1]][2]
  prior_a_2<-prior[[stage2]][1]
  prior_b_2<-prior[[stage2]][2]

  if(prior_a_2>0 & prior_b_2>0){

    y1<-sample_sum[[stage1]][1]
    y2<-sample_sum[[stage2]][1]

    t1<-sample_sum[[stage1]][2]
    t2<-sample_sum[[stage2]][2]

    post_a_1<-prior_a_1+y1
    post_b_1<-prior_b_1+t1
    post_a_2<-prior_a_2+y2
    post_b_2<-prior_b_2+t2

    new_prior_a<-prior_a_1+prior_a_2
    new_prior_b<-prior_b_1+prior_b_2

    new_y<-y1+y2
    new_t<-t1+t2

    new_post_a<-new_prior_a+new_y
    new_post_b<-new_prior_b+new_t

    diff<--(prior_a_1*log(prior_b_1)+prior_a_2*log(prior_b_2)-(post_a_1*log(post_b_1)+post_a_2*log(post_b_2))+
              lgamma(post_a_1)+lgamma(post_a_2)-lgamma(prior_a_1)-lgamma(prior_a_2)-
              new_prior_a*log(new_prior_b)+new_post_a*log(new_post_b)-lgamma(new_post_a)+lgamma(new_prior_a)) #note the minus sign
  }else{diff<-0}
  return(diff)
}

#' The PCEG function
#'
#' This function fits a PCEG model to a data set.
#'
#' This function takes a data set, a prior structure, and various other inputs to fit a PCEG model. It is flexible enough to fit a vanilla CEG, and perform variable discretisation methods.
#'
#' @param data A data set where the observed response vector and time vector (if applicable and variable) are the last two columns.
#' @param equivsize A numeric value specifying the equivalent sample size for the prior, a measure of confidence in the prior.
#' @param poisson_response A logical value indicating whether the response variable is Poisson (TRUE) or categorical (FALSE).
#' @param variable_time A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE), if applicable.
#' @param zip A logical value indicating whether the model specified is zero-inflated (TRUE) or not (FALSE).
#' @param gamma_alpha A numeric value for the shape hyperparameter of the Gamma prior for the Poisson rate, if applicable.
#' @param gamma_beta A numeric value for the rate hyperparameter of the Gamma prior for the Poisson rate.
#' @param structural_zero A logical value indicating whether zero counts in the data set should be considered as structural (TRUE) or sampling (FALSE) for the setting of the prior.
#' @param var_disc An integer value specifying which variable to discretise. If 0, no discretisation is necessary.
#' @param disc_length An integer value specifying how many neighbours can be searched over for the purposes of variable discetisation. If 0, all other possible stages may be merged over.
#' @param restrict A logical value indicating whether variable discretisation should be restricted to stages with the same unfolding of the process (TRUE) or not (FALSE).
#' @param mirror A logical value indicating whether variable discretisation should be equivalent across each unfolding of the process (TRUE) or not (FALSE).
#' @param cat_limit An integer value specifying the minimum number of categories to the variable can be discretised to. If 0, there is no minimum number of categories.
#' @param collapse A logical value indicating whether, when a restricted discretisation occurs, the final output should be compact (TRUE) or not (FALSE).
#'
#' @return A list specifying a PCEG model. The list contains: the prior for the final model, the data for the final model, the number of nodes at each level of the tree, the stage numbers for the final model,
#' the stage structure for the final model, the vector of likelihoods after each merging, the details of the stages merged at each step,
#' the comparison set of stages left to be merged, and the log marginal likelihood of the final model.
#' @export
#'
#' @examples
#' mod<-pceg(knee_pain_obs,2,TRUE,TRUE)
#' mod$result
pceg<-function(data ,equivsize=3,  poisson_response = FALSE, variable_time = FALSE, zip=FALSE,
                        gamma_alpha =1, gamma_beta = 2,structural_zero = FALSE, var_disc = 0, disc_length = 0,
                        restrict = FALSE, mirror = FALSE, cat_limit=0, collapse = FALSE){

  exampledata<-data

  no_disc_length_spec<-(disc_length == 0) #if discretisation is chosen but no interval length given
  cat_limit_ind<-(cat_limit>0)

  if(!poisson_response & variable_time){
    stop("Variable time requires a Poisson response")
  }

  if(!poisson_response & zip){
    stop("Zero Inflated Poisson requires Poisson response")
  }

  if(var_disc == 0 & (disc_length > 0 | restrict | mirror)){
    stop("Variable discretisation requires a variable to discretise")
  }

  if(equivsize <= 0){
    stop("Equivalent Sample size should be greater than 0")
  }

  if((var_disc > 0) & collapse & !(restrict | mirror )){
    stop("Collapsed results only possible when merging is restricted")
  }


  numbvariables<-dim(exampledata)[2] - 1*variable_time #this means if there is a variable time, it must come last
  numbcat <-sapply(exampledata[,1:numbvariables],FUN=nlevels) #number of categories at each level
  numb<-c(1,cumprod(numbcat[1:(numbvariables-1)])) #number of nodes at each level
  numb_out<-numb

  if(var_disc > numbvariables){
    stop("Variable to discretise does not exist")
  }

  situation_data<-c(list(rbind(table(exampledata[,1]))))
  for (i in 2:(numbvariables-1*poisson_response)){
    for (j in 1:numb[i]){
      situation_data<-c(situation_data ,list(rbind(ftable(exampledata[,1:i])[j,])))
    }
  }

  nv<-numbvariables-1*poisson_response

  if(any(numbcat[1:nv]==0)){
    stop("Covariate and/or non Poisson response should be a factor - check again")
  }

  prior_output<-prior_set(situation_data,structural_zero,nv,numbcat[1:nv],numb[1:nv],equivsize)
  prior<-prior_output$prior
  no_merge_ind<-prior_output$no_merge

  prior_response<-c()
  if(poisson_response ){

    if(length(gamma_alpha) != length(gamma_beta)){
      warning("Alpha and Beta hyperparameter inputs not same length")
    }


    if(length(gamma_alpha) == 1){
      prior_alpha=c(rep(gamma_alpha,numb[numbvariables]))
    }else{
      prior_alpha=gamma_alpha
      if(length(gamma_alpha)!=numb[numbvariables]){
        warning("Alpha hyperparameter incorrect length")
      }
    }
    if(length(gamma_beta) == 1){
      prior_beta=c(rep(gamma_beta,numb[numbvariables]))
    }else{
      prior_beta=gamma_beta
      if(length(gamma_alpha)!=numb[numbvariables]){
        warning("Beta hyperparameter incorrect length")
      }}

    for(j in 1:numb[numbvariables]){
      prior_response<-c(prior_response ,list(cbind(prior_alpha[j],prior_beta[j])))
    }
    prior<-c(prior,prior_response)
  }

  #moved labelling up here

  labelling <-c()
  for (k in 1:(numbvariables -1)){
    label <-c(1,rep("NA",sum(numb[1:k]) -1))
    label<-c(label ,rep(levels(exampledata[,k]),numb[k]))
    if (k<(numbvariables -1)){
      for (i in (k+1):(numbvariables -1)){
        label<-c(label ,rep(levels(exampledata[,k]),each=numb[i+1]/numb[k+1],numb[k+1]
                            /numbcat[k]))
      }
    }
    labelling<-cbind(labelling ,label)
  }

  #for poisson response, need to compute summary statistics
  #total sum of y and time (which will be just count, when time is not variable) for each path
  #the two for loops are really slow
  if(poisson_response){
    data_sum<-data.frame(cbind(labelling[(sum(numb[-numbvariables])+1):(sum(numb[-numbvariables])+numb[numbvariables]),],y=0,t=0))
    data_counts<-vector(mode = "list",length = numb[numbvariables])
    for(i in 1:(numbvariables-1)){
      data_sum[,i]=factor(data_sum[,i],levels=levels(exampledata[,i]))
    }
    data_sum$y=as.numeric(data_sum$y)
    data_sum$t=as.numeric(data_sum$t)
    colnames(data_sum)[1:(numbvariables-1)]<-colnames(exampledata)[1:(numbvariables-1)]
    for(k in 1:numb[numbvariables]){
      y_temp=0
      t_temp=0
      v<-data_sum[k,1:(numbvariables-1)]
      ind<-which(row.match(exampledata[,1:(numbvariables-1)],v)==1 )
      data_counts[[k]]<-exampledata[ind,-(1:(numbvariables-1))]
      if(!variable_time){
        len<-length(data_counts[[k]])
        data_counts[[k]]<-cbind(data_counts[[k]],rep(1,len))
      }
      data_sum$y[k]=sum(data_counts[[k]][,1])
      data_sum$t[k]=sum(data_counts[[k]][,2])
    }
  }

  data<-situation_data

  if(poisson_response){
    for (j in 1:numb[numbvariables]){
      data<-c(data,list(cbind(data_sum$y[j],data_sum$t[j])))
    }
  }

  # List of the stages that can be merged in the first step
  comparisonset<-c()
  for (i in 2:numbvariables){
    ind_temp<-c((sum(numb[1:(i-1)])+1):(sum(numb[1:i])))
    check_temp<-which(ind_temp %in% no_merge_ind)
    if(length(check_temp)>0 & structural_zero){
      ind_temp<-ind_temp[-check_temp]
    }
    comparisonset <-c(comparisonset ,list(ind_temp))
  }
  labelling <-c()
  for (k in 1:(numbvariables -1)){
    label <-c(1,rep("NA",sum(numb[1:k]) -1))
    label<-c(label ,rep(levels(exampledata[,k]),numb[k]))
    if (k<(numbvariables -1)){
      for (i in (k+1):(numbvariables -1)){
        label<-c(label ,rep(levels(exampledata[,k]),each=numb[i+1]/numb[k+1],numb[k+1]
                            /numbcat[k]))
      }
    }
    labelling<-cbind(labelling ,label)
  }
  mergedlist <-c()
  for (i in 1:sum(numb)){
    mergedlist<-c(mergedlist ,list(labelling[i,]))
  }
  merged1<-c() #this will be the recording of all mergings
  merged2<-c() #this will be the temporary recording to add to merged1
  merged_out<-c() #this will be the final output to avoid repetition

  lik<-0
  if(poisson_response){numb=numb[-numbvariables]}
  for( i in 1: sum(numb)){
    alpha<-unlist(prior[i])
    N<-unlist(data[i])
    lik<-lik+sum(lgamma(alpha+N)-lgamma(alpha))+sum(lgamma(sum(alpha))-lgamma(
      sum(alpha+N)))
  }

  if(poisson_response){
    lik=lik+poisson_score(data_counts,prior_response)
  }

  score<-c(lik)
  #At each step we calculate the difference between the current CEG and the CEG
  #in which two stages in the current comparison set have been merged.
  #We go through every possible combination of stages that can be merged . k is
  # an index for the comparisonset we are in ,
  # and i and j the position of the stages within the comparison set .

  n_iter<-length(comparisonset)

  for (k in 1:n_iter){

    diff.end<-1 #to start the algorithm

    final_ind<-k == n_iter #final level indicator

    if(poisson_response & final_ind){
      poisson<-TRUE #poisson response indicator
    }else{poisson<-FALSE}

    disc<-(k == var_disc)*1 #discretisation indicator

    num_var<-numbcat[k]
    groups<-length(comparisonset[[k]])/num_var #the number of covariate combinations excluding the most recent covariate
    intervals<-seq(1,num_var*(groups+1),by=num_var) #the index for start of each previous covariate combination
    intervals_init<-comparisonset[[k]][intervals[1:numb[k]]] #the nodes itself
    intervals_init<-c(intervals_init,max(comparisonset[[k]])+1) #include an extra element for purpose of functions
    total_cats<-rep(num_var,numb[k]) #the total number of categories

    if(mirror & disc){
      comparisonset1=comparisonset[[k]][c(1:num_var)] #only look at the first set of situations if mirroring
      restrict = TRUE #automatically restrict if we're gonna mirror
      if(no_disc_length_spec){
        disc_length = num_var #if mirroring and no search length given, then can search over all nodes with the same covariate combo
      }
    }else{comparisonset1 = comparisonset[[k]]}


    if(zip & final_ind){
      len = length(comparisonset1)/2 #only half of the leaves are at risk
      ind<-seq(1,length.out = len, by = 2)
      no_risk_ind<-comparisonset1[ind] #these are the risk free leaf indices

      for(i in 2:len){

        result2 = bayes_factor(data, prior,comparisonset1[1], no_risk_ind[i]) #note just regular data and prior, due to what the stages are

        prior[[comparisonset1[1]]]<-prior[[comparisonset1[1]]]+ prior[[no_risk_ind[i]]] #combine priors
        prior[[ no_risk_ind[i]]] <-cbind(NA ,NA) #overwrite priors
        data[[comparisonset1[1]]]<-data[[comparisonset1[1]]]+data[[no_risk_ind[i]]] #combine data
        data[[ no_risk_ind[i]]] <-cbind(NA,NA) #overwrite data
        comparisonset[[k]]<-comparisonset[[k]][-(which(comparisonset
                                                       [[k]]== no_risk_ind[i]))]
        comparisonset1<-comparisonset1[-which(comparisonset1==no_risk_ind[i])]
        mergedlist[[comparisonset1[1]]]<-cbind(mergedlist[[comparisonset1[1]]],mergedlist[[no_risk_ind[i]]])
        mergedlist [[ no_risk_ind[i]]] <-cbind(NA ,NA)
        lik<-lik+result2
        score<-c(score ,lik)
        merged2<-c(comparisonset1[1] ,no_risk_ind[i] ,k)
        merged1<-cbind(merged1 ,merged2)
      }

      #remove the need to compare no risk and risk
      no_risk_stage<-comparisonset1[1]
      comparisonset1<-comparisonset1[-1]
      comparisonset[[k]]<-comparisonset[[k]][-1]
    }

    comparisonset_split<-split(comparisonset[[k]],ceiling(seq_along(comparisonset[[k]])/num_var))

    stop_searching<-c()

    if(restrict & cat_limit_ind & disc){
      min_category<-rep(cat_limit,numb[[k]])
    }else{min_category<-numeric(numb[[k]])}


    while(diff.end >0){ #We stop when no positive difference is obtained by merging two stages

      difference<-0 #keep an eye on this, need to reset each run

      num_nodes<-length(comparisonset1)

      if(num_nodes >1){ # can only merge if more than one stage in the comparisonset, might need to change this

        #first need to decide what function can search over, based on discretisation restrictions

        if(cat_limit_ind & restrict & disc){
          no_include<-c()
          no_include_ind<-c()
          if(length(stop_searching)>0){ #find the nodes to exclude from search based on restrictions
            for(j in 1:length(stop_searching)){
              no_include<-c(no_include,comparisonset_split[[stop_searching[j]]])
            }
            no_include_ind=which(comparisonset1 %in% no_include)
          }
          if(length(no_include_ind)>0){ #If there are nodes to exclude, do it here
            comparisonset2<-comparisonset1[-no_include_ind]
            comparisonset3<-comparisonset[[k]][-no_include_ind]
            possible_i<-c(1:(length(comparisonset2)-1))
          }else{ #if there are no nodes to exclude
            comparisonset2<-comparisonset1
            comparisonset3<-comparisonset[[k]]
            possible_i<-c(1:(length(comparisonset2)-1))
          }
        }else{ #if there's no restrictions at all
          comparisonset2<-comparisonset1
          comparisonset3<-comparisonset[[k]]
          possible_i<-c(1:(length(comparisonset2)-1))
        }

        #comparisonset2 will be comparisonset1 after the necessary nodes are excluded
        #comparisonset3 will be comparisonset[[k]] after the necessary nodes are excluded

        if(disc & restrict){
          for(j in 2:(numb[[k]]+1)){
            intervals[j] = intervals[j-1]+total_cats[j-1]
          }
        }

        diff_cats<-numeric(numb[[k]]+1)

        if(length(stop_searching)>0){
          for(j in stop_searching){
            if( j <= numb[[k]]){
              diff_cats[c((j+1):(numb[[k]]+1))]=diff_cats[c((j+1):(numb[[k]]+1))]+total_cats[j]
            }
          }

          intervals_temp<-intervals-diff_cats
          intervals_temp<-intervals_temp[-stop_searching]
          intervals_init_temp<-intervals_init[-stop_searching]
        }else{
          intervals_temp<-intervals
          intervals_init_temp<-intervals_init
        }

        for (i in possible_i){
          compare1_temp <-comparisonset2[i]
          num_nodes = length(comparisonset2) #might not need this but I checked and its necessary, unless more changes are made
          disc_length_temp<-disc_length
          if(no_disc_length_spec | (disc_length > num_nodes)){disc_length_temp = num_nodes} #might not need this
          #they may not be required due to "restrict" being in place anyway
          if(restrict & disc){
            interval_num<-findInterval(compare1_temp,intervals_init_temp) #added temp
            max_val<-intervals_temp[interval_num + 1]-1 #added temp
          }else{max_val<-length(comparisonset[[k]])}

          if(i >= max_val){
            int_set<-c()
          }else if(disc & (disc_length_temp > 0) ){
            interval_num<-findInterval(compare1_temp,intervals_init) #this is a different interval_num to use, no temp for interval_init
            max_search<-min(disc_length_temp,total_cats[interval_num]-min_category[interval_num])
            if(max_search > 0){
              int_set<-seq(i+1,min(i+max_search,max_val),by=1)
            }else{int_set<-c()}
          }else{int_set <- seq(i+1,max_val,by=1)}

          if(disc&mirror){
            i<-seq(i,i+num_nodes*(groups-1),by = num_nodes)
          }

          for (j in int_set){
            #to compare
            if(mirror){
              h = seq(j,j+num_nodes*(groups-1),by = num_nodes)
            }else{h = j}
            compare1 <-comparisonset3[i]
            compare2 <-comparisonset3[h]
            #we calculate the difference between
            # the CEG where two stages are merged
            result<-0
            merged_temp<-c()
            L<-length(compare1)
            for(l in 1:L){
              if(poisson){
                result = result + bayes_factor(data, prior,compare1[l], compare2[l] ) #note just regular data and prior, due to what the stages are
              }
              else{ result=result+lgamma(sum(prior[[compare1[l]]]+prior[[compare2[l]]]))-lgamma(sum(prior[[
                compare1[l]]]+data[[compare1[l]]]+prior[[compare2[l]]]+data[[compare2[l]]]))+
                sum(lgamma(prior[[compare1[l]]]+data[[compare1[l]]]+prior[[compare2[l]]]+data[[
                  compare2[l]]]))-sum(lgamma(prior[[compare1[l]]]+prior[[compare2[l]]]))-
                # and the CEG where the two stages are not merged
                (lgamma(sum(prior[[compare1[l]]]))-lgamma(sum(prior[[compare1[l]]]+data[[compare1[l]
                ]]))+sum(lgamma(prior[[compare1[l]]]+data[[compare1[l]]]))-
                  sum(lgamma(prior[[compare1[l]]]))+lgamma(sum(prior[[compare2[l]]]))-lgamma(sum(
                    prior[[compare2[l]]]+data[[compare2[l]]]))+
                  sum(lgamma(prior[[compare2[l]]]+data[[compare2[l]]]))-sum(lgamma(prior[[compare2[l]]])
                  ) )
              }

              merged_temp<-cbind(merged_temp,c(compare1[l],compare2[l],k))
            }
            #   print(result) use for checking
            #if the resulting difference is greater than the current difference then we replace it
            if (result > difference){
              difference<-result
              merged<-merged_temp
            }
          }
        }
      }
      diff.end<-difference
      #We update our priorlist , datalist and comparisonset to obtain the priorlist ,
      #datalist and comparisonlist for C_ {1}
      if(diff.end >0){

        l1<-which(comparisonset1==merged[1,1])
        l2<-which(comparisonset1==merged[2,1])
        cat_minus<-length(merged[1,])*(l2-l1) #how many nodes will be in between
        if(mirror & disc){
          total_cats=total_cats-cat_minus/numb[[k]]
        }else if(restrict & disc) {
          cat_interval<-findInterval(merged[1],intervals_init)
          total_cats[cat_interval]=total_cats[cat_interval]-cat_minus
          if(disc & (total_cats[cat_interval]==min_category[cat_interval])){
            stop_searching = c(stop_searching,cat_interval) #when number of categories for a group hits limit, stop searching that group
          }
        }

        for(l in 1:L){
          if(disc){
            cluster_index<-seq(from=merged[1,l],to=merged[2,l],by = 1)
          }else{cluster_index=c(merged[1,l],merged[2,l])}

          n=length(cluster_index)

          for(i in 2:n){
            if(!(cluster_index[i] %in% merged1[2,])){
              if(poisson){
                result2 = bayes_factor(data, prior,cluster_index[1], cluster_index[i] ) #note just regular data and prior, due to what the stages are
              }
              else{
                result2<-lgamma(sum(prior[[cluster_index[1]]]+prior[[cluster_index[i]]]))-lgamma(sum(prior[[
                  cluster_index[1]]]+data[[cluster_index[1]]]+prior[[cluster_index[i]]]+data[[cluster_index[i]]]))+
                  sum(lgamma(prior[[cluster_index[1]]]+data[[cluster_index[1]]]+prior[[cluster_index[i]]]+data[[
                    cluster_index[i]]]))-sum(lgamma(prior[[cluster_index[1]]]+prior[[cluster_index[i]]]))-
                  # and the CEG where the two stages are not merged
                  (lgamma(sum(prior[[cluster_index[1]]]))-lgamma(sum(prior[[cluster_index[1]]]+data[[cluster_index[1]
                  ]]))+sum(lgamma(prior[[cluster_index[1]]]+data[[cluster_index[1]]]))-
                    sum(lgamma(prior[[cluster_index[1]]]))+lgamma(sum(prior[[cluster_index[i]]]))-lgamma(sum(
                      prior[[cluster_index[i]]]+data[[cluster_index[i]]]))+
                    sum(lgamma(prior[[cluster_index[i]]]+data[[cluster_index[i]]]))-sum(lgamma(prior[[cluster_index[i]]])
                    ) )
              }
              prior[[cluster_index[1]]]<-prior[[cluster_index [1]]]+ prior[[cluster_index [i]]]
              prior[[ cluster_index [i]]] <-cbind(NA ,NA)
              data[[cluster_index[1]]]<-data[[cluster_index[1]]]+data[[cluster_index[i]]]
              data[[ cluster_index[i]]] <-cbind(NA,NA)
              comparisonset[[k]]<-comparisonset[[k]][-(which(comparisonset
                                                             [[k]]== cluster_index[i]))]
              if(l==1){ #not sure if necessary but i think it is
                comparisonset1<-comparisonset1[-(which(comparisonset1== cluster_index[i]))] #not sure if necessary but i think it is
              }
              mergedlist[[cluster_index[1]]]<-cbind(mergedlist[[cluster_index[1]]],mergedlist[[cluster_index
                                                                                               [i]]])
              mergedlist [[ cluster_index[i]]] <-cbind(NA ,NA)
              lik<-lik+result2
              score<-c(score ,lik)
              merged2<-c(cluster_index[1] ,cluster_index[i] ,k)
              merged1<-cbind(merged1 ,merged2)
            }else{merged2<-c(cluster_index[1] ,cluster_index[i] ,k)
            merged1<-cbind(merged1 ,merged2)}
          }
          merged_out<-cbind(merged_out,c(cluster_index[1],cluster_index[n],k)) #added this to remove double counting in output
        }
      }
      if(disc & cat_limit_ind){
        if(sum(abs(min_category-total_cats))==0){ break} #break when minimum categories have been reached
      }
    }

    if(zip & k == length(comparisonset)){
      comparisonset1<-c(no_risk_stage,comparisonset1)
      comparisonset[[k]]<-c(no_risk_stage,comparisonset[[k]])
    }

  }
  # Output : stages of the finest partition to be combined to obtain the most probable CEG structure
  stages<-c(1)
  for (i in 2:numbvariables){
    stages<-c(stages ,comparisonset[[i-1]])
  }

  if(restrict & (var_disc > 1) & collapse){
    stages_disc<-comparisonset[[var_disc]]
    mergedlist_disc<-mergedlist[stages_disc]
    for(j in 1:length(mergedlist_disc)){
      v<-as.matrix(mergedlist_disc[[j]])
      n<-dim(v)[2]
      if(n>1){
      start<-mergedlist_disc[[j]][var_disc,1]
      end<-mergedlist_disc[[j]][var_disc,n]
      vec.out<-mergedlist_disc[[j]][,1]
      state.out<-paste0(start," - ",end)
      }else{
      start<-mergedlist_disc[[j]][var_disc]
      vec.out<-mergedlist_disc[[j]]
      state.out<-paste0(start)
      }
      vec.out[var_disc]<-state.out
      mergedlist_disc[[j]]<-vec.out
    }
    mergedlist[stages_disc]<-mergedlist_disc
  }

  result<-mergedlist[stages]
  newlist <-list(prior=prior ,data=data ,numb=numb_out,stages=stages ,result=result ,score=score ,
                 merged=merged_out ,comparisonset=comparisonset ,lik=lik)
  return(newlist)
}
