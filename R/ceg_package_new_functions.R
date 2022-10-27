#' The OAHC Definer
#'
#' This function creates an object of class "OAHC" from a pceg model.
#'
#' @param cluster A list detailing the stage structure for a level
#' @param score A numeric value of the score contribution from this level
#'
#' @return An object of class "OAHC" with a stage structure and score from a level of the tree
#'
#' @examples
oahc_definer<-function(cluster,score=0){
  return(new("OAHC",score=score,cluster=cluster))
}

#' The Stage Structure Function
#'
#' This function takes a pceg model and creates a stage structure that is compatible with the "Stratified.staged.tree" S4 class.
#'
#' @param mod A list detailing the output from [pceg()]
#' @param score A numeric vector detailing the change in score for each level of the tree
#' @param zip A logical value indicating whether the model specified is zero-inflated (TRUE) or not (FALSE).
#'
#' @return A list specifying the stage structure from the output of [pceg()] in a way that is compatible with the "Stratified.staged.tree" S4 class.
#'
#' @examples
#' mod1<-pceg(knee_pain_obs,2,TRUE,TRUE)
#' stage_structure(mod1)
stage_structure<-function(mod,score=0,zip=FALSE){
  output<-list()
  M<-mod$merged
  comparisonset<-mod$comparisonset
  numb<-mod$numb


  for(i in 1:length(numb)){

    if(score!=0){
      lik<-score[i]
    }else{lik<-0}

    cluster<-vector(mode="list",length=numb[i])
    m<-as.matrix(M[,which(M[3,]==(i-1))])

    if(i>1){
      comp<-comparisonset[[i-1]]
      start<-sum(numb[1:(i-1)])
      comp<-comp-start
    }else{
      start<-1
      comp<-1
    }

    m[1:2,]<-m[1:2,]-start

    if(length(m)>0){
      check<-merged_list_extractor(m)
    }else{
      check<-NA
    }

    zip_ind<-(i==length(numb))& zip

    if(zip_ind){
      no_risk_ind<-seq(from=1,to=numb[i]-1,by=2)
    }else{no_risk_ind<-NA}

    k<-1

    for(j in 1:numb[[i]]){

      if((j==1)&zip_ind ){
        cluster[[j]]<-no_risk_ind
      } else if((j %in% no_risk_ind[-1])&zip_ind){
        cluster[[j]]<-NA
      }else if(all(is.na(check))){
        cluster[[j]]<-j
      }else if(!(j %in% comp)){
        cluster[[j]]<-NA
      }else if(!(j%in%m[1:2,])){
        cluster[[j]]<-j
      }else{
        cluster[[j]]<-as.vector(unlist(check[k]))
        k<-k+1
      }
     # if((j==1)&zip ){
     #   cluster[[j]]<-no_risk_ind
    #  }

    # if(all(is.na(check))){
     #   cluster[[j]]=j
    #  }else{
    #    if(!(j %in% comp)){
     #     cluster[[j]]<-NA
     #   } else if(!(j%in%m[1:2,])){
    #      if(zip & (j %in% no_risk_ind[-1])){
      #      cluster[[j]]<-NA
      #    }else{cluster[[j]]<-j}
      #  }else{
    #      cluster[[j]]<-as.vector(unlist(check[k]))
     #     k<-k+1
    #    }
   #   }
    }
    output[[i]]=oahc_definer(cluster,score=lik)
  }

  return(output)

}

#' The Output List Converter Function
#'
#' A function to transform the list of prior and data outputs into the S4 class versions of the prior and posterior.
#'
#' @param mod A list from the output of the [pcegr()] function.
#' @param poisson_response A logical value indicating whether the response variable is Poisson (TRUE) or categorical (FALSE).
#'
#' @return A list for the prior distribution and a list for the posterior distribution.
#'
#' @examples
#' pmod<-pceg(knee_pain_obs,2,TRUE,TRUE)
#' output_list_converter(pmod,TRUE)
output_list_converter<-function(mod,poisson_response=TRUE){
  prior<-mod$prior
  data<-mod$data
  numb<-mod$numb
  n<-length(numb)

  numb_temp<-c(1,numb)

  nodes_start<-cumsum(numb_temp)[1:n]
  nodes_end<-cumsum(numb)

  nodes<-nodes_end-nodes_start + 1

  prior_out<-vector(mode="list",length=n)
  data_out<-vector(mode="list",length=n)
  post_out<-vector(mode="list",length=n)

  for(i in 1:n){
    prior_mat_temp<-matrix(nrow=nodes[i],ncol=dim(prior[[nodes_start[i]]])[2])
    data_mat_temp<-matrix(nrow=nodes[i],ncol=dim(prior[[nodes_start[i]]])[2])

    nodes_temp<-seq(nodes_start[i],nodes_end[i],by=1)
    for(j in 1:nodes[i]){
      prior_mat_temp[j,]<-prior[[nodes_temp[j]]]
      data_mat_temp[j,]<-data[[nodes_temp[j]]]
    }
    prior_out[[i]]<-prior_mat_temp
    data_out[[i]]<-data_mat_temp

    sum_out<-prior_mat_temp+data_mat_temp

    if((i==n)&poisson_response){
      post_out[[i]]<-sum_out[,1]/sum_out[,2]
    }else{
      post_out[[i]]<-sum_out/rowSums(sum_out)
      }
    }
  return(list(prior=prior_out,posterior=post_out))
}



#' The Event Tree Creator
#'
#' A function to create an object of the "Stratified.event.tree" S4 class that is compatible with Poisson response variables.
#'
#' @param data A data set where the observed response vector and time vector (if applicable and variable) are the last two columns.
#' @param poisson_response A logical value indicating whether the response variable is Poisson (TRUE) or categorical (FALSE).
#' @param variable_time A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE), if applicable.
#'
#' @return An object of the class "Stratified.event.tree"
#' @export
#'
#' @examples
#' tree<-event.tree.creator(knee_pain_obs,TRUE,TRUE)
#' plot(tree)
event.tree.creator<-function(data,poisson_response=TRUE,variable_time=TRUE){

  if(!poisson_response & variable_time){
    stop("Poisson response needed for variable time")
  }

  n<-dim(data)[2]-1*poisson_response - 1*variable_time

  if(poisson_response){
    empty.resp<-factor(rep("y",length(data[,1])),levels=c("y"))
    data.final<-data.frame(data[,1:n],resp=empty.resp)
  }else{
    data.final<-data
  }

  tree<-set(data.final)
  return(tree)
}

#' The Staged Tree Creator
#'
#'#' A function to create an object of the "Stratified.event.tree" S4 class that is compatible with the outputs of the [pceg()] function.
#'
#' @param data A data set where the observed response vector and time vector (if applicable and variable) are the last two columns.
#' @param mod A list from the output of the [pceg()] function.
#' @param poisson_response A logical value indicating whether the response variable is Poisson (TRUE) or categorical (FALSE).
#' @param variable_time A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE), if applicable.
#' @param zip A logical value indicating whether the model specified is zero-inflated (TRUE) or not (FALSE).
#' @param score A numeric vector of the score contributions from each level
#'
#' @return An object of the class "Stratified.staged.tree"
#' @export
#'
#' @examples
#'  mod1<-pceg(knee_pain_obs,2,TRUE,TRUE)
#'  staged.tree<-staged.tree.creator(mod1)
#'  plot(staged.tree)
staged.tree.creator<-function(data,mod,poisson_response=TRUE,variable_time=TRUE,zip=FALSE,score=0){

  if(!poisson_response & variable_time){
    stop("Variable time requires a Poisson response")
  }

  if(!poisson_response & zip){
    stop("Zero Inflated Poisson requires Poisson response")
  }

  if(zip){
    n<-dim(data)[2]-1*poisson_response -1*variable_time
    state<-factor(rep("No Risk",length(data[,1])),levels<-c("No Risk", "Risk"))
    data.final<-data.frame(data[,1:n], State = state, data[,-(1:n)])
  }else{
    data.final<-data
  }
  event.tree<-event.tree.creator(data.final,poisson_response,variable_time)
  temp<-output_list_converter(mod,poisson_response)
  prior.struc<-temp$prior
  post.struc<-temp$posterior
  stage.struc<-stage_structure(mod,score,zip=zip)
  staged.tree<-new("Stratified.staged.tree", event.tree,
                   situation = list(), contingency.table = list(), stage.structure = stage.struc,
                   stage.probability = list(), prior.distribution=prior.struc, posterior.distribution=post.struc,
                   model.score=mod$lik)
  return(staged.tree)
}
