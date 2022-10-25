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
#'
#' @return A list specifying the stage structure from the output of [pceg()] in a way that is compatible with the "Stratified.staged.tree" S4 class.
#'
#' @examples
#' mod1<-pceg(knee_pain_obs,2,TRUE,TRUE)
#' stage_structure(mod1)
stage_structure<-function(mod,score=0){
  output<-list()
  M<-mod$merged
  comparisonset<-mod$comparisonset
  numb<-mod$numb


  for(i in 1:length(numb)){

    if(score!=0){
      lik<-score[i]
    }else{lik<-0}

    cluster<-vector(mode="list",length=numb[i])
    m<-M[,which(M[3,]==(i-1))]

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

    k<-1

    for(j in 1:numb[[i]]){

      if(all(is.na(check))){
        cluster[[j]]=j
      }else{
        if(!(j %in% comp)){
          cluster[[j]]<-NA
        } else if(!(j%in%m[1:2,])){
          cluster[[j]]<-j
        }else{
          cluster[[j]]<-as.vector(unlist(check[k]))
          k<-k+1
        }
      }
    }
    output[[i]]=oahc_definer(cluster,score=lik)
  }

  return(output)

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
#' @param mod
#' @param poisson_response A logical value indicating whether the response variable is Poisson (TRUE) or categorical (FALSE).
#' @param variable_time A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE), if applicable.
#' @param score A numeric vector of the score contributions from each level
#'
#' @return An object of the class "Stratified.staged.tree"
#' @export
#'
#' @examples
#'  mod1<-pceg(knee_pain_obs,2,TRUE,TRUE)
#'  staged.tree<-staged.tree.creator(mod1)
#'  plot(staged.tree)
staged.tree.creator<-function(data,mod,poisson_response=TRUE,variable_time=TRUE,score=0){
  event.tree<-event.tree.creator(data,poisson_response,variable_time)
  stage.struc<-stage_structure(mod,score)
  staged.tree<-new("Stratified.staged.tree", event.tree,
                   situation = list(), contingency.table = list(), stage.structure = stage.struc,
                   stage.probability = list(), prior.distribution=list(), posterior.distribution=list(),
                   model.score=mod$lik)
  return(staged.tree)
}
