#' The Situation Stage monitor
#'
#' A function to check if any situations in a stage are different to the rest.
#'
#' @param data A data set where the observed response vector and time vector (if applicable and variable) are the last two columns.
#' @param mod A StagedTree model fit to the data set, as produced by [pceg()].
#' @param var An integer value detailing which variable, if applicable.
#' @param stage An integer value detailing either which stage in the tree to investigate, or which stage for a chosen variable.
#' @param cat An integer value detailing which category of the applicable variable should be used as reference.
#' @param ci A numeric value between 0 and 1 specifying the percentage confidence for the highest density interval.
#' @param N An integer specifying the number of iterations for the empirical functions used.
#' @param zip A logical value indicating whether the model specified is zero-inflated (TRUE) or not (FALSE).
#'
#' @return A plot of confidence intervals versus observed when each situation in the stage is removed.
#' @export
#'
#' @examples
#' mod<-pceg(knee_pain_obs,2)
#' situation_stage_monitor(knee_pain_obs,mod,var=4,stage=3) #the third stage for variable 4
#' situation_stage_monitor(knee_pain_obs,mod,stage=10) #the 10th stage in the tree, which is the same as the above
situation_stage_monitor<-function(data,mod,var = NULL,stage = NULL,cat=1,N=10000,ci=0.95,zip=FALSE){

  poisson_response<-mod$event.tree$poisson.response
  remove_risk_free<-mod$remove.risk.free.edges
  variable_time<-mod$event.tree$variable.time
  prior<-mod$prior.distribution
  stage_struct<-mod$stage.structure
  levels<-mod$event.tree$num.category
  labels<-mod$event.tree$label.category
  stages<-mod$stages
  num_var<-mod$event.tree$num.variable
  situations_level<-mod$event.tree$num.situation

  num_stages<-length(stages)

  situation_starts<-c(0,cumsum(situations_level[1:(num_var-1)]))+1
  situation_ends<-cumsum(situations_level)

  if(length(situation_starts)!=length(situation_ends)){
    stop("Situation start and end points don't match")
  }

  situation_matrix<-cbind(situation_starts,situation_ends)

  stages_level<-sapply(stages,findInterval,vec=situation_starts) #the level for each stage
  num_stages_level<-sapply(1:num_var,counter,stages_level) #the number of stages by level


  if(is.null(var) & is.null(stage)){
    stop("At least one of variable or stage should be provided")
  }

  #function currently only works if stage is a single value, can try make it more functional with vectors

  if(is.null(var)){
    if(stage>length(stages)){
      stop("Stage specified does not exist")
    }
    var<-stages_level[stage]
    stage_index<-stage
  }else{
    if(is.null(stage)){
      stage<-1 #default just look at first stage for a variable
    }else{
      if(stage>num_stages_level[var]){
        stop("Stage chosen greater than number of stages for chosen level - check that var hasn't been specified also")
      }
    }

    stage_starts<-c(0,cumsum(num_stages_level[1:(num_var-1)]))

    stage_index=stage_starts[var]+stage

  }

  if(var>num_var){
    stop("Choice of var must be less than  or equal to number of variables in tree")
  }

  if(var<=0){
    stop("Choice of var must be positive")
  }

  if(var==1){
    stop("Situation stage monitor only well defined for variables after first")
  }

  if(cat<=0){
    stop("Choice of category must be positive")
  }

  if(cat>levels[var]){
    stop("Choice of category must be less than or equal to number of categories for chosen variable")
  }

  #create a model which is totally saturated

  alpha=sum(prior[[1]])
  sat_mod<-pceg(data,alpha,poisson_response,variable_time,zip,remove_risk_free,saturated = c(2:num_var))

  stage_situation_index<-stages[stage_index]

  relative_stage_index<-stage_situation_index-situation_starts[var]+1

  situation_index<-as.vector(unlist(stage_struct[[var]][relative_stage_index]))
  overall_situation_index<-situation_starts[var]+situation_index-1

  data_stage<-sat_mod$data.summary[[var]][situation_index,]
  prior_stage<-sat_mod$prior.distribution[[var]][situation_index,]

  if(levels[var]>2){
    data_summ<-cbind(data_stage[,cat],rowSums(data_stage[,-cat]))
    prior_summ<-cbind(prior_stage[,cat],rowSums(prior_stage[,-cat]))
  }else{
    data_summ<-cbind(data_stage[,cat],data_stage[,-cat])
    prior_summ<-cbind(prior_stage[,cat],prior_stage[,-cat])
  }

  m<-length(situation_index)

  output<-data.frame(situation=overall_situation_index,post=numeric(m),ci_low=numeric(m),ci_high=numeric(m),obs=numeric(m))

  if(m<=1){
    stop("Chosen stage must not be singleton")
  }

  for(j in 1:m){

    if(poisson_response & var==num_var){

    output$obs[j]<-data_summ[j,1]/data_summ[j,2]

    data_excl<-matrix(data_summ[-j,],ncol=2)
    prior_excl<-matrix(prior_summ[-j,],ncol=2)

    post_excl<-colSums(data_excl+prior_excl)
    shape=post_excl[1]
    scale=1/post_excl[2]

    marg_post<-distribution_gamma(N,shape=shape,scale=scale)

    output$post[j]<-shape*scale
    output$ci_low[j]<-hdi(marg_post,ci=ci)$CI_low
    output$ci_high[j]<-hdi(marg_post,ci=ci)$CI_high

    }else{

   output$obs[j]<-data_summ[j,1]/sum(data_summ[j,])

   data_excl<-matrix(data_summ[-j,],ncol=2)
   prior_excl<-matrix(prior_summ[-j,],ncol=2)

   post_excl<-data_excl+prior_excl
   marg_post<-distribution_beta(N,shape1=post_excl[1],shape2=post_excl[2])

   output$post[j]<-post_excl[1]/sum(post_excl)
   output$ci_low[j]<-hdi(marg_post,ci=ci)$CI_low
   output$ci_high[j]<-hdi(marg_post,ci=ci)$CI_high

  }

  }

  output$situation=as.factor(output$situation)

  ggplot(data=output,mapping=aes(x=situation))+geom_point(mapping=aes(y=post),col="red")+
    ggtitle(paste0("Situation Stage Monitor for Stage ",stage_index))+
    geom_point(mapping=aes(y=obs),col="blue",pch=17)+
    geom_linerange(mapping=aes(ymin=ci_low,ymax=ci_high),col="red")
}
