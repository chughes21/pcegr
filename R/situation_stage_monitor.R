situation_stage_monitor<-function(data,mod,var,stage = NULL,zip=FALSE){

  poisson_response<-mod$event.tree$poisson.response
  remove_risk_free<-mod$remove.risk.free.edges
  variable_time<-mod$event.tree$variable.time
  data_summ<-mod$data.summary
  prior<-mod$prior.distribution
  post<-mod$posterior.expectation

  num_var<-mod$event.tree$num.variable

  alpha=sum(prior[[1]])

  if(var>num_var){
    stop("Choice of var must be less than  or equal to number of variables in tree")
  }

  if(var<=0){
    stop("Choice of var must be positive")
  }

  if(var==1){
    stop("Situation stage monitor only well defined for variables after first")
  }

  #create a model which is totally saturated

  sat_mod<-pceg(data,alpha,poisson_response,variable_time,zip,remove_risk_free,saturated = c(2:num_var))

}
