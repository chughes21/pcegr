#' Knee Pain Data
#'
#' A data set of the instances of knee pain for individuals.over a variable observation time.
#'
#' @format A data frame with 4 covariates, a response and an observed time:
#' \describe{
#' \item{Age}{A binary variable indicating whether an individual is young or old.}
#' \item{Weight}{A binary variable indicating if an individual is overweight or normal weight.}
#' \item{History}{A binary variable indicating if an individual has a history of knee injury or not.}
#' \item{State}{A binary variable indicating if an individual is at risk of knee pain or not.}
#' \item{y}{An integer recording the number of instances of knee pain suffered by the individual over the observation time}
#' \item{t}{The observation time in years.}
#'
#' }

set.seed(1234)

N<-10000

p_weight<-c(0.3,0.8)
p_history<-c(0.2,0.9)
lambda_pain_overweight<-c(3,5)
p_pain_overweight<-c(0.7,0.9)
lambda_pain_no_overweight<-matrix(data=c(0.5,2.5,1.5,4),nrow=2,ncol=2,byrow=TRUE)
p_pain_no_overweight<-matrix(data=c(0.2,0.6,0.45,0.8),nrow=2,ncol=2,byrow=TRUE)

data<-data.frame(age=factor(rep("Young",N),levels=c("Young","Old")),weight = factor(rep("Normal",N),levels=c("Normal","Over")),history = factor(rep("No",N),levels=c("No","Yes")),state=factor(rep("No Risk",N),levels=c("No Risk","Risk")),y=numeric(N),t=rep(1,N))

poss_interval_max<-2

min_t<-1/365

data$t = runif(N)*poss_interval_max
data$t[which(data$t < min_t)]<-min_t

for(i in 1:N){

  randoms<-runif(4,0,1)

  if(randoms[1]<=0.5){
    data$age[i]<-"Young"
    age_ind <- 1
  }else{
    data$age[i]<-"Old"
    age_ind <- 2}

  if(randoms[2]<=p_weight[age_ind]){
    data$weight[i]<-"Over"
    weight_ind<-2
  }else{
    data$weight[i]<-"Normal"
    weight_ind<-1
  }

  if(randoms[3]<=p_history[age_ind]){
    data$history[i]<-"Yes"
    history_ind<-2
  }else{
    data$history[i]<-"No"
    history_ind<-1
  }

  if(weight_ind == 2){
    lambda = lambda_pain_overweight[history_ind]
    state_check<-(randoms[4]<=p_pain_overweight[history_ind])*1 # 1 means at risk, 0 means not
    data$y[i]<-rpois(1,lambda*state_check*data$t[i])
    if(state_check){
      data$state[i]<-"Risk"
    }
  }else{
    lambda = lambda_pain_no_overweight[age_ind,history_ind]
    state_check<-(randoms[4]<=p_pain_no_overweight[age_ind,history_ind])*1
    data$y[i]<-rpois(1,lambda*state_check*data$t[i])
    if(state_check){
      data$state[i]<-"Risk"
    }
  }
}

knee_pain<-data

knee_pain

usethis::use_data(knee_pain, overwrite = TRUE)


