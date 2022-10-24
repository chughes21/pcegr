#' Knee Pain Age Data
#'
#' A data set of the instances of knee pain for individuals over a variable observation time.  Age is now an integer.
#'
#' @format A data frame with 4 covariates, a response and an observed time:
#' \describe{
#' \item{Age}{An integer variable indicating the age of an individual between 1 and 80.}
#' \item{Weight}{A binary variable indicating if an individual is overweight or normal weight.}
#' \item{History}{A binary variable indicating if an individual has a history of knee injury or not.}
#' \item{State}{A binary variable indicating if an individual is at risk of knee pain or not.}
#' \item{y}{An integer recording the number of instances of knee pain suffered by the individual over the observation time}
#' \item{t}{The observation time in years.}
#'
#' }
knee_pain_age<-knee_pain

young_limit<-40
old_limit<-80

ind_old<-which(knee_pain_age$age == "Old")
ind_young<-which(knee_pain_age$age == "Young")

knee_pain_age$age<-numeric(length(knee_pain_age$age))

knee_pain_age$age[ind_old]<-sample(c((young_limit+1):old_limit),length(ind_old),replace=TRUE)
knee_pain_age$age[ind_young]<-sample(c(1:young_limit),length(ind_young),replace=TRUE)

knee_pain_age

usethis::use_data(knee_pain_age, overwrite = TRUE)
