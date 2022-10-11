#' Knee Pain Observed Data
#'
#' A data set of the instances of knee pain for individuals.over a variable observation time without the risk state included.
#'
#' @format A data frame with 3 covariates, a response and an observed time:
#' \describe{
#' \item{Age}{A binary variable indicating whether an individual is young or old.}
#' \item{Weight}{A binary variable indicating if an individual is overweight or normal weight.}
#' \item{History}{A binary variable indicating if an individual has a history of knee injury or not.}
#' \item{y}{An integer recording the number of instances of knee pain suffered by the individual over the observation time}
#' \item{t}{The observation time in years.}
#'
#' }

knee_pain_obs<-knee_pain[,-4]
knee_pain_obs

usethis::use_data(knee_pain_obs, overwrite = TRUE)
