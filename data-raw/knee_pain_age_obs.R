#' Knee Pain Age Observed Data
#'
#' A data set of the instances of knee pain for individuals over a variable observation time without the risk state included. Age is now an integer.
#'
#' @format A data frame with 3 covariates, a response and an observed time:
#' \describe{
#' \item{Age}{An integer variable indicating the age of an individual between 1 and 80.}
#' \item{Weight}{A binary variable indicating if an individual is overweight or normal weight.}
#' \item{History}{A binary variable indicating if an individual has a history of knee injury or not.}
#' \item{y}{An integer recording the number of instances of knee pain suffered by the individual over the observation time}
#' \item{t}{The observation time in years.}
#'
#' }

knee_pain_age_obs<-knee_pain_age[,-4]
knee_pain_age_obs

usethis::use_data(knee_pain_age_obs, overwrite = TRUE)
