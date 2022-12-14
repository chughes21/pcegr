#' Knee Pain Parameters
#'
#' 2 data sets, one with the conditional transition probabilities for the covariate situations (cov_probs), and one for the risk probabilities and rates for the leaves (leaf_params).
#'
#' @format 2 data frames with several variables
#' \describe{
#' \item{Statement}{A string detailing a specific conditional probability or dependence.}
#' \item{Stage}{A string detailing which stage the statement belongs to, for the covariates only.}
#' \item{Situation}{A string detailing which situations are in the stage, for the covariates only.}
#' \item{Prob}{A numeric specifying the conditional transition probability for the stage, for the covariates only.}
#' \item{Leaves}{A string detailing which leaves are associated with the specific conditional dependence statement.}
#' \item{Risk_Prob}{A numeric specifying the risk probability probability for the leaf stage.}
#' \item{Rate}{A numeric specifying the rate for the leaf stage.}
#' }


cov_statement<-c("P(A = 1)","P(W = 0 | A = 1 )","P(W = 0 | A = 0 )","P(H = 0 | A = 1 )","P(H = 0 | A = 0 )")
cov_stages<-c("u1","u2","u3","u4","u5")
sit_stages<-c("s0","s1","s2","s3,s4","s5, s6")

probs<-c(0.5, 0.3, 0.8, 0.2, 0.9)

cov_probs<-data.frame(Statement=cov_statement, Stage=cov_stages, Situations=sit_stages, Prob = probs)

cov_combo<-c("A = 1, W = 1, H = 1","A = 1, W = 1, H = 0","W = 0, H = 1","W = 0, H = 0","A = 0, W = 1, H = 1","A = 0, W = 1, H = 0" )
leaf_statement<-paste0("K | ",cov_combo)
leaves<-c("l1","l2","l3, l7","l4, l8","l5","l6")
prob<-c(0.2,0.6, 0.7, 0.9, 0.45, 0.8)
lambda<-c(0.5,2.5,3,5,1.5,4)

leaf_params<-data.frame(Statement = leaf_statement, Leaves = leaves, Risk_Prob = prob, Rate = lambda)

usethis::use_data(cov_probs, overwrite = TRUE)
usethis::use_data(leaf_params, overwrite = TRUE)
