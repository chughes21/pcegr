% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zip_functions.R
\name{zipceg}
\alias{zipceg}
\title{The ZIPCEG function}
\usage{
zipceg(
  data,
  method = "Gibbs",
  iter = 10000,
  equivsize = 2,
  poisson_response = TRUE,
  variable_time = TRUE,
  remove_risk_free = FALSE,
  stoch_imputation = TRUE,
  gamma_alpha = 1,
  gamma_beta = 2,
  beta_c = 1,
  beta_d = 1,
  p_0 = 0.5,
  l_0 = 1,
  tol = 1e-10,
  var_disc = 0,
  disc_length = 0,
  restrict = FALSE,
  mirror = FALSE,
  cat_limit = 0
)
}
\arguments{
\item{data}{A data set where the observed count vector and time vector (if variable) are the last two columns.}

\item{method}{A character string indicating the method for parameter estimation. The character string can be an element of c("Gibbs","nlm","EM","mle","mm").}

\item{iter}{The number of iterations for the Gibbs sampler or Expectation-Maximisation algorithm.}

\item{equivsize}{A numeric value specifying the equivalent sample size for the prior, a measure of confidence in the prior.}

\item{poisson_response}{A logical value indicating whether the response variable is Poisson (TRUE) or categorical (FALSE).}

\item{variable_time}{A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE), if applicable.}

\item{remove_risk_free}{A logical value indicating whether the risk free leaves and edges should be removed (TRUE) or not (FALSE).}

\item{stoch_imputation}{A logical value indicating whether stochastic (TRUE) or deterministic (FALSE) imputation should be used.}

\item{gamma_alpha}{A numeric value for the shape hyperparameter of the Gamma prior for the Poisson rate, if applicable.}

\item{gamma_beta}{A numeric value for the rate hyperparameter of the Gamma prior for the Poisson rate.}

\item{beta_c}{A numeric value for the alpha hyperparameter of the Beta prior for the risk probability.}

\item{beta_d}{A numeric value for the beta hyperparameter of the Beta prior for the risk probability.}

\item{p_0}{A numeric initial value for the risk probability.}

\item{l_0}{A numeric initial value for the rate.}

\item{tol}{A numeric which represents the minimum change in the complete data log likelihood needed to continue the Expectation-Maximisation algorithm.}

\item{var_disc}{An integer value specifying which variable to discretise. If 0, no discretisation is necessary.}

\item{disc_length}{An integer value specifying how many neighbours can be searched over for the purposes of variable discetisation. If 0, all other possible stages may be merged over.}

\item{restrict}{A logical value indicating whether variable discretisation should be restricted to stages with the same unfolding of the process (TRUE) or not (FALSE).}

\item{mirror}{A logical value indicating whether variable discretisation should be equivalent across each unfolding of the process (TRUE) or not (FALSE).}

\item{cat_limit}{An integer value specifying the minimum number of categories to the variable can be discretised to. If 0, there is no minimum number of categories.}
}
\value{
An object of the S3 class StagedTree.
}
\description{
This function fits a ZIPCEG model to the chosen data set, based on chosen methods of parameter estimation and state imputation.
}
\details{
This function takes the same inputs as the \code{\link[=pceg]{pceg()}}, along with chosen methods of parameter estimation and state imputation, to fit a ZIPCEG model. As with \code{\link[=pceg]{pceg()}} it can perform variable discretisation methods, but is not to be used for vanilla CEGs.
}
\examples{
mod<-zipceg(knee_pain_obs,"nlm")
summary(mod)
}
