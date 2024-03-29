% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uniform_time_functions.R
\name{mle_variable_time_zip}
\alias{mle_variable_time_zip}
\title{The Maximum Likelihood Estimator under Variable Time}
\usage{
mle_variable_time_zip(data)
}
\arguments{
\item{data}{A data set, where the observed count vector and time vector (if included) are the last two columns.}
}
\value{
A list containing mles for at-risk rate and risk proportion for each leaf
}
\description{
Calculates the maximum likelihood estimator (mle) for a ZIP distribution under variable time, by dividing the counts by observation time.
}
\examples{
mle_variable_time_zip(knee_pain_obs)

}
