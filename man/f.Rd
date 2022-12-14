% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/numerical_optimisation_functions.R
\name{f}
\alias{f}
\title{The zero-inflated Poisson likelihood for a single observation.}
\usage{
f(p, lambda, y, t)
}
\arguments{
\item{p}{A numeric value between 0 and 1 specifying the risk probability.}

\item{lambda}{A numeric value specifying the rate of the Poisson process.}

\item{y}{An integer value for the observed event count.}

\item{t}{A positive numeric value for the observation time.}
}
\value{
A numeric value specifying the likelihood.
}
\description{
This function calculates the likelihood for a single observation assuming a zero-inflated Poisson distribution
}
