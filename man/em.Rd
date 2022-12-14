% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/em_algorithm_functions.R
\name{em}
\alias{em}
\title{The EM algorithm for a leaf}
\usage{
em(p_0, l_0, y, t, max_iter = 10000, tol = 1e-10)
}
\arguments{
\item{p_0}{A numeric initial value for the risk probability}

\item{l_0}{A numeric initial value for the rate.}

\item{y}{An integer vector of observed event counts.}

\item{t}{A numeric vector of observed times.}

\item{max_iter}{An integer for the maximum number of iterations for the algorithm.}

\item{tol}{A numeric value which represents the minimum change in the complete data log likelihood needed to continue the algorithm.}
}
\value{
A list containing the estimates of the risk probability and the rate, and the number of iterations necessary
}
\description{
This function carries out the Expectation-Maximisation algorithm for a leaf, given initial values of the parameters and the count and time vectors.
}
\details{
For a specific leaf, which is a unique evolution of the process, the event counts and observed times are used, alongside the \code{\link[=cdll]{cdll()}} function and some initial parameter values,
to carry out the Expectation-Maximisation algorithm.
}
\examples{
y<-c(1,0,2,0)
t<-c(0.5,0.7,0.8,1.2)
params<-c(0.5,1.5)
em(params[1],params[2],y,t)
}
