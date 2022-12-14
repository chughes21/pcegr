% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/numerical_optimisation_functions.R
\name{glogit}
\alias{glogit}
\title{The logit function}
\usage{
glogit(p)
}
\arguments{
\item{p}{A numeric value between 0 and 1}
}
\value{
A numeric value on the real line.
}
\description{
This function maps a probability to a value on the real line, for the purposes of the \code{\link[=nlm_zip]{nlm_zip()}} function.
}
