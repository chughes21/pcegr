% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kl_divergence.R
\name{kl_divergence}
\alias{kl_divergence}
\title{The Kullback-Leibler Divergence function}
\usage{
kl_divergence(new_mod, old_mod = NULL, old_prob = NULL, precision = 5)
}
\arguments{
\item{new_mod}{A StagedTree object.}

\item{old_mod}{A StagedTree object. Must be included if old_prob is not.}

\item{old_prob}{A numeric vector of path probabilities that sums to one. Must be included if old_mod is not.}

\item{precision}{An integer value specifying how many decimal places numerical checks should be done to.}
}
\value{
A numeric value for the KL divergence between the StagedTree object and reference model.
}
\description{
Calculate the Kullback-Leibler (KL) divergence between a StagedTree object and a reference model, either another StagedTree object or a vector of path probabilities that sum to one.
}
\examples{
mod<-pceg(knee_pain_obs,2,TRUE,TRUE)
mod2<-pceg(knee_pain_obs,2,TRUE,TRUE,indep=c(2:3))
kl_divergence(mod,old_mod=mod2)
equiv_prob<-rep(1/8,8)
kl_divergence(mod,old_prob=equiv_prob)
}
