% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_synthetic.R
\name{create_synthetic}
\alias{create_synthetic}
\title{create_synthetic}
\usage{
create_synthetic(K, X.tar, n.src, r, B)
}
\arguments{
\item{K}{Number of source sites that individual data cannot be shared}

\item{X.tar}{X from target data}

\item{n.src}{Source sample size}

\item{r}{A positive integer}

\item{B}{A matrix of beta and w, beta is in the 1st column}
}
\value{
A list of synthetic X and Y
}
\description{
Create synthetic data under federated settings
}
