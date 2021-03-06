% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lar.R
\name{lar}
\alias{lar}
\title{Least angle regression for compressed sensing}
\usage{
lar(need, train, seed = 1218)
}
\arguments{
\item{need}{vector, numeric vector of gene that is needed}

\item{train}{matrix, numeric expression matrix of training data}

\item{seed}{numeric, random seed}
}
\value{
list with coefficients and predictive R2
}
\description{
The function runs least angle regression to
train compression for a single gene
}
