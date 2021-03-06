% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findNumberCells.R
\name{findNumberCells}
\alias{findNumberCells}
\title{Select the number of cell types using SVD methods}
\usage{
findNumberCells(yref, scree = c("drop", "cumvar", "residual"))
}
\arguments{
\item{yref}{matrix, numeric expression matrix}

\item{scree}{character, method to estimate n.types}
}
\value{
number of cell types
}
\description{
The function estimates the number of cell-types using one
of three SVD methods on the input expression matrix
}
