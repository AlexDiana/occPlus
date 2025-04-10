% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output.R
\name{plotSigElementsCovMatrix}
\alias{plotSigElementsCovMatrix}
\title{plotSigElementsCovMatrix}
\usage{
plotSigElementsCovMatrix(fitmodel, idx_species = NULL)
}
\arguments{
\item{fitmodel}{Output from the function runOccPlus}
}
\value{
A ggplot object
}
\description{
Plot the significant element of the covariance matrix
}
\details{
Plots the 95% credible interval of density of true positives and false positives
}
\examples{
\dontrun{
plotSigElementsCovMatrix(fitmodel)
}

}
