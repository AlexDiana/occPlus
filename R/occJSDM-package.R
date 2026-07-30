#' occJSDM: Occupancy and Joint Species Distribution Models for eDNA Metabarcoding Data
#'
#' Fits occupancy and joint species distribution models (JSDMs) to
#' environmental DNA (eDNA) metabarcoding read data within a hierarchical
#' Bayesian framework, explicitly modeling the two-stage detection process
#' of field sample collection followed by PCR amplification and
#' sequencing, while optionally accounting for species traits and spatial
#' autocorrelation via a Gaussian process. The main entry point for
#' fitting a model is \code{\link{runOccJSDM}}; \code{\link{simulateOccJSDMData}}
#' simulates data in the format it expects. See Ji et al. (2025)
#' \doi{10.1111/ele.70302} for methodological details.
#'
#' @return This man page documents the package as a whole; it has no
#'   return value of its own. See \code{\link{runOccJSDM}} for the
#'   function that fits a model and returns results.
#'
#' @keywords internal
#' @useDynLib occJSDM, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats binomial cov cov2cor dbinom dgamma dnorm dpois glm kmeans logLik median predict quantile rbeta rbinom reorder rgamma rnbinom rnorm rpois runif sd setNames var
#' @importFrom tidyr pivot_longer
"_PACKAGE"
