.onLoad <- function(libname, pkgname) {
  # Set mc.cores (but limit to 2 for CRAN compliance)
  options(mc.cores = parallel::detectCores())
  rstan::rstan_options(auto_write = TRUE)
}


