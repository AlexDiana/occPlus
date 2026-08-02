if(F){

rm(list = ls())

library(Rcpp); library(RcppArmadillo)
library(here); library(ggplot2); library(beepr); library(ggtern)
library(matrixNormal); library(splines); library(tidyverse)

sourceCpp(here("src/jsdm.cpp"))
source(here("R/jsdmfun.R"))
source(here("R/diagnostics.R"))

# SETTINGS -----

n <- 1000 # sites
ns <- 1000 # unique spatial locations
S <- 5 # species

p <- 3 # environmental covariates
g <- 0 # observed traits
gt <- 0 # unobserved traits

d <- 0 # number of factors
ds <- 2 # latent spatial factors

model <- "continuous" # "binary" #  "count"

# variances
{
  if(model == "continuous"){
    tau <- rep(0.001, S)# rgamma(S, 5, 5)
    rnb <- NULL
  } else if(model == "count"){
    rnb <- rpois(S, lambda = 50)
    tau <- NULL
  } else {
    tau <- NULL
    rnb <- NULL
  }

  # variation of residual environmental covariates
  sigma_b <- .01

  # variation of spatial traits
  sigma_ts <- .0001

  # variation of residual spatial field
  sigma_bs <- .0001

  # variation of factor scores
  sigma_h <- 1

  # variation of spatial field
  sigma_s <- .5

  # spatial field scale
  length_grid_ls <- 20
  l_s_grid <- seq(0.01, 0.3, length.out = length_grid_ls)
  idx_ls <- 5
  l_s <- l_s_grid[idx_ls]

}

list_simData <- simulateData(
  n, S, p,
  g, gt, d, tau, rnb, ds, ns,
  sigma_b, sigma_bs, sigma_ts, sigma_h, sigma_s, l_s,
  useSpatField = T, usingSplines = F, model)

# data
{
  z <- list_simData$data$z
  X <- list_simData$data$X
  Tr <- list_simData$data$Tr
  Xs <- list_simData$data$Xs
}

# true params
{
  B0_true <- list_simData$trueParams$B0
  B_true <- list_simData$trueParams$B
  G_true <- list_simData$trueParams$G
  A_true <- list_simData$trueParams$A
  C_true <- list_simData$trueParams$C
  U_true <- list_simData$trueParams$U
  Bs_true <- list_simData$trueParams$Bs
  Gs_true <- list_simData$trueParams$Gs
  As_true <- list_simData$trueParams$As
  Cs_true <- list_simData$trueParams$Cs
  U_true <- list_simData$trueParams$U
  L_true <- list_simData$trueParams$L
  SE_true <- list_simData$trueParams$SE
  spatField_true <- list_simData$trueParams$spatField
  sigma_b_true <- list_simData$trueParams$sigma_b
  sigma_bs_true <- list_simData$trueParams$sigma_bs
  sigma_s_true <- list_simData$trueParams$sigma_s
  l_s_true <- list_simData$trueParams$l_s
  tau_true <- list_simData$trueParams$tau
  rnb_true <- list_simData$trueParams$rnb
}

# PRIORS ---------

a_sigmab <- .01; b_sigmab <- .01
a_sigmabs <- .01; b_sigmabs <- .01
a_tau <- 5; b_tau <- 5
a_l_s <- 1; b_l_s <- 1

# MCMC --------

nchain <- 1
nburn <- 1000
niter <- 1000
nthin <- 1

updateBBsL <- T; trueBBsL <- !updateBBsL
updateGC <- F; trueGC <- !updateGC
updateA <- F; trueA <- !updateA
updateGCs <- F; trueGCs <- !updateGCs
updateAs <- F; trueAs <- !updateAs
updateU <- F; trueU <- !updateU
updateSigma <- F; trueSigma <- !updateSigma
update_ls <- F; truels <- !update_ls
update_taurnb <- T; truetaurnb <- !update_taurnb

# precompute projection and spatial quantities
{
  if(p > 0){
    P_X <- diag(1, nrow = n) - X %*% solve(t(X) %*% X) %*% t(X)
  }

  ps <- getDefaultSupportPoints(n)

  # Spatial covariates matrix
  list_Xs <- computeSpatialSummaries(Xs, ps, maxPoints = 5)
  Xs_centers <- list_Xs$Xs_centers
  Xs_index <- list_Xs$Xs_index
  X_s_index <- list_Xs$X_s_index
  X_s_centers <- list_Xs$X_s_centers
  X_tilde <- list_Xs$X_tilde
  X_s <- list_Xs$X_s

  list_SoRSummaries <- precomputeSORmatrices(l_s_grid, list_Xs)
}

# output
{
  B0_output <- array(NA, dim = c(S, niter, nchain))
  B_output <- array(NA, dim = c(p, S, niter, nchain))
  G_output <- array(NA, dim = c(g, p, niter, nchain))
  A_output <- array(NA, dim = c(S, gt, niter, nchain))
  C_output <- array(NA, dim = c(gt, p, niter, nchain))
  Bs_output <- array(NA, dim = c(ps, S, niter, nchain))
  Gs_output <- array(NA, dim = c(g, ps, niter, nchain))
  As_output <- array(NA, dim = c(S, gt, niter, nchain))
  Cs_output <- array(NA, dim = c(gt, ps, niter, nchain))
  U_output  <- array(NA, dim = c(n, d, niter, nchain))
  L_output  <- array(NA, dim = c(d, S, niter, nchain))
  tau_output <- array(NA, dim = c(S, niter, nchain))
  idx_ls_output <- array(NA, dim = c(niter, nchain))
  sigmab_output <- array(NA, dim = c(niter, nchain))
  sigmabs_output <- array(NA, dim = c(niter, nchain))
  varPart_output <- array(NA, dim = c(S, 4, niter, nchain))
}

}
