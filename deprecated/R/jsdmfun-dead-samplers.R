# Superseded MCMC samplers, moved out of R/jsdmfun.R on 30 July 2026.
#
# Scope was set by Doug: only the four sample_* functions from TODO.md group D
# item 2. The other four functions that item lists (computePredictiveProbs,
# partition_r2, returnSpatialEffectMean, plotSpatialEffect) are dead by the
# same test but are NOT sampler code, and two of them look like unfinished
# groundwork rather than abandoned work, so they stay in R/ pending a decision.
#
# All four below were unreachable: no callers in R/, tests/ or vignettes/, not
# in NAMESPACE, and not referenced as strings anywhere, so nothing could reach
# them via do.call(), get() or match.fun() either. Each also referenced
# undefined globals, so they could not have run as written:
#
#   sample_BCsL      undefined U, gt, gts
#   sample_U         globals S, U, n
#   sample_Br        globals S, U, n
#   sample_BL_fixed  globals S, U, n
#
# Superseded by the _cpp implementations in src/. Kept for reference, per
# "DEPRECATE, BUT KEEP FOR REFERENCE". This directory is excluded from the
# build by .Rbuildignore, so nothing here ships.

# sample the fixed effects and the factor loadings
sample_BCsL <- function(k, X, H, G, Tr,
                      A, C, sigma_b,
                      Ks, Xs_centers, Gs, As,
                      Omega) {

  p <- ncol(X)
  d <- ncol(U)
  S <- ncol(Omega)
  ps <- nrow(As)

  M_B <- t(computeBtcoef(G, Tr, A, C, matrix(0, S, p)))

  B <- matrix(0, p, S)
  Cs <- matrix(0, gt, S)
  L <- matrix(0, d, S)

  # no spatial factor loadings if there is no spatial effect
  ps <- ifelse(ps > 0, ps, 0)

  if(p + gts + d > 0){

    for (s in 1:S) {

      k_current <- k[,s]

      XU <- cbind(X, U)

      b_current <- c(M_B[,s], rep(0, ps), rep(0, d))
      B_current <- diag(1, nrow = p + gts + d)
      diag(B_current)[seq_len(p)] <- sigma_b^2

      # if(ps == 0){
      #   BL <- sampleB(XU, B_current, b_current, Omega[,s], k_current)
      # } else {
        invB <- diag(1 / diag(B_current), nrow = p + ps + d)

        BCsL <- sampleB_SoR(XU, B_current, b_current, k_current,
                            Omega[,s], Xs_centers, Ks, ps)
        if(gts > 0){
          Cs[seq_len(gt),s] <- BCsL[p + d + seq_len(gts)]
        }

      # }
      B[seq_len(p),s] <- BCsL[seq_len(p)]
      L[seq_len(d),s] <- BCsL[p + seq_len(d)]


    }

  }

  Bt <- t(B) - Tr %*% G - A %*% C

  Bs <- computeBscoef(Tr, Gs, As, Cs)

  list("B" = B,
       "L" = L,
       "Bt" = Bt,
       "Bs" = Bs,
       "Cs" = Cs)
}

# sample factor scores
sample_U <- function(k, L, X, B, SE, B0, Omega, model){

  d <- nrow(L)
  n <- nrow(k)

  U <- matrix(NA, n, d)

  B0_mat <- matrix(B0, n, S, byrow = T)

  if(d > 0){

    if(model == "continuous"){
      k_new <- k - (B0_mat + X %*% B + SE)
    } else if(model == "binary"){
      k_new <- k - Omega * (B0_mat + X %*% B + SE)
    }

    B_current <- diag(1, nrow = d)
    b_current <- rep(0, d)

    for (i in 1:n) {

      U[i,] <- sampleB(t(L), B_current, b_current, Omega[i,], k_new[i,])

    }
  }

  U

}

# sample the ordination covariate coefficients
sample_Br <- function(k, L, Xr, Utilde, eta, Omega){

  ncov_ord <- ncol(Xr)
  d <- nrow(L)
  S <- ncol(L)

  Omega_vec <- as.vector(Omega)

  ktilde <- k - (eta + Utilde %*% L) * Omega

  X <- kronecker(t(L), Xr)

  if(ncov_ord > 0 & d > 0){
    Br <- sample_beta_cpp(X,
                          diag(100, nrow = ncol(X)),
                          rep(0, ncol(X)),
                          Omega_vec,
                          ktilde)

  }

  matrix(Br, ncov_ord, d, byrow = F)
}

# sample the fixed effects and the factor loadings
sample_BL_fixed <- function(k, X, H, G, Tr,
                            A, C, sigma_b,
                            Omega) {

  p <- ncol(X)
  d <- ncol(U)
  S <- ncol(Omega)

  M_B <- t(computeBtcoef(G, Tr, A, C, matrix(0, S, p)))

  B <- matrix(0, p, S)
  L <- matrix(0, d, S)

  for (s in 1:S) {

    # elemZero <- max(d - s, 0)
    # elem1 <- ifelse(s <= d, 1, 0)
    # elemNonZero <- d - elemZero - elem1
    #
    # if(elem1 > 0){
    #
    #   k_current <- k[,s] - Omega[,s] * U[,s]
    #
    # } else {
    #
    k_current <- k[,s]
    #
    # }

    elemNonZero <- d

    if(p + elemNonZero > 0){

      XU <- cbind(X, matrix(U[,seq_len(elemNonZero)], n, elemNonZero))

      b_current <- c(M_B[,s], rep(0, elemNonZero))
      B_current <- diag(1, nrow = p + elemNonZero)
      diag(B_current)[seq_len(p)] <- sigma_b^2

      BL <- sampleB(XU, B_current, b_current, Omega[,s], k_current)
      B[seq_len(p),s] <- BL[seq_len(p)]
      L[seq_len(elemNonZero),s] <- BL[p + seq_len(elemNonZero)]

    }

    # if(elem1 > 0) L[s,s] <- 1

  }

  Bt <- t(B) - Tr %*% G - A %*% C

  list("B" = B,
       "L" = L,
       "Bt" = Bt)
}
