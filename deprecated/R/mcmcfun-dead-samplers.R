# Superseded MCMC samplers, moved out of R/mcmcfun.R on 30 July 2026.
#
# TODO.md group D item 1, "DEPRECATE, BUT KEEP FOR REFERENCE". All four were
# unreachable: no callers in R/, tests/ or vignettes/, absent from NAMESPACE,
# and not referenced as string literals anywhere, so nothing could reach them
# via do.call(), get() or match.fun().
#
#   sample_z       references undeclared M, sumM, n, z, w, y
#   sample_w       same
#   sample_cimk    same
#   loglik_sigma1  an unfinished stub whose entire body is p[primerIdx[]]
#
# The three samplers are superseded by the _cpp implementations in src/.
# loglik_sigma1 was never finished. This directory is excluded from the build
# by .Rbuildignore, so nothing here ships.

# OK
sample_z <- function( w, psi, theta, theta0){

  S <- ncol(psi)
  n <- nrow(psi)

  z <- matrix(NA, n, S)

  for (s in 1:S) {

    for (i in 1:n) {

      p_zsequal1 <-
        sum(
          sapply(1:M[i], function(m){
            dbinom(w[sumM[i] + m,s], 1, theta[sumM[i] + m,s], log = T)
          })
        ) + dbinom(1, 1, psi[i,s], log = T)

      p_zsequal0 <-
        sum(
          sapply(1:M[i], function(m){
            dbinom(w[sumM[i] + m,s], 1, theta0[s], log = T)
          })
        ) + dbinom(0, 1, psi[i,s], log = T)

      p_1 <- exp(p_zsequal1) / (exp(p_zsequal1) + exp(p_zsequal0))

      z[i,s] <- rbinom(1, 1, p_1)

    }

  }

  z
}

# OK
sample_w <- function(y, theta, theta0, p, q,
                     M, K, sumL, sumM, sumK, maxL){

  S <- ncol(theta)
  N <- nrow(theta)
  w <- matrix(NA, nrow = N, ncol = S)

  for (s in 1:S) {

    for (i in 1:n) {

      for (m in 1:M[i]) {

        p_wsequal1 <-
          sum(
            sapply(1:maxL, function(l){
              sapply(1:K[sumL[sumM[i] + m] + l], function(k){
                dbinom(y[sumK[sumL[sumM[i] + m] + l] + k,s], 1,
                       p[l,s], log = T)
              })
            })
          )

        p_wsequal0 <-
          sum(
            sapply(1:maxL, function(l){
              sapply(1:K[sumL[sumM[i] + m] + l], function(k){
                dbinom(y[sumK[sumL[sumM[i] + m] + l] + k,s], 1,
                       q[l,s], log = T)
              })
            })
          )

        if(z[i, s] == 1){

          p_wsequal1 <- p_wsequal1 + dbinom(1, 1, theta[sumM[i] + m,s], log = T)
          p_wsequal0 <- p_wsequal0 + dbinom(0, 1, theta[sumM[i] + m,s], log = T)

        } else {

          p_wsequal1 <- p_wsequal1 + dbinom(1, 1, theta0[s], log = T)
          p_wsequal0 <- p_wsequal0 + dbinom(0, 1, theta0[s], log = T)

        }

        p_ws1 <- exp(p_wsequal1) / (exp(p_wsequal1) + exp(p_wsequal0))

        w[sumM[i] + m,s] <- rbinom(1, 1, p_ws1)

      }

    }

  }

  w
}

sample_cimk <- function(logy1, mu1, sigma1, pi0, sigma0,
                        p, q, idx_k, primerIdx){

  N3 <- length(idx_k)
  S <- ncol(p)

  c_imk <- matrix(NA, N3, S)

  for (i in 1:N3) {
    for (s in 1:S) {

      term1_loglik <- dnorm(logy1[i,s], mu1, sigma1, log = T)
      term2_loglik <- ifelse(logy1[i,s] == 0, log(pi0),
                             dnorm(y[i,s], 0, sigma0, log = T) - log(.5))
      # term2_loglik <- ifelse(logy1[i,s] == 0, log(pi0),
      #                        dlaplace(logy1[i,s], 0, sigma0, log = T) - log(.5))
                             # dnorm(logy1[i,s], 0, sigma0, log = T) - log(.5))

        # log(pi0 * dnorm(logy1[i,s], 0, sigma0, log = T) +
        #                     (1 - pi0))

      if(w[idx_k[i],s] == 1){
        term1_prior <- dbinom(1, 1, p[primerIdx[i],s], log = T)
        term2_prior <- dbinom(0, 1, p[primerIdx[i],s], log = T)
      } else {
        term1_prior <- dbinom(1, 1, q[primerIdx[i],s], log = T)
        term2_prior <- dbinom(0, 1, q[primerIdx[i],s], log = T)
      }

      term12_diff <- (term1_loglik + term1_prior) - (term2_loglik + term2_prior)

      # p_cimk1 <- exp(term12_diff) / (exp(term12_diff) + 1)
      p_cimk1 <- 1 / (exp(-term12_diff) + 1)


      c_imk[i,s] <- rbinom(1, 1, p_cimk1)
    }
  }

  c_imk
}

loglik_sigma1 <- function(w, logy1){

  sum(
    sapply(1:S, function(s){
      sapply(1:N3, function(i){

        p[primerIdx[]]

      })
    })
  )

}
