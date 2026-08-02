
buildGrid <- function(XY_sp, gridStep){

  x_grid <- seq(min(XY_sp[,1]) - (1.5) * gridStep,
                max(XY_sp[,1]) + (1.5) * gridStep, by = gridStep)
  y_grid <- seq(min(XY_sp[,2]) - (1.5) * gridStep,
                max(XY_sp[,2]) + (1.5) * gridStep, by = gridStep)

  pointInGrid <- matrix(T, nrow = length(x_grid), ncol = length(y_grid))

  for (i in 2:(length(x_grid) - 1)) {

    for (j in 2:(length(y_grid) - 1)) {

      isAnyPointInBandRight <- isPointInBandRight(XY_sp, x_grid, y_grid, i - 1, j - 1)

      isAnyPointInBandLeft <- isPointInBandLeft(XY_sp, x_grid, y_grid, i - 1, j - 1)

      isAnyPointInBandUp <- isPointInBandUp(XY_sp, x_grid, y_grid, i - 1, j - 1)

      isAnyPointInBandDown <- isPointInBandDown(XY_sp, x_grid, y_grid, i - 1, j - 1)

      if(!isAnyPointInBandRight | !isAnyPointInBandLeft | !isAnyPointInBandUp | !isAnyPointInBandDown){
        pointInGrid[i,j] <- F
      }

    }

  }

  pointInGrid <- pointInGrid[-c(1,nrow(pointInGrid)),]
  pointInGrid <- pointInGrid[,-c(1,ncol(pointInGrid))]
  x_grid <- x_grid[-c(1,length(x_grid))]
  y_grid <- y_grid[-c(1,length(y_grid))]

  allPoints <- cbind(expand.grid(x_grid, y_grid), as.vector((pointInGrid)))
  allPoints <- allPoints[allPoints[,3],-3]

  allPoints
}

logistic <- function(x){
  1 / (1 + exp(-x))
}

computeU <- function(X_ord, beta_ord, E){
  X_ord %*% beta_ord + E
}

computePsiE <- function(X_psi, beta_psi, X_ord, beta_ord, LL){

  logit_psi = X_psi %*% beta_psi + (X_ord %*% beta_ord) %*% LL

  logistic(logit_psi)

}

computeTheta <- function(X_theta, beta_theta){
  logistic(X_theta %*% beta_theta)
}

# OK
sample_pq <- function(c_imk, w, primerIdx, idx_k, maxL, a_p,
                      b_p, a_q, b_q){

  p <- matrix(NA, maxL, S)
  q <- matrix(NA, maxL, S)

  w_all <- w[idx_k,,drop=F]

  for (s in 1:S) {

    for (l in 1:maxL) {

      w1_primerl_cases_1 <- sum(primerIdx == l & w_all[,s] == 1 & c_imk[,s] == 1)
      w1_primerl_cases_0 <- sum(primerIdx == l & w_all[,s] == 1 & c_imk[,s] == 0)

      w0_primerl_cases_1 <- sum(primerIdx == l & w_all[,s] == 0 & c_imk[,s] == 2)
      w0_primerl_cases_0 <- sum(primerIdx == l & w_all[,s] == 0 & c_imk[,s] == 0)

      p[l,s] <- rbeta(1, a_p + w1_primerl_cases_1, b_p + w1_primerl_cases_0)
      q[l,s] <- rbeta(1, a_q + w0_primerl_cases_1, b_q + w0_primerl_cases_0)


    }

  }

  list("p" = p,
       "q" = q)

}

# OK
sample_betatheta <- function(w, z, beta_theta, idx_z, X_theta,
                             b_betatheta, B_betatheta){

  ncov_theta <- nrow(beta_theta)
  S <- ncol(beta_theta)

  z_all <- z[idx_z,]

  for (s in 1:S) {

    k <- as.vector(t(w[z_all[,s]==1,s])) - .5

    n <- rep(1, length(k))

    X_thetasubset <- X_theta[z_all[,s]==1,,drop=F]

    beta_theta[,s] <- sample_beta_nocov_cpp(beta_theta[,s], X_thetasubset,
                                            b_betatheta, B_betatheta, n, k)

  }

  beta_theta
}

# OK
sample_theta0 <- function(z, w, idx_z, a_theta0, b_theta0){

  S <- ncol(z)

  z_all <- z[idx_z,]

  theta0 <- rep(NA, S)

  for (s in 1:S) {

    z0_1 <- sum(z_all[,s] == 0 & w[,s] == 1)
    z0_0 <- sum(z_all[,s] == 0 & w[,s] == 0)

    theta0[s] <- rbeta(1, a_theta0 + z0_1, b_theta0 + z0_0)

  }

  theta0
}

# OK
sample_betapsiLL <- function(beta_psi, LL, Omega, X_psi, k, U,
                             prior_beta_psi, prior_beta_psi_sd){

  ncov_psi <- ncol(X_psi)
  d <- nrow(LL)
  S <- ncol(LL)

  beta_all <- matrix(NA, ncov_psi + d, S)

  for (s in 1:S) {

    elemZero <- max(d - s, 0)
    elem1 <- ifelse(s <= d, 1, 0)
    # elemZero <- elem1 <- 0
    elemNonZero <- d - elemZero - elem1

    if(elem1 > 0){

      k_current <- k[,s] - Omega[,s] * U[,s]

    } else {

      k_current <- k[,s]

    }

    X_all <- cbind(X_psi, U[,seq_len(elemNonZero)])

    b_current <- rep(prior_beta_psi, ncov_psi + elemNonZero)
    B_current <- diag(prior_beta_psi_sd, nrow = ncov_psi + elemNonZero)

    # betapntsiLL_current <- c(beta_psi[,s], LL[,s])
    if(ncol(X_all) > 0){
      betapsiLL <- sample_beta_cpp(X_all, B_current, b_current, Omega[,s], k_current)
      beta_all[seq_len(ncov_psi),s] <- betapsiLL[seq_len(ncov_psi)]
    }

    if(elemNonZero > 0){

      beta_all[ncov_psi + seq_len(elemNonZero),s] <-
        betapsiLL[ncov_psi + seq_len(elemNonZero)]

    }


  }

  if(T){

    # set lower diagonal to 0
    beta_all[ncov_psi + seq_len(d),seq_len(d)][!upper.tri(beta_all[ncov_psi + seq_len(d),seq_len(d)])] <- 0

    # set diagonal to 1
    if(d > 1){

      diag(beta_all[ncov_psi + seq_len(d),seq_len(d)]) <- 1
    } else {
      beta_all[ncov_psi + 1,1] <- 1
    }

  }


  beta_psi <- beta_all[seq_len(ncov_psi),,drop=F]
  LL <- beta_all[ncov_psi + seq_len(d),,drop=F]

  list("beta_psi" = beta_psi,
       "LL" = LL)

}

# OK
sample_LL <- function(beta_psi, LL, Omega, X_psi, k, U,
                             prior_beta_psi, prior_beta_psi_sd){

  ncov_psi <- ncol(X_psi)
  d <- nrow(LL)
  S <- ncol(LL)
  n <- nrow(U)

  LL <- matrix(NA, d, S)

  k <- k - Omega * X_psi %*% beta_psi

  for (s in 1:S) {

    elemZero <- max(d - s, 0)
    elem1 <- ifelse(s <= d, 1, 0)
    # elemZero <- elem1 <- 0
    elemNonZero <- d - elemZero - elem1

    if(elem1 > 0){

      k_current <- k[,s] - Omega[,s] * U[,s]

    } else {

      k_current <- k[,s]

    }

    # betapntsiLL_current <- c(beta_psi[,s], LL[,s])

    if(elemNonZero > 0){

      X_all <- matrix(U[,seq_len(elemNonZero)], n, elemNonZero)

      b_current <- rep(prior_beta_psi, elemNonZero)
      B_current <- diag(prior_beta_psi_sd, nrow = elemNonZero)

      LL_s <- sample_beta_cpp(X_all, B_current, b_current, Omega[,s], k_current)
      LL[seq_len(elemNonZero),s] <- LL_s

    }

  }

  if(T){

    # set lower diagonal to 0
    LL[seq_len(d),seq_len(d)][!upper.tri(LL[seq_len(d),seq_len(d)])] <- 0

    # set diagonal to 1
    if(d > 1){
      for (i in 1:d) {
        LL[i,i] <- 1
      }
      # diag(LL[ncov_psi + seq_len(d),seq_len(d)]) <- 1
    } else {
      LL[1,1] <- 1
    }

  }

  LL

}

sample_betapsis <- function(beta_psi, beta_s, LL, Omega, X_psi,
                            k, U, Ks, As,
                            prior_beta_psi, prior_beta_psi_sd){

  ncov_psi <- ncol(X_psi)
  d <- nrow(LL)
  S <- ncol(LL)
  ncov_spat <- ncol(Ks)
  n_spatfactors <- ncol(As)

  beta_psi <- matrix(NA, ncov_psi + n_spatfactors, S)

  for (s in 1:S) {

    # if(elem1 > 0){

    k_current <- k[,s] - Omega[,s] * U %*% LL[,s]

    # } else {
    #
    #   k_current <- k[,s]
    #
    # }

    # X_psi <- cbind()

    b_current <- rep(prior_beta_psi, ncov_psi)
    B_current <- diag(prior_beta_psi_sd, nrow = ncov_psi)

    beta_psi[,s] <- sample_beta_cpp(X_psi, B_current, b_current, Omega[,s], k_current)

  }

  beta_psi

}

sample_musigma <- function(sigma1,
                            logy1, c_imk,
                            mu_mu1, sd_mu1,
                            a_sigma1, b_sigma1, tpfp){

  idx_tpfp <- ifelse(tpfp, 1, 2)

  logy1_cimk1 <- logy1[c_imk == idx_tpfp]
  n_samples <- length(logy1_cimk1)

  if(n_samples > 0){

    sample_mean <- mean(logy1_cimk1)

    ww <- (1/sd_mu1^2) / ((1/sd_mu1^2) + (n_samples / sigma1^2))
    mu_n <- ww * mu_mu1 + (1-ww) * sample_mean
    sigma2_n <- 1 / ((1 / sd_mu1^2) + (n_samples / sigma1^2))


  } else {

    mu_n <- mu_mu1
    sigma2_n <- sd_mu1

  }

  mu1 <- rnorm(1, mean = mu_n, sd = sqrt(sigma2_n))

  if(n_samples > 0){

    alpha_n <- a_sigma1 + n_samples / 2
    beta_n <- b_sigma1 + 0.5 * sum((logy1_cimk1 - mu1)^2)

  } else {

    alpha_n <- a_sigma1
    beta_n <- b_sigma1

  }

  sigma1 <- sqrt(1 / rgamma(1, shape = alpha_n, rate = beta_n))

  list("mu1" = mu1,
       "sigma1" = sigma1)
}

sample_sigma0 <- function(logy1, c_imk,
                          a_sigma1, b_sigma1){

  logy1_cimk1 <- logy1[c_imk == 2]
  n_samples <- length(logy1_cimk1)

  if(n_samples > 0){

    alpha_n <- a_sigma1 + n_samples / 2
    beta_n <- b_sigma1 + 0.5 * sum((logy1_cimk1)^2)

  } else {

    alpha_n <- a_sigma1
    beta_n <- b_sigma1

  }

  sigma0 <- sqrt(1 / rgamma(1, shape = alpha_n, rate = beta_n))

  sigma0
}

sample_mu0sigma0 <- function(y, c_imk,
                            a_pi0, b_pi0,
                            a_sigma0, b_sigma0){

  logy1_cimk0 <- y[c_imk == 0]
  logy1_cimk0_pos <- logy1_cimk0[logy1_cimk0 > 0]

  n_samples <- length(logy1_cimk0)

  num0 <- sum(logy1_cimk0 == 0)
  num1 <- sum(logy1_cimk0 > 0)

  pi0 <- rbeta(1, a_pi0 + num0, b_pi0 + (n_samples - num0))

  alpha_n <- a_sigma0 + num1 / 2
  beta_n <- b_sigma0 + 0.5 * sum(logy1_cimk0_pos^2)

  sigma0 <- sqrt(1 / rgamma(1, shape = alpha_n, rate = beta_n))

  list("pi0" = pi0,
       "sigma0" = sigma0)
}



# WAIC CALCULATIONS -----

update_waic_summary <- function(logliks, list_waic, iter){

  M2 <- list_waic$M2
  mean_log <- list_waic$mean_log
  mean_lik <- list_waic$mean_lik

  delta_log <- logliks - mean_log
  mean_log <- mean_log + delta_log / iter

  delta2_log <- logliks - mean_log
  M2 <- M2 + delta_log * delta2_log

  delta_lik <- exp(logliks) - mean_lik
  mean_lik <- mean_lik + delta_lik / iter

  list("M2" = M2,
       "mean_log" = mean_log,
       "mean_lik" = mean_lik)

}

compute_waic <- function(list_waic, numIters){

  mean_lik <- list_waic$mean_lik
  var_loglik <- list_waic$M2 / (numIters - 1)

  lppd <- sum(log(mean_lik), na.rm = T)
  p_waic <- sum(var_loglik , na.rm = T)

  WAIC <- - 2 * (lppd - p_waic)

  WAIC
}

computeModelLoglikJSDM <- function(z, eta, model, tau = NULL){

  n <- nrow(z)
  S <- ncol(z)

  if(model == "continuous"){

    logliks <- as.vector(
      sapply(1:n, function(i){
        sapply(1:S, function(j){
          dnorm(z[i,j], eta[i,j], tau[j], log = T)
        })
      })
    )

  } else if(model == "binary"){

    psi <- logistic(eta)

    logliks <- as.vector(
      sapply(1:n, function(i){
        sapply(1:S, function(j){
          dbinom(z[i,j], 1, psi[i,j], log = T)
        })
      })
    )

  } else if (model == "counts"){

    mu <- exp(eta)

    logliks <- as.vector(
      sapply(1:n, function(i){
        sapply(1:S, function(j){
          dpois(z[i,j], lambda = mu[i,j], log = T)
        })
      })
    )

  }

  logliks

}

computeModelLoglikFirstStage <- function(w, z, theta, theta0, idx_z_w){

  z_rep <- z[idx_z_w,]

  logliks <- as.vector(
    sapply(1:nrow(w), function(i){
      sapply(1:ncol(w), function(s){
        if(z_rep[i,s] == 1){
          dbinom(w[i,s], 1, theta[i,s], log = T)
        } else {
          dbinom(w[i,s], 1, theta0[s], log = T)
        }

      })
    })
  )

  logliks
}

computeModelLoglikSecondStage <- function(y, w, p, q, idx_w_k, idx_p_k){

    w_rep <- w[idx_w_k,]

    logliks <- as.vector(
      sapply(1:ncol(y), function(s){
        sapply(1:nrow(y), function(i){
          if(w_rep[i,s] == 1){
            dbinom(y[i,s], 1, p[idx_p_k[i],s], log = T)
          } else {
            dbinom(y[i,s], 1, q[idx_p_k[i],s], log = T)
          }
        })
      })
    )


  logliks
}
