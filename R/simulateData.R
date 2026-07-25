

logit <- function(x) log(x / (1-x))

logistic <- function(x) 1 / (1 + exp(-x))


#' simulateOccJSDMData
#'
#' Simulate data for any supported model type
#'
#' @param list_datasettings List of data dimension settings (n, S, g, M, P, K,
#'   ncov_psi, ncov_theta).
#' @param list_params List of detection/amplification parameters (p, q, theta0,
#'   theta_baseline). Optionally also `mu1`, `sigma1`, `mu0`, `sigma0`
#'   controlling the simulated read-count intensities (see Details); if
#'   omitted these default to `mu1 = 5`, `sigma1 = 1`, `mu0 = 1.5`,
#'   `sigma0 = 1`.
#' @param list_jsdmParams List of JSDM parameters (gt, d, ds, sigma_b, sigma_bs,
#'   sigma_ts, sigma_h, l_s).
#' @param model Supported data formats are "binary","continuous","occupancy"
#'  and "two_stage"
#'
#' @details
#' General-purpose simulation function that supports multiple model types.
#' For two-stage eDNA data specifically, use \code{\link{simulateOccJSDMData}}.
#'
#' @return A list with simulated data and true parameters.
#'
#' @export
#'
simulateOccJSDMData <- function(list_datasettings,
                                list_params,
                                list_jsdmParams,
                                model){

  if(!(model %in% c("binary","occupancy","continuous","two_stage"))){
    stop("The model used is not supported")
  }

  # read data settings
  {
    n <- list_datasettings$n
    S <- list_datasettings$S
    g <- list_datasettings$g
    M <- list_datasettings$M
    P <- list_datasettings$P
    K <- list_datasettings$K
    ncov_psi <- list_datasettings$ncov_psi
    ncov_theta <- list_datasettings$ncov_theta
  }

  # read jsdm params
  {
    gt <- list_jsdmParams$gt
    d <- list_jsdmParams$d
    ds <- list_jsdmParams$ds
    sigma_b <- list_jsdmParams$sigma_b
    sigma_bs <- list_jsdmParams$sigma_bs
    sigma_ts <- list_jsdmParams$sigma_ts
    sigma_h <- list_jsdmParams$sigma_h
    sigma_s <- list_jsdmParams$sigma_s
    l_s <- list_jsdmParams$l_s
    tau <- list_jsdmParams$tau
    rnb <- list_jsdmParams$rnb
    useSpatField <- list_jsdmParams$useSpatField
  }

  # read param settings
  {
    p <- list_params$p
    q <- list_params$q
    theta0 <- list_params$theta0
    theta_baseline <- list_params$theta_baseline
  }

  if(model %in% c("binary","occupancy","two_stage")){
    jSDMsimModel <- "binary"
  } else if(model == "continuous"){
    jSDMsimModel <- "continuous"
  }

  list_simjSDMData <- simulateData(
    n, S, ncov_psi,
    g, gt, d, tau, rnb, ds, n,
    sigma_b, sigma_bs, sigma_ts, sigma_h, sigma_s, l_s,
    useSpatField = useSpatField, usingSplines = F, model = jSDMsimModel)

  # data
  {
    z <- list_simjSDMData$data$z
    X_psi <- list_simjSDMData$data$X
    Tr <- list_simjSDMData$data$Tr
    Xs <- list_simjSDMData$data$Xs
  }

  # create data structure (for occupancy or two-stage model)
  {
    N <- sum(M)

    if(model %in% c("occupancy","two_stage")){
      N2 <- P * N
      N3 <- sum(K)

      sumM <- c(0, cumsum(M)[-n])
      sumP <- c(0, cumsum(rep(P, N))[-N])
      sumK <- c(0, cumsum(K)[-N2])

      list_idx <- createDataIdx(n, M, P, K, model == "two_stage")
      idx_z_w <- list_idx$idx_z_w
      idx_z_k <- list_idx$idx_z_k
      idx_w_p <- list_idx$idx_w_p
      idx_z_p <- list_idx$idx_z_p
      idx_w_k <- list_idx$idx_w_k
      idx_p_k <- list_idx$idx_p_k
    }

  }

  if(model %in% c("occupancy","two_stage")){
    X_theta <- cbind(1, matrix(rnorm(N * ncov_theta), N, ncov_theta))

    beta_theta_true <- matrix(sample(c(-1,1,0), (ncov_theta + 1) * S, replace = T), ncov_theta + 1, S)
    beta_theta_true[1,] <- logit(theta_baseline)
    theta_true <- logistic(X_theta %*% beta_theta_true)

    theta0_true <- theta0

    z_rep <- z[idx_z_w,]
    w <- sapply(1:S, function(s){
      sapply(1:N, function(i){
        if(z_rep[i,s] == 1){
          rbinom(1, 1, theta_true[i,s])
        } else {
          rbinom(1, 1, theta0_true[s])
        }
      })
    })

    if(model == "two_stage"){
      p_true <- p
      q_true <- q

      y <- matrix(NA, N3, S)
      w_rep <- w[idx_w_k,]

      cimk_true <- sapply(1:S, function(s){
        sapply(1:N3, function(i){
          if(w_rep[i,s] == 1){
            rbinom(1, 1, p_true[idx_p_k[i],s])
          } else {
            2 * rbinom(1, 1, q_true[idx_p_k[i],s])
          }
        })
      })

      # Read-count model (matches the threshold = 0 mode of runOccJSDM()): given
      # a true detection (cimk_true == 1), log(y + 1) ~ Normal(mu1, sigma1); given
      # a false-positive/contamination detection (cimk_true == 2),
      # log(y + 1) ~ Normal(mu0, sigma0). No reads (cimk_true == 0) gives y = 0.
      mu1_true <- get_param(list_params, "mu1", 5)
      sigma1_true <- get_param(list_params, "sigma1", 1)
      mu0_true <- get_param(list_params, "mu0", 1.5)
      sigma0_true <- get_param(list_params, "sigma0", 1)

      logy1 <- matrix(0, N3, S)
      logy1[cimk_true == 1] <- rnorm(sum(cimk_true == 1), mu1_true, sigma1_true)
      logy1[cimk_true == 2] <- rnorm(sum(cimk_true == 2), mu0_true, sigma0_true)

      y <- round(exp(logy1) - 1)
      y[y < 0] <- 0
    }

  }

  if(model %in% c("continuous","binary")){
    y <- z
  } else if(model == "occupancy"){
    y <- w
  } else if(model == "two_stage"){
    y <- y
  }

  speciesNames <- paste0("OTU_", 1:S)

  colnames(y) <- speciesNames
  rownames(Tr) <- speciesNames
  colnames(Tr) <- paste0("Trait_", seq_len(ncol(Tr)))

  if(model %in% c("binary","continuous")){
    data_info <- data.frame(
      X_psi = X_psi,
      Xs = Xs
    )
  } else if(model %in% c("occupancy")){
    data_info <- data.frame(
      Site = idx_z_w,
      X_psi = X_psi[idx_z_w,],
      Xs = Xs[idx_z_w,],
      X_theta = X_theta[,-1]
    )
  } else if (model == "two_stage"){
    data_info <- data.frame(
      Site = idx_z_k,
      Sample = idx_w_k,
      Primer = idx_p_k,
      X_psi = X_psi[idx_z_k,],
      Xs = Xs[idx_z_k,],
      X_theta = X_theta[idx_w_k,-1]
    )
  }

  data <- list(info = data_info,
               OTU = y,
               traits = Tr)

  true_params <- list(
    "jsdmParams_true" = list_simjSDMData$trueParams
  )

  if(model %in% c("occupancy","two_stage")){
    true_params$beta_theta_true <- beta_theta_true
    true_params$z_true <- z
  }

  if(model == "two_stage"){
    true_params$w_true <- w
    true_params$p_true <- p_true
    true_params$q_true <- q_true
  }

  list(true_params = true_params,
       data_list = data)
}
