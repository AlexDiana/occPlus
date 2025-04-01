#' runOccPlus
#'
#' Run the OccPlus model.
#'
#' @details
#' Fit the model described in Diana et al..
#'
#' @param data A list with two elements:
#' \describe{
#'   \item{info}{A data.frame of N rows (where N is the total number of samples
#'   analyzed). Further information in the vignette.}
#'   \item{OTU}{A matrix of dimension (N x S), where S is the number of species,
#'   with the number of reads of each species in each sample.}
#' }
#' @param d Number of factors.
#' @param occCovariates vector of the name of the covariates for the occupancy probabilities.
#' Names should match the column name in data$info.
#' @param ordCovariates vector of the name of the covariates for the occupancy probabilities.
#' Names should match the column name in data$info.
#' @param detCovariates vector of the name of the covariates for the detection probabilities.
#' Names should match the column name in data$info.
#'
#' @return Description of the return value (e.g., a list, data frame, or numeric output).
#'
#' @examples
#' \dontrun{
#' # Example usage
#' fitmodel  <- runOccPlus(data,
#' d = 2,
#' occCovariates = c("X_psi.1","X_psi.2"),
#' ordCovariates = c("X_ord.1","X_ord.2"),
#' detCovariates = c("X_theta"))
#' }
#'
runOccPlus <- function(data,
                          d,
                          occCovariates = c(),
                          ordCovariates = c(),
                          detCovariates = c(),
                       numSamples = 1000){

  data_info <- as.data.frame(data$info)
  OTU <- data$OTU

  # sort the data
  {
    data_info <- data_info %>%
      arrange(Site, Sample, Primer)
  }

  # samples per site
  {
    M_df <- data_info %>% group_by(Site, Sample) %>%
      slice(1) %>%
      group_by(Site) %>%
      summarise(M = n())

    M <- M_df$M
    names(M) <- M_df$Site

    n <- length(M)

    sumM <- c(0, cumsum(M)[-n])

    siteNames <- unique(data_info$Site)

  }

  # marker per samples
  {
    # number of markers

    L_df <- data_info %>% group_by(Site, Sample, Primer) %>%
      slice(1) %>%
      group_by(Site, Sample) %>%
      summarise(L = n())

    L <- L_df$L
    names(L) <- L_df$Sample

    # number of observations
    {
      L_all <- data_info %>% group_by(Site, Sample, Primer) %>%
        slice(1) %>%
        group_by(Site, Sample) %>%
        summarise(L_m = n())
      L_all <- L_all$L_m

      sumL <- c(0, cumsum(L_all)[-length(L_all)])

    }

    primerNames <- unique(data_info$Primer)
    # sumL <- c(0, cumsum(L)[-length(L)])

  }

  # pcr per marker
  {
    data_K <- data_info %>%
      group_by(Site, Sample, Primer) %>%
      summarise(K = n(),
                across(contains("Species"),function(x){sum(x > 0)})
      ) %>%
      ungroup()#%>%
    # ungroup() %>%
    # mutate(across(contains("Species"), function(x){sum(x > 0)}))

    K <- data_K$K

    numL <- length(K)

    sumK <- c(0, cumsum(K)[-numL])
  }

  logy1 <- log(OTU + 1)

  # data infos
  {
    n <- length(M)
    S <- ncol(logy1)
    N <- sum(M)
    N2 <- numL
    N3 <- nrow(logy1)

    speciesNames <- colnames(data$OTU)
    if(is.null(speciesNames)){
      speciesNames <- 1:S
    }

  }

  # create delta
  {
    {
      M_marker_df <- data_info %>% group_by(Site, Sample) %>%
        summarise(M = n())

      M_marker <- M_marker_df$M
      # names(M) <- M_df$Site

      sumM_marker <- c(0, cumsum(M_marker)[-N])

    }

    delta <- matrix(NA, N, S)

    for (s in 1:S) {

      for (i in 1:n) {

        for (m in 1:M[i]) {

          delta[sumM[i] + m,s] <-
            as.numeric(all(OTU[sumM_marker[sumM[i] + m] + 1:M_marker[sumM[i] + m], s] == 0))

        }

      }

    }

    delta[is.na(delta)] <- 0

  }

  # create covariates matrix
  {

    if(is.null(dim(OTU))){
      S <- 1
    } else {
      S <- ncol(OTU)
    }

    if(length(ordCovariates) > 0){

      X_ord <- data_info %>%
        group_by(Site) %>%
        summarise(across(all_of(ordCovariates),
                         function(x) {x[1]}))

      sitesNames <- X_ord$Site

      X_ord <- X_ord %>%
        select(-Site) %>%
        mutate_if(is.numeric, scale) %>%
        mutate(across(where(~ !is.numeric(.x)), as.factor)) %>%
        mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>%
        model.matrix(~., .)

      X_ord <- X_ord[,-1,drop=F]

      rownames(X_ord) <- sitesNames

    } else {
      X_ord <- matrix(0, n, 0)
    }

    if(length(occCovariates) > 0){

      X_psi <- data_info %>%
        group_by(Site) %>%
        summarise(across(all_of(occCovariates),
                         function(x) {x[1]})) %>%
        select(-Site) %>%
        mutate_if(is.numeric, scale) %>%
        mutate(across(where(~ !is.numeric(.x)), as.factor)) %>%
        model.matrix(~., .)

      X_psi <- X_psi[,-1,drop=F]

    } else {

      X_psi <- matrix(0, n, 0)

    }

    if(length(detCovariates) > 0){

      X_theta <- data_info %>%
        group_by(Sample) %>%
        summarise(across(all_of(detCovariates),
                         function(x) {x[1]})) %>%
        select(-Sample) %>%
        mutate_if(is.numeric, scale) %>%
        mutate(across(where(~ !is.numeric(.x)), as.factor)) %>%
        model.matrix(~., .)

    } else {

      X_theta <- matrix(1, N, 1)

    }

    ncov_psi <- ncol(X_psi)
    ncov_theta <- ncol(X_theta)
    ncov_ord <- ncol(X_ord)


  }

  # check for nas
  {
    logy_na <- is.na(logy1)
    mode(logy_na) <- "integer"
    logy1[is.na(logy1)] <- -1
  }

  # priors
  {
    prior_beta_psi <- 0
    prior_beta_psi_sd <- 1
    prior_beta_theta <- 0
    prior_beta_theta_sd <- 1
    prior_atheta0 <- 1
    prior_btheta0 <- 20
    prior_ap <- 20
    prior_bp <- 1
    prior_aq <- 1
    prior_bq <- 20
    a_sigma0 <- 1
    b_sigma0 <- 1
    a_sigma1 <- 1
    b_sigma1 <- 1
  }

  edna_dat <- list(n = n,
                   N = N,
                   N2 = N2,
                   N3 = N3,
                   M = M,
                   sumM = sumM,
                   L = L,
                   d = d,
                   sumL = sumL,
                   sumK = sumK,
                   maxL = max(L),
                   S = S,
                   X_psi = X_psi,
                   ncov_psi = ncov_psi,
                   X_theta = X_theta,
                   ncov_theta = ncov_theta,
                   X_ord = X_ord,
                   ncov_ord = ncov_ord,
                   K = K,
                   logy1 = logy1,
                   logy_na = logy_na,
                   delta = delta,
                   prior_beta_psi = prior_beta_psi,
                   prior_beta_psi_sd = prior_beta_psi_sd,
                   prior_beta_theta = prior_beta_theta,
                   prior_beta_theta_sd = prior_beta_theta_sd,
                   a_sigma0 = a_sigma0,
                   b_sigma0 = b_sigma0,
                   a_sigma1 = a_sigma1,
                   b_sigma1 = b_sigma1,
                   a_p = prior_ap,
                   b_p = prior_bp,
                   a_q = prior_aq,
                   b_q = prior_bq,
                   a_theta0 = prior_atheta0,
                   b_theta0 = prior_btheta0)

  init_fun <- function(...) list(
    theta0 = rep(0.05, S),
    p = matrix(.95, L, S),
    q = matrix(.05, L, S),
    sigma0 = .5,
    sigma1 = 1,
    mu1 = 3
  )

  print("Loading STAN")

  model <- rstan::stan_model(file = system.file("stan/code.stan", package = "occPlus"))

  vb_fit <- rstan::vb(model, data = edna_dat,
               algorithm = "meanfield",
               pars = c("beta_psi","beta_ord","beta_theta",
                        "mu1","sigma0", "sigma1","log_psi","beta0_psi",
                        "U", "LL","E","p", "pi0","q","theta0"
               ),
               init = init_fun,
               elbo_samples = 500,
               tol_rel_obj = 0.0005,
               output_samples = numSamples)

  matrix_of_draws <- as.matrix(vb_fit)

  infos <- list(
    "S" = S,
    "speciesNames" = speciesNames,
    "primerNames" = primerNames,
    "siteNames" = siteNames,
    "d" = d,
    "ncov_theta" = ncov_theta,
    "ncov_psi" = ncov_psi,
    "ncov_ord" = ncov_ord,
    "maxexplogy1" = max(exp(logy1), na.rm = T)
  )

  list(
      "vb_fit" = vb_fit,
      "matrix_of_draws" = matrix_of_draws,
      "infos" = infos,
      "X_ord" = X_ord,
      "X_theta" = X_theta,
      "X_psi" = X_psi)

}
