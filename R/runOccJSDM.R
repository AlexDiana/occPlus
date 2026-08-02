get_param <- function(params, key, default = 0) {
  if (!is.null(params[[key]])) {
    return(params[[key]])
  } else {
    return(default)
  }
}

process_covariates <- function(data_info, covariates, group_by_col, n_obs,
                               remove_intercept = FALSE,
                               spline_vars = F) {

  if (length(covariates) > 0) {

    # Process the covariates
    df <- data_info %>%
      dplyr::group_by(!!rlang::sym(group_by_col)) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(covariates), ~ dplyr::first(.x))) %>%
      dplyr::select(-dplyr::all_of(group_by_col)) %>%
      dplyr::mutate(dplyr::across(dplyr::where(~ !is.numeric(.x)), as.factor))

    if (any(is.infinite(as.matrix(df)))) stop("Infinite values (Inf or -Inf) detected in covariates.")
    if (any(is.nan(as.matrix(df)))) stop("NaN values detected in covariates.")

    # old code
    {

      #     is_numeric <- sapply(df, is.numeric)
      #
      #     names_df <- colnames(df)
      #
      #     means_df <- sapply(df, function(x) if(is.numeric(x)) mean(x, na.rm = TRUE) else NA)
      #     sd_df   <- sapply(df, function(x) if(is.numeric(x)) sd(x, na.rm = TRUE) else NA)
      #
      #     if (any(sd_df == 0, na.rm = TRUE)) {
      #       zero_var_cols <- names_df[which(sd_df == 0)]
      #       stop(paste("The following covariates have constant values:",
      #                  paste(zero_var_cols, collapse = ", ")))
      #     }
      #
      #     cat_levels <- list()
      #     for (col in 1:ncol(df)) {
      #       if(is_numeric[col]){
      #         cat_levels[[col]] <- NA
      #       }else{
      #         cat_levels[[col]] <- levels(as.factor(df[[col]]))
      #     }
      #   }
      #
      #   list_matrix <- list(
      #     "names_df" = names_df,
      #     "mean_df" = means_df,
      #     "sd_df" = sd_df,
      #     "cat_levels" = cat_levels,
      #     "is_numeric" = is_numeric
      #   )
      #
      # out_matrix <- transformCovariatesMatrix(df, list_matrix, remove_intercept)
    }

    list_X <- create_covariates_matrix(
      df,
      spline_vars = spline_vars,
      remove_intercept = remove_intercept)
    out_matrix <- list_X$X
    list_matrix <- list_X$list_matrix
    X0 <- list_X$X0

  } else {

    if (!remove_intercept) {
      out_matrix <- (matrix(1, n_obs, 1))
    } else {
      out_matrix <- (matrix(0, n_obs, 0))
    }

    list_matrix <- NULL
    X0 <- NULL
  }

  list("df" = out_matrix,
       "list_matrix" = list_matrix,
       "X0" = X0)

}

createDataIdx <- function(n, M, P, K, twostage, primerId_p = NULL){
  # P is now a vector with one entry per sample (length N = sum(M)), giving
  # the number of primers actually used for that sample -- it no longer has
  # to be constant across samples. primerId_p (length N2 = sum(P)) gives the
  # *global* primer identity (an index into 1:length(primerNames)) for each
  # (sample, primer) block, in the same order samples/primers are visited
  # below; this is what lets p/q be looked up by actual primer identity
  # rather than by a primer's position within a given sample's block.

  N <- sum(M)
  sumM <- c(0, cumsum(M)[-n])
  idx_z_w <- rep(NA, N)

  if(twostage){
    N2 <- sum(P)
    N3 <- sum(K)
    sumP <- c(0, cumsum(P)[-N])
    sumK <- c(0, cumsum(K)[-N2])

    idx_z_p <- rep(NA, N2)
    idx_w_p <- rep(NA, N2)
    idx_z_k <- rep(NA, N3)
    idx_p_k <- rep(NA, N3)
    idx_w_k <- rep(NA, N3)
  } else {
    idx_z_p <- NULL
    idx_w_p <- NULL
    idx_z_k <- NULL
    idx_p_k <- NULL
    idx_w_k <- NULL
  }

  idx_z <- 1
  idx_w <- 1
  idx_p <- 1
  idx_k <- 1

  for (i in 1:n) {
    for (m in 1:M[i]) {
      idx_z_w[idx_w] <- i

      if(twostage){

        sample_idx <- sumM[i] + m

        for (pp in seq_len(P[sample_idx])) {
          idx_z_p[idx_p] <- i

          primer_id <- if (!is.null(primerId_p)) primerId_p[idx_p] else pp

          for (k in 1:K[sumP[sample_idx] + pp]) {

            idx_z_k[idx_k] <- idx_z
            idx_w_k[idx_k] <- idx_w
            idx_p_k[idx_k] <- primer_id

            idx_k <- idx_k + 1
          }
          idx_z_p[idx_p] <- idx_z
          idx_w_p[idx_p] <- idx_w

          idx_p <- idx_p + 1
        }

      }


      idx_w <- idx_w + 1
    }
    idx_z <- idx_z + 1
  }

  list("idx_z_w" = idx_z_w,
       "idx_z_p" = idx_z_p,
       "idx_w_p" = idx_w_p,
       "idx_p_k" = idx_p_k,
       "idx_w_k" = idx_w_k,
       "idx_z_k" = idx_z_k)

}

inferDataModel <- function(data){

  stages_type <- NULL

  data_info <- data$info
  OTU <- data$OTU

  site_present <- !is.null(data_info$Site)
  sample_present <- !is.null(data_info$Sample)
  primer_present <- !is.null(data_info$Primer)

  if(site_present){ # site column present

    if(all(!duplicated(data_info$Site))){ # no repeated observation at each site

      stages_type <- "no_stage"
      message("No repeated observation for any site")

    } else { # repeated observation at each site

      if(!sample_present){ # sample column missing

        message("Sample column missing, but detected repeated observations for some sites")
        stages_type <- "one_stage"

      } else { # sample column present

        if(all(!duplicated(data_info$Sample))){ # no repeated observation in each sample

          stages_type <- "one_stage"
          message("No repeated observation for any sample")

        } else {

          stages_type <- "two_stage"
          message("Detected repeated observation for some sample")

        }

      }

    }
  } else { # site column missing
    stages_type <- "no_stage"
    message("Site column missing")
  }

  # infer the data type
  {
    # if(stages_t)
    if(mode(OTU) != "numeric"){
      stop("Data in OTU need to be numeric")
    } else if(all(round(OTU) == OTU, na.rm = T)) { ## TODO: write that the data are integer
      if(all(unique(OTU) %in% c(0,1,NA), na.rm = T)){
        message("Only 0 and 1 detected")
        data_type <- "binary"
      } else {
        message("Count data detected")
        data_type <- "counts"
      }
    } else {
      message("Continuous observations detected")
      data_type <- "continuous"
    }
  }

  if(stages_type == "no_stage"){
    if(data_type == "binary"){
      model <- "binary"
    } else if(data_type == "continuous"){
      model <- "continuous"
    } else if(data_type == "counts"){
      model <- "counts"
    }
  } else if(stages_type == "one_stage"){
    model <- "occupancy"
  } else if(stages_type == "two_stage"){
    model <- "two_stage"
  } else {
    stop("No model recognised")
  }

  if(model == "binary") message(paste("occJSDM has inferred one-stage binary data"))
  if(model == "continuous") message(paste("occJSDM has inferred one-stage continuous data"))
  if(model == "counts") message(paste("occJSDM has inferred one-stage counts data"))
  if(model == "occupancy") message(paste("occJSDM has inferred occupancy data"))
  if(model == "two_stage") message(paste("occJSDM has inferred two stage (eDNA style) data"))


  message("Check your data if this is not what you were expecting!")

  if(model == "counts") stop("Counts model not supported yet :(")

  model

}

create_waic_quantities <- function(n_obs){

  M2 <- rep(0, n_obs)
  mean_log <- rep(0, n_obs)
  mean_lik <- rep(0, n_obs)

  list("M2" = M2,
       "mean_log" = mean_log,
       "mean_lik" = mean_lik)
}

#' runOccJSDM
#'
#' Fit an occupancy joint species distribution model (occJSDM), optionally
#' accounting for environmental and detection covariates, species traits,
#' spatial autocorrelation, and for eDNA-style data, a two-stage observation
#' process (false-negative and false-positive detection errors in the field
#' and in the lab).
#'
#' @details
#' The model actually fit is inferred automatically from the shape of
#' \code{data} (via the internal \code{inferDataModel()}): whether
#' \code{data$info} has repeated observations per \code{Site} (and, within a
#' \code{Site}, per \code{Sample}) determines whether the data are treated
#' as single-visit JSDM data, classical (one-stage) occupancy data, or
#' two-stage (eDNA-style) occupancy/detection data. A message is printed
#' describing which model was inferred -- check that it matches expectations
#' for your data.
#'
#' @param data A list with two elements:
#' \describe{
#'   \item{info}{A data.frame with one row per observation (i.e. per PCR
#'   replicate for two-stage eDNA data), containing \code{Site}, \code{Sample},
#'   and \code{Primer} id columns (as applicable to the inferred model) plus
#'   the occupancy/collection/spatial covariate columns named in
#'   \code{occCovariates}/\code{collCovariates}/\code{spatCovariates}. See the
#'   vignette for the full expected shape.}
#'   \item{OTU}{A matrix of dimension (N x S), where N is
#'   \code{nrow(data$info)} and S is the number of species, containing the
#'   number of reads of each species in each observation.}
#'   \item{traits}{(Optional) species (rows) by trait (columns) matrix of
#' species traits, matched to \code{colnames(data$OTU)} by row name. }
#' }
#' @param listParams (Optional) list of model-size hyperparameters:
#' \describe{
#'   \item{n_factors}{Number of latent JSDM factors used for the residual
#'   species covariance. Capped to the number of species if larger.}
#'   \item{n_lattrait}{Number of latent trait factors (\code{gt}). Default
#'   \code{2}.}
#'   \item{n_supportpoints}{Number of spatial support points used to
#'   approximate the Gaussian process over site coordinates when
#'   \code{spatCovariates} is non-empty. Defaults to
#'   \code{getDefaultSupportPoints(n)}.}
#' }
#' @param threshold Threshold used to truncate the reads to binary detections
#' for occupancy/two-stage models. Reads greater than or equal to the
#' threshold are considered a detection (default \code{1}). Must be
#' \code{>= 1} -- continuous detection modeling (\code{threshold = 0}) is not
#' supported; any value less than \code{1} raises an error.
#' @param occCovariates (Optional) vector of the names of the occupancy
#' covariates (site-level). Names should match column names in
#' \code{data$info}.
#' @param collCovariates (Optional) vector of the names of the collection
#' covariates (sample-level, relevant for occupancy/two-stage models only).
#' Names should match column names in \code{data$info}.
#' @param spatCovariates (Optional) vector of the names of the site
#' coordinate columns used to fit a spatially autocorrelated random field
#' over occupancy. Names should match column names in \code{data$info}. If
#' omitted, no spatial field is fit.
#' @param MCMCparams (Optional) list of MCMC settings: \code{nchain} (number
#' of chains), \code{nburn} (number of burn-in iterations to discard),
#' \code{niter} (number of post-burn-in iterations to keep), and
#' \code{nthin} (thinning interval, default \code{1}). Defaults to
#' \code{list(nchain = 2, nburn = 5000, niter = 5000, nthin = 1)}.
#' @param summarisedLatentPresences (Optional) logical, default \code{TRUE}.
#' If \code{TRUE}, the returned latent-presence/probability arrays
#' (\code{z_output}, \code{w_output}, \code{psi_output}, \code{theta_output})
#' are collapsed to posterior means, which keeps the fitted object small. If
#' \code{FALSE}, the full posterior samples are kept instead, which allows
#' computing credible intervals for conditional/predictive occupancy
#' probabilities but produces a much larger object.
#' @param listPriors (Optional) list of prior hyperparameters. Currently
#' supported: \code{a_p}/\code{b_p} (Beta prior shape parameters for the true
#' positive lab stage rate,
#' default \code{5}/\code{1}), \code{a_q}/\code{b_q} (Beta prior shape parameters
#' for the false positive
#' lab stage rate, default \code{1}/\code{20})  and
#' \code{a_theta0}/\code{b_theta0} (Beta prior shape parameters for the
#' false positive field stage rate, default \code{1}/\code{20}), and
#' \code{b_betatheta_slope_var} (prior variance on the collection-covariate
#' slopes, default \code{2}). The last of these is provisional: it exists so
#' the slope prior can be varied without editing the source, and the default
#' reproduces the previous hard-coded behaviour exactly. Its value is still
#' under review, so treat a non-default setting as a diagnostic rather than a
#' recommended configuration.
#'
#' @return A list with:
#' \describe{
#'   \item{results_output}{Posterior samples/summaries, including
#'   \code{jsdm_output} (JSDM coefficient posteriors), \code{beta_theta_output},
#'   \code{p_output}, \code{q_output}, \code{theta0_output} (detection-process
#'   posteriors), \code{WAIC}, and (per \code{summarisedLatentPresences})
#'   either posterior means or full samples of \code{z_output} (latent
#'   occupancy), \code{w_output} (latent collection), \code{psi_output}
#'   (occupancy probabilities), and \code{theta_output} (collection
#'   probabilities).}
#'   \item{infos}{Metadata about the fitted model (e.g. \code{S}, \code{n},
#'   \code{M}, \code{K}, \code{speciesNames}, \code{primerNames},
#'   \code{siteNames}, \code{ncov_psi}, \code{ncov_theta}, the \code{OTU}
#'   matrix used for fitting, and standardisation info for \code{X_psi}/
#'   \code{Xs}), used internally by the plotting/output helpers (e.g.
#'   \code{\link{plotOccupancyRates}}, \code{\link{plotDetectionRates}}).}
#'   \item{Tr}{The species trait matrix actually used to fit the model.}
#'   \item{X_theta}{The collection covariate design matrix (including
#'   intercept) used to fit the model.}
#'   \item{X_s}{The site coordinates used for the spatial field.}
#'   \item{X_psi}{The (standardised) occupancy covariate design matrix used
#'   to fit the model.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage
#' fitmodel <- runOccJSDM(data,
#' listParams = list(n_factors = 2),
#' occCovariates = c("X_psi.1", "X_psi.2"),
#' collCovariates = c("X_theta"))
#' }
#'
#' @export
#' @import dplyr
#'
runOccJSDM <- function(data,
                       listParams = list(),
                       threshold = 1,
                       occCovariates = c(),
                       collCovariates = c(),
                       spatCovariates = c(),
                       MCMCparams = list(nchain = 2,
                                         nburn = 5000,
                                         niter = 5000,
                                         nthin = 1),
                       summarisedLatentPresences = T,
                       listPriors = list()){

  {
    # Debugging scaffold: a set of arguments to paste in when stepping through
    # this function by hand. Everything here must stay commented out. In 8f9f315
    # part of it was left live, which overwrote every argument with a dataset
    # that does not exist in the package.
    #
    data = occ_data_effort
    listParams = list(n_factors = 3)
    threshold = 1
    occCovariates = c("season_year","z_prop_closed")
    collCovariates = c("predator_season_year", "z_log_effort_m_total")
    spatCovariates = c()
    MCMCparams = list(
    nchain = 2,
    nburn = 5000,
    niter = 5000,
    nthin = 1 )
    listPriors = list(
    a_theta0 = 1,
    b_theta0 = 100)# listParams = list(n_factors = 3, splineVars = F)
    summarisedLatentPresences = T
    # threshold = 1
    # occCovariates = c( "X_psi.EnvCov.1", "X_psi.EnvCov.2", "X_psi.EnvCov.3",
    #                    "X_psi.EnvCov.4", "X_psi.EnvCov.5" ,"X_psi.EnvCov.6")
    # collCovariates = c("X_theta.1","X_theta.2")
    # spatCovariates = c("Xs.1","Xs.2")
    # MCMCparams = list(
    #   nchain = 2,
    #   nburn = 200,
    #   niter = 200,
    #   nthin = 1 )
    # listPriors = list(
    #   a_theta0 = 1,
    #   b_theta0 = 100)
    # listPriors = list()
  }

  # data structure infer
  {
    if(is.null(data$info) || is.null(data$OTU)){
      stop("data_info or OTU missing")
    } else {
      data_info <- as.data.frame(data$info)
      OTU <- data$OTU
    }

    if(nrow(data_info) != nrow(OTU)){
      stop("OTU and data_info cannot have different number of rows")
    }

    model <- inferDataModel(data)

    if(model %in% c("binary","occupancy","two_stage")){
      jsdmModel <- "binary"
    } else if(model == "continuous") {
      jsdmModel <- "continuous"
    }  else if(model == "counts") {
      jsdmModel <- "counts"
    }

    if(is.null(data_info$Site)){
      data_info$Site <- 1:nrow(data_info)
    }

    if(model == "occupancy" & is.null(data_info$Sample)){
      data_info$Sample <- 1:nrow(data_info)
    }

    data_info <- data_info %>%
      mutate(rownum = row_number())

    if(model == "occupancy"){
      data_info <- data_info %>%
        dplyr::arrange(Site, Sample)
    } else if(model == "two_stage"){
      data_info <- data_info %>%
        dplyr::arrange(Site, Sample, Primer)

    }

    OTU <- OTU[data_info$rownum, , drop = FALSE]

    data_info <- data_info %>% select(-rownum)
  }

  # read OTU data
  {
    y <- OTU
    if(is.null(dim(y))){
      S <- 1
    } else {
      S <- ncol(y)
    }

    if(model %in% c("occupancy","two_stage")){

      # truncate data

      if (!is.numeric(threshold) || length(threshold) != 1) {
        stop("'threshold' must be a single numeric value.")
      }
      if (threshold <= 0) { # You already note threshold > 0 in your code, enforce it cleanly!
        stop("'threshold' must be strictly greater than 0.")
      }

      y[y >= threshold] <- 1
      y[y < threshold] <- 0

    }

    if(model %in% c("binary","continuous","counts")){
      z <- y
    } else if(model == "occupancy"){
      w <- y
    }

    # check for nas
    {
      y_NA <- is.na(y)
      mode(y_NA) <- "integer"

      if(sum(y_NA) > 0 & model != "two_stage"){
        stop("NAs are allowed only in the two-stage model", .call = F)
      }
    }

    speciesNames <- colnames(data$OTU)
    if(is.null(speciesNames)){
      speciesNames <- 1:S
    }

    if (any(duplicated(speciesNames))) stop("Duplicated species names found in OTU matrix.")

    n_obs <- nrow(y)


  }

  # data checks
  {

    if(any(is.na(data_info$Site)) |
       any(is.na(data_info$Sample)) |
       any(is.na(data_info$Primer))){
      stop("NA in Site, Sample or Primer columns")
    }

    if(is.null(occCovariates)) occCovariates <- c()
    if(is.null(collCovariates)) collCovariates <- c()
    if(is.null(spatCovariates)) spatCovariates <- c()

    if(!(length(spatCovariates) %in% c(0,2))) stop("Two columns needed for the spatial covariates")

    if(!all(sapply(data_info[,spatCovariates], is.numeric))) {
      stop("Non numeric columns in spatial covariates")
    }

    if(!all(c(occCovariates, collCovariates, spatCovariates) %in% colnames(data$info))){
      stop("Covariate names provided not in data$info")
    }

  }

  # clean the data
  {

    # samples per site
    {
      if(model %in% c("occupancy","two_stage")){
        M_df <- data_info %>%
          dplyr::group_by(Site, Sample) %>%
          slice(1) %>%
          dplyr::group_by(Site) %>%
          dplyr::summarise(M = n(),
                           .groups = "keep")

        M <- M_df$M
        names(M) <- M_df$Site
        siteNames <- M_df$Site

        n <- length(M)
        N <- sum(M)

        sumM <- c(0, cumsum(M)[-n])

        data_info$SiteSample <- paste(data_info$Site,data_info$Sample, sep = "-")

      } else {
        n <- nrow(data_info)
        siteNames <- 1:n
        M <- NULL
      }

    }

    # marker per samples
    {
      if(model == "two_stage"){
        # number of markers

        P_df <- data_info %>%
          dplyr::group_by(Site, Sample, Primer) %>%
          dplyr::slice(1) %>%
          dplyr::group_by(Site, Sample) %>%
          dplyr::summarise(P = n(),
                           .groups = "keep")

        # P is now allowed to differ across samples: each sample can use a
        # different subset (and number) of primers.
        P <- P_df$P
        names(P) <- P_df$Sample

        # primerNames is the *global* set of distinct primers across the
        # whole dataset; p/q are indexed by primer identity (row into this
        # vector), not by a primer's position within a given sample.
        primerNames <- sort(unique(data_info$Primer))
        maxP <- length(primerNames)

        N2 <- sum(P)
        sumP <- c(0, cumsum(P)[-N])

      } else {
        primerNames <- NULL
      }
    }

    # pcr per marker
    {
      if(model == "two_stage"){
        data_K <- data_info %>%
          dplyr::group_by(Site, Sample, Primer) %>%
          dplyr::summarise(K = n(),
                           dplyr::across(contains("Species"),function(x){sum(x > 0)}),
                           .groups = "keep"
          ) %>%
          dplyr::ungroup()

        K <- data_K$K
        sumK <- c(0, cumsum(K)[-N2])

        N3 <- sum(K)

        # global primer identity (index into primerNames) for each
        # (Site, Sample, Primer) block -- data_K is at the same
        # (Site, Sample, Primer) granularity and in the same order (both
        # derived from the already Site/Sample/Primer-sorted data_info), so
        # this aligns row-for-row with K/sumP.
        primerId_p <- match(data_K$Primer, primerNames)

      } else {
        K <- NULL
        primerId_p <- NULL
      }

    }

    if(model %in% c("occupancy","two_stage")){
      list_idx <- createDataIdx(n, M, P, K, model == "two_stage", primerId_p = primerId_p)
      idx_z_w <- list_idx$idx_z_w
      idx_z_k <- list_idx$idx_z_k
      idx_w_p <- list_idx$idx_w_p
      idx_z_p <- list_idx$idx_z_p
      idx_w_k <- list_idx$idx_w_k
      idx_p_k <- list_idx$idx_p_k
    } else {
      list_idx <- NULL
    }

  }

  # create covariates matrix
  {

    splineVars <- get_param(listParams, "splineVars", F)

    # For occupancy covariates (group by Site, includes intercept)
    {
      list_X_psi <- process_covariates(data_info, occCovariates, "Site", n,
                                  remove_intercept = TRUE,
                                  spline_vars = splineVars)
      X_psi <- list_X_psi$df
      list_Xpsi_mat <- list_X_psi$list_matrix
      X0_psi <- list_X_psi$X0
    }

    # For the spatial field
    {
      list_Xs <- process_covariates(data_info, spatCovariates, "Site", n,
                               remove_intercept = TRUE,
                               spline_vars = F)
      Xs <- list_Xs$df
      list_Xs_mat <- list_Xs$list_matrix
    }

    # For collection covariates (group by Sample, includes intercept)
    {
      if(model %in% c("occupancy","two_stage")){
        list_X_theta <- process_covariates(data_info, collCovariates, "SiteSample",
                                           N, remove_intercept = FALSE,
                                           spline_vars = F)
        X_theta <- list_X_theta$df
        list_X_theta_mat <- list_X_theta$list_matrix
      } else {
        X_theta <- NULL
        list_X_theta_mat <- NULL
      }
    }

    # For the traits
    {
      if (!is.null(data$traits) && any(duplicated(rownames(data$traits)))) {
        stop("Duplicated species names found in traits matrix rownames.")
      }

      if(!is.null(data$traits)){
        speciesNamesInTraitsMatrix <- rownames(data$traits)
        if(!all(speciesNames %in% speciesNamesInTraitsMatrix)){
          stop("Species names in OTU not present in traits matrix")
        }

        idx_speciesNames <- match(speciesNames, speciesNamesInTraitsMatrix)
        Tr <- data$traits
        Tr <- Tr[idx_speciesNames,]
        Tr <- as.matrix(Tr)
        traitsNames <- colnames(Tr)
      } else {
        Tr <- matrix(NA, S, 0)
      }

    }

    ncov_psi <- ncol(X_psi)
    ncov_theta <- ncol(X_theta)
    g <- ncol(Tr)

  }

  # set unknonwn parameters
  {
    d <- get_param(listParams, "n_factors")

    if(d > ncol(OTU)){
      message("More species than factors. The number of factors will be capped to the
            number of species")
      d <- ncol(OTU)
    }

    gt_default <- floor(sqrt(min(S, ncov_psi)))
    gt <- get_param(listParams, "n_lattrait", gt_default)

    if(ncol(Xs) > 0){
      ns <- nrow(unique(Xs))
      ps <- get_param(listParams, "n_supportpoints", getDefaultSupportPoints(ns))

    } else {
      ps <- 0
    }


  }

  # priors
  {
    a_theta0 <- ifelse(is.null(listPriors$a_theta0), 1, listPriors$a_theta0)
    b_theta0 <- ifelse(is.null(listPriors$b_theta0), 20, listPriors$b_theta0)
    a_p <- ifelse(is.null(listPriors$a_p), 5, listPriors$a_p)
    b_p <- ifelse(is.null(listPriors$b_p), 1, listPriors$b_p)
    a_q <- ifelse(is.null(listPriors$a_q), 1, listPriors$a_q)
    b_q <- ifelse(is.null(listPriors$b_q), 20, listPriors$b_q)

    if(model %in% c("occupancy","two_stage")){
      b_betatheta <- rep(0, ncov_theta)

      # Slope prior variance, previously hard-coded. Exposed here so it can be
      # tuned/tested (see the b_betatheta variance-decision item in
      # TODO.md group B, and PLAN.md 13.8) without a
      # source edit each time. Default of 2 is unchanged from before this
      # override existed. ALEX TO REVIEW before treating a non-default value
      # as a real fix rather than a diagnostic.
      b_betatheta_slope_var <- ifelse(is.null(listPriors$b_betatheta_slope_var),
                                      2, listPriors$b_betatheta_slope_var)
      B_betatheta <- diag(b_betatheta_slope_var, nrow = ncov_theta)

      prior_beta_theta <- 0
      prior_beta_theta_sd <- 1
      b_betatheta[1] <- prior_beta_theta
      B_betatheta[1,1] <- prior_beta_theta_sd^2
    }

  }

  # chain parameters
  {
    nchain <- get_param(MCMCparams, "nchain", 2)
    nburn <- get_param(MCMCparams, "nburn", 1000)
    niter <- get_param(MCMCparams, "niter", 1000)
    nthin <- get_param(MCMCparams, "nthin", 1)
  }

  # precompute spatial quantities
  {
    # Spatial covariates matrix
    list_Xs <- computeSpatialSummaries(Xs, ps, maxPoints = 5)
    Xs_centers <- list_Xs$Xs_centers
    Xs_index <- list_Xs$Xs_index
    X_s_centers <- list_Xs$X_s_centers
    X_tilde <- list_Xs$X_tilde
    X_s <- list_Xs$X_s
    ps <- list_Xs$ps

    length_grid_ls <- 10
    l_s_grid <- seq(0.01, 0.3, length.out = length_grid_ls)
    list_SoRSummaries <- precomputeSORmatrices(l_s_grid, list_Xs)
  }

  # precompute jsdm params
  {
    list_data <- list(
      "X" = X_psi,
      "Tr" = Tr)

    a_sigmab <- 10; b_sigmab <- 1
    a_sigmabs <- 10; b_sigmabs <- 1
    a_sigmah <- 10; b_sigmah <- 1
    a_tau <- 5; b_tau <- 5
    a_l_s <- 1; b_l_s <- 1

    list_priors <- list(
      "a_sigmab" = a_sigmab,
      "b_sigmab" = b_sigmab,
      "a_sigmabs" = a_sigmabs,
      "b_sigmabs" = b_sigmabs,
      "a_sigmah" = a_sigmah,
      "b_sigmah" = b_sigmah,
      "a_tau" = a_tau,
      "b_tau" = b_tau,
      "a_l_s" = a_l_s,
      "b_l_s" = b_l_s
    )
  }

  # WAIC calculation
  {

    list_waic_jsdm <- create_waic_quantities(n * S)

    if(model %in% c("occupancy","two_stage")){
      list_waic_w <- create_waic_quantities(N * S)
    }

    if(model %in% "two_stage"){
      list_waic_y <- create_waic_quantities(N3 * S)
    }

    currentWAICiter <- 1

  }

  # chain output
  {
    if(model %in% c("occupancy","two_stage")){
      beta_theta_output <- array(NA, dim = c(ncov_theta, S, niter, nchain))
      theta0_output <- array(NA, dim = c(S, niter, nchain))
      if(summarisedLatentPresences){
        z_output_mean <- matrix(0, n, S)
        psi_output_mean <- matrix(0, n, S)
      } else {
        z_output <- array(NA, dim = c(n, S, niter, nchain))
        psi_output <- array(NA, dim = c(n, S, niter, nchain))
      }
      theta_output_mean <- matrix(0, N, S)

    } else {
      beta_theta_output <- NULL
      theta0_output <- NULL
      if(summarisedLatentPresences){
        z_output_mean <- NULL
        psi_output_mean <- NULL
      } else {
        z_output <- NULL
        psi_output <- NULL
      }
      theta_output_mean <- NULL

    }

    if(model == "two_stage"){
      p_output <- array(NA, dim = c(maxP, S, niter, nchain))
      q_output <- array(NA, dim = c(maxP, S, niter, nchain))
      w_output_mean <- matrix(0, N, S)
    } else {
      p_output <- NULL
      q_output <- NULL
      w_output_mean <- NULL
    }

    # jsdm params

    {
      B0_output <- array(NA, dim = c(S, niter, nchain))
      B_output <- array(NA, dim = c(ncov_psi, S, niter, nchain))
      G_output <- array(NA, dim = c(g, ncov_psi, niter, nchain))
      A_output <- array(NA, dim = c(S, gt, niter, nchain))
      C_output <- array(NA, dim = c(gt, ncov_psi, niter, nchain))
      Bs_output <- array(NA, dim = c(ps, S, niter, nchain))
      Gs_output <- array(NA, dim = c(g, ps, niter, nchain))
      As_output <- array(NA, dim = c(S, gt, niter, nchain))
      Cs_output <- array(NA, dim = c(gt, ps, niter, nchain))
      U_output  <- array(NA, dim = c(n, d, niter, nchain))
      L_output  <- array(NA, dim = c(d, S, niter, nchain))
      tau_output <- array(NA, dim = c(S, niter, nchain))
      sigmab_output <- array(NA, dim = c(niter, nchain))
      sigmabs_output <- array(NA, dim = c(niter, nchain))
      sigmah_output <- array(NA, dim = c(niter, nchain))
      idx_ls_output <- array(NA, dim = c(niter, nchain))
      varPart_output <- array(NA, dim = c(S, 4, niter, nchain))
    }

  }

  # run MCMC
  message("Running MCMC")

  # Seed the C++ samplers from R's RNG. The samplers use per-thread mt19937
  # engines rather than R's RNG (which is a single global and unsafe to call
  # from inside the OpenMP regions), so they have to be seeded explicitly --
  # without this, set.seed() has no effect on the fit. Drawing the seed here
  # rather than using a constant also makes two consecutive fits in the same
  # session independent. Chains deliberately share the draw: each continues the
  # stream where the previous one stopped, so they differ from one another
  # while the fit as a whole stays reproducible.
  setOccJSDMSeed(sample.int(.Machine$integer.max, 1L))

  for (chain in 1:nchain) {

    # chain output
    {
      if(model %in% c("occupancy","two_stage")){
        beta_theta_output_chain <- array(NA, dim = c(ncov_theta, S, niter))
        theta0_output_chain <- array(NA, dim = c(S, niter))
        if(!summarisedLatentPresences){
          z_output_chain <- array(NA, dim = c(n, S, niter))
          psi_output_chain <- array(NA, dim = c(n, S, niter))
          theta_output_chain <- array(NA, dim = c(N, S, niter))
        }
      }

      if(model == "two_stage"){
        p_output_chain <- array(NA, dim = c(maxP, S, niter))
        q_output_chain <- array(NA, dim = c(maxP, S, niter))
      }

      # jsdm params
      {
        B0_output_chain <- array(NA, dim = c(S, niter))
        B_output_chain <- array(NA, dim = c(ncov_psi, S, niter))
        G_output_chain <- array(NA, dim = c(g, ncov_psi, niter))
        A_output_chain <- array(NA, dim = c(S, gt, niter))
        C_output_chain <- array(NA, dim = c(gt, ncov_psi, niter))
        Bs_output_chain <- array(NA, dim = c(ps, S, niter))
        Gs_output_chain <- array(NA, dim = c(g, ps, niter))
        As_output_chain <- array(NA, dim = c(S, gt, niter))
        Cs_output_chain <- array(NA, dim = c(gt, ps, niter))
        U_output_chain <- array(NA, dim = c(n, d, niter))
        L_output_chain <- array(NA, dim = c(d, S, niter))
        sigmab_output_chain <- rep(NA, niter)
        sigmabs_output_chain <- rep(NA, niter)
        sigmah_output_chain <- rep(NA, niter)
        idx_ls_output_chain <- rep(NA, niter)
        tau_output_chain <- matrix(NA, S, niter)
        varPart_output_chain <-  array(NA, dim = c(S, 4, niter))
      }

    }

    # starting values
    {
      if(model == "two_stage"){

        w <- matrix(NA, N, S)

        for (s in 1:S) {
          for (i in 1:n) {
            for (m in 1:M[i]) {
              idx_im1 <- sumK[sumP[sumM[i] + m] + 1] + 1
              idx_im2 <- sumK[sumP[sumM[i] + m] + P[sumM[i] + m]] +
                K[sumP[sumM[i] + m] + P[sumM[i] + m]]
              y_subset <- y[idx_im1:idx_im2,s]
              y_subset <- y_subset[!is.na(y_subset)]
              if(length(y_subset) > 0){
                w[sumM[i] + m,s] <- as.numeric(any(y_subset > 0))
              } else {
                w[sumM[i] + m,s] <- 1
              }
            }
          }
        }

        c_imk <- array(as.numeric(y > 0), dim = dim(y), dimnames = dimnames(y))

        p <- matrix(.9, maxP, S)
        q <- matrix(.05, maxP, S)

        y_pos <- (y > 0)

      }

      if(model %in% c("occupancy","two_stage")){

        z <- matrix(NA, n, S)

        for (s in 1:S) {
          for (i in 1:n) {
            idx_i1 <- sumM[i] + 1
            idx_i2 <- sumM[i] + M[i]

            z[i,s] <- as.numeric(any(w[idx_i1:idx_i2, s] > 0))
          }
        }

        beta_theta <- matrix(0, ncov_theta, S)
        theta <- computeTheta(X_theta, beta_theta)
        theta0 <- rep(.05, S)
      }

      # jsdmParams
      {
        B0 <- rep(0, S)
        B <- matrix(0, ncov_psi, S)
        Bs <- matrix(0, ps, S)
        L <- matrix(1, d, S)
        diag(L) <- 1
        L[lower.tri(L)] <- 0
        G <- matrix(0, g, ncov_psi)
        C <- matrix(1, gt, ncov_psi)
        diag(C) <- 1
        C[lower.tri(C)] <- 0
        A <- matrix(0, S, gt)
        Gs <- matrix(0, g, ps)
        Cs <- matrix(1, gt, ps)
        diag(Cs) <- 1
        Cs[lower.tri(Cs)] <- 0
        As <- matrix(0, S, gt)
        U <- matrix(0, n, d)
        sigma_b <- 1
        sigma_bs <- .001
        sigma_h <- 1
        idx_ls <- 3 # dim(list_SoRSummaries$Ks_all)[3][5]
        tau <- rep(1, S)

        Bt <- t(B) - computeBtcoef(G, Tr, A, C, matrix(0, S, ncov_psi))
        Bst <- t(Bs) - computeBtcoef(Gs, Tr, As, Cs, matrix(0, S, ps))

        Ks <- list_SoRSummaries$Ks_all[,,idx_ls]

        list_jSDMparams <- list(
          "B0" = B0,
          "B" = B,
          "G" = G,
          "A" = A,
          "C" = C,
          "Bt" = Bt,
          "Bs" = Bs,
          "Gs" = Gs,
          "As" = As,
          "Cs" = Cs,
          "Bst" = Bst,
          "U" = U,
          "L" = L,
          "sigma_b" = sigma_b,
          "sigma_bs" = sigma_bs,
          "sigma_h" = sigma_h,
          "idx_ls" = idx_ls,
          "tau" = tau
        )

        list_psiCoef <- computePsiCoef(
          X_psi, Ks, list_Xs$Xs_centers, Tr,
          B0, G, A, C, Bt,
          Gs, As, Cs, Bst,
          U, L)

      }

      psi <- logistic(list_psiCoef$eta)

    }

    for (iter in 1:(nburn + niter * nthin)) {

      if(iter > nburn){
        currentIter <- (iter - nburn) / nthin
        if(currentIter %% 100 == 0){
          message(paste0("Chain ", chain, " - Iteration ",currentIter))
        }
      } else {
        if(iter %% 100 == 0){
          message(paste0("Chain ", chain, " - Burn in Iteration ",iter))
        }
      }

      # sample z
      if(model %in% c("occupancy","two_stage")){
        z <- sample_z_cpp_parallel(w, psi, theta, theta0, M, sumM)
      }

      # sample jsdm coef
      {
        list_data$z <- z

        list_jSDMparams <- update_jSDMcoef(
          list_data,
          list_jSDMparams,
          list_priors,
          list_Xs,
          list_SoRSummaries,
          model = jsdmModel,
          computeVarPart = (iter > nburn & (iter - nburn) %% nthin == 0)
        )
        psi <- list_jSDMparams$psi
        eta <- list_jSDMparams$eta
      }

      # sample w
      if(model == "two_stage"){
        if(threshold > 0){

          w <- sample_w_cim_cipp_parallel(y, y_NA, theta, theta0, p, q,
                                 M, K, sumP, sumM, sumK, P, primerId_p, maxP, z)

        }
      }

      # sample cimk
      if(model == "two_stage"){
        if (threshold > 0) {

          w_all <- w[idx_w_k,,drop=F]

          # faster way to assign c_imk to 1 if logy1 > 0 for w_all = 1 and to 2
          # if log1 > 0 when w_all = 0
          c_imk <- y_pos * (2 - (w_all == 1))

        }
      }

      # sample theta
      if(model %in% c("occupancy","two_stage")){
        beta_theta <- sample_betatheta_cpp_parallel(w, z, beta_theta, idx_z_w, X_theta,
                                           b_betatheta, B_betatheta)
        theta <- computeTheta(X_theta, beta_theta)
      }

      # sample pq
      if(model == "two_stage"){
        list_pq <- sample_pq_cpp(c_imk, y_NA, w, idx_p_k, idx_w_k, maxP, a_p, b_p, a_q, b_q)
        p <- list_pq$p
        q <- list_pq$q
      }

      # sample theta0
      if(model %in% c("occupancy","two_stage")){
        theta0 <- sample_theta0(z, w, idx_z_w, a_theta0, b_theta0)
      }

      {

        if(iter > nburn & (iter - nburn) %% nthin == 0){
          currentIter <- (iter - nburn) / nthin

          if(model %in% c("occupancy","two_stage")){
            beta_theta_output_chain[,,currentIter] <- beta_theta
            theta0_output_chain[,currentIter] <- theta0
            if(summarisedLatentPresences){
              z_output_mean <- z_output_mean +
                (1 / (niter * nchain)) * z
              psi_output_mean <- psi_output_mean +
                (1 / (niter * nchain)) * psi
            } else {
              z_output_chain[,,currentIter] <- z
              psi_output_chain[,,currentIter] <- psi
            }
            theta_output_mean <- theta_output_mean +
              (1 / (niter * nchain)) * theta
          }

          if(model == "two_stage"){
            p_output_chain[,,currentIter] <- p
            q_output_chain[,,currentIter] <- q
            w_output_mean <- w_output_mean +
              (1 / (niter * nchain)) * w
          }

          # save jsdm params
          {
            B0_output_chain[,currentIter] <- list_jSDMparams$B0
            B_output_chain[,,currentIter] <- list_jSDMparams$B
            G_output_chain[,,currentIter] <- list_jSDMparams$G
            A_output_chain[,,currentIter] <- list_jSDMparams$A
            C_output_chain[,,currentIter] <- list_jSDMparams$C
            Bs_output_chain[,,currentIter] <- list_jSDMparams$Bs
            Gs_output_chain[,,currentIter] <- list_jSDMparams$Gs
            As_output_chain[,,currentIter] <- list_jSDMparams$As
            Cs_output_chain[,,currentIter] <- list_jSDMparams$Cs
            U_output_chain[,,currentIter] <- list_jSDMparams$U
            L_output_chain[,,currentIter] <- list_jSDMparams$L
            sigmab_output_chain[currentIter] <- list_jSDMparams$sigma_b
            sigmabs_output_chain[currentIter] <- list_jSDMparams$sigma_bs
            sigmah_output_chain[currentIter] <- list_jSDMparams$sigma_h
            idx_ls_output_chain[currentIter] <- list_jSDMparams$idx_ls
            if(model == "continuous") tau_output_chain[,currentIter] <- list_jSDMparams$tau

            varPart_output_chain[,,currentIter] <-
              as.matrix(list_jSDMparams$variancePartitioning)
          }

          # update WIAC
          {

            loglik_jsdm <- computeModelLoglikJSDM_cpp(z, eta, jsdmModel, tau)
            list_waic_jsdm <- update_waic_summary(loglik_jsdm, list_waic_jsdm, currentWAICiter)

            if(model %in% c("occupancy","two_stage")){
              logliks_firststage <- computeModelLoglikFirstStage_cpp(w, z, theta, theta0,
                                                                 list_idx$idx_z_w)
              list_waic_w <- update_waic_summary(logliks_firststage, list_waic_w, currentWAICiter)
            }

            if(model == "two_stage"){
              logliks_secondstage <- computeModelLoglikSecondStage_cpp(y, w, p, q,
                                                                   list_idx$idx_w_k,
                                                                   list_idx$idx_p_k)
              list_waic_y <- update_waic_summary(logliks_secondstage, list_waic_y, currentWAICiter)
            }

            currentWAICiter <- currentWAICiter + 1

          }

        }

      }

    }

    if(model %in% c("occupancy","two_stage")){
      beta_theta_output[,,,chain] <- beta_theta_output_chain
      theta0_output[,,chain] <- theta0_output_chain
      if(!summarisedLatentPresences){
        z_output[,,,chain] <- z_output_chain
        psi_output[,,,chain] <- psi_output_chain
      }
    }

    if(model == "two_stage"){
      p_output[,,,chain] <- p_output_chain
      q_output[,,,chain] <- q_output_chain
    }

    # save jsdm params
    {
      B0_output[,,chain] <- B0_output_chain
      B_output[,,,chain] <- B_output_chain
      G_output[,,,chain] <- G_output_chain
      A_output[,,,chain] <- A_output_chain
      C_output[,,,chain] <- C_output_chain
      Bs_output[,,,chain] <- Bs_output_chain
      Gs_output[,,,chain] <- Gs_output_chain
      As_output[,,,chain] <- As_output_chain
      Cs_output[,,,chain] <- Cs_output_chain
      U_output[,,,chain] <- U_output_chain
      L_output[,,,chain] <- L_output_chain
      sigmab_output[,chain] <- sigmab_output_chain
      sigmabs_output[,chain] <- sigmabs_output_chain
      sigmah_output[,chain] <- sigmah_output_chain
      idx_ls_output[,chain] <- idx_ls_output_chain
      varPart_output[,,,chain] <- varPart_output_chain
      tau_output[,,chain] <- tau_output_chain

    }

  }

  # compute WAIC
  {
    numIters <- niter * nchain
    if(numIters != (currentWAICiter-1)) stop("Current wAIC iter wrong")
    WAIC_jsdm <- compute_waic(list_waic_jsdm, numIters)

    if(model %in% c("occupancy","two_stage")){
      WAIC_w <- compute_waic(list_waic_w, numIters)
    } else {
      WAIC_w <- 0
    }
    if(model %in% "two_stage"){
      WAIC_y <- compute_waic(list_waic_y, numIters)
    } else {
      WAIC_y <- 0
    }

    WAIC <- WAIC_jsdm + WAIC_w + WAIC_y
  }

  # reparametrise factor coefficients
  {
    if(d > 0){

      list_UL_output_reparm <- reparamFactorModel(U_output, L_output)
      U_output <- list_UL_output_reparm$U_output
      L_output <- list_UL_output_reparm$L_output

    }

    if(ncov_psi > 0 & gt > 0){

      list_AC_output_reparm <- reparamFactorModel(A_output, C_output)
      A_output <- list_AC_output_reparm$U_output
      C_output <- list_AC_output_reparm$L_output

    }

  }

  jsdm_results_output <- list(
    "B0_output" = B0_output,
    "B_output" = B_output,
    "G_output" = G_output,
    "A_output" =  A_output,
    "C_output" = C_output,
    "Bs_output" = Bs_output,
    "Gs_output" = Gs_output,
    "As_output" = As_output,
    "Cs_output" = Cs_output,
    "U_output" = U_output,
    "L_output" = L_output,
    "sigmab_output" = sigmab_output,
    "sigmabs_output" = sigmabs_output,
    "sigmah_output" = sigmah_output,
    "idx_ls_output" = idx_ls_output,
    "varPart_output" = varPart_output,
    "tau_output" = tau_output
  )

  results_output <- list(
    "jsdm_output" = jsdm_results_output,
    "beta_theta_output" = beta_theta_output,
    "p_output" = p_output,
    "q_output" = q_output,
    "theta0_output" = theta0_output,
    "WAIC" = WAIC
  )

  if(summarisedLatentPresences){
    results_output$z_output <- z_output_mean
    results_output$psi_output <- psi_output_mean
  } else {
    results_output$z_output <- z_output
    results_output$psi_output <- psi_output
  }
  results_output$w_output <- w_output_mean
  results_output$theta_output <- theta_output_mean

  # print diagnostics
  computeDiagnostics(results_output)

  infos <- list(
    "S" = S,
    "g" = g,
    "M" = M,
    "n" = n,
    "K" = K,
    "ps" = ps,
    "n_factors" = d,
    "list_idx" = list_idx,
    "data_info" = data_info,
    "speciesNames" = speciesNames,
    "primerNames" = primerNames,
    "siteNames" = siteNames,
    "ncov_theta" = ncov_theta,
    "ncov_psi" = ncov_psi,
    "OTU" = OTU,
    "X0_psi" = X0_psi,
    "list_Xs" = list_Xs,
    "list_X_psi_mat" = list_Xpsi_mat,
    "list_Xs_mat" = list_Xs_mat,
    "list_X_theta_mat" = list_X_theta_mat,
    "l_s_grid" = l_s_grid,
    "model" = model,
    "jsdmModel" = jsdmModel
  )

  list(
    "results_output" = results_output,
    "infos" = infos,
    "Tr" = Tr,
    "X_theta" = X_theta,
    "Xs" = Xs,
    "X_psi" = X_psi)

}
