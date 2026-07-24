#' Pad a posterior array to 4 dimensions
#'
#' Internal helper. MCMC output arrays produced by [runOccJSDM()] vary in
#' shape depending on whether a parameter is indexed by nothing (scalar,
#' `[niter, nchain]`), one index (e.g. species-only, `[dim1, niter, nchain]`),
#' or two indices (e.g. covariate x species, `[dim1, dim2, niter, nchain]`).
#' This pads the first two cases up to the common 4-dimensional shape so
#' downstream diagnostics can treat every parameter uniformly.
#'
#' @param param_output A numeric array with 2, 3, or 4 dimensions, with the
#'   last two dimensions always being `niter` and `nchain`.
#' @return A 4-dimensional array `[dim1, dim2, niter, nchain]`.
#' @noRd
as4d <- function(param_output) {
  d <- dim(param_output)
  if (is.null(d)) {
    stop("param_output must be an array with 2, 3, or 4 dimensions")
  }

  nd <- length(d)

  if (nd == 4) {
    return(param_output)
  } else if (nd == 3) {
    return(array(param_output, dim = c(d[1], 1, d[2], d[3])))
  } else if (nd == 2) {
    return(array(param_output, dim = c(1, 1, d[1], d[2])))
  } else {
    stop("param_output must be an array with 2, 3, or 4 dimensions")
  }
}

#' Compute effective sample size for an MCMC output array
#'
#' @param param_output A numeric array with 2, 3, or 4 dimensions
#'   (`[dim1, dim2, niter, nchain]`, or a 2- or 3-dimensional array missing
#'   one or both of the leading index dimensions).
#'
#' @return A `dim1` x `dim2` matrix of effective sample sizes (pooled across
#'   chains).
#' @noRd
computeESSparams <- function(param_output) {
  param_output <- as4d(param_output)

  dim1 <- dim(param_output)[1]
  dim2 <- dim(param_output)[2]
  nchain <- dim(param_output)[4]

  ESS_param <- matrix(NA_real_, dim1, dim2)
  if (dim1 > 0 & dim2 > 0) {
    for (x in 1:dim1) {
      for (y in 1:dim2) {
        chains <- lapply(1:nchain, function(i) coda::mcmc(param_output[x, y, , i]))
        mcmc_list <- coda::mcmc.list(chains)
        ESS_param[x, y] <- coda::effectiveSize(mcmc_list)
      }
    }
  }

  ESS_param
}

computeMinESS <- function(results_output) {

  beta0_psi_output <- results_output$jsdm_output$B0_output
  beta_psi_output <- results_output$jsdm_output$B_output
  beta_theta_output <- results_output$beta_theta_output
  L_output <- results_output$jsdm_output$L_output

  ESS_beta0psi <- rep(NA, dim(beta0_psi_output)[1])
  for (x in 1:length(ESS_beta0psi)) {
    chains <- lapply(1:dim(beta0_psi_output)[3],
                     function(i) coda::mcmc(beta0_psi_output[x, , i]))
    mcmc_list <- coda::mcmc.list(chains)
    ESS_beta0psi[x] <- coda::effectiveSize(mcmc_list)
  }

  ESS_betapsi <- matrix(NA, dim(beta_psi_output)[1],dim(beta_psi_output)[2])
  if(dim(beta_psi_output)[1] > 1){
    for (x in 1:nrow(ESS_betapsi)) {
      for (y in 1:ncol(ESS_betapsi)) {
        chains <- lapply(1:dim(beta_psi_output)[4], function(i) coda::mcmc(beta_psi_output[x, y, , i]))
        mcmc_list <- coda::mcmc.list(chains)
        ESS_betapsi[x,y] <- coda::effectiveSize(mcmc_list)
      }
    }
  }

  if(!is.null(beta_theta_output)){
    ESS_betatheta <- matrix(NA, dim(beta_theta_output)[1],dim(beta_theta_output)[2])
    for (x in 1:nrow(ESS_betatheta)) {
      for (y in 1:ncol(ESS_betatheta)) {
        chains <- lapply(1:dim(beta_theta_output)[4], function(i) coda::mcmc(beta_theta_output[x, y, , i]))
        mcmc_list <- coda::mcmc.list(chains)
        ESS_betatheta[x,y] <- coda::effectiveSize(mcmc_list)
      }
    }
  } else {
    ESS_betatheta <- NULL
  }

  ESS_L <- matrix(NA, dim(L_output)[1],dim(L_output)[2])
  if(dim(L_output)[1] > 0){
    for (x in 1:nrow(ESS_L)) {
      for (y in 1:ncol(ESS_L)) {
        chains <- lapply(1:dim(L_output)[4], function(i) coda::mcmc(L_output[x, y, , i]))
        mcmc_list <- coda::mcmc.list(chains)
        ESS_L[x,y] <- coda::effectiveSize(mcmc_list)
      }
    }
  }


  return(min(ESS_beta0psi,
             ESS_betatheta,
             ESS_betapsi[!is.na(ESS_betapsi)],
             ESS_L[ESS_L != 0]))

}

#' Compute Gelman-Rubin Rhat for an MCMC output array
#'
#' Per-element potential scale reduction factor (Rhat), via
#' [coda::gelman.diag()]. Elements with fewer than two chains, or where any
#' chain never moved (zero variance, which makes [coda::gelman.diag()]
#' error), are returned as `NA`.
#'
#' @param param_output A numeric array with 2, 3, or 4 dimensions
#'   (`[dim1, dim2, niter, nchain]`, or a 2- or 3-dimensional array missing
#'   one or both of the leading index dimensions).
#'
#' @return A `dim1` x `dim2` matrix of Rhat values.
#' @noRd
computeRhat <- function(param_output) {
  param_output <- as4d(param_output)

  dim1 <- dim(param_output)[1]
  dim2 <- dim(param_output)[2]
  nchain <- dim(param_output)[4]

  Rhat <- matrix(NA_real_, dim1, dim2)

  if (nchain < 2) {
    return(Rhat)
  }

  for (x in seq_len(dim1)) {
    for (y in seq_len(dim2)) {
      chains <- lapply(seq_len(nchain), function(i) coda::mcmc(param_output[x, y, , i]))
      variances <- vapply(chains, function(ch) stats::var(as.numeric(ch)), numeric(1))
      if (any(variances == 0)) {
        Rhat[x, y] <- NA_real_
      } else {
        mcmc_list <- coda::mcmc.list(chains)
        Rhat[x, y] <- tryCatch(
          coda::gelman.diag(mcmc_list, autoburnin = FALSE)$psrf[1, 1],
          error = function(e) NA_real_
        )
      }
    }
  }

  Rhat
}

#' Reshape a posterior array into a long-format tibble
#'
#' Internal helper feeding [plotTraceplot()]. Produces one row per draw per
#' chain per parameter element.
#'
#' @param param_output A numeric array with 2, 3, or 4 dimensions.
#' @param param_name A label for the parameter, stored in the `param` column.
#' @param dimnames1 Optional labels for the first index dimension (e.g.
#'   covariate names). Defaults to `1:dim1`.
#' @param dimnames2 Optional labels for the second index dimension (e.g.
#'   species names). Defaults to `1:dim2`.
#'
#' @return A tibble with columns `param`, `label1`, `label2`, `chain`, `iter`,
#'   `value`.
#' @noRd
paramOutputToLong <- function(param_output, param_name = "parameter",
                               dimnames1 = NULL, dimnames2 = NULL) {
  param_output <- as4d(param_output)

  dim1 <- dim(param_output)[1]
  dim2 <- dim(param_output)[2]
  niter <- dim(param_output)[3]
  nchain <- dim(param_output)[4]

  label1 <- if (!is.null(dimnames1)) dimnames1 else as.character(seq_len(dim1))
  label2 <- if (!is.null(dimnames2)) dimnames2 else as.character(seq_len(dim2))

  purrr::map_dfr(seq_len(dim1), function(x) {
    purrr::map_dfr(seq_len(dim2), function(y) {
      purrr::map_dfr(seq_len(nchain), function(ch) {
        tibble::tibble(
          param = param_name,
          label1 = label1[x],
          label2 = label2[y],
          chain = factor(ch),
          iter = seq_len(niter),
          value = param_output[x, y, , ch]
        )
      })
    })
  })
}

#' Summarise a posterior array into a tidy diagnostics table
#'
#' Computes posterior mean, standard deviation, 95% credible interval,
#' Rhat, and effective sample size for every element of a posterior array,
#' pooling draws across chains for the summary statistics.
#'
#' @param param_output A numeric array with 2, 3, or 4 dimensions
#'   (`[dim1, dim2, niter, nchain]`, or a 2- or 3-dimensional array missing
#'   one or both of the leading index dimensions).
#' @param param_name A label for the parameter, stored in the `param` column.
#' @param dimnames1 Optional labels for the first index dimension (e.g.
#'   covariate names). Defaults to `1:dim1`.
#' @param dimnames2 Optional labels for the second index dimension (e.g.
#'   species names). Defaults to `1:dim2`.
#'
#' @return A tibble with columns `param`, `idx1`, `idx2`, `label1`, `label2`,
#'   `mean`, `sd`, `q2.5`, `q97.5`, `rhat`, `ess`.
#' @noRd
summarisePosterior <- function(param_output, param_name = "parameter",
                                dimnames1 = NULL, dimnames2 = NULL) {
  param_output <- as4d(param_output)

  dim1 <- dim(param_output)[1]
  dim2 <- dim(param_output)[2]

  ess <- computeESSparams(param_output)
  rhat <- computeRhat(param_output)

  label1 <- if (!is.null(dimnames1)) dimnames1 else as.character(seq_len(dim1))
  label2 <- if (!is.null(dimnames2)) dimnames2 else as.character(seq_len(dim2))

  purrr::map_dfr(seq_len(dim1), function(x) {
    purrr::map_dfr(seq_len(dim2), function(y) {
      draws <- as.vector(param_output[x, y, , ])
      tibble::tibble(
        param = param_name,
        idx1 = x,
        idx2 = y,
        label1 = label1[x],
        label2 = label2[y],
        mean = mean(draws),
        sd = stats::sd(draws),
        q2.5 = as.numeric(stats::quantile(draws, .025)),
        q97.5 = as.numeric(stats::quantile(draws, .975)),
        rhat = rhat[x, y],
        ess = ess[x, y]
      )
    })
  })
}

#' Plot MCMC trace plots for a posterior array
#'
#' Faceted per-chain trace plots for every element of a posterior array.
#'
#' @param param_output A numeric array with 2, 3, or 4 dimensions
#'   (`[dim1, dim2, niter, nchain]`, or a 2- or 3-dimensional array missing
#'   one or both of the leading index dimensions).
#' @param param_name A label for the parameter, used for the y-axis title.
#' @param dimnames1 Optional labels for the first index dimension (e.g.
#'   covariate names). Defaults to `1:dim1`.
#' @param dimnames2 Optional labels for the second index dimension (e.g.
#'   species names). Defaults to `1:dim2`.
#'
#' @return A `ggplot` object.
#' @export
plotTraceplot <- function(param_output, param_name = "parameter",
                           dimnames1 = NULL, dimnames2 = NULL) {
  df <- paramOutputToLong(param_output, param_name, dimnames1, dimnames2)

  facet_formula <- if (dplyr::n_distinct(df$label2) > 1) {
    ggplot2::vars(label1, label2)
  } else {
    ggplot2::vars(label1)
  }

  ggplot2::ggplot(df, ggplot2::aes(x = iter, y = value, color = chain)) +
    ggplot2::geom_line(alpha = .7) +
    ggplot2::facet_wrap(facet_formula, scales = "free_y") +
    ggplot2::labs(x = "Iteration", y = param_name, color = "Chain")
}

#' Assemble a tidy MCMC convergence diagnostics table for a fitted model
#'
#' Convenience wrapper around [summarisePosterior()] that assembles one
#' tidy table across the main occupancy/detection coefficients of a fitted
#' [runOccJSDM()] model: `beta0_psi`, `beta_psi`, `beta_theta`, `p`, `q`, and
#' `theta0`. Components not present in `fitmodel$results_output` (e.g.
#' `beta_theta`/`p`/`q` for a JSDM-only fit) are silently skipped.
#'
#' @param fitmodel A model object returned by [runOccJSDM()].
#'
#' @return A tibble with columns `param`, `idx1`, `idx2`, `label1`, `label2`,
#'   `mean`, `sd`, `q2.5`, `q97.5`, `rhat`, `ess`, one row per parameter
#'   element.
#' @export
returnConvergenceDiagnostics <- function(fitmodel) {
  results_output <- fitmodel$results_output

  speciesNames <- fitmodel$infos$speciesNames
  psiCovNames <- colnames(fitmodel$X_psi)
  thetaCovNames <- colnames(fitmodel$X_theta)
  primerNames <- if (!is.null(fitmodel$infos$primerNames)) {
    as.character(fitmodel$infos$primerNames)
  } else {
    NULL
  }

  pieces <- list()

  if (!is.null(results_output$jsdm_output$B0_output)) {
    pieces$beta0_psi <- summarisePosterior(
      results_output$jsdm_output$B0_output,
      param_name = "beta0_psi",
      dimnames1 = speciesNames
    )
  }

  if (!is.null(results_output$jsdm_output$B_output)) {
    pieces$beta_psi <- summarisePosterior(
      results_output$jsdm_output$B_output,
      param_name = "beta_psi",
      dimnames1 = psiCovNames,
      dimnames2 = speciesNames
    )
  }

  if (!is.null(results_output$beta_theta_output)) {
    pieces$beta_theta <- summarisePosterior(
      results_output$beta_theta_output,
      param_name = "beta_theta",
      dimnames1 = thetaCovNames,
      dimnames2 = speciesNames
    )
  }

  if (!is.null(results_output$p_output)) {
    pieces$p <- summarisePosterior(
      results_output$p_output,
      param_name = "p",
      dimnames1 = primerNames,
      dimnames2 = speciesNames
    )
  }

  if (!is.null(results_output$q_output)) {
    pieces$q <- summarisePosterior(
      results_output$q_output,
      param_name = "q",
      dimnames1 = primerNames,
      dimnames2 = speciesNames
    )
  }

  if (!is.null(results_output$theta0_output)) {
    pieces$theta0 <- summarisePosterior(
      results_output$theta0_output,
      param_name = "theta0",
      dimnames1 = speciesNames
    )
  }

  dplyr::bind_rows(pieces)
}

#' Compute diagnostics
#'
#' @details
#' Prints MCMC convergence diagnostics (R-hat and effective sample size) to
#' the console for each block of parameters in a fitted model's
#' `results_output` (e.g. `beta0_psi`, `beta_psi`, `beta_theta`, `p`, `q`,
#' `theta0`), skipping latent-state components (`z_output`, `w_output`,
#' `theta_output`, `idx_ls_output`, `varPart_output`) and any block that is
#' entirely fixed/constant. Emits a warning for any parameter block whose
#' max R-hat exceeds 1.1 or minimum ESS falls below 50. Called automatically
#' at the end of \code{\link{runOccJSDM}}.
#'
#' @param results_output The `results_output` element of a fitted model, as
#'   returned by \code{\link{runOccJSDM}} (`fitModel$results_output`).
#'
#' @return `NULL`, invisibly. Called for its side effect of printing
#'   diagnostics to the console.
#' @export
computeDiagnostics <- function(results_output){

  nchain <- dim(results_output$jsdm_output$B0_output)[3]
  niter <- dim(results_output$jsdm_output$B0_output)[2]

  results_output_all <- c(
    results_output[names(results_output) != "jsdm_output"],
    results_output$jsdm_output)

  for (param_name in names(results_output_all)) {

    if(!(param_name %in% c("z_output","w_output",
                           "theta_output","idx_ls_output","varPart_output"))){

      arr <- results_output_all[[param_name]]
      dims <- dim(arr)

      if (is.null(dims) || length(dims) < 2 || any(dims== 0 ) || any(is.na(arr)) ) {
        next
      }

      ndim <- length(dims)

      if (ndim == 2) {
        n_params <- 1
      } else {
        n_params <- prod(dims[1:(ndim - 2)])
      }

      # Reshape uniformly to 3D: (n_params, niter, nchain)
      arr_3d <- array(arr, dim = c(n_params, niter, nchain))

      arr_flat <- matrix(arr_3d, nrow = n_params, ncol = niter * nchain)
      param_vars <- apply(arr_flat, 1, var, na.rm = TRUE)

      # Keep only parameters with variance greater than a tiny threshold
      is_variable <- !is.na(param_vars) & param_vars > 1e-12
      valid_params <- sum(is_variable)

      if (valid_params == 0) {
        message(sprintf("\nParameter Set: %s", param_name))
        message("--------------------------------------------------")
        message("  [Skipped: All parameters in this set are fixed/constant]")
        next
      }

      # Subset to remove the fixed parameters
      arr_3d <- arr_3d[is_variable, , , drop = FALSE]

      # Convert to coda mcmc.list
      chain_list <- lapply(1:nchain, function(c) {
        # Use matrix() to prevent R from dropping dimensions on single parameters
        chain_slice <- matrix(arr_3d[, , c], nrow = valid_params, ncol = niter)
        coda::mcmc(t(chain_slice)) # Transpose so rows = iterations, cols = parameters
      })

      mcmc_obj <- coda::mcmc.list(chain_list)

      # Diagnostics
      ess <- coda::effectiveSize(mcmc_obj)

      if (nchain > 1) {
        rhat <- tryCatch({
          diag_res <- coda::gelman.diag(mcmc_obj, autoburnin = FALSE, multivariate = FALSE)
          diag_res$psrf[, 1]
        }, error = function(e) rep(NA, n_params))
      } else {
        rhat <- rep(NA, n_params)
      }

      # Print Summary
      param_name_to_print <- gsub("_output","",param_name)
      message(sprintf("\nParameter Set: %s (Flattened N = %d)", param_name_to_print, n_params))
      message("--------------------------------------------------")

      # Use suppressWarnings for min/max to prevent warnings if all values are NA (e.g., fixed params)
      max_rhat <- suppressWarnings(max(rhat, na.rm = TRUE))
      mean_rhat <- mean(rhat, na.rm = TRUE)
      min_ess <- suppressWarnings(min(ess, na.rm = TRUE))
      mean_ess <- mean(ess, na.rm = TRUE)

      if (nchain > 1) {
        max_rhat <- suppressWarnings(max(rhat, na.rm = TRUE))
        mean_rhat <- mean(rhat, na.rm = TRUE)
        message(sprintf("  R-hat | Max:  %7.3f | Mean: %7.3f", max_rhat, mean_rhat))

        if (is.finite(max_rhat) && max_rhat > 1.1) {
          warning(sprintf("Convergence issue in '%s'. Max R-hat > 1.1", param_name), call. = FALSE)
        }
      } else {
        message("  R-hat | [Requires multiple chains to compute]")
      }

      message(sprintf("  ESS   | Min:  %7.1f | Mean: %7.1f", min_ess, mean_ess))

      if (is.finite(min_ess) && min_ess < 50) {
        warning(sprintf("Low ESS in '%s'. Min ESS < 50", param_name_to_print), call. = FALSE)
      }

    }

  }

  message("\n==================================================")

  invisible(NULL)

}


