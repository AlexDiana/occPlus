#' thinOutput
#'
#' Thins the MCMC output.
#'
#' @details
#' Returns the same fitModel object, but with the MCMC iterations in
#' `results_output` thinned, keeping only every `thin`-th iteration.
#'
#' @param fitModel Output from the function runOccJSDM
#' @param thin Thinning interval; every `thin`-th iteration is kept (default 5)
#'
#' @return A fitModel object with the same structure as the input, but with
#' the MCMC output in `results_output` thinned
#'
#' @import dplyr
#' @import ggplot2
#'
thinOutput <- function(fitModel, thin = 5){

  niter <- dim(fitModel$results_output$jsdm_output$B0_output)[2]
  idx_thinned <- seq(1, niter, by = thin)

  results_output <- fitModel$results_output

  results_output_thinned <- lapply(1:length(results_output), function(i){

    nameElem <- names(results_output)[i]

    x <- results_output[[i]]

    if(nameElem == "z_output"){

      if(length(dim(x))== 4){

        x[,,idx_thinned,,drop=F]

      } else if(length(dim(x))== 2){

        x

      } else {

        print("Dimension not recognised")

      }

    } else {

      if(length(dim(x)) == 4){

        x[,,idx_thinned,,drop=F]

      } else if (length(dim(x)) == 3) {

        x[,idx_thinned,,drop=F]

      } else if (length(dim(x)) == 2) {

        x[idx_thinned,,drop=F]

      } else {

        print("Dimension not recognised")

      }

    }

  })

  names(results_output_thinned) <- names(results_output)

  fitModel$results_output <- results_output_thinned

  fitModel

}

#' logit
#'
#' Logit transformation.
#'
#' @param x A numeric value or vector with values in (0, 1)
#'
#' @return The logit of `x`
#'
#' @noRd
logit <- function(x){
  log(x / (1 - x))
}

#' logistic
#'
#' Logistic (inverse logit) transformation.
#'
#' @param x A numeric value or vector
#'
#' @return The logistic transform of `x`, mapped to (0, 1)
#'
#' @noRd
logistic <- function(x) 1 / (1 + exp(-x))

# TRAITS ---------

#' returnTraitsCoeff
#'
#' Traits covariate coefficients.
#'
#' @details
#' Returns the traits covariates coefficients posterior sample
#'
#' @param fitModel Output from the function runOccJSDM
#'
#' @return An array of posterior samples of size
#' (iterations x occupancy covariates x traits)
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
returnTraitsCoeff <- function(fitModel){

  # S <- fitModel$infos$S
  # g <-
  # ncov_psi <- fitModel$infos$ncov_psi
  # # speciesNames <- fitModel$infos$speciesNames
  occCovNames <- colnames(fitModel$X_psi)
  traitNames <- colnames(fitModel$Tr)

  g <- fitModel$infos$g

  if(g > 0){
    traitsCoeffOutput <- fitModel$results_output$jsdm_output$G_output
    traitsCoeffOutput <- apply(traitsCoeffOutput, c(1,2), c)

    niter <- dim(traitsCoeffOutput)[1]
    dimnames(traitsCoeffOutput)[[2]] <- traitNames
    dimnames(traitsCoeffOutput)[[3]] <- occCovNames

    traitsCoeffOutput <- aperm(traitsCoeffOutput, c(1,3,2))

    return(traitsCoeffOutput)
  } else {

    stop("No Traits Present")
  }

}

#' plotTraitsCoefficients
#'
#' Traits covariate coefficients.
#'
#' @details
#' Plots the 95% credible interval of the occupancy covariates coefficients
#'
#' @param fitModel Output from the function runOccJSDM
#' @param covName Name of the covariate to be plotted (same name as in data$info)
#' @param idx_traits Indexes of the traits to be plotted (leave out to plot all the traits).
#'
#' @return A ggplot object
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotTraitsCoefficients <- function(fitModel,
                                   covName = NULL,
                                   idx_traits = NULL
){

  traits_output <- returnTraitsCoeff(fitModel)

  if(is.null(covName)){
    stop("No name provided")
  }

  plotCoefficient(traits_output, covName, idx_traits) + xlab("Traits")

}

# OCCUPANCY COVARIATES --------

#' returnOccupancyCovariates
#'
#' Occupancy covariate coefficients.
#'
#' @details
#' Returns the occupancy covariates coefficients posterior sample
#'
#' @param fitModel Output from the function runOccJSDM
#'
#' @return An array of posterior samples of size
#' (iterations x occupancy covariates x species)
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
returnOccupancyCovariates <- function(fitModel){

  occCovNames <- colnames(fitModel$X_psi)
  speciesNames <- fitModel$infos$speciesNames

  occCoeffOutput <- fitModel$results_output$jsdm_output$B_output
  occCoeffOutput <- apply(occCoeffOutput, c(1,2), c)

  niter <- dim(occCoeffOutput)[1]
  dimnames(occCoeffOutput)[[2]] <- occCovNames
  dimnames(occCoeffOutput)[[3]] <- speciesNames

  occCoeffOutput

}

#' plotOccupancyCovariates
#'
#' Occupancy covariate coefficients.
#'
#' @details
#' Plots the 95% credible interval of the occupancy covariates coefficients
#'
#' @param fitModel Output from the function runOccJSDM
#' @param covName Name of the covariate to be plotted (same name as in data$info)
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotOccupancyCovariates <- function(fitModel,
                                    covName = NULL,
                                    idx_species = NULL
){

  occCov_output <- returnOccupancyCovariates(fitModel)

  if(is.null(covName)){
    stop("No name provided")
  }

  plotCoefficient(occCov_output, covName, idx_species) + xlab("Species")

}

#' returnOccupancyGradient
#'
#' Predicted occupancy probability across a covariate gradient.
#'
#' @details
#' Returns, for a single occupancy covariate, the posterior credible
#' interval of predicted occupancy probability across a grid of values
#' of that covariate, for each species. All other occupancy covariates
#' are held fixed at their median (observed) value. Predictions are
#' computed from row-matched posterior draws -- i.e. the same MCMC
#' iteration is used for every grid point of a given species -- so the
#' resulting interval reflects genuine posterior uncertainty in the
#' species' intercept and coefficients rather than treating each grid
#' point as an independent estimate.
#'
#' This complements \code{\link{plotOccupancyCovariates}}, which shows
#' the credible interval of the covariate's raw coefficient (on the
#' logit scale, relative to zero) -- useful for assessing whether an
#' effect is credibly nonzero and in which direction. Because the logit
#' link is nonlinear, a given coefficient can correspond to a large or
#' small change in occupancy probability depending on where the species'
#' baseline occupancy sits; \code{returnOccupancyGradient} (and
#' \code{\link{plotOccupancyGradient}}) show that effect directly, on
#' the natural probability scale.
#'
#' @param fitModel Output from the function runOccJSDM
#' @param covName Name of the covariate to vary (same name as in data$info)
#' @param idx_species Indexes of the species to include (leave out for all species)
#' @param n_grid Number of grid points spanning the covariate's range (default 40)
#' @param confidence Confidence level of the credible interval, default .95
#' @param quantile_range Quantiles of the observed covariate values used to
#' set the lower/upper bounds of the grid, default c(.02, .98) (avoids
#' extrapolating too far beyond the bulk of the observed data)
#'
#' @return A tibble with columns \code{species}, \code{x} (covariate value),
#' and the lower/median/upper quantiles of predicted occupancy probability
#'
#' @examples
#' \dontrun{
#' returnOccupancyGradient(fitModel, covName = "X_psi.1")
#' }
#'
#' @export
#' @import dplyr
#' @importFrom purrr map_dfr
#' @importFrom tibble tibble
#'
returnOccupancyGradient <- function(fitModel,
                                     covName = NULL,
                                     idx_species = NULL,
                                     n_grid = 40,
                                     confidence = .95,
                                     quantile_range = c(.02, .98)){

  if(is.null(covName)){
    stop("No name provided")
  }

  occCovNames <- colnames(fitModel$X_psi)
  speciesNames <- fitModel$infos$speciesNames
  S <- fitModel$infos$S

  if(!(covName %in% occCovNames)){
    stop("Covariate name not found. If you are using a categorical covariates,
         the name might have changed to code the level. Use
         colnames(fitModel$X_psi) to find the new names")
  }

  if(is.null(idx_species)){
    idx_species <- seq_len(S)
  }

  otherCovNames <- setdiff(occCovNames, covName)
  otherCovMedian <- apply(fitModel$X_psi, 2, median)

  B0_output <- fitModel$results_output$jsdm_output$B0_output
  B_output <- fitModel$results_output$jsdm_output$B_output

  beta0_draws <- apply(B0_output, 1, c)
  dimnames(beta0_draws) <- list(NULL, speciesNames)

  betaAll_draws <- apply(B_output, c(1,2), c)
  niter <- dim(betaAll_draws)[1]
  dimnames(betaAll_draws)[[2]] <- occCovNames
  dimnames(betaAll_draws)[[3]] <- speciesNames

  grid_vals <- seq(
    quantile(fitModel$X_psi[, covName], quantile_range[1]),
    quantile(fitModel$X_psi[, covName], quantile_range[2]),
    length.out = n_grid
  )

  conflevels <- c((1 - confidence) / 2, .5, (1 + confidence) / 2)

  purrr::map_dfr(speciesNames[idx_species], function(sp){

    otherContribution <- if(length(otherCovNames) > 0){
      B_other <- matrix(betaAll_draws[, otherCovNames, sp],
                        nrow = niter, ncol = length(otherCovNames))
      as.vector(B_other %*% otherCovMedian[otherCovNames])
    } else {
      0
    }

    # Same posterior draw reused across the whole grid for a given
    # species, so uncertainty is not double-counted across grid points
    eta <- outer(betaAll_draws[, covName, sp], grid_vals) +
      (beta0_draws[, sp] + otherContribution)

    psi <- logistic(eta)
    q <- apply(psi, 2, quantile, probs = conflevels)

    tibble::tibble(
      species = sp,
      x = grid_vals,
      low = q[1, ],
      med = q[2, ],
      high = q[3, ]
    )
  })
}


#' plotOccupancyGradient
#'
#' Plot predicted occupancy probability across a covariate gradient.
#'
#' @details
#' Plots, for a single occupancy covariate, a smooth curve of predicted
#' occupancy probability (with a credible band) across the observed
#' range of that covariate, faceted by species, with a rug of the
#' observed covariate values. See \code{\link{returnOccupancyGradient}}
#' for details on how the curve and interval are computed, and how this
#' plot complements \code{\link{plotOccupancyCovariates}}.
#'
#' @param fitModel Output from the function runOccJSDM
#' @param covName Name of the covariate to vary (same name as in data$info)
#' @param idx_species Indexes of the species to include (leave out for all species)
#' @param n_grid Number of grid points spanning the covariate's range (default 40)
#' @param confidence Confidence level of the credible interval, default .95
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotOccupancyGradient(fitModel, covName = "X_psi.1")
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotOccupancyGradient <- function(fitModel,
                                   covName = NULL,
                                   idx_species = NULL,
                                   n_grid = 40,
                                   confidence = .95){

  speciesNames <- fitModel$infos$speciesNames

  if(is.null(idx_species)){
    idx_species <- seq_len(fitModel$infos$S)
  }

  gradient_df <- returnOccupancyGradient(fitModel, covName, idx_species,
                                          n_grid, confidence)

  rug_df <- data.frame(x = fitModel$X_psi[, covName])

  gradient_df %>%
    mutate(species = factor(species, levels = speciesNames[idx_species])) %>%
    ggplot(aes(x = x, y = med)) +
    geom_ribbon(aes(ymin = low, ymax = high), fill = "steelblue", alpha = 0.25) +
    geom_line(color = "steelblue") +
    geom_rug(data = rug_df, aes(x = x, y = NULL), sides = "b",
             alpha = 0.3, inherit.aes = FALSE) +
    facet_wrap(~ species) +
    ylim(c(0,1)) +
    labs(x = covName, y = "Occupancy probability",
         title = paste0("Predicted occupancy vs ", covName)) +
    theme_bw()

}

#' plotSpeciesResponseCurve
#'
#' Plot predicted species response curve
#'
#' @details
#'
#' @param fitModel Output from the function runOccJSDM
#' @param covName Name of the covariate to vary (same name as in data$info)
#' @param idx_species Indexes of the species to include (leave out for all species)
#' @param n_grid Number of grid points spanning the covariate's range (default 40)
#' @param confidence Confidence level of the credible interval, default .95
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotSpeciesResponseCurve(fitModel, covName = "X_psi.1")
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotSpeciesResponseCurve <- function(species_name,
                               target_cov,
                               beta_mcmc_j,
                               list_matrix,
                               raw_df,
                               n_points = 200) {

  # 1. Create a grid for the target covariate across its original range
  cov_seq <- seq(min(raw_df[[target_cov]], na.rm = TRUE),
                 max(raw_df[[target_cov]], na.rm = TRUE),
                 length.out = n_points)

  # 2. Build a synthetic grid dataframe holding other covariates constant at their mean/mode
  grid_df <- data.frame(matrix(ncol = length(list_matrix$names_df), nrow = n_points))
  colnames(grid_df) <- list_matrix$names_df

  for (col_name in list_matrix$names_df) {
    if (col_name == target_cov) {
      grid_df[[col_name]] <- cov_seq
    } else if (list_matrix$is_numeric[[col_name]]) {
      grid_df[[col_name]] <- list_matrix$mean_df[[col_name]] # Hold continuous at original mean
    } else {
      grid_df[[col_name]] <- list_matrix$cat_levels[[col_name]][1] # First factor level
    }
  }

  # 3. Transform synthetic grid into X_pred using your saved function
  # Note: set remove_intercept = FALSE if your MCMC beta vector includes the intercept
  X_pred <- transform_new_covariates(grid_df, list_matrix, remove_intercept = FALSE)

  # 4. Multiply X_pred by MCMC posterior draws for species j
  # beta_mcmc_j is a matrix of size (N_draws x P_cols)
  # posterior_fit will be (n_points x N_draws)
  posterior_fit <- X_pred %*% t(beta_mcmc_j)

  # If you are fitting Presence-Absence (Probit link), apply normal CDF:
  # posterior_fit <- pnorm(posterior_fit)

  # 5. Summarize mean response and 95% Credible Intervals across MCMC draws
  plot_data <- data.frame(
    x     = cov_seq,
    mean  = apply(posterior_fit, 1, mean),
    lower = apply(posterior_fit, 1, quantile, probs = 0.025),
    upper = apply(posterior_fit, 1, quantile, probs = 0.975)
  )

  # 6. Plot using ggplot2
  p <- ggplot(plot_data, aes(x = x, y = mean)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "skyblue", alpha = 0.4) +
    geom_line(color = "blue", linewidth = 1) +
    labs(
      title = paste("Response of", species_name, "to", target_cov),
      x = target_cov,
      y = "Predicted Response (Linear Predictor / Probability)"
    ) +
    theme_minimal()

  return(p)
}


#' returnCovariateEffect
#'
#' Return predicted species response curve
#'
#' @details
#'
#'
#' @param fitModel Output from the function runOccJSDM
#' @param covName Character vector. Name(s) of the covariate(s) to evaluate
#' @param idx_species Indexes of the species to include (leave out for all species)
#' @param confidence Numeric scalar between 0 and 1. The width of the Bayesian
#'   credible interval to compute (default is \code{0.95} for a 95\% interval).
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
returnCovariateEffect <- function(fitModel,
                                  covName,
                                  idx_species = NULL,
                                  confidence = .95){

  B_output <- fitModel$results_output$jsdm_output$B_output
  B0_output <- fitModel$results_output$jsdm_output$B0_output
  speciesNames <- fitModel$infos$speciesNames

  if (is.null(idx_species)) {
    idx_species <- seq_along(speciesNames)
  }

  sp_name <- speciesNames[idx_species]

  B_output_vec <- apply(B_output, c(1,2), c)
  B0_output_vec <- apply(B0_output, 1, c)

  X_psi <- fitModel$X_psi
  X0_psi <- fitModel$infos$X0_psi

  link_model <- ifelse(fitModel$infos$jsdmModel == "continuous","identity","logit")

  list_matrix <- fitModel$infos$list_X_psi_mat

  returnCovariateEffect_base(
    covName,
    idx_species,
    sp_name,
    B0_output_vec,
    B_output_vec,
    list_matrix,
    speciesNames,
    X0 = X0_psi,
    X = X_psi,
    link = link_model
  )

}

#' plotCovariateEffect
#'
#' Plot predicted species response curve
#'
#' @details
#'
#'
#' @param fitModel Output from the function runOccJSDM
#' @param covNames Character vector. Name(s) of the covariate(s) to evaluate
#' @param idx_species Indexes of the species to include (leave out for all species)
#' @param confidence Numeric scalar between 0 and 1. The width of the Bayesian
#'   credible interval to compute (default is \code{0.95} for a 95\% interval).
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotCovariateEffect <- function(fitModel,
                                covNames,
                                idx_species = NULL,
                                confidence = .95){

  B_output <- fitModel$results_output$jsdm_output$B_output
  B0_output <- fitModel$results_output$jsdm_output$B0_output
  speciesNames <- fitModel$infos$speciesNames

  if (is.null(idx_species)) {
    idx_species <- seq_along(speciesNames)
  }

  sp_name <- speciesNames[idx_species]

  B_output_vec <- apply(B_output, c(1,2), c)
  B0_output_vec <- apply(B0_output, 1, c)

  X_psi <- fitModel$X_psi
  X0_psi <- fitModel$infos$X0_psi

  link_model <- ifelse(fitModel$infos$jsdmModel == "continuous","identity","logit")

  list_matrix <- fitModel$infos$list_X_psi_mat

  plot_list <- plotCovariateEffect_base(
    idx_species,
    covNames,
    B0_output_vec,
    B_output_vec,
    list_matrix,
    speciesNames,
    X0 = X0_psi,
    X = X_psi,
    link = link_model
  )

  plot_list
}

#' returnBaselineOccupancyRates
#'
#' Baseline occupancy rate for each species.
#'
#' @details
#' Returns the posterior samples of the baseline occupancy probability
#' (intercept-only, on the probability scale) for each species.
#'
#' @param fitModel Output from the function runOccJSDM
#'
#' @return A matrix of size (species x iterations) with the posterior
#' samples of the baseline occupancy probability for each species
#'
#' @examples
#' \dontrun{
#' returnBaselineOccupancyRates(fitModel)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
returnOccupancyRates <- function(fitModel){

  speciesNames <- fitModel$infos$speciesNames

  logisticB0_output <- logistic(fitModel$results_output$jsdm_output$B0_output)
  logisticB0_output <- apply(logisticB0_output, 1, c)

  dimnames(logisticB0_output)[[2]] <- speciesNames

  logisticB0_output
}

#' plotOccupancyRates
#'
#' Baseline occupancy rate for each species.
#'
#' @details
#' Plots the 95% credible interval of the baseline occupancy rates
#'
#' @param fitModel Output from the function runOccJSDM
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#' @param confidence Confidence level of the estimate,  default to .95
#'
#' @return The credible interval plot
#'
#' @examples
#' \dontrun{
#' plotOccupancyRates(fitModel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotOccupancyRates <- function(fitModel,
                               idx_species = NULL,
                               confidence = .95){

  confInt <- c((1 - confidence) / 2, (1 + confidence) / 2)

  psi0_output <- returnOccupancyRates(fitModel)
  S <- fitModel$infos$S
  speciesNames <- fitModel$infos$speciesNames

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  data_plot <- apply(psi0_output, 2, function(x) {
    quantile(x, probs = confInt)
  }) %>%
    t %>%
    as.data.frame

  names(data_plot) <- c("Min", "Max")

  plotSpeciesRates(data_plot, idx_species, speciesNames) +
    ggtitle("Baseline Occupancy rates")

}


# COLLECTION COVARIATES --------

#' returnCollectionCovariates
#'
#' Collection covariate coefficients.
#'
#' @details
#' Returns the collection covariates coefficients posterior sample
#'
#' @param fitModel Output from the function runOccJSDM
#'
#' @return An array of posterior samples of size
#' (iterations x collection covariates x species)
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
returnCollectionCovariates <- function(fitModel){

  matrix_of_draws <- fitModel$results_output

  S <- fitModel$infos$S
  ncov_theta <- fitModel$infos$ncov_theta
  speciesNames <- fitModel$infos$speciesNames
  collcovNames <- colnames(fitModel$X_theta)
  # idxcov <- which(occCovNames == covName)

  # samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]
  # samples_subset <- samples_subset[,idxcov + 0:(S - 1)*ncov_psi]

  beta_theta_output <- matrix_of_draws$beta_theta_output

  beta_theta_output <- apply(beta_theta_output, c(1,2), c)
  # beta_theta_output <- aperm(beta_theta_output, c(2,3,1))

  dimnames(beta_theta_output)[[2]] <- collcovNames
  dimnames(beta_theta_output)[[3]] <- speciesNames

  beta_theta_output

}

#' plotCollectionCovariates
#'
#' Collection covariate coefficients.
#'
#' @details
#' Plots the 95% credible interval of the collection covariates coefficients
#'
#' @param fitModel Output from the function runOccJSDM
#' @param covName Name of the covariate to be plotted (same name as in data$info)
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotCollectionCovariates <- function(fitModel,
                                     covName = NULL,
                                     idx_species = NULL
){

  collCov_output <- returnCollectionCovariates(fitModel)

  if(is.null(covName)){
    stop("No name provided")
  }

  plotCoefficient(collCov_output, covName, idx_species) + xlab("Species")

}

#' plotSpeciesRates
#'
#' Plot the 95% credible interval of a per-species rate for a subset of species.
#'
#' @details
#' Internal helper that plots error bars for a per-species rate. `data_plot`
#' must have one row per species, in species order, with numeric columns `Min`
#' and `Max` giving the interval bounds. Subsetting and labelling are done here
#' rather than by the caller, so that the ordering is computed on exactly the
#' rows that get plotted.
#'
#' @param data_plot A data frame with one row per species and columns `Min`,
#'   `Max`
#' @param idx_species Indexes of the species to display
#' @param speciesNames Character vector of all species names, in species order
#'   (i.e. not pre-subset by `idx_species`)
#'
#' @return A ggplot object
#'
#' @noRd
plotSpeciesRates <- function(data_plot,
                             idx_species,
                             speciesNames){

  data_plot <- data_plot[idx_species, , drop = FALSE]
  data_plot$Species <- speciesNames[idx_species]

  # Order on the subset that is actually plotted. Ordering the full set and
  # then indexing the subset (or the reverse) mismatches labels to bars -- the
  # defect in TODO.md group B item 1, which the other three rate plots still
  # have.
  speciesNameOrdered <- data_plot$Species[order(data_plot$Min)]

  plot_speciesrates <- data_plot %>%
    ggplot(aes(x =  factor(Species, levels = speciesNameOrdered),
               ymin = Min,
               ymax = Max)) + geom_errorbar() +
    xlab("Species") +
    theme_bw() +
    ylim(c(0,1)) +
    ylab("") +
    theme(
      axis.text = element_text(angle = 0,
                               size = 8),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = .5,
                                size = 15)
    ) + coord_flip()

  plot_speciesrates

}



#' plotCollectionRates
#'
#' Baseline collection rate for each species.
#'
#' @details
#' Plots the 95% credible interval of the baseline collection rates
#'
#' @param fitModel Output from the function runOccJSDM
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return The credible interval plot
#'
#' @examples
#' \dontrun{
#' plotCollectionRates(fitModel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotCollectionRates <- function(fitModel,
                                idx_species = NULL){

  S <- fitModel$infos$S
  speciesNames <- fitModel$infos$speciesNames

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  beta_theta_output <- fitModel$results_output$beta_theta_output
  if(is.null(beta_theta_output)){
    stop("plotCollectionRates() needs collection-stage parameters, which model '",
         fitModel$infos$model, "' does not have. It applies to 'two_stage' and ",
         "'occupancy' fits.")
  }

  # Row 1 of beta_theta is the collection-rate intercept, so this is the
  # baseline rate with covariates at zero. drop = FALSE keeps the array 4-D so
  # the apply() below behaves for S == 1.
  samples_subset <- beta_theta_output[1, , , , drop = FALSE]
  samples_subset <- apply(samples_subset, 2, c)

  data_plot <- apply(samples_subset, 2, function(x) {
    quantile(logistic(x), probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame

  # plotSpeciesRates() works in Min/Max; quantile() names these "2.5%"/"97.5%".
  names(data_plot) <- c("Min", "Max")

  plotSpeciesRates(data_plot, idx_species, speciesNames) +
    ggtitle("Collection rates")


  # %>%
  #   mutate(Species = speciesNames) %>%
  # mutate(speciesOrder = order(`2.5%`)) %>%
  # filter(Species %in% speciesNames[idx_species])
  #
  # orderSpecies <- order(data_plot$`2.5%`)
  #
  # plot_collectionrates <- data_plot %>%
  #   ggplot(aes(x =  factor(Species, level = speciesNames[orderSpecies]),
  #              ymin = `2.5%`,
  #              ymax = `97.5%`)) + geom_errorbar() +
  #   xlab("Species") +
  #   ggtitle("Collection rates") +
  #   theme_bw() +
  #   ylim(c(0,1)) +
  #   theme(
  #     axis.text = element_text(angle = 0,
  #                              size = 8),
  #     plot.title = element_text(hjust = .5,
  #                               size = 15)
  #   ) + coord_flip()
  #
  # plot_collectionrates

}

# SECOND STAGE RATES ----

#' plotFPTPStage2Rates
#'
#' True and false positives rates at the lab stage for each species.
#'
#' @details
#' Plots the 95% credible interval of the true and false positives at the lab stage
#'
#' @param fitModel Output from the function runOccJSDM
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#' @param primerName Name of the primer to plot (defaults to the first primer in `fitModel$infos$primerNames`)
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotFPTPStage2Rates(fitModel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import stringr
#' @import ggplot2
#'
plotFPTPStage2Rates <- function(fitModel,
                                idx_species = NULL,
                                primerName = NULL){

  S <- fitModel$infos$S
  speciesNames <- fitModel$infos$speciesNames
  primerNames <- fitModel$infos$primerNames

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  if(is.null(primerName)){
    primerName <- primerNames[1]
  }

  idx_primer <- match(as.character(primerName), as.character(primerNames))

  if(is.na(idx_primer)){
    stop(paste0("primerName '", primerName, "' not found. Available primers: ",
                paste(primerNames, collapse = ", ")))
  }

  # Subset to the requested primer before summarising, so that the credible
  # intervals describe that primer alone. p_output/q_output are 4-D arrays
  # indexed [primer, species, iteration, chain]; dropping only the primer
  # margin leaves the species margin (2) in place for the apply() below.
  p_output <- fitModel$results_output$p_output[idx_primer, , , , drop = FALSE]
  q_output <- fitModel$results_output$q_output[idx_primer, , , , drop = FALSE]

  data_plot_p <- apply(p_output, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame %>%
    rename(p1 = `2.5%`,
           p2 = `97.5%`)

  data_plot_q <- apply(q_output, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame %>%
    rename(q1 = `2.5%`,
           q2 = `97.5%`)

  data_plot <- cbind(data_plot_p, data_plot_q) %>%
    mutate(Species = speciesNames) %>%
    filter(Species %in% speciesNames[idx_species])

  # Order and label using the filtered subset that is actually plotted --
  # computing order() on the filtered rows and then indexing the unfiltered
  # speciesNames (as before) mismatches labels to bars for any idx_species
  # other than a prefix 1:k (TODO.md group C item 1).
  speciesNameOrdered <- data_plot$Species[order(data_plot$p1)]

  detectionRates <- data_plot %>%
    ggplot()  +
    geom_errorbar(aes(x = factor(Species, level = speciesNameOrdered),
                      ymin = p1,
                      ymax = p2,
                      color = "TP rate"), position = position_dodge(width = .15), # Use the SAME width as geom_col
                  width = .5) +
    geom_errorbar(aes(x = factor(Species, level = speciesNameOrdered),
                      ymin = q1,
                      ymax = q2,
                      color = "FP rate"), position = position_dodge(width = .15), # Use the SAME width as geom_col
                  width = .5) +
    xlab("Species") +
    ggtitle(paste0("Detection rates (primer ", primerName, ")")) +
    theme_bw() +
    # ylim(c(0,1)) +
    ylab("Detection probability") +
    scale_color_manual(
      name = "Colour",
      values = c("TP rate" = "blue", "FP rate" = "red")
    ) +
    theme(
      axis.text = element_text(angle = 0,
                               size = 8),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = .5,
                                size = 15)
    ) + coord_flip()

  detectionRates

}

#' plotDetectionRates
#'
#' True positives rates at the lab stage for each species.
#'
#' @details
#' Plots the 95% credible interval of the true positives at the lab stage
#'
#' @param fitModel Output from the function runOccJSDM
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotDetectionRates(fitModel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import stringr
#' @import ggplot2
#'
plotDetectionRates <- function(fitModel,
                               idx_species = NULL){

  if(fitModel$infos$model == "two_stage"){


    S <- fitModel$infos$S
    speciesNames <- fitModel$infos$speciesNames
    primerNames <- fitModel$infos$primerNames

    if(is.null(idx_species)){
      idx_species <- 1:S
    }

    p_output <- fitModel$results_output$p_output

    data_plot_array <- apply(p_output, c(1,2), function(x) {
      quantile(x, probs = c(0.025, 0.975))
    })

    dimnames(data_plot_array)[[2]] <- primerNames

    data_plot <- data.frame(
      lower = as.vector(data_plot_array[1,,]),
      upper = as.vector(data_plot_array[2,,]),
      Primer = factor(rep(primerNames, times = dim(data_plot_array)[3]),
                      levels = primerNames),
      Species = rep(speciesNames, each = dim(data_plot_array)[2])
    ) %>%
      filter(Species %in% speciesNames[idx_species]) %>%
      group_by(Species) %>%
      mutate(mean_lower = mean(lower)) %>%
      ungroup() %>%
      mutate(Species = reorder(factor(Species), mean_lower))

    plotDetectionRates <- ggplot() +
      geom_errorbar(data = data_plot, aes(x = Species,
                                          ymin = lower, ymax = upper, color = Primer),
                    position = position_dodge(width = 0.6)) +
      labs(
        title = "Detection rates",
        x = "Species",
        y = "p",
        color = "Primer"
      ) +
      # position_dodge() orders groups by factor level left-to-right, and
      # coord_flip() inverts that into top-to-bottom; reverse the legend so
      # it matches the on-screen vertical order of the errorbars
      guides(color = guide_legend(reverse = TRUE)) +
      theme_bw() + coord_flip()

    plotDetectionRates
  } else {
    stop("Cannot generate detection rates if the model is not two-stage")
  }
}

#' plotStage1FPRates
#'
#' False positives rates at the field stage for each species.
#'
#' @details
#' Plots the 95% credible interval of the false positives rate at the field stage
#'
#' @param fitModel Output from the function runOccJSDM
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotStage1FPRates(fitModel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotStage1FPRates <- function(fitModel,
                              idx_species = NULL){

  S <- fitModel$infos$S
  # ncov_theta <- fitModel$infos$ncov_psi
  speciesNames <- fitModel$infos$speciesNames
  primerNames <- fitModel$infos$primerNames

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  samples_subset <- fitModel$results_output$theta0_output
  samples_subset <- apply(samples_subset, 1, c)

  data_plot <- apply(samples_subset, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame

  names(data_plot) <- c("Min", "Max")

  plotSpeciesRates(data_plot, idx_species, speciesNames) +
    ggtitle("Stage 1 FP rates")

}

#' plotStage2FPRates
#'
#' False positives rates at the lab stage for each species.
#'
#' @details
#' Plots the 95% credible interval of the false positives at the lab stage
#'
#' @param fitModel Output from the function runOccJSDM
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotStage2FPRates(fitModel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import stringr
#' @import ggplot2
#'
plotStage2FPRates <- function(fitModel,
                              idx_species = NULL){

  S <- fitModel$infos$S
  speciesNames <- fitModel$infos$speciesNames
  primerNames <- fitModel$infos$primerNames
  maxP <- length(primerNames)

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  # q_output: (maxP, S, niter, nchain)
  q_output <- fitModel$results_output$q_output

  # Build a data frame of quantiles per (Primer, Species)
  data_plot <- expand.grid(Primer = seq_len(maxP), Species = seq_len(S)) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      `2.5%`  = quantile(q_output[Primer, Species, , ], probs = 0.025, na.rm = TRUE),
      `97.5%` = quantile(q_output[Primer, Species, , ], probs = 0.975, na.rm = TRUE)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(Species = speciesNames[Species],
                  Primer  = primerNames[Primer]) |>
    dplyr::filter(Species %in% speciesNames[idx_species])

  # Order species by lower bound of first primer
  first_primer <- data_plot |>
    dplyr::filter(Primer == primerNames[1])
  species_order <- first_primer$Species[order(first_primer$`2.5%`)]

  detectionRates <- data_plot |>
    ggplot(aes(x = factor(Species, levels = species_order),
               ymin = `2.5%`,
               ymax = `97.5%`,
               color = factor(Primer))) +
    geom_errorbar(position = position_dodge(width = .15),
                  width = .5) +
    xlab("Species") +
    ggtitle("Stage 2 FP rates") +
    theme_bw() +
    ylab("q") +
    labs(color = "Primer") +
    theme(
      axis.text = element_text(angle = 0, size = 8),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = .5, size = 15)
    ) +
    coord_flip()

  detectionRates

}

# CORRELATION MATRIX -----

#' plotResidualCorrelationMatrix
#'
#' Plot the residual correlation matrix after accounting for the covariates
#'
#' @details
#' Plots the posterior median of the correlation matrix,
#' with nonsignificant correlations marked with an X
#'
#' @param fitModel Output from the function runOccJSDM
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#' @param showSignificance Should an X be shown for non significant elements?
#' @param confidence Confidence level used to assess significance, default to .95
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotResidualCorrelationMatrix(fitModel)
#' }
#'
#' @export
#' @import dplyr
#' @importFrom ggcorrplot ggcorrplot
#'
plotResidualCorrelationMatrix <- function(fitModel,
                                          idx_species = NULL,
                                          showSignificance = T,
                                          confidence = .95){

  L_output <- fitModel$results_output$jsdm_output$L_output
  speciesNames <- fitModel$infos$speciesNames

  plotCorrelationMatrix(L_output,
                        idx_species,
                        speciesNames,
                        showSignificance,
                        confidence)

}

#' returnResidualCorrelationMatrix
#'
#' Return the quantiles of the residual correlation matrix after accounting for the covariates
#'
#' @details
#' Return the quantiles of the residual correlation matrix after accounting for the covariates
#'
#' @param fitModel Output from the function runOccJSDM
#' @param confidence Confidence level used to assess significance, default to .95
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' returnResidualCorrelationMatrix(fitModel)
#' }
#'
#' @export
#' @import dplyr
#'
returnResidualCorrelationMatrix <- function(fitModel,
                                            confidence = .95){

  L_output <- fitModel$results_output$jsdm_output$L_output
  speciesNames <- fitModel$infos$speciesNames

  Sigma_output <- returnCorrelationMatrixOutput(L_output,
                                                seq_along(speciesNames),
                                                speciesNames)

  conflevels <- c((1 - confidence)/2, .5, (1 + confidence)/2)

  Sigma_quantiles <- apply(Sigma_output, c(2,3),
                           function(x){quantile(x, probs = conflevels)})

  Sigma_quantiles
}


# ORDINATON -----

#' returnOrdinationScores
#'
#' Return the quantiles of the factors scores for each observation.
#'
#' @details
#'  Return the quantiles of the factors scores for each observation.
#'
#' @param fitModel Output from the function runOccJSDM
#' @param confidence Confidence level
#'
#' @return An object of size (quantiles x number of sites x number of factors)
#'
#' @export
#' @import dplyr
#'
returnOrdinationScores <- function(fitModel,
                                  confidence = .95){

  if(fitModel$infos$n_factors == 0) stop("No factor used")

  siteNames <- fitModel$infos$siteNames
  factorScoresOutput <- returnFactorScores(fitModel$results_output$jsdm_output,
                                           confidence)

  dimnames(factorScoresOutput)[[2]] <- fitModel$infos$siteNames

  factorScoresOutput
}

#' plotOrdinationScores
#'
#' Plot the ordination scores with their credible interval
#'
#' @details
#' Plot the ordination scores with their credible interval
#'
#' @param fitModel Output from the function runOccJSDM
#' @param idx_factors Which factors to plot (2 should be selected)
#' @param confidence Confidence of the quantiles, default to
#'
#' @return A ggplot object
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotOrdinationScores <- function(fitModel,
                                 idx_factors = c(1,2),
                                 confidence = .95){

  n_factors <- fitModel$infos$n_factors
  siteNames <- fitModel$infos$siteNames

  plotFactorScores(fitModel$results_output$jsdm_output,
                   idx_factors,
                   siteNames)

}

#' returnFactorLoadings
#'
#' Return the quantiles of the factors loadings for each species
#'
#' @details
#' Return the quantiles of the factors loadings for each species
#'
#' @param fitModel Output from the function runOccJSDM
#' @param confidence Confidence level
#'
#' @return An object of size (quantiles x number of factors x number of species)
#'
#' @export
#' @import dplyr
#'
returnFactorLoadings <- function(fitModel,
                                  confidence = .95){

  if(fitModel$infos$n_factors == 0) stop("No factor used")

  speciesNames <- fitModel$infos$speciesNames
  factorLoadingsOutput <- returnFactorLoadings_jsdm(fitModel$results_output$jsdm_output,
                                           confidence)

  dimnames(factorLoadingsOutput)[[3]] <- speciesNames

  factorLoadingsOutput
}

#' plotFactorLoadings
#'
#' Plot the ordination loadings with their credible interval
#'
#' @details
#' Plot the ordination loadings with their credible interval
#'
#' @param fitModel Output from the function runOccJSDM
#' @param idx_factors Which factors to plot (2 should be selected)
#' @param confidence Confidence of the quantiles, default to 95
#'
#' @return A ggplot object
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotFactorLoadings <- function(fitModel,
                                 idx_factors = c(1,2),
                                 confidence = .95){

  n_factors <- fitModel$infos$n_factors
  speciesNames <- fitModel$infos$speciesNames

  plotFactorLoadings_jsdm(fitModel$results_output$jsdm_output,
                          idx_factors,
                          speciesNames)

}

#' plotBiplot
#'
#' Plot a biplot of site ordination scores and species factor loadings,
#' using posterior medians only (no credible intervals).
#'
#' @details
#' Sites are plotted as points at their posterior median factor scores
#' (see \code{\link{returnOrdinationScores}}), and species are plotted as
#' arrows at their posterior median factor loadings (see
#' \code{\link{returnFactorLoadings}}). Species arrows are rescaled so
#' their extent matches \code{arrow_scale} times the site score radius,
#' purely for visual balance -- this does not affect relative arrow
#' lengths or directions.
#'
#' Under the model's identifiability constraint, species \code{i} anchors
#' factor \code{i} (its loading is fixed to 1 rather than estimated; see
#' \code{\link{returnFactorLoadings}}). When one of the fixed
#' species/factor pairs falls within \code{idx_factors}, a caption flags
#' this so the fixed arrow(s) are not mistaken for an estimated pattern.
#'
#' @param fitModel Output from the function runOccJSDM
#' @param idx_factors Which factors to plot (2 should be selected)
#' @param arrow_scale Target extent of the species arrows, as a fraction
#' of the maximum site score radius (default 0.8)
#' @param label_size Species label text size (default 3.5)
#'
#' @return A ggplot object
#'
#' @export
#' @import dplyr
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#'
plotBiplot <- function(fitModel,
                        idx_factors = c(1, 2),
                        arrow_scale = 0.8,
                        label_size = 3.5){

  n_factors <- fitModel$infos$n_factors

  if(length(idx_factors) != 2) stop("idx_factors must have length 2")
  if(any(idx_factors > n_factors)){
    stop("idx_factors exceeds the number of factors fit (", n_factors, ")")
  }

  siteScoresArr <- returnOrdinationScores(fitModel)
  speciesLoadingsArr <- returnFactorLoadings(fitModel)

  siteNames <- dimnames(siteScoresArr)[[2]]
  speciesNames <- dimnames(speciesLoadingsArr)[[3]]

  site_df <- data.frame(
    site = siteNames,
    x = siteScoresArr["50%", , idx_factors[1]],
    y = siteScoresArr["50%", , idx_factors[2]]
  )

  species_df <- data.frame(
    species = speciesNames,
    x = speciesLoadingsArr["50%", idx_factors[1], ],
    y = speciesLoadingsArr["50%", idx_factors[2], ]
  )

  # Rescale species arrows to match the site score spread (standard
  # biplot scaling); guard against a degenerate all-zero loadings case
  site_radius <- max(sqrt(site_df$x^2 + site_df$y^2))
  species_radius <- max(sqrt(species_df$x^2 + species_df$y^2))
  scale_factor <- if(species_radius > 0){
    arrow_scale * site_radius / species_radius
  } else {
    1
  }
  species_df$x <- species_df$x * scale_factor
  species_df$y <- species_df$y * scale_factor

  # Species i anchors factor i under the model's identifiability
  # constraint (reparamFactorModel()); only flag anchors whose factor is
  # actually being plotted
  fixed_factors <- idx_factors[idx_factors <= n_factors & idx_factors <= length(speciesNames)]
  fixed_species <- speciesNames[fixed_factors]

  caption_text <- NULL
  if(length(fixed_species) > 0){
    pairs_txt <- paste0(fixed_species, " (Factor ", fixed_factors, ")")
    verb <- if(length(pairs_txt) > 1) "loadings are" else "loading is"
    caption_text <- paste0(
      paste(pairs_txt, collapse = " and "), " ", verb,
      "\nfixed to 1 for identifiability and ",
      if(length(pairs_txt) > 1) "are" else "is",
      " not estimated."
    )
  }

  ggplot() +
    geom_point(data = site_df, aes(x, y), color = "grey50", alpha = 0.6) +
    geom_segment(
      data = species_df,
      aes(x = 0, y = 0, xend = x, yend = y),
      arrow = arrow(length = unit(0.2, "cm")),
      color = "steelblue"
    ) +
    ggrepel::geom_text_repel(
      data = species_df,
      aes(x, y, label = species),
      color = "steelblue",
      size = label_size
    ) +
    labs(
      x = paste0("Factor ", idx_factors[1]),
      y = paste0("Factor ", idx_factors[2]),
      caption = caption_text
    ) +
    theme_minimal() +
    theme(plot.caption = element_text(hjust = 0, face = "italic", size = 8))
}

# PREDICTIONS --------

createSpatialPredMatrix <- function(Xs, l_s_grid, X_tilde, list_Xs_mat){

  ps <- nrow(X_tilde)

  Xs <- transformCovariatesMatrix(Xs, list_Xs_mat, remove_intercept = T)
  n <- nrow(Xs)

  length_grid_ls <- length(l_s_grid)
  Ks_all <- array(NA, dim = c(n, ps, length_grid_ls))

  for (j in 1:length_grid_ls) {

    l_s <- l_s_grid[j]

    K_uu <- K2(X_tilde, X_tilde, 1, l_s) + diag(0.0001, nrow = nrow(X_tilde))
    L_Kmm <- FastGP::rcppeigen_get_chol(K_uu)
    invL_Kmm <- FastGP::rcppeigen_invert_matrix(L_Kmm)
    K_staru <- K2(Xs, X_tilde, 1, l_s)
    KnmLmt <- K_staru %*% t(invL_Kmm)

    Ks_all[,,j] <- KnmLmt

  }

  Ks_all

}

#' predictNewSites
#'
#' Computes the quantiles of the predictive occupancy probability at new sites
#'
#' @details
#' Compute the credible interval of the occupancy probability at new sites
#'
#' @param fitModel Output from the function runOccJSDM
#' @param X_psi Occupancy covariates matrix for the new locations. Required
#' only when the occupancy-covariate term is used; defaults to \code{NULL}.
#' @param X_s Spatial/ordination covariates matrix for the new locations.
#' Required only when the spatial term is used; defaults to \code{NULL}.
#' @param useEnvCov Logical, or \code{NULL} (default). Whether to include the
#' effect of the occupancy covariates (\code{X_psi}) in the prediction. If
#' \code{NULL}, the term is used when occupancy covariates were estimated in
#' \code{fitModel} and skipped otherwise. Setting \code{TRUE} explicitly on a
#' fit that estimated none is an error rather than a silent no-op.
#' @param useSpatial Logical, or \code{NULL} (default). Whether to include the
#' spatially autocorrelated random field (based on \code{X_s}) in the
#' prediction. Resolved the same way as \code{useEnvCov}.
#' @param useBiotic Logical, or \code{NULL} (default). Whether to include
#' the latent-factor (residual species covariance) term in the prediction.
#' Resolved the same way as \code{useEnvCov}.
#' @param summarised Should the output be return in the form of quantiles? Set to TRUE if the number of sites is very large
#' @param confidence If quantiles are returned, the confidence level of the quantiles.
#'
#' @return An array of size (,sites,species) with either the quantiles or the
#' iterations in the first dimension
#'
#'
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
predictNewSites <- function(fitModel,
                            X_psi = NULL,
                            X_s = NULL,
                            useEnvCov = NULL,
                            useSpatial = NULL,
                            useBiotic = NULL,
                            summarised = T,
                            confidence = .95
){

  S <- fitModel$infos$S
  speciesNames <- fitModel$infos$speciesNames

  areEnvCovEstimated <- fitModel$infos$ncov_psi > 0
  isSpatFieldEstimated <- fitModel$infos$ps > 0
  areFactorsEstimated <- fitModel$infos$n_factors > 0

  # Resolve the three component switches. NULL means "use it if the fit
  # estimated it", TRUE means "use it, and complain if it was not estimated",
  # FALSE means "skip it". useBiotic already worked this way; useEnvCov and
  # useSpatial now match, which is what their own documentation always
  # claimed. Previously both defaulted to TRUE and useSpatial hard-stopped on
  # a fit with no spatial field, so the documented "ignored if not estimated"
  # behaviour was unreachable.
  resolveTerm <- function(flag, estimated, what){
    if(is.null(flag)) return(estimated)
    if(isTRUE(flag) && !estimated){
      stop("Cannot use ", what, ": it was not estimated in this fit.",
           call. = FALSE)
    }
    isTRUE(flag)
  }

  useEnvCov  <- resolveTerm(useEnvCov,  areEnvCovEstimated,   "occupancy covariates")
  useSpatial <- resolveTerm(useSpatial, isSpatFieldEstimated, "the spatial field")
  useBiotic  <- resolveTerm(useBiotic,  areFactorsEstimated,  "latent factors")

  # `&&` rather than `&` throughout: `&` evaluates both sides, which forced
  # the X_psi/X_s promises even when the caller had asked for neither term.
  # With no defaults on those arguments that made every such call fail with
  # R's generic "argument is missing" rather than the message below.
  if(useEnvCov && is.null(X_psi)) {
    stop("X_psi is required to predict with occupancy covariates. Supply a ",
         "matrix of covariates for the new sites, or set useEnvCov = FALSE.",
         call. = FALSE)
  }

  if(useSpatial && is.null(X_s)) {
    stop("X_s is required to predict with the spatial field. Supply the ",
         "coordinates of the new sites, or set useSpatial = FALSE.",
         call. = FALSE)
  }

  # create env cov matrix
  if(useEnvCov) {

    X_psi <- transform_new_covariates(
      X_psi,
      fitModel$infos$list_X_psi_mat,
      remove_intercept = TRUE
    )

  } else {

    X_psi <- matrix(NA, 0, 0)

  }

  # create spatial matrix
  if(useSpatial){
    Ks <-
      createSpatialPredMatrix(X_s,
                              fitModel$infos$l_s_grid,
                              fitModel$infos$list_Xs$X_tilde,
                              fitModel$infos$list_Xs_mat)
  } else {
    Ks <- array(NA, dim = c(0,0,0))
  }

  B0_output <- fitModel$results_output$jsdm_output$B0_output
  B_output <- fitModel$results_output$jsdm_output$B_output
  Bs_output <- fitModel$results_output$jsdm_output$Bs_output
  L_output <- fitModel$results_output$jsdm_output$L_output
  sigmah_output <- fitModel$results_output$jsdm_output$sigmah_output
  idx_ls_output <- fitModel$results_output$jsdm_output$idx_ls_output

  if(!summarised){

    stop("Only summarised version for now")

  } else {

    conflevels <- c((1 - confidence)/2, .5, (1 + confidence)/2)

    B0_output_vec <- aperm(apply(B0_output, 1, c), c(2,1))
    B_output_vec <- aperm(apply(B_output, c(1,2), c), c(2,3,1))
    L_output_vec <- aperm(apply(L_output, c(1,2), c), c(2,3,1))

    # A fit with no spatial field still returns Bs_output, but with a
    # zero-length first dimension (0 x S x niter x nchain). apply() over
    # margins (1,2) then has no cells to work on and drops the dimensions,
    # so the aperm() below fails with "'perm' is of wrong length". Only
    # reshape it when the spatial term is actually in play; otherwise hand
    # the C++ an empty cube, exactly as Ks already is.
    if(useSpatial){
      Bs_output_vec <- aperm(apply(Bs_output, c(1,2), c), c(2,3,1))
    } else {
      Bs_output_vec <- array(NA_real_, dim = c(0,0,0))
    }

    sigmah_output_vec <- as.vector(sigmah_output)
    idx_ls_output_vec <- as.vector(idx_ls_output)

    pred_output <- computeNewOutputs(
      X_psi,
      B0_output_vec,
      B_output_vec,
      Ks,
      Bs_output_vec,
      L_output_vec,
      sigmah_output_vec,
      idx_ls_output_vec,
      conflevels,
      useEnvCov, useSpatial, useBiotic,
      fitModel$infos$jsdmModel)

  }

  pred_output

}

# SITE-SAMPLE SUMMARIES ----------

#' computePredictiveOccupancyProbs
#'
#' Computes the average predictive occupancy probabilities at the existing site
#'
#' @details
#' Computes the average predictive occupancy probabilities
#'
#' @param fitModel Output from the function runOccJSDM
#'
#' @return An matrix of size (sites,species) with the posterior means
#'
#'
#' @examples
#' \dontrun{
#' computePredictiveOccupancyProbs(fitModel)
#' }
#'
#' @export
#' @import dplyr
#'
computePredictiveOccupancyProbs <- function(fitModel){

  psi_output <- fitModel$results_output$psi_output

  if(length(dim(psi_output))== 4){
    psi_mean <- apply(psi_output, c(1,2), mean)
  } else {
    psi_mean <- psi_output
  }

  rownames(psi_mean) <- fitModel$infos$siteNames
  colnames(psi_mean) <- fitModel$infos$speciesNames

  psi_mean

}

#' computeAverageCollectionProbs
#'
#' Computes the average collection probabilities
#'
#' @details
#' Computes the average predictive occupancy probabilities
#'
#' @param fitModel Output from the function runOccJSDM
#'
#' @return An array of size (samples,species)
#'
#'
#' @examples
#' \dontrun{
#' computeAverageCollectionProbs(fitModel)
#' }
#'
#' @export
#' @import dplyr
#'
computeAverageCollectionProbs <- function(fitModel){

  theta_mean <- fitModel$results_output$theta_output

  colnames(theta_mean) <- fitModel$infos$speciesNames

  theta_mean

}

#' computeConditionalOccupancyProbs
#'
#' Computes the posterior mean of the conditional occupancy probability
#'
#' @details
#' Computes the posterior mean of the conditional occupancy probability
#'
#' @param fitModel Output from the function runOccJSDM
#'
#' @return A matrix of size (site X species) with the posterior mean occupancy at
#' each site for each species.
#'
#' @examples
#' \dontrun{
#' computeConditionalOccupancyProbs(fitModel)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
computeConditionalOccupancyProbs <- function(fitModel){

  z_output <- fitModel$results_output$z_output

  if(length(dim(z_output))== 4){
    z_mean <- apply(z_output, c(1,2), mean)
  } else {
    z_mean <- z_output
  }

  rownames(z_mean) <- fitModel$infos$siteNames
  colnames(z_mean) <- fitModel$infos$speciesNames

  z_mean

}

#' computeConditionalSamplePresenceProbs
#'
#' Computes the posterior mean of a sample being occupied
#'
#' @details
#' Computes the posterior mean of a sample being occupied
#'
#' @param fitModel Output from the function runOccJSDM
#'
#' @return A matrix of size (samples X species) with the posterior mean presence
#' at each sample for each species.
#'
#' @examples
#' \dontrun{
#' computeConditionalSamplePresenceProbs(fitModel)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
computeConditionalSamplePresenceProbs <- function(fitModel){

  w_mean <- fitModel$results_output$w_output

  colnames(w_mean) <- fitModel$infos$speciesNames

  w_mean

}

#' returnLatentPresences
#'
#' Compute the latent presences
#'
#' @details
#' Computes, for a single species, the posterior mean of the latent
#' occupancy state `z` and the latent detection state `w` across sites.
#'
#' @param fitModel Output from the function runOccJSDM
#' @param idx_species Index of the species to compute latent presences for (default 1)
#'
#' @return A matrix
#'
#' @note As currently written, the function computes `z_mean` and `w_mean`
#' but then returns `returnVariancePartitioningMatrix(varPart_output, speciesNames)`,
#' where `varPart_output` is not defined anywhere in the function body. This
#' looks like a pre-existing bug (likely leftover from a copy-paste) that will
#' error at runtime and should be fixed separately.
#'
#' @examples
#' \dontrun{
#' returnLatentPresences(fitModel)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
returnLatentPresences <- function(fitModel, idx_species = 1){

  if(fitModel$infos$model == "two_stage"){
    z_mean <- computeConditionalOccupancyProbs(fitModel)
    w_mean <- computeConditionalSamplePresenceProbs(fitModel)
    psi_mean <- computePredictiveOccupancyProbs(fitModel)
    theta_mean <- computeAverageCollectionProbs(fitModel)

    p_mean <- apply(fitModel$results_output$p_output, c(1,2), mean)

    siteNames <- fitModel$infos$siteNames

    df <- fitModel$infos$data_info[,c("Site","Sample","Primer")]

    list_idx <- fitModel$infos$list_idx

    df$OTU <- fitModel$infos$OTU[,idx_species]

    df$CondOccProb <- z_mean[list_idx$idx_z_k,idx_species]
    df$CondSampleProb <- w_mean[list_idx$idx_w_k,idx_species]
    df$PredOccProb<- psi_mean[list_idx$idx_z_k,idx_species]
    df$CollectionProb <- theta_mean[list_idx$idx_w_k,idx_species]
    df$DetectionProb <- p_mean[list_idx$idx_p_k,idx_species]

    data_info <- fitModel$infos$data_info
    idx_columns <-  which(!(names(data_info) %in% c("Site", "Sample","Primer")))
    data_info_covariates <- data_info[,idx_columns]

    df <- cbind(df, data_info_covariates)
    df
  } else {
    stop("Only two stage model supported")
  }

}


#' plotLatentPresences
#'
#' Display latent presences as a colour-banded gt table
#'
#' @details
#' Takes the data frame returned by \code{\link{returnLatentPresences}} and
#' renders a \pkg{gt} table with a nested banding scheme: fill colour changes
#' with each new Site, and alternates between lighter and darker shades for
#' each Sample and Primer within a Site. This makes the hierarchical
#' site/sample/primer structure visually clear.
#'
#' Requires the \pkg{gt} and \pkg{colorspace} packages.
#'
#' @param latentPresences Data frame returned by \code{\link{returnLatentPresences}}
#' @param sticky_header Logical; if \code{TRUE}, the table header stays fixed
#'   while scrolling in the RStudio Viewer (default \code{TRUE})
#' @param title Character string for the table title
#'   (default \code{"Latent presence summary"})
#' @param subtitle Character string for the subtitle (default describes the
#'   colour scheme)
#' @param container_height Height of the scrollable container in pixels
#'   (default 1200)
#' @param base_colors Character vector of base colours to cycle through for
#'   successive sites (default: a 4-colour Tableau-inspired palette)
#' @param shade_amounts Numeric vector of length 4 giving the lighten amounts
#'   for the four sub-shades within each site colour (ordered lightest to
#'   darkest; default \code{c(0.78, 0.60, 0.45, 0.28)})
#' @param species_name Optional character string identifying the species. If
#'   provided, it is appended to \code{title} in parentheses
#'   (e.g. \code{"Latent presence summary (Lepisosteus osseus)"})
#' @param decimals Number of decimal places for numeric columns (default 2)
#' @param columns Character vector of column names to display. If \code{NULL}
#'   (default), all columns are shown
#'
#' @return A \code{gt_tbl} object
#'
#' @examples
#' \dontrun{
#' lp <- returnLatentPresences(fitModel, idx_species = 1)
#' plotLatentPresences(lp, title = "Species 1")
#' }
#'
#' @export
#'
plotLatentPresences <- function(latentPresences,
                                sticky_header = TRUE,
                                title = "Latent presence summary",
                                subtitle = paste0(
                                  "Colour: Site | ",
                                  "Shade: Sample (coarse) x Primer (fine)"
                                ),
                                container_height = 1200,
                                base_colors = c("#4E79A7", "#F28E2B",
                                                "#59A14F", "#E15759"),
                                shade_amounts = c(0.78, 0.60, 0.45, 0.28),
                                species_name = NULL,
                                decimals = 2,
                                columns = NULL) {

  if (!is.null(species_name)) {
    title <- paste0(title, " (", species_name, ")")
  }

  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' is required. Install it with install.packages('gt').")
  }
  if (!requireNamespace("colorspace", quietly = TRUE)) {
    stop("Package 'colorspace' is required. Install it with install.packages('colorspace').")
  }

  stopifnot(length(shade_amounts) == 4)

  # Build shade lookup: rows = base colours, cols = 4 sub-shades

  shade_lookup <- sapply(shade_amounts, function(amt) {
    colorspace::lighten(base_colors, amt)
  })

  df_banded <- latentPresences |>
    dplyr::mutate(
      site_group   = cumsum(Site   != dplyr::lag(Site,   default = dplyr::first(Site))),
      sample_group = cumsum(Sample != dplyr::lag(Sample, default = dplyr::first(Sample))),
      primer_group = cumsum(Primer != dplyr::lag(Primer, default = dplyr::first(Primer))),
      site_color_idx  = (site_group %% length(base_colors)) + 1,
      sample_parity   = sample_group %% 2,
      primer_parity   = primer_group %% 2,
      shade_level     = sample_parity * 2 + primer_parity + 1,
      row_color       = shade_lookup[cbind(site_color_idx, shade_level)]
    )

  # Select columns

  if (!is.null(columns)) {
    display_cols <- intersect(columns, names(df_banded))
    if (length(display_cols) == 0) {
      stop("None of the requested columns found in latentPresences.")
    }
    df_display <- df_banded[, c(display_cols, "row_color"), drop = FALSE]
  } else {
    # Drop the helper columns used for banding (except row_color, hidden later)
    helper_cols <- c("site_group", "sample_group", "primer_group",
                     "site_color_idx", "sample_parity", "primer_parity",
                     "shade_level")
    df_display <- df_banded[, !(names(df_banded) %in% helper_cols), drop = FALSE]
  }

  # Identify numeric columns to format (exclude row_color and id columns)
  id_cols <- c("Site", "Sample", "Primer", "OTU", "row_color")
  numeric_cols <- setdiff(
    names(df_display)[vapply(df_display, is.numeric, logical(1))],
    id_cols
  )

  tbl <- df_display |>
    gt::gt() |>
    gt::fmt_number(columns = dplyr::all_of(numeric_cols), decimals = decimals) |>
    gt::cols_hide("row_color") |>
    gt::tab_header(title = title, subtitle = subtitle) |>
    gt::tab_options(
      container.height = gt::px(container_height),
      container.overflow.y = "auto"
    ) |>
    gt::tab_source_note(
      source_note = gt::html(paste0(
        "<b>Colour</b> = Site &middot; <b>shade</b> = Sample (coarse) ",
        "within a Site &middot; <b>fine shade</b> = Primer within a Sample."
      ))
    )

  # Apply row-level fill colours
  for (col in unique(df_banded$row_color)) {
    tbl <- tbl |>
      gt::tab_style(
        style = gt::cell_fill(color = col),
        locations = gt::cells_body(rows = df_banded$row_color == col)
      )
  }

  # Sticky header CSS
  if (sticky_header) {
    tbl <- tbl |>
      gt::opt_css(
        css = "
        .gt_table thead {
          position: sticky;
          top: 0;
          z-index: 2;
        }
        "
      )
  }

  tbl
}


# OTHER -------

#' extractWAIC
#'
#' Compute the WAIC for model comparison
#'
#' @details
#' Compute the WAIC for model comparison
#'
#' @param fitModel Output from the function runOccJSDM
#'
#' @return The WAIC
#'
#' @examples
#' \dontrun{
#' extractWAIC(fitModel)
#' }
#'
#' @export
#'
extractWAIC <- function(fitModel){

  fitModel$results_output$WAIC

}

#' returnVariancePartitioning
#'
#' Compute the variance partitioning for each species
#'
#' @details
#' Compute the variance partitioning for each species
#'
#' @param fitModel Output from the function runOccJSDM
#'
#' @return A matrix of size (species x 4)
#'
#' @examples
#' \dontrun{
#' returnVariancePartitioning(fitModel)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
returnVariancePartitioning <- function(fitModel){

  varPart_output <- fitModel$results_output$jsdm_output$varPart_output

  speciesNames <- fitModel$infos$speciesNames

  returnVariancePartitioningMatrix(varPart_output, speciesNames)

}

#' plotVariancePartitioning
#'
#' Plot the variance partitioning for each species
#'
#' @details
#' Plot the variance partitioning for each species
#'
#' @param fitModel Output from the function runOccJSDM
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotVariancePartitioning(fitModel)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#' @importFrom ggtern ggtern theme_showarrows
#'
plotVariancePartitioning <- function(fitModel){

  varPart_output <- fitModel$results_output$jsdm_output$varPart_output

  speciesNames <- fitModel$infos$speciesNames

  plotVarPart(varPart_output, speciesNames)

}

#' plotReadIntensity
#'
#' Plot the reads distribution under the true positives and false positives
#'
#' @details
#' Plots the 95% credible interval of density of true positives and false positives
#'
#' @param fitModel Output from the function runOccJSDM
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotReadIntensity(fitModel)
#' }
#'
#' @import dplyr
#' @import ggplot2
#'
plotReadIntensity <- function(fitModel){

  mu1_output <- fitModel$results_output$mu1_output
  mu0_output <- fitModel$results_output$mu0_output
  sigma1_output <- fitModel$results_output$sigma1_output
  sigma0_output <- fitModel$results_output$sigma0_output

  niter <- length(mu1_output)

  # x_grid <- seq(1, fitModel$infos$maxexplogy1, by = 5)
  x_grid <- exp(seq(log(1), log(fitModel$infos$maxexplogy1), length.out = 250))

  # seq(1, fitModel$infos$maxexplogy1, by = 5)

  x <- log(x_grid + 1)

  densities_plot_pos <- matrix(NA, length(x_grid), niter)
  densities_plot_neg <- matrix(NA, length(x_grid), niter)

  for (iter in 1:niter) {

    mu1 <- mu1_output[iter]
    mu0 <- mu0_output[iter]
    sigma1 <- sigma1_output[iter]
    sigma0 <- sigma0_output[iter]

    densities_plot_pos[,iter] <- dnorm(x, mean = mu1, sd = sigma1)
    densities_plot_neg[,iter] <- dnorm(x, mean = mu0, sd = sigma0)
  }

  densities_plot_pos_quantiles <-
    apply(densities_plot_pos, 1,
          function(x) {quantile(x, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))}) %>% t %>%
    as.data.frame %>%
    mutate(x = x_grid,
           Type = "True Positive")

  densities_plot_neg_quantiles <-
    apply(densities_plot_neg, 1,
          function(x) {quantile(x, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))}) %>% t %>%
    as.data.frame %>%
    mutate(x = x_grid,
           Type = "False Positives")

  densities_plot_quantiles <-
    rbind(densities_plot_pos_quantiles,
          densities_plot_neg_quantiles)

  # ggplot() +
  #   geom_ribbon(data = densities_plot_pos_quantiles,
  #               aes(x = x_grid,
  #                   ymax = `97.5%`,
  #                   ymin = `2.5%`)) +
  #   geom_line(data = densities_plot_pos_quantiles,
  #             aes(x = x_grid,
  #                 y = `50%`))
  #
  # df0 <- as_tibble(densities_plot_neg) %>%
  #   mutate(x = x_grid) %>%
  #   pivot_longer(cols = -x, names_to = "iter", values_to = "density") %>%
  #   mutate(Type = "False positives")
  #
  # df1 <- as_tibble(densities_plot_pos) %>%
  #   mutate(x = x_grid) %>%
  #   pivot_longer(cols = -x, names_to = "iter", values_to = "density") %>%
  #   mutate(Type = "True positives")
  #
  # # Combine into one data frame
  # df_combined <- bind_rows(df0, df1) %>%
  #   mutate(iter = as.numeric(gsub("V", "", iter)))  # Clean iteration labels

  x_grid_breaks <- round( c(0, 10, 20,
                            x_grid[seq(1, length(x_grid), by = 10)] - 1))

  # ggplot(df_combined, aes(x = x, y = density, group = interaction(iter, Type), color = Type)) +
  #   geom_line(alpha = 0.1, aes(color = Type)) +
  #   scale_color_manual(values = c("False positives" = "blue", "True positives" = "red")) +
  #   labs(title = "Reads distributions",
  #        x = "x", y = "Density") +
  #   theme_minimal() +
  #   scale_x_continuous(
  #     name = "Number of reads",
  #     breaks = x_grid_breaks,  # Custom breaks (e.g., e⁻², e⁻¹, e⁰, e¹, ...)
  #     # labels = function(x) sprintf("%.2f", exp(x) - 1)  # Format labels
  #     trans = "log"
  #   ) +
  #   theme(axis.text.x = element_text(angle = 90),
  #         plot.title = element_text(hjust = 0.5,
  #                                   size = 16,
  #                                   face = "bold"))


  ggplot() +
    geom_ribbon(data = densities_plot_quantiles,
                aes(x = x,
                    ymax = `97.5%`,
                    ymin = `2.5%`,
                    fill = Type),
                alpha = .5) +
    geom_ribbon(data = densities_plot_quantiles,
                aes(x = x,
                    ymax = `90%`,
                    ymin = `10%`,
                    fill = Type),
                alpha = .75) +
    theme_minimal() +
    scale_x_continuous(
      name = "Number of reads",
      breaks = x_grid_breaks,  # Custom breaks (e.g., e⁻², e⁻¹, e⁰, e¹, ...)
      # labels = function(x) sprintf("%.2f", exp(x) - 1)  # Format labels
      trans = "log"
    ) +
    theme(axis.text.x = element_text(angle = 90,
                                     size = 12),
          axis.text.y = element_text(angle = 90,
                                     size = 12),
          plot.title = element_text(hjust = 0.5,
                                    size = 16,
                                    face = "bold"))

}


#' plotCumulativeSpeciesDetections
#'
#' Cumulative species detection plot
#'
#' @details
#' Plot the credible interval of the overall species detections for several
#' number of PCR, conditional on species presence in the sample
#'
#' @param fitModel Output from the function runOccJSDM
#' @param M Maximum number of environmental samples to consider
#' @param K Maximum number of technical replicates (PCRs) to consider
#' @param primer Index of the primer to use; 0 (default) pools across all primers
#' @param alpha Confidence level of the credible interval, default to .95
#' @param byK Should K (technical replicates) be on the x-axis (instead of M)?
#'
#' @return A plot with the credible interval of cumulative detections
#'
#' @examples
#' \dontrun{
#' plotCumulativeSpeciesDetections(fitModel, M = 3, K = 5)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotCumulativeSpeciesDetections <- function(fitModel, M, K, primer = 0, alpha = .95,
                                            byK = T){

  if(fitModel$infos$model != "two_stage") stop("Only allowed for a two-stage model")

  beta_theta_output <- fitModel$results_output$beta_theta_output
  p_output <- fitModel$results_output$p_output

  speciesDetected <- computeSpeciesDetected(beta_theta_output, p_output,
                                               M, K, primer, alpha)

  df <- expand.grid(M = 1:M, K = 1:K)
  df$lower <- as.vector(speciesDetected[1, , ])
  df$upper <- as.vector(speciesDetected[2, , ])

  if(byK){
    p <- ggplot(df, aes(x = K)) +
      geom_errorbar(aes(ymin = lower, ymax = upper)) +
      facet_grid(rows =  vars(M), labeller = as_labeller(function(x) paste("M = ", x))) +
      labs(
        x = "K",
        y = "Value"
      ) +
      scale_x_continuous(breaks = 1:K, name = "Number of technical replicates")

  } else {
    p <- ggplot(df, aes(x = M)) +
      geom_errorbar(aes(ymin = lower, ymax = upper)) +
      facet_grid(rows =  vars(K), labeller = as_labeller(function(x) paste("K = ", x))) +
      labs(
        x = "M",
        y = "Value"
      ) +
      scale_x_continuous(breaks = 1:M, name = "Number of environmental samples")
  }

  p + scale_y_continuous(name = "Species Detected") + theme_minimal()

}

#' computeSpeciesDetected
#'
#' Simulate the number of species detected as a function of the number of
#' technical replicates (PCRs) and environmental samples.
#'
#' @details
#' Internal helper for `plotCumulativeSpeciesDetections`.
#'
#' @noRd
computeSpeciesDetected <- function(beta_theta_output, p_output, M, K, primer, alpha){


  beta_theta_output <- apply(beta_theta_output[1,,,,drop=F], 2, c)
  p_output <- apply(p_output, c(1,2), c)

  S <- dim(beta_theta_output)[2]
  numPrimers <- dim(p_output)[2]
  niter <- dim(beta_theta_output)[1]

  idxObs <- unique(floor(seq(1, niter, length.out = 500)))
  B <- length(idxObs)

  if(primer == 0){
    idx_primer <- 1:numPrimers
  } else {
    idx_primer <- primer
  }

  P <- length(idx_primer)

  detectionsPerSamplePCR <- array(NA, dim = c(B, M, K))

  for (b in 1:B) {

    # detections for each species, sample and PCR
    detectionsPerSamplePCRSpecies <- array(0, dim = c(S, M, K))

    w <- sapply(1:M, function(m){
      sapply(1:S, function(s){
        rbinom(1, 1, logistic(beta_theta_output[idxObs[b],s]))
      })
    })

    # simulate detections
    for (s in 1:S) {

      for (m in 1:M) {


        if(w[s,m] == 1){

          for (k in 1:K) {

            detectionsPerSamplePCRSpecies[s,m,k] <-
              sum(sapply(1:P, function(p){
                rbinom(1, 1, p_output[idxObs[b],p,s])
              }) > 0)

          }

        }

      }

    }

    # cumulative detections for each species, sample and PCR
    detectionsPerSamplePCRSpeciesSummarised <- array(0, dim = c(S, M, K))

    # summarise across the different samples and PCR
    for (s in 1:S) {
      for (m in 1:M) {
        for (k in 1:K) {
          detectionsPerSamplePCRSpeciesSummarised[s,m,k] <-
            any(detectionsPerSamplePCRSpecies[s,1:m,1:k] > 0)
        }
      }
    }


    detectionsPerSamplePCR[b,,] <-
      apply(detectionsPerSamplePCRSpeciesSummarised, c(2,3), sum)
  }

  output <- apply(detectionsPerSamplePCR, c(2,3), function(x){
    quantile(x, probs = c((1 - alpha)/2, (1 + alpha)/2))
  })

  output

}


#' plotOccupancyStates
#'
#' Plot a heatmap comparing estimated latent occupancy against observed
#' detection frequencies, for each site and species.
#'
#' @details
#' For each site and species, plots a tile coloured by the posterior mean
#' latent occupancy probability, labelled with the observed frequency of
#' detection across samples at that site.
#'
#' @param fitModel Output from the function runOccJSDM
#'
#' @return A ggplot object
#'
#' @note This function references `data_info` and `OTU` directly rather than
#' via `fitModel` or function arguments; these objects are not defined in the
#' function body and are not present in `fitModel`, so as written this will
#' error unless matching objects happen to exist in the calling environment.
#' This looks like a pre-existing bug worth fixing separately.
#'
#' @noRd
plotOccupancyStates <- function(fitModel){

  speciesNames <- fitModel$infos$speciesNames
  siteNames <- fitModel$infos$siteNames

  z_output <- fitModel$results_output$z_output
  z_mean <- apply(z_output, c(1,2), mean)
  dimnames(z_mean) <- list(siteNames, speciesNames)

  observedOccupancies <- cbind(data_info, OTU) %>%
    group_by(Site, Sample) %>%
    summarise(
      across(starts_with("OTU_"), function(x) sum(x > 0))
    ) %>%
    group_by(Site) %>%
    summarise(
      across(starts_with("OTU_"), function(x) round(mean(x > 0), 2))
    )  %>% as.data.frame %>% dplyr::select(-Site) %>% as.matrix
  dimnames(observedOccupancies) <- list(siteNames, speciesNames)

  df_colors  <- as.data.frame.table(z_mean) %>%
    rename(Site = Var1, Species = Var2, Occupancy = Freq)

  df_numbers <- as.data.frame.table(observedOccupancies) %>%
    rename(Site = Var1, Species = Var2, Frequency = Freq)

  plot_data <- left_join(df_colors, df_numbers, by = c("Site", "Species")) #%>%

  ggplot(plot_data, aes(x = Species, y = Site)) +
    geom_tile(aes(fill = Occupancy), color = "white", linewidth = 0.5) +
    geom_text(aes(label = Frequency), color = "black", fontface = "bold", size = 4) +
    scale_fill_viridis_c(option = "plasma", name = "Occupancy") +
    # scale_y_reverse() +
    # scale_x_continuous(breaks = 1:5, position = "top") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 12, face = "bold")
    )

}
