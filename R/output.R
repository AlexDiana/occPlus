
logit <- function(x){
  log(x / (1 - x))
}

logistic <- function(x) 1 / (1 + exp(-x))

plotSpeciesRates <- function(data_plot,
                             orderSpecies,
                             subset){

  data_plot %>%
    filter(Species %in% speciesNames[orderSpecies[subset]]) %>%
    ggplot(aes(x =  factor(Species, level = speciesNames[orderSpecies]),
               ymin = `2.5%`,
               ymax = `97.5%`)) + geom_errorbar() +
    xlab("Species") +
    # ylim(c(0,1)) +
    ggtitle("Collection rates") +
    theme_bw() +
    ylim(c(0,1)) +
    theme(
      axis.text = element_text(angle = 90,
                               size = 8),
      plot.title = element_text(hjust = .5,
                                size = 15)
    )

}



#' plotOccupancyCovariates
#'
#' Occupancy covariate coefficients.
#'
#' @details
#' Plots the 95% credible interval of the occupancy covariates coefficients
#'
#' @param fitmodel Output from the function runOccPlus
#' @param covName Name of the covariate to be plotted (same name as in data$info)
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotOccupancyCovariates <- function(fitmodel,
                                    covName = NULL,
                                    idx_species = NULL
){

  if(is.null(covName)){
    stop("No name provided")
  }

  matrix_of_draws <- fitmodel$matrix_of_draws

  S <- fitmodel$infos$S
  ncov_psi <- fitmodel$infos$ncov_psi
  speciesNames <- fitmodel$infos$speciesNames
  occCovNames <- colnames(fitmodel$X_psi)
  idxcov <- which(occCovNames == covName)

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  param <- "beta_psi"

  samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]
  samples_subset <- samples_subset[,idxcov + 0:(S - 1)*ncov_psi]

  data_plot <- apply(samples_subset, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame %>%
    mutate(Species = speciesNames) %>%
    mutate(speciesOrder = order(`2.5%`)) %>%
    filter(Species %in% speciesNames[idx_species])

  orderSpecies <- order(data_plot$`2.5%`)

  plot_occcovs <- data_plot %>%
    ggplot(aes(x =  factor(Species, level = speciesNames[orderSpecies]),
               ymin = `2.5%`,
               ymax = `97.5%`)) + geom_errorbar() +
    xlab("Species") +
    # ylim(c(0,1)) +
    ggtitle(covName) +
    theme_bw() +
    theme(
      axis.text = element_text(angle = 90,
                               size = 8),
      plot.title = element_text(hjust = .5,
                                size = 15)
    )

  plot_occcovs

}


#' plotDetectionCovariates
#'
#' Detection covariate coefficients.
#'
#' @details
#' Plots the 95% credible interval of the occupancy covariates coefficients
#'
#' @param fitmodel Output from the function runOccPlus
#' @param covName Name of the covariate to be plotted (same name as in data$info)
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotDetectionCovariates <- function(fitmodel,
                                    covName = NULL,
                                    idx_species = NULL
){

  if(is.null(covName)){
    stop("No name provided")
  }

  matrix_of_draws <- fitmodel$matrix_of_draws

  S <- fitmodel$infos$S
  ncov_theta <- fitmodel$infos$ncov_theta
  speciesNames <- fitmodel$infos$speciesNames
  occCovNames <- colnames(fitmodel$X_theta)
  idxcov <- which(occCovNames == covName)

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  param <- "beta_theta"

  samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]
  samples_subset <- samples_subset[,idxcov + 0:(S - 1)*ncov_theta]

  data_plot <- apply(samples_subset, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame %>%
    mutate(Species = speciesNames) %>%
    mutate(speciesOrder = order(`2.5%`)) %>%
    filter(Species %in% speciesNames[idx_species])

  orderSpecies <- order(data_plot$`2.5%`)

  plot_occcovs <- data_plot %>%
    ggplot(aes(x =  factor(Species, level = speciesNames[orderSpecies]),
               ymin = `2.5%`,
               ymax = `97.5%`)) + geom_errorbar() +
    xlab("Species") +
    # ylim(c(0,1)) +
    ggtitle(covName) +
    theme_bw() +
    theme(
      axis.text = element_text(angle = 90,
                               size = 8),
      plot.title = element_text(hjust = .5,
                                size = 15)
    )

  plot_occcovs

}

#' plotOrdinationCovariates
#'
#' Ordination covariate coefficients.
#'
#' @details
#' Plots the 95% credible interval of the ordination covariates coefficients
#'
#' @param fitmodel Output from the function runOccPlus
#' @param covName Name of the covariate to be plotted (same name as in data$info)
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotOrdinationCovariates <- function(fitmodel,
                                    covName = NULL,
                                    idx_factors = NULL
){

  if(is.null(covName)){
    stop("No name provided")
  }

  matrix_of_draws <- fitmodel$matrix_of_draws

  d <- fitmodel$infos$d
  ncov_ord <- fitmodel$infos$ncov_ord
  speciesNames <- fitmodel$infos$speciesNames
  ordCovNames <- colnames(fitmodel$X_ord)
  idxcov <- which(ordCovNames == covName)

  if(is.null(idx_factors)){
    idx_factors <- 1:d
  }

  param <- "beta_ord"

  samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]
  samples_subset <- samples_subset[,idxcov + 0:(d - 1)*ncov_ord]

  data_plot <- apply(samples_subset, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame %>%
    mutate(Factor = 1:d) %>%
    filter(Factor %in% idx_factors)

  plot_covs <- data_plot %>%
    ggplot(aes(x = Factor,
               ymin = `2.5%`,
               ymax = `97.5%`)) + geom_errorbar() +
    xlab("Factors") +
    # ylim(c(0,1)) +
    ggtitle(covName) +
    theme_bw() +
    theme(
      axis.text = element_text(angle = 90,
                               size = 8),
      plot.title = element_text(hjust = .5,
                                size = 15)
    )

  plot_covs

}

#' plotOccupancyRates
#'
#' Baseline occupancy rate for each species.
#'
#' @details
#' Plots the 95% credible interval of the baseline occupancy rates
#'
#' @param fitmodel Output from the function runOccPlus
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return The credible interval plot
#'
#' @examples
#' \dontrun{
#' plotSpeciesRates(fitmodel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotOccupancyRates <- function(fitmodel,
                               idx_species = NULL){

  matrix_of_draws <- fitmodel$matrix_of_draws

  samples_subset <-
    matrix_of_draws[,grepl("beta0_psi\\[", colnames(matrix_of_draws))]

  niter <- nrow(samples_subset)
  S <- fitmodel$infos$S
  speciesNames <- fitmodel$infos$speciesNames

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  data_plot <- apply(samples_subset, 2, function(x) {
    quantile(logistic(x), probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame %>%
    mutate(Species = speciesNames) %>%
    mutate(speciesOrder = order(`2.5%`)) %>%
    filter(Species %in% speciesNames[idx_species])

  orderSpecies <- order(data_plot$`2.5%`)


  plot_occupancyrates <- data_plot %>%
    ggplot(aes(x =  factor(Species, level = speciesNames[orderSpecies]),
               ymin = `2.5%`,
               ymax = `97.5%`)) + geom_errorbar() +
    xlab("Species") +
    # ylim(c(0,1)) +
    ggtitle("Baseline Occupancy rates") +
    theme_bw() +
    # ylim(c(0,1)) +
    ylab("") +
    theme(
      axis.text = element_text(angle = 90,
                               size = 8),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = .5,
                                size = 15)
    )

  plot_occupancyrates



}

#' plotCollectionRates
#'
#' Baseline detection rate for each species.
#'
#' @details
#' Plots the 95% credible interval of the baseline detection rates
#'
#' @param fitmodel Output from the function runOccPlus
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return The credible interval plot
#'
#' @examples
#' \dontrun{
#' plotSpeciesRates(fitmodel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotCollectionRates <- function(fitmodel,
                                idx_species = NULL){

  matrix_of_draws <- fitmodel$matrix_of_draws

  S <- fitmodel$infos$S
  ncov_theta <- fitmodel$infos$ncov_theta
  speciesNames <- fitmodel$infos$speciesNames

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  param <- "beta_theta"

  samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]
  samples_subset <- samples_subset[,1 + 0:(S-1)*ncov_theta]


  data_plot <- apply(samples_subset, 2, function(x) {
    quantile(logistic(x), probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame %>%
    mutate(Species = speciesNames) %>%
    mutate(speciesOrder = order(`2.5%`)) %>%
    filter(Species %in% speciesNames[idx_species])

  orderSpecies <- order(data_plot$`2.5%`)

  plot_collectionrates <- data_plot %>%
    ggplot(aes(x =  factor(Species, level = speciesNames[orderSpecies]),
               ymin = `2.5%`,
               ymax = `97.5%`)) + geom_errorbar() +
    xlab("Species") +
    # ylim(c(0,1)) +
    ggtitle("Collection rates") +
    theme_bw() +
    ylim(c(0,1)) +
    theme(
      axis.text = element_text(angle = 90,
                               size = 8),
      plot.title = element_text(hjust = .5,
                                size = 15)
    )

  plot_collectionrates

}

#' plotDetectionRates
#'
#' True positives rates at the lab stage for each species.
#'
#' @details
#' Plots the 95% credible interval of the true positives at the lab stage
#'
#' @param fitmodel Output from the function runOccPlus
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotDetectionRates(fitmodel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotDetectionRates <- function(fitmodel,
                               idx_species = NULL){

  matrix_of_draws <- fitmodel$matrix_of_draws

  S <- fitmodel$infos$S
  # ncov_theta <- fitmodel$infos$ncov_psi
  speciesNames <- fitmodel$infos$speciesNames
  primerNames <- fitmodel$infos$primerNames

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  param <- "p\\["

  samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]

  data_plot <- apply(samples_subset, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame

  texts <- rownames(data_plot)
  idx_speciesprimer <- str_match(texts, "\\[(\\d+),(\\d+)\\]")

  data_plot <- data_plot %>%
    mutate(Species = as.numeric(idx_speciesprimer[,3]),
           Primer = as.numeric(idx_speciesprimer[,2])) %>%
    mutate(Species = speciesNames[Species],
           Primer = primerNames[Primer]) %>%
    mutate(speciesOrder = order(`2.5%`)) %>%
    filter(Species %in% speciesNames[idx_species])

  orderSpecies <- order(data_plot$`2.5%`[data_plot$Primer == data_plot$Primer[1]])

  detectionRates <- data_plot %>%
    ggplot(aes(x =  factor(Species, level = speciesNames[orderSpecies]),
               ymin = `2.5%`,
               ymax = `97.5%`,
               color = factor(Primer))) + geom_errorbar() +
    xlab("Species") +
    # ylim(c(0,1)) +
    ggtitle("Detection rates") +
    theme_bw() +
    # ylim(c(0,1)) +
    ylab("p") +
    theme(
      axis.text = element_text(angle = 90,
                               size = 8),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = .5,
                                size = 15)
    )

  detectionRates

}

#' plotFPDetectionRates
#'
#' False positives rates at the lab stage for each species.
#'
#' @details
#' Plots the 95% credible interval of the false positives at the lab stage
#'
#' @param fitmodel Output from the function runOccPlus
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotFPDetectionRates(fitmodel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotFPDetectionRates <- function(fitmodel,
                                 idx_species = NULL){

  matrix_of_draws <- fitmodel$matrix_of_draws

  S <- fitmodel$infos$S
  # ncov_theta <- fitmodel$infos$ncov_psi
  speciesNames <- fitmodel$infos$speciesNames
  primerNames <- fitmodel$infos$primerNames

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  param <- "q\\["

  samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]

  data_plot <- apply(samples_subset, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame

  texts <- rownames(data_plot)
  idx_speciesprimer <- str_match(texts, "\\[(\\d+),(\\d+)\\]")

  data_plot <- data_plot %>%
    mutate(Species = as.numeric(idx_speciesprimer[,3]),
           Primer = as.numeric(idx_speciesprimer[,2])) %>%
    mutate(Species = speciesNames[Species],
           Primer = primerNames[Primer]) %>%
    mutate(speciesOrder = order(`2.5%`)) %>%
    filter(Species %in% speciesNames[idx_species])

  orderSpecies <- order(data_plot$`2.5%`[data_plot$Primer == data_plot$Primer[1]])

  detectionRates <- data_plot %>%
    ggplot(aes(x =  factor(Species, level = speciesNames[orderSpecies]),
               ymin = `2.5%`,
               ymax = `97.5%`,
               color = factor(Primer))) + geom_errorbar() +
    xlab("Species") +
    # ylim(c(0,1)) +
    ggtitle("Detection rates") +
    theme_bw() +
    # ylim(c(0,1)) +
    ylab("p") +
    theme(
      axis.text = element_text(angle = 90,
                               size = 8),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = .5,
                                size = 15)
    )

  detectionRates

}

#' plotReadIntensity
#'
#' Plot the reads distribution under the true positives and false positives
#'
#' @details
#' Plots the 95% credible interval of density of true positives and false positives
#'
#' @param fitmodel Output from the function runOccPlus
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotReadIntensity(fitmodel)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotReadIntensity <- function(fitmodel){

  matrix_of_draws <- fitmodel$matrix_of_draws

  mu1_output <-
    matrix_of_draws[,grepl("mu1", colnames(matrix_of_draws))]
  sigma0_output <-
    matrix_of_draws[,grepl("sigma0", colnames(matrix_of_draws))]
  sigma1_output <-
    matrix_of_draws[,grepl("sigma1", colnames(matrix_of_draws))]

  iter <- 1

  niter <- length(mu1_output)

  # x_grid <- seq(1, fitmodel$infos$maxexplogy1, by = 5)
  x_grid <- exp(seq(log(1), log(fitmodel$infos$maxexplogy1), length.out = 1000))

  # seq(1, fitmodel$infos$maxexplogy1, by = 5)

  x <- log(x_grid + 1)

  densities_plot_pos <- matrix(NA, length(x_grid), niter)
  densities_plot_neg <- matrix(NA, length(x_grid), niter)

  for (iter in 1:niter) {
 # print(iter)
    mu1 <- mu1_output[iter]
    sigma1 <- sigma1_output[iter]
    sigma0 <- sigma0_output[iter]

    densities_plot_pos[,iter] <- dnorm(x, mean = mu1, sd = sigma1)
    densities_plot_neg[,iter] <- dnorm(x, mean = 0, sd = sigma0)

    # Create data frame with both densities
    # df <- data.frame(
    #   x = x,
    #   Density0 = ,
    #   Density1 = dnorm(x, mean = mu1, sd = sigma1)
    # )

    # Plot using ggplot2
    # ggplot(df) +
    #   geom_line(aes(x = x, y = Density0, color = "N(0, σ₀²)"), linewidth = 1) +
    #   geom_line(aes(x = x, y = Density1, color = "N(μ₁, σ₁²)"), linewidth = 1) +
    #   scale_color_manual(values = c("N(0, σ₀²)" = "blue", "N(μ₁, σ₁²)" = "red"),
    #                      name = "Distribution") +
    #   labs(title = "Comparison of Two Normal Distributions",
    #        x = "Value",
    #        y = "Density") +
    #   theme_minimal() +
    #   scale_x_continuous(
    #     name = "x",
    #     breaks = x,  # Custom breaks (e.g., e⁻², e⁻¹, e⁰, e¹, ...)
    #     labels = function(x) sprintf("%.2f", exp(x - 1))  # Format labels
    #   )


  }

  df0 <- as_tibble(densities_plot_neg) %>%
    mutate(x = x_grid) %>%
    pivot_longer(cols = -x, names_to = "iter", values_to = "density") %>%
    mutate(Type = "False positives")

  df1 <- as_tibble(densities_plot_pos) %>%
    mutate(x = x_grid) %>%
    pivot_longer(cols = -x, names_to = "iter", values_to = "density") %>%
    mutate(Type = "True positives")

  # Combine into one data frame
  df_combined <- bind_rows(df0, df1) %>%
    mutate(iter = as.numeric(gsub("V", "", iter)))  # Clean iteration labels

  x_grid_breaks <- round( c(0, 10, 20,
                     x_grid[seq(1, length(x_grid), by = 20)] - 1))

  ggplot(df_combined, aes(x = x, y = density, group = interaction(iter, Type), color = Type)) +
    geom_line(alpha = 0.3, aes(color = Type)) +
    scale_color_manual(values = c("False positives" = "blue", "True positives" = "red")) +
    labs(title = "Reads distributions",
         x = "x", y = "Density") +
    theme_minimal() +
    scale_x_continuous(
      name = "x",
      breaks = x_grid_breaks,  # Custom breaks (e.g., e⁻², e⁻¹, e⁰, e¹, ...)
      # labels = function(x) sprintf("%.2f", exp(x) - 1)  # Format labels
      trans = "log"
    ) +
    theme(axis.text.x = element_text(angle = 90))

}

outputCovarianceMatrix <- function(fitmodel,
                                   idx_species = NULL){

  matrix_of_draws <- fitmodel$matrix_of_draws

  L_output <-
    matrix_of_draws[,grepl("LL\\[", colnames(matrix_of_draws))]
  S <- fitmodel$infos$S

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  niter <- nrow(L_output)

  Lambda_output <- array(NA, dim = c(niter, S, S))

  for (iter in 1:niter) {
    L_output_current <- matrix(L_output[iter,], S, d, byrow = T)
    Lambda_output[iter,,] <- cov2cor(L_output_current %*% t(L_output_current))
  }

  Lambda_output <- Lambda_output[,idx_species, idx_species]

  # Lambda_mean <- apply(Lambda_output, c(2,3), mean)

  Lambda_quantiles <- apply(Lambda_output, c(2,3),
                            function(x){quantile(x, probs = c(0.025, 0.5, 0.975))})



  ggcorrplot(Lambda_quantiles[1,,], method = "square", type = "lower",
             lab = F, lab_size = 3,
             colors = c("blue", "white", "red"),
             title = "Covariance Matrix (as Correlation)")

  ggcorrplot(Lambda_quantiles[3,,], method = "square", type = "lower",
             lab = F, lab_size = 3,
             colors = c("blue", "white", "red"),
             title = "Covariance Matrix (as Correlation)")


  Lambda_CI_lower <- Lambda_quantiles[1,,]
  Lambda_CI_lower <- Lambda_CI_lower[lower.tri(Lambda_CI_lower)]

  Lambda_CI_upper <- Lambda_quantiles[3,,]
  Lambda_CI_upper <- Lambda_CI_upper[lower.tri(Lambda_CI_upper)]

  Lambda_true <- Lambda_true[lower.tri(Lambda_true)]

  data.frame(x = 1:length(Lambda_CI_lower),
             y = Lambda_true,
             ymin = Lambda_CI_lower,
             ymax = Lambda_CI_upper) %>%
    ggplot(aes(x = x, y = y, ymin = ymin, ymax = ymax)) +
    geom_point() + geom_errorbar()

  ggcorrplot(Lambda_mean, method = "square", type = "lower",
             lab = F, lab_size = 3,
             colors = c("blue", "white", "red"),
             title = "Covariance Matrix (as Correlation)")



  Lambda_true <- t(L_true) %*% L_true

  Lambda_true <- cov2cor(Lambda_true)
  #
  # ggcorrplot(Lambda_true, method = "square", type = "lower",
  #            lab = F, lab_size = 3,
  #            colors = c("blue", "white", "red"),
  #            title = "Covariance Matrix (as Correlation)")

  list("Median" = Lambda_quantiles[2,,],
       "Lower" = Lambda_quantiles[1,,],
       "Upper" = Lambda_quantiles[3,,])

}

plotCovarianceMatrix <- function(fitmodel,
                                 idx_species = NULL){

  matrix_of_draws <- fitmodel$matrix_of_draws

  L_output <-
    matrix_of_draws[,grepl("LL\\[", colnames(matrix_of_draws))]
  S <- fitmodel$infos$S

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  niter <- nrow(L_output)

  Lambda_output <- array(NA, dim = c(niter, S, S))

  for (iter in 1:niter) {
    L_output_current <- matrix(L_output[iter,], S, d, byrow = T)
    Lambda_output[iter,,] <- cov2cor(L_output_current %*% t(L_output_current))
  }

  Lambda_output <- Lambda_output[,idx_species, idx_species]

  Lambda_quantiles <- apply(Lambda_output, c(2,3),
                            function(x){quantile(x, probs = c(0.025, 0.5, 0.975))})

  ggcorrplot(Lambda_quantiles[2,,], method = "square", type = "lower",
             lab = F, lab_size = 3,
             colors = c("blue", "white", "red"),
             title = "Covariance Matrix (as Correlation)")


}

# to write
plotFactorScores <- function(fitmodel,
                             idx_factor = c(1,2)){

  d <- fitmodel$d

  if(d == 0){

  }

}
