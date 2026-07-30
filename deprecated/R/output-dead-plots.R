# Plotting functions that can no longer work, moved out of R/output.R on
# 30 July 2026.
#
# TODO.md group D item 3, "DEPRECATE BOTH". Neither was reachable: no callers
# in R/, tests/ or vignettes/, absent from NAMESPACE, no string references.
#
#   plotReadIntensity    reads results_output$mu1_output / mu0_output /
#                        sigma1_output / sigma0_output and infos$maxexplogy1,
#                        none of which runOccJSDM() produces any more. It was
#                        de-exported earlier rather than fixed.
#   plotOccupancyStates  references undefined data_info and OTU.
#
# Their roxygen blocks moved with them, deliberately: leaving a block behind
# would re-attach it to whatever function followed. This directory is excluded
# from the build by .Rbuildignore, so nothing here ships.

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
