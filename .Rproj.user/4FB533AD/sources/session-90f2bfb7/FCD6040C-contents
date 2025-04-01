# library(here); library(stringr);
# library(tidyverse); library(mixtools)
#
# library("rstan")
# options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = TRUE)
# rm(list = ls())

# LOAD DATA ------

load(here("Stan/Count model","wangdata.rda"))

data$OTU[data$info$elution_volume == 0,] <- NA

# fit_ednaplus <- runOccPlusPro(data, d = 1,
#                               ordCovariates = c("altitude"),
#                               detCovariates = c("input_volume"))

# SIMULATE DATA ------

# param count distr
{
  lambda0 <- exp(2.5)
  lambda1 <- exp(5)

  size0 <- 3
  size1 <- 3

  pi0 <- .95

}

r_f0y <- function(n, pi0, lambda0, size0){
  rbinom(n, size = 1, prob = 1 - pi0) * rnbinom(n, mu = lambda0, size = size0)
}

r_f1y <- function(n, lambda1, size1){
  rnbinom(n, size = size1, mu = lambda1)
}

# settings
{
  n <- 20
  S <- 5
  M <- rep(2, n)
  N <- sum(M)
  L <- rep(2, N)
  N2 <- sum(L)
  K <- rep(2, N2)
  N3 <- sum(K)

  sumM <- c(0, cumsum(M)[-n])
  sumL <- c(0, cumsum(L)[-N])
  sumK <- c(0, cumsum(K)[-N2])

  rep(rep(rep(1:n, M), L), K)

  ncov_psi <- 2
  ncov_theta <- 2
  ncov_ord <- 2

  X_psi <- matrix(rnorm(n * ncov_psi), n, ncov_psi)
  X_theta <- cbind(1, matrix(rnorm(N * (ncov_theta - 1)), N, ncov_theta - 1))
  X_ord <- matrix(rnorm(n * ncov_ord), n, ncov_ord)

  # na samples
  idx_nasample <- sample(1:N, N * .01, replace = T)

}

d <- 2

logistic <- function(x) 1 / (1 + exp(-x))

beta0_psi_true <- rnorm(S, sd = .1)
beta_psi_true <- matrix(sample(c(-1,1), ncov_psi * S, replace = T), ncov_psi, S)
# beta_psi_true <- matrix(rnorm(ncov_psi * S), ncov_psi, S)
# beta_theta_true <- matrix(rnorm(ncov_theta * S), ncov_theta, S)
beta_theta_true <- rbind(1,
                         matrix(sample(c(-1,1), (ncov_theta - 1) * S, replace = T),
                                (ncov_theta - 1), S))
beta_ord_true <- matrix(sample(c(-1,1,0), ncov_ord * d, replace = T), ncov_ord, d)

E_true <- matrix(rnorm(n * d, sd = .5), n, d)
U_true <- X_ord %*% beta_ord_true + E_true
L_true <- matrix(rnorm(d * S), d, S)

psi_true <- logistic(matrix(beta0_psi_true, n, S, byrow = T) +
                       X_psi %*% beta_psi_true + U_true %*% L_true)
theta_true <- logistic(X_theta %*% beta_theta_true)
theta0_true <- rep(0.01, S)
p_true <- matrix(.9, max(L), S)
q_true <- matrix(.02, max(L), S)

z <- matrix(0, n, S)
delta <- matrix(0, N, S)
cimk <- matrix(0, N3, S)
y <- matrix(NA, N3, S)
for (s in 1:S) {

  idx_sample <- 1

  for (i in 1:n) {

    z[i,s] <- rbinom(1, 1, prob = psi_true[i,s])

    if(z[i,s] == 1){
      for (m in 1:M[i]) {
        delta[sumM[i] + m, s] <- rbinom(1, 1, theta_true[sumM[i] + m, s])
      }
    } else {
      for (m in 1:M[i]) {
        delta[sumM[i] + m, s] <- rbinom(1, 1, theta0_true[s])
      }
    }

    for (m in 1:M[i]) {
      if(delta[sumM[i] + m, s] == 1){
        for(l in 1:L[sumM[i] + m]){
          for(k in 1:K[sumL[sumM[i] + m] + l]){
            cimk[sumL[sumM[i] + m] + l, s] <- rbinom(1, 1, p_true[l, s])
          }
        }
      } else {
        for(l in 1:L[sumM[i] + m]){
          for(k in 1:K[sumL[sumM[i] + m] + l]){
            cimk[sumL[sumM[i] + m] + l, s] <- rbinom(1, 1, q_true[l, s])
          }
        }
      }
    }

    for (m in 1:M[i]) {
      for(l in 1:L[sumM[i] + m]){
        for(k in 1:K[sumL[sumM[i] + m] + l]){

          if(!(idx_sample %in% idx_nasample)){
            if(cimk[sumL[sumM[i] + m] + l, s] == 1){
              y[sumK[sumL[sumM[i] + m] + l] + k, s] <-
                r_f1y(1, lambda1, size = size1)
            } else {
              y[sumK[sumL[sumM[i] + m] + l] + k, s] <-
                r_f0y(1, pi0, lambda0, size0 = size0)
            }
          } else {
            y[sumK[sumL[sumM[i] + m] + l] + k, s] <-
              NA
          }

        }
      }
      idx_sample <- idx_sample + 1
    }

  }

}

OTU <- y

result <- unlist(mapply(function(k, x) rep(x, k), K, lapply(K, function(x) 1:max(L))))

data_info <- data.frame(Site = rep(rep(rep(1:n, M), L), K),
                        Sample = rep(rep(1:N, L), K),
                        Primer = rep(  rep(1:max(L), each = K[1]), times = N),
                        X_psi = X_psi[rep(rep(rep(1:n, M), L), K),],
                        X_ord = X_ord[rep(rep(rep(1:n, M), L), K),],
                        X_theta = X_theta[rep(rep(1:N, L), K),-1])

data <- list(info = data_info,
             OTU = OTU)

# SIMULATE BASED ON REAL ------

# param count distr
{
  lambda0 <- exp(1)
  lambda1 <- exp(5)

  size0 <- 2
  size1 <- 2

  pi0 <- .95

}

r_f0y <- function(n, pi0, lambda0, size0){
  rbinom(n, size = 1, prob = 1 - pi0) * rnbinom(n, mu = lambda0, size = size0)
}

r_f1y <- function(n, lambda1, size1){
  rnbinom(n, size = size1, mu = lambda1)
}

# settings
{
  n <- 200
  S <- 50
  M <- rep(3, n)
  N <- sum(M)
  L <- rep(2, N)
  N2 <- sum(L)
  K <- rep(2, N2)
  N3 <- sum(K)

  sumM <- c(0, cumsum(M)[-n])
  sumL <- c(0, cumsum(L)[-N])
  sumK <- c(0, cumsum(K)[-N2])

  ncov_psi <- 0
  ncov_theta <- 1
  ncov_ord <- 2

  X_psi <- matrix(rnorm(n * ncov_psi), n, ncov_psi)
  X_theta <- cbind(1, matrix(rnorm(N * (ncov_theta - 1)), N, ncov_theta - 1))
  X_ord <- matrix(rnorm(n * ncov_ord), n, ncov_ord)
}

d <- 1

logistic <- function(x) 1 / (1 + exp(-x))

beta_psi_true <- matrix(sample(c(-1,1), ncov_psi * S, replace = T), ncov_psi, S)
# beta_psi_true <- matrix(rnorm(ncov_psi * S), ncov_psi, S)
# beta_theta_true <- matrix(rnorm(ncov_theta * S), ncov_theta, S)
beta_theta_true <- matrix(-1, ncov_theta, S)
beta_ord_true <- matrix(sample(c(-1,1), ncov_ord * d, replace = T), ncov_ord, d)

E_true <- matrix(rnorm(n * d, sd = .025), n, d)
U_true <- X_ord %*% beta_ord_true + E_true
L_true <- matrix(rnorm(d * S), d, S)

# psi_true <- logistic(X_psi %*% beta_psi_true)
psi_true <- logistic(X_psi %*% beta_psi_true + U_true %*% L_true)
theta_true <- logistic(X_theta %*% beta_theta_true)
p_true <- matrix(.3, max(L), S)

z <- matrix(0, n, S)
delta <- matrix(0, N, S)
cimk <- matrix(0, N3, S)
y <- matrix(NA, N3, S)
for (s in 1:S) {
  for (i in 1:n) {
    z[i,s] <- rbinom(1, 1, prob = psi_true[i,s])
    if(z[i,s] == 1){
      for (m in 1:M[i]) {
        delta[sumM[i] + m, s] <- rbinom(1, 1, theta_true[sumM[i] + m, s])
        if(delta[sumM[i] + m, s] == 1){
          for(l in 1:L[sumM[i] + m]){

            for(k in 1:K[sumL[sumM[i] + m] + l]){
              cimk[sumL[sumM[i] + m] + l, s] <- rbinom(1, 1, p_true[l, s])
              if(cimk[sumL[sumM[i] + m] + l, s] == 1){

                y[sumK[sumL[sumM[i] + m] + l] + k, s] <-
                  r_f1y(1, lambda1, size = size1)

              } else {

                y[sumK[sumL[sumM[i] + m] + l] + k, s] <-
                  r_f0y(1, pi0, lambda0, size0 = size0)

              }

            }
          }
        } else {

          for(l in 1:L[sumM[i] + m]){
            for(k in 1:K[sumL[sumM[i] + m] + l]){
              y[sumK[sumL[sumM[i] + m] + l] + k, s] <-
                r_f0y(1, pi0, lambda0, size0 = size0)
            }

          }

        }
      }
    } else {

      for (m in 1:M[i]) {
        for(l in 1:L[sumM[i] + m]){
          for(k in 1:K[sumL[sumM[i] + m] + l]){
            y[sumK[sumL[sumM[i] + m] + l] + k, s] <-
              r_f0y(1, pi0, lambda0, size0 = size0)
          }

        }
      }

    }
  }
}

logy1 <- log(y + 1)

# MODEL -----

fitmodel  <- runOccPlus(data,
                        d = 2,
                        occCovariates = c("X_psi.1","X_psi.2"),
                        ordCovariates = c("X_ord.1","X_ord.2"),
                        detCovariates = c("X_theta"))

# OUTPUT USING FUNCTIONS ---------

plotOccupancyCovariates(fitmodel,
                        covName = "X_psi.1")

plotDetectionCovariates(fitmodel,
                        covName = "X_theta")

plotOrdinationCovariates(fitmodel,
                        covName = "X_ord.2")

plotOccupancyRates(fitmodel,
                    idx_species = 1:5)

plotCollectionRates(fitmodel,
                    idx_species = 1:5)

plotDetectionRates(fitmodel,
                   idx_species = 1:5)

plotFPDetectionRates(fitmodel,
                     idx_species = 1:5)

plotReadIntensity(fitmodel)

# EXTRACT AND REPARAM OUTPUT ------

matrix_of_draws <- fit_ednaplus$matrix_of_draws

U_output <-
  matrix_of_draws[,grepl("U\\[", colnames(matrix_of_draws))]
L_output <-
  matrix_of_draws[,grepl("LL\\[", colnames(matrix_of_draws))]
E_output <-
  matrix_of_draws[,grepl("E\\[", colnames(matrix_of_draws))]
beta_psi_output <-
  matrix_of_draws[,grepl("beta_psi\\[", colnames(matrix_of_draws))]
beta0_psi_output <-
  matrix_of_draws[,grepl("beta0_psi\\[", colnames(matrix_of_draws))]
beta_ord_output <-
  matrix_of_draws[,grepl("beta_ord\\[", colnames(matrix_of_draws)),drop=F]
beta_theta_output <-
  matrix_of_draws[,grepl("beta_theta\\[", colnames(matrix_of_draws))]
theta0_output <-
  matrix_of_draws[,grepl("theta0\\[", colnames(matrix_of_draws))]
p_output <-
  matrix_of_draws[,grepl("p\\[", colnames(matrix_of_draws))]
q_output <-
  matrix_of_draws[,grepl("q\\[", colnames(matrix_of_draws))]
mu1_output <-
  matrix_of_draws[,grepl("mu1", colnames(matrix_of_draws))]
pi0_output <-
  matrix_of_draws[,grepl("pi0", colnames(matrix_of_draws))]
sigma0_output <-
  matrix_of_draws[,grepl("sigma0", colnames(matrix_of_draws))]
sigma1_output <-
  matrix_of_draws[,grepl("sigma1", colnames(matrix_of_draws))]

# reparam the true
{
  L1_true <- L_true[1]
  L_true_reparam <- L_true / L1_true
  E_true_reparam <- E_true * L1_true
  U_true_reparam <- U_true * L1_true
  beta_ord_true_reparam <- beta_ord_true * L1_true
}

niter <- nrow(L_output)

L_output_reparam <- L_output
U_output_reparam <- U_output
E_output_reparam <- E_output
beta_ord_output_reparam <- beta_ord_output

for (iter in 1:niter) {
  print(iter)

  L1 <- L_output[iter,1]
  L_output_reparam[iter,] <- L_output[iter,] / L1
  E_output_reparam[iter,] <- E_output[iter,] * L1
  U_output_reparam[iter,] <- U_output[iter,] * L1
  beta_ord_output_reparam[iter,] <- beta_ord_output[iter,] * L1

}

plotParam(beta_ord_output_reparam, beta_ord_true_reparam, c(-3,3))
#
# iter <- 2
# term1 <- (X_ord %*% beta_ord_output_reparam[iter,] + E_output_reparam[iter,]) %*% L_output_reparam[iter,]
# term2 <- (X_ord %*% beta_ord_output[iter,] + E_output[iter,]) %*% L_output[iter,]
# apply(beta_ord_output_reparam, 2, function(x){
#   quantile(x, probs = c(0.025, 0.975))
# })
#
#
#
# qplot(X_ord[,4], apply(U_output_reparam, 2, mean)) + geom_smooth(method = "lm")
#

# OUTPUT SIMS --------

plotParam <- function(samples_subset, trueParams_mat, lims){

  params_CI <- apply(samples_subset, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  })

  trueParams <- as.vector(trueParams_mat)

  data_plot <- cbind(
    t(params_CI), trueParams
  ) %>% as.data.frame

  trueVars <- (data_plot$trueParams < data_plot$`97.5%` &
                 data_plot$trueParams > data_plot$`2.5%`)

  data_plot$trueVars <- trueVars

  data_plot$x <- rownames(data_plot)

  ggplot(data_plot, aes(x = x,
                        y = trueParams,
                        ymin = `2.5%`,
                        ymax = `97.5%`)) + geom_errorbar() +
    geom_point() +
    theme(
      axis.text = element_text(angle = 45)
    ) + ylim(lims)

}

# p
{
  plotParam(p_output, p_true, c(0,1))
}

# q
{
  plotParam(q_output, q_true, c(0,1))
}

# theta0
{
  plotParam(theta0_output, theta0_true, c(0,1))
}

# beta_theta0
{
  plotParam(beta_theta_output, beta_theta_true, c(-2,2))
}

# beta_ord
{
  plotParam(beta_ord_output, beta_ord_true, c(-3,3))
  plotParam(beta_ord_output_reparam, beta_ord_true_reparam, c(-2.5,2.5))
}

# beta_psi
{
  plotParam(beta0_psi_output, beta0_psi_true, c(-4,4))
}

# beta_psi
{
  plotParam(beta_psi_output, beta_psi_true, c(-3,3))
}

# U
{
  plotParam(U_output_reparam, U_true_reparam, c(-1,1))
}

# L
{
  # plotParam(L_output, L_true, c(-2,2))
  plotParam(L_output_reparam, L_true_reparam, c(-20,20))
}

quantile(sigma0_output, probs = c(0.025, 0.975))
log(lambda0)

quantile(mu1_output, probs = c(0.025, 0.975))
log(lambda1)

quantile(mu1m_mu0_output + mu0_output, probs = c(0.025, 0.975))
log(lambda0 + lambda1)

quantile(pi0_output, probs = c(0.025, 0.975))
pi0

# OUTPUT SIMS BETA PSI --------

library(tidyverse)

param <- "beta_psi"

samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]

params_CI <- apply(samples_subset, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
})

trueParams <- as.vector(beta_psi_true)

data_plot <- cbind(
  t(params_CI), trueParams
) %>% as.data.frame

trueVars <- (data_plot$trueParams < data_plot$`97.5%` &
               data_plot$trueParams > data_plot$`2.5%`)

data_plot$trueVars <- trueVars

data_plot$x <- rownames(data_plot)

ggplot(data_plot, aes(x = x,
                      y = trueParams,
                      ymin = `2.5%`,
                      ymax = `97.5%`)) + geom_errorbar() +
  geom_point() +
  theme(
    axis.text = element_text(angle = 45)
  ) + ylim(c(-2,2))

# OUTPUT SIMS BETA THETA --------

library(tidyverse)

param <- "beta_theta"

samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]

params_CI <- apply(samples_subset, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
})

trueParams <- as.vector(beta_theta_true)

data_plot <- cbind(
  t(params_CI), trueParams
) %>% as.data.frame

trueVars <- (data_plot$trueParams < data_plot$`97.5%` &
               data_plot$trueParams > data_plot$`2.5%`)

data_plot$trueVars <- trueVars

data_plot$x <- rownames(data_plot)

ggplot(data_plot, aes(x = x,
                      y = trueParams,
                      ymin = `2.5%`,
                      ymax = `97.5%`)) + geom_errorbar() +
  geom_point() +
  theme(
    axis.text = element_text(angle = 45)
  ) + ylim(c(-2,2))

# OUTPUT SIMS BETA ORD --------

library(tidyverse)

param <- "beta_ord"

samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws)), drop= F]

params_CI <- apply(samples_subset, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
})

trueParams <- as.vector(beta_ord_true)

data_plot <- cbind(
  t(params_CI), trueParams
) %>% as.data.frame

trueVars <- (data_plot$trueParams < data_plot$`97.5%` &
               data_plot$trueParams > data_plot$`2.5%`)

data_plot$trueVars <- trueVars

data_plot$x <- rownames(data_plot)

ggplot(data_plot, aes(x = x,
                      y = trueParams,
                      ymin = `2.5%`,
                      ymax = `97.5%`)) + geom_errorbar() +
  geom_point() +
  theme(
    axis.text = element_text(angle = 45)
  ) + ylim(c(-2,2))

# OUTPUT - COVARIATE ------

library(tidyverse)

param <- "beta_psi"

samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]

data_plot <- apply(samples_subset, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
}) %>% t %>% as.data.frame


covNames <- colnames(X_psi)

# species names
{
  speciesNames <- colnames(y)
  speciesNames <- gsub("Species_","",speciesNames)
}


data_plot$x <- rownames(data_plot)

data_plot$cov <- rep(covNames, S)
data_plot$species <- rep(speciesNames, each = ncov_psi)

for (cov in 2:ncov_psi) {

  # speciesOrder <-

  # speciesSubset <- speciesNames
  # speciesOrder <- order(data_plot$)

  data_plot_cov <- data_plot %>%
    filter(str_detect(cov, (!!covNames[cov])))

  orderSpecies <- order(data_plot_cov$`2.5%`)

  plot1 <- data_plot_cov %>%
    filter(species %in% speciesNames[orderSpecies[1:72]]) %>%
    # filter(str_detect(cov, (!!covNames[cov]))) %>%
    # arrange(`2.5%`) %>%
    mutate(speciesOrder = order(`2.5%`)) %>%
    ggplot(aes(x = factor(species, level = species[speciesOrder]),
               ymin = `2.5%`,
               ymax = `97.5%`)) + geom_errorbar() +
    # facet_grid(rows = vars(cov)) +
    geom_hline(aes(yintercept = 0), color = "red") +
    ggtitle(covNames[cov]) +
    xlab("Species") +
    ylim(c(-3,3)) +
    theme_bw() +
    theme(
      axis.text = element_text(angle = 90,
                               size = 8),
      plot.title = element_text(hjust = .5,
                                size = 15)
    )

  plot2 <- data_plot_cov %>%
    filter(species %in% speciesNames[orderSpecies[73:144]]) %>%
    # filter(str_detect(cov, (!!covNames[cov]))) %>%
    # arrange(`2.5%`) %>%
    mutate(speciesOrder = order(`2.5%`)) %>%
    ggplot(aes(x = factor(species, level = species[speciesOrder]),
               ymin = `2.5%`,
               ymax = `97.5%`)) + geom_errorbar() +
    # facet_grid(rows = vars(cov)) +
    geom_hline(aes(yintercept = 0), color = "red") +
    ggtitle(covNames[cov]) +
    xlab("Species") +
    ylim(c(-3,3)) +
    theme_bw() +
    theme(
      axis.text = element_text(angle = 90,
                               size = 8),
      plot.title = element_text(hjust = .5,
                                size = 15)
    )

  plot3 <- data_plot_cov %>%
    filter(species %in% speciesNames[orderSpecies[145:218]]) %>%
    # filter(str_detect(cov, (!!covNames[cov]))) %>%
    # arrange(`2.5%`) %>%
    mutate(speciesOrder = order(`2.5%`)) %>%
    ggplot(aes(x = factor(species, level = species[speciesOrder]),
               ymin = `2.5%`,
               ymax = `97.5%`)) + geom_errorbar() +
    # facet_grid(rows = vars(cov)) +
    geom_hline(aes(yintercept = 0), color = "red") +
    ggtitle(covNames[cov]) +
    xlab("Species") +
    ylim(c(-3,3)) +
    theme_bw() +
    theme(
      axis.text = element_text(angle = 90,
                               size = 8),
      plot.title = element_text(hjust = .5,
                                size = 15)
    )

  setwd(here("Code","Results"))

  ggsave(plot1, file = paste0(covNames[cov], "1.jpeg"))
  ggsave(plot2, file = paste0(covNames[cov], "2.jpeg"))
  ggsave(plot3, file = paste0(covNames[cov], "3.jpeg"))


}

# OUTPUT - COVARIATE ORD ------

library(tidyverse)

param <- "beta_ord"

samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]
# samples_subset <- beta_ord_output_reparam#matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]

data_plot <- apply(samples_subset, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
}) %>% t %>% as.data.frame

covNames <- colnames(X_ord)

# species names
{
  # speciesNames <- colnames(y)
  # speciesNames <- gsub("Species_","",speciesNames)
}


data_plot$x <- rownames(data_plot)

data_plot$cov <- rep(covNames, d)
data_plot$factor <- rep(1, each = ncov_ord)


plotCov <- data_plot %>%
  ggplot(aes(x = cov,
             ymin = `2.5%`,
             ymax = `97.5%`)) + geom_errorbar() +
  # facet_grid(rows = vars(cov)) +
  geom_hline(aes(yintercept = 0), color = "red") +
  # ggtitle(covNames[cov]) +
  xlab("Covariates") +
  # ylim(c(-,1)) +
  theme_bw() +
  theme(
    axis.text = element_text(angle = 90,
                             size = 8),
    plot.title = element_text(hjust = .5,
                              size = 15)
  )

setwd(here("Code Counts","Results"))
ggsave("CovariatesEffect.jpeg", plotCov)

for (cov in 1:ncov_ord) {

  # speciesOrder <-

  # speciesSubset <- speciesNames
  # speciesOrder <- order(data_plot$)

  data_plot_cov <- data_plot %>%
    filter(str_detect(cov, (!!covNames[cov])))

  orderSpecies <- order(data_plot_cov$`2.5%`)

  plot1 <- data_plot_cov %>%
    filter(species %in% speciesNames[orderSpecies[1:72]]) %>%
    # filter(str_detect(cov, (!!covNames[cov]))) %>%
    # arrange(`2.5%`) %>%
    mutate(speciesOrder = order(`2.5%`)) %>%
    ggplot(aes(x = factor(species, level = species[speciesOrder]),
               ymin = `2.5%`,
               ymax = `97.5%`)) + geom_errorbar() +
    # facet_grid(rows = vars(cov)) +
    geom_hline(aes(yintercept = 0), color = "red") +
    ggtitle(covNames[cov]) +
    xlab("Species") +
    ylim(c(-3,3)) +
    theme_bw() +
    theme(
      axis.text = element_text(angle = 90,
                               size = 8),
      plot.title = element_text(hjust = .5,
                                size = 15)
    )

  plot2 <- data_plot_cov %>%
    filter(species %in% speciesNames[orderSpecies[73:144]]) %>%
    # filter(str_detect(cov, (!!covNames[cov]))) %>%
    # arrange(`2.5%`) %>%
    mutate(speciesOrder = order(`2.5%`)) %>%
    ggplot(aes(x = factor(species, level = species[speciesOrder]),
               ymin = `2.5%`,
               ymax = `97.5%`)) + geom_errorbar() +
    # facet_grid(rows = vars(cov)) +
    geom_hline(aes(yintercept = 0), color = "red") +
    ggtitle(covNames[cov]) +
    xlab("Species") +
    ylim(c(-3,3)) +
    theme_bw() +
    theme(
      axis.text = element_text(angle = 90,
                               size = 8),
      plot.title = element_text(hjust = .5,
                                size = 15)
    )

  plot3 <- data_plot_cov %>%
    filter(species %in% speciesNames[orderSpecies[145:218]]) %>%
    # filter(str_detect(cov, (!!covNames[cov]))) %>%
    # arrange(`2.5%`) %>%
    mutate(speciesOrder = order(`2.5%`)) %>%
    ggplot(aes(x = factor(species, level = species[speciesOrder]),
               ymin = `2.5%`,
               ymax = `97.5%`)) + geom_errorbar() +
    # facet_grid(rows = vars(cov)) +
    geom_hline(aes(yintercept = 0), color = "red") +
    ggtitle(covNames[cov]) +
    xlab("Species") +
    ylim(c(-3,3)) +
    theme_bw() +
    theme(
      axis.text = element_text(angle = 90,
                               size = 8),
      plot.title = element_text(hjust = .5,
                                size = 15)
    )

  setwd(here("Code","Results"))

  ggsave(plot1, file = paste0(covNames[cov], "1.jpeg"))
  ggsave(plot2, file = paste0(covNames[cov], "2.jpeg"))
  ggsave(plot3, file = paste0(covNames[cov], "3.jpeg"))


}

# OUTPUT - READ FP RATES ------

samples_subset_mu0 <- matrix_of_draws[,grepl("mu0", colnames(matrix_of_draws)) &
                                        !grepl("mu1", colnames(matrix_of_draws))]

samples_subset_mu1 <- matrix_of_draws[,grepl("mu1_m", colnames(matrix_of_draws))]

samples_subset_sigma0 <- matrix_of_draws[,grepl("sigma0", colnames(matrix_of_draws))]
samples_subset_sigma1 <- matrix_of_draws[,grepl("sigma1", colnames(matrix_of_draws))]

samples_subset_pi0 <- matrix_of_draws[,grepl("pi0", colnames(matrix_of_draws))]

mean(samples_subset_mu0)
mean(samples_subset_mu1)

data_plot <- apply(samples_subset, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
}) %>% t %>% as.data.frame

#

samples_subset <- matrix_of_draws[,grepl("sigma0", colnames(matrix_of_draws))]

data_plot <- apply(samples_subset, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
}) %>% t %>% as.data.frame



samples_subset <- matrix_of_draws[,grepl("sigma1", colnames(matrix_of_draws))]

data_plot <- apply(samples_subset, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
}) %>% t %>% as.data.frame




# species names
{
  speciesNames <- colnames(y)
  speciesNames <- gsub("Species_","",speciesNames)
}

data_plot$x <- rownames(data_plot)

data_plot %>%
  ggplot(aes(x = x,
             ymin = `2.5%`,
             ymax = `97.5%`)) + geom_errorbar() +
  # facet_grid(rows = vars(cov)) +
  # geom_hline(aes(yintercept = 0), color = "red") +
  # ggtitle(covNames[cov]) +
  # xlab("Covariates") +
  # ylim(c(-,1)) +
  theme_bw() +
  theme(
    axis.text = element_text(angle = 90,
                             size = 8),
    plot.title = element_text(hjust = .5,
                              size = 15)
  )


# OUTPUT - READ FP SIGMA RATES ------

library(tidyverse)

param <- "sigma0"

samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]

data_plot <- apply(samples_subset, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
}) %>% t %>% as.data.frame


# species names
{
  # speciesNames <- colnames(y)
  # speciesNames <- gsub("Species_","",speciesNames)
}

data_plot$x <- rownames(data_plot)

data_plot %>%
  ggplot(aes(x = x,
             ymin = `2.5%`,
             ymax = `97.5%`)) + geom_errorbar() +
  # facet_grid(rows = vars(cov)) +
  # geom_hline(aes(yintercept = 0), color = "red") +
  # ggtitle(covNames[cov]) +
  # xlab("Covariates") +
  # ylim(c(-,1)) +
  theme_bw() +
  theme(
    axis.text = element_text(angle = 90,
                             size = 8),
    plot.title = element_text(hjust = .5,
                              size = 15)
  )


# OUTPUT - COUNT RATES ------

samples_mu1_subset <- matrix_of_draws[,grepl("mu1", colnames(matrix_of_draws))]
samples_sigma1_subset <- matrix_of_draws[,grepl("sigma1", colnames(matrix_of_draws))]
samples_sigma0_subset <- matrix_of_draws[,grepl("sigma0", colnames(matrix_of_draws))]
samples_pi0_subset <- matrix_of_draws[,grepl("pi0", colnames(matrix_of_draws))]

niter <- length(samples_sigma1_subset)

x_grid <- seq(0, 10, length.out = 200)

dens1_val <- matrix(NA, niter, 200)
dens2_val <- matrix(NA, niter, 200)

for (iter in 1:niter) {

  mu1_current <- samples_mu1_subset[iter]
  sigma1_current <- samples_sigma1_subset[iter]
  sigma0_current <- samples_sigma0_subset[iter]
  pi0_current <- samples_pi0_subset[iter]

  fdens1 <- function(x){
    dnorm(x, 0, sigma0_current) #+
    # .5 * dnorm(x, mean_mu0 + mean_mu1, mean_sigma1)
  }

  fdens2 <- function(x){
    # dnorm(x, mean_mu0, mean_sigma0) #+
    dnorm(x, mu1_current, sigma1_current)
  }

  dens1_val[iter,] <- fdens1(x_grid)
  dens2_val[iter,] <- fdens2(x_grid)

}

dens1_mean <- apply(dens1_val, 2, function(x){
  quantile(x, probs = c(0.025, 0.5, 0.975))
})
dens2_mean <- apply(dens2_val, 2, function(x){
  quantile(x, probs = c(0.025, 0.5, 0.975))
})

plotrates <- ggplot() +
  geom_line(data = NULL, aes(x = x_grid,
                             y = dens1_mean[2,]),
            color = "blue", size = 2) +
  geom_line(data = NULL, aes(x = x_grid,
                             y = dens2_mean[2,]),
            color = "red", fill = "red") +
  geom_ribbon(data = NULL, aes(x = x_grid,
                               ymin = dens1_mean[1,],
                               ymax = dens1_mean[3,]),
              color = "blue", fill = "blue", alpha = .3) +
  geom_ribbon(data = NULL, aes(x = x_grid,
                               ymin = dens2_mean[1,],
                               ymax = dens2_mean[3,]),
              color = "red", fill = "red", alpha = 1) +
  scale_x_continuous(breaks = seq(0, 10, length.out = 10),
                     labels = function(x) round(exp(x)),
                     limits = c(0, 10)) + theme_bw() + xlab("Reads") + ylab("")


setwd(here("Code Counts","Results"))
ggsave(plotrates, file = "CountRates.jpeg")


# OUTPUT - READ TP RATES ------

param <- "mu1"

samples_mu1_subset <- matrix_of_draws[,grepl("mu1", colnames(matrix_of_draws))]

samples_subset <- samples_mu1_subset

data_plot <- apply(samples_subset, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
}) %>% t %>% as.data.frame

# species names
{
  # speciesNames <- colnames(y)
  # speciesNames <- gsub("Species_","",speciesNames)
}

data_plot$x <- rownames(data_plot)

data_plot %>%
  ggplot(aes(x = x,
             ymin = `2.5%`,
             ymax = `97.5%`)) + geom_errorbar() +
  # facet_grid(rows = vars(cov)) +
  # geom_hline(aes(yintercept = 0), color = "red") +
  # ggtitle(covNames[cov]) +
  # xlab("Covariates") +
  # ylim(c(-,1)) +
  theme_bw() +
  theme(
    axis.text = element_text(angle = 90,
                             size = 8),
    plot.title = element_text(hjust = .5,
                              size = 15)
  )


# OUTPUT - OCCUPANCY RATES ------

param <- "beta_psi"

logistic <- function(x) 1 / (1 + exp(-x))

samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]

samples_subset <- samples_subset[,1 + 0:(S-1)*ncov_psi]

{

  plot1 <- data_plot_cov %>%
    filter(species %in% speciesNames[orderSpecies[1:72]]) %>%
    # filter(str_detect(cov, (!!covNames[cov]))) %>%
    # arrange(`2.5%`) %>%
    mutate(speciesOrder = order(`2.5%`)) %>%
    ggplot(aes(x = factor(species, level = species[speciesOrder]),
               ymin = `2.5%`,
               ymax = `97.5%`)) + geom_errorbar() +
    # facet_grid(rows = vars(cov)) +
    geom_hline(aes(yintercept = 0), color = "red") +
    ggtitle(covNames[cov]) +
    xlab("Species") +
    ylim(c(-3,3)) +
    theme_bw() +
    theme(
      axis.text = element_text(angle = 90,
                               size = 8),
      plot.title = element_text(hjust = .5,
                                size = 15)
    )

  plot2 <- data_plot_cov %>%
    filter(species %in% speciesNames[orderSpecies[73:144]]) %>%
    # filter(str_detect(cov, (!!covNames[cov]))) %>%
    # arrange(`2.5%`) %>%
    mutate(speciesOrder = order(`2.5%`)) %>%
    ggplot(aes(x = factor(species, level = species[speciesOrder]),
               ymin = `2.5%`,
               ymax = `97.5%`)) + geom_errorbar() +
    # facet_grid(rows = vars(cov)) +
    geom_hline(aes(yintercept = 0), color = "red") +
    ggtitle(covNames[cov]) +
    xlab("Species") +
    ylim(c(-3,3)) +
    theme_bw() +
    theme(
      axis.text = element_text(angle = 90,
                               size = 8),
      plot.title = element_text(hjust = .5,
                                size = 15)
    )

  plot3 <- data_plot_cov %>%
    filter(species %in% speciesNames[orderSpecies[145:218]]) %>%
    # filter(str_detect(cov, (!!covNames[cov]))) %>%
    # arrange(`2.5%`) %>%
    mutate(speciesOrder = order(`2.5%`)) %>%
    ggplot(aes(x = factor(species, level = species[speciesOrder]),
               ymin = `2.5%`,
               ymax = `97.5%`)) + geom_errorbar() +
    # facet_grid(rows = vars(cov)) +
    geom_hline(aes(yintercept = 0), color = "red") +
    ggtitle(covNames[cov]) +
    xlab("Species") +
    ylim(c(-3,3)) +
    theme_bw() +
    theme(
      axis.text = element_text(angle = 90,
                               size = 8),
      plot.title = element_text(hjust = .5,
                                size = 15)
    )

}

data_plot <- apply(samples_subset, 2, function(x) {
  quantile(logistic(x), probs = c(0.025, 0.975))
}) %>%
  t %>%
  as.data.frame %>%
  mutate(Species = speciesNames) %>%
  mutate(speciesOrder = order(`2.5%`))

orderSpecies <- order(data_plot$`2.5%`)

plotSpeciesRates <- function(subset){

  data_plot %>%
    filter(Species %in% speciesNames[orderSpecies[subset]]) %>%
    ggplot(aes(x =  factor(Species, level = speciesNames[orderSpecies]),
               ymin = `2.5%`,
               ymax = `97.5%`)) + geom_errorbar() +
    xlab("Species") +
    # ylim(c(0,1)) +
    ggtitle("Occupancy rates") +
    theme_bw() +
    ylim(c(0,1)) +
    theme(
      axis.text = element_text(angle = 90,
                               size = 8),
      plot.title = element_text(hjust = .5,
                                size = 15)
    )

}

(plot1 <- plotSpeciesRates(1:72))
plot2 <- plotSpeciesRates(73:144)
plot3 <- plotSpeciesRates(145:218)

setwd(here("Code","Results"))
ggsave(plot1, file = "OccupancyRates1.jpeg")
ggsave(plot2, file = "OccupancyRates2.jpeg")
ggsave(plot3, file = "OccupancyRates3.jpeg")

# OUTPUT - COLLECTION RATES ------

S <- ncol(data$OTU)

ncov_theta <- ncol(fit_ednaplus$X_theta)
speciesNames <- ncol(data$OTU)

param <- "beta_theta"

logistic <- function(x) 1 / (1 + exp(-x))

samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]

samples_subset <- samples_subset[,1 + 0:(S-1)*ncov_theta]

# data_plot <- apply(samples_subset, 2, function(x) {
#   quantile(logistic(x), probs = c(0.025, 0.975))
# }) %>%
#   t %>%
#   as.data.frame %>%
#   mutate(Species = colnames(y)) %>%
#   mutate(speciesOrder = order(`2.5%`)) %>%
#   ggplot(aes(x =  factor(Species, level = Species[speciesOrder]),
#              ymin = `2.5%`,
#              ymax = `97.5%`)) + geom_errorbar() +
#   xlab("Species") +
#   # ylim(c(0,1)) +
#   ggtitle("Collection rates") +
#   theme_bw() +
#   theme(
#     axis.text = element_text(angle = 90,
#                              size = 8),
#     plot.title = element_text(hjust = .5,
#                               size = 15)
#   )


data_plot <- apply(samples_subset, 2, function(x) {
  quantile(logistic(x), probs = c(0.025, 0.975))
}) %>%
  t %>%
  as.data.frame %>%
  mutate(Species = speciesNames) %>%
  mutate(speciesOrder = order(`2.5%`))

orderSpecies <- order(data_plot$`2.5%`)

plotSpeciesRates <- function(subset){

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


(plot1 <- plotSpeciesRates(1:50))
plot2 <- plotSpeciesRates(51:100)

setwd(here("Code Counts","Results"))
ggsave(plot1, file = "CollectionRates1.jpeg")
ggsave(plot2, file = "CollectionRates2.jpeg")
# ggsave(plot3, file = "CollectionRates3.jpeg")

# ggsave(data_plot, file = "CollectionRates.jpeg")


# OUTPUT - DETECTION RATES ------

param <- "p\\["

samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]

#
#
# (data_plot <- apply(samples_subset, 2, function(x) {
#   quantile(logistic(x), probs = c(0.025, 0.975))
# }) %>%
#   t %>%
#   as.data.frame %>%
#   mutate(Marker = factor(rep(c("16SFish","MiFish"), S)),
#          Species = rep(colnames(y), each = 2)) %>%
#   # mutate(speciesOrder = order(`2.5%`)) %>%
#   ggplot(aes(x =  Species,
#              color = Marker,
#              ymin = `2.5%`,
#              ymax = `97.5%`)) + geom_errorbar(linewidth = 1) +
#   xlab("Species") +
#   # ylim(c(0,1)) +
#   ggtitle("Collection rates") +
#   theme_bw() +
#   theme(
#     axis.text = element_text(angle = 90,
#                              size = 8),
#     plot.title = element_text(hjust = .5,
#                               size = 15)
#   )
# )

{

  data_plot <- apply(samples_subset, 2, function(x) {
    quantile(logistic(x), probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame %>%
    mutate(Marker = factor(rep(c("16SFish","MiFish"), S)),
           Species = rep(speciesNames, each = 2)) %>%
    mutate(speciesOrder = order(`2.5%`))

  orderSpecies <- order(data_plot$`2.5%`[data_plot$Marker == "16SFish"])

  minDetRates <- min(data_plot$`2.5%`)
  maxDetRates <- max(data_plot$`97.5%`)

  plotSpeciesRates <- function(subset){

    data_plot %>%
      filter(Species %in% speciesNames[orderSpecies[subset]]) %>%
      # ggplot(aes(x =  factor(Species, level = speciesNames[orderSpecies]),
      #            ymin = `2.5%`,
      #            ymax = `97.5%`)) +
      ggplot(aes(x =  factor(Species, level = speciesNames[orderSpecies]),
                 color = Marker,
                 ymin = `2.5%`,
                 ymax = `97.5%`)) + geom_errorbar() +
      xlab("Species") +
      ggtitle("Detection rates") +
      theme_bw() +
      ylim(c(minDetRates,maxDetRates)) +
      theme(
        axis.text = element_text(angle = 90,
                                 size = 8),
        plot.title = element_text(hjust = .5,
                                  size = 15)
      )

  }

  (plot1 <- plotSpeciesRates(1:50))
  (plot2 <- plotSpeciesRates(51:100))

  setwd(here("Code Counts","Results"))
  ggsave(plot1, file = "DetectionRates1.jpeg")
  ggsave(plot2, file = "DetectionRates2.jpeg")
  # ggsave(plot3, file = "DetectionRates3.jpeg")
}


# OUTPUT - FACTOR SCORES ------

param <- "U\\["

samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]

#
#
# (data_plot <- apply(samples_subset, 2, function(x) {
#   quantile(logistic(x), probs = c(0.025, 0.975))
# }) %>%
#   t %>%
#   as.data.frame %>%
#   mutate(Marker = factor(rep(c("16SFish","MiFish"), S)),
#          Species = rep(colnames(y), each = 2)) %>%
#   # mutate(speciesOrder = order(`2.5%`)) %>%
#   ggplot(aes(x =  Species,
#              color = Marker,
#              ymin = `2.5%`,
#              ymax = `97.5%`)) + geom_errorbar(linewidth = 1) +
#   xlab("Species") +
#   # ylim(c(0,1)) +
#   ggtitle("Collection rates") +
#   theme_bw() +
#   theme(
#     axis.text = element_text(angle = 90,
#                              size = 8),
#     plot.title = element_text(hjust = .5,
#                               size = 15)
#   )
# )

{
  apply(samples_subset, 2, function(x){
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
    t %>%
    as.data.frame %>%
    arrange(`50%`) %>%
    ggplot(aes(x = sites,
               y = `50%`,
               ymin = `2.5%`,
               ymax = `97.5%`)) +
    geom_point() +
    geom_errorbar() +
    coord_flip() + theme_bw() +
    ylab("Factor score")
}

# OUTPUT - FACTOR LOADINGS ------

param <- "LL\\["

samples_subset <- L_output_reparam# matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]

#
#
# (data_plot <- apply(samples_subset, 2, function(x) {
#   quantile(logistic(x), probs = c(0.025, 0.975))
# }) %>%
#   t %>%
#   as.data.frame %>%
#   mutate(Marker = factor(rep(c("16SFish","MiFish"), S)),
#          Species = rep(colnames(y), each = 2)) %>%
#   # mutate(speciesOrder = order(`2.5%`)) %>%
#   ggplot(aes(x =  Species,
#              color = Marker,
#              ymin = `2.5%`,
#              ymax = `97.5%`)) + geom_errorbar(linewidth = 1) +
#   xlab("Species") +
#   # ylim(c(0,1)) +
#   ggtitle("Collection rates") +
#   theme_bw() +
#   theme(
#     axis.text = element_text(angle = 90,
#                              size = 8),
#     plot.title = element_text(hjust = .5,
#                               size = 15)
#   )
# )

{
  apply(samples_subset, 2, function(x){
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
    t %>%
    as.data.frame %>%
    arrange(`50%`) %>%
    ggplot(aes(x = 1:S,
               y = `50%`,
               ymin = `2.5%`,
               ymax = `97.5%`)) +
    geom_point() +
    geom_errorbar() +
    coord_flip() + theme_bw() +
    ylab("Factor score")
}
