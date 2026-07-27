# Tier-1 fixtures. See dev/simstudy/PLAN.md sections 5 and 6.2.
#
# Two rules govern everything in this file and the tier-1 tests that use it:
#
#   1. Speed. The whole tier is budgeted under 30 s because it ships to CRAN
#      and runs on every check.
#
#   2. No numeric reproducibility. While TODO.Rmd group A item 2 is open,
#      randinvg() draws from R's global RNG inside an OpenMP loop, so a fixed
#      seed does NOT reproduce across platforms -- deterministic on stock
#      macOS clang (no OpenMP), racy on Linux/Windows. A test asserting
#      expect_equal() on fitted values would pass locally and fail
#      intermittently on CRAN, which is the worst available failure mode.
#      Tier-1 assertions are therefore structural: does it run, are the
#      dimensions right, is a block populated, does a value vary.
#
#      Simulated *data* is reproducible (simulateOccJSDMData() uses only R's
#      RNG, no OpenMP), so tests on design matrices are safe.

# Small by design. n = 40 clears the 31-site spatial floor (TODO.Rmd group B
# item 3) with room to spare.
FIXTURE_N <- 40L
FIXTURE_S <- 4L

# Pinned rather than left to getDefaultSupportPoints(), which would request 30
# knots for 40 sites -- nearly one per site. See PLAN.md 5.4.
FIXTURE_KNOTS <- 8L

FIXTURE_MCMC <- list(nchain = 2, nburn = 20, niter = 20, nthin = 1)

fixture_occ_covariates  <- function() c("X_psi.EnvCov.1", "X_psi.EnvCov.2")
fixture_spat_covariates <- function() c("Xs.1", "Xs.2")

#' Simulate a small dataset for tier-1 use.
#'
#' Returns the full simulateOccJSDMData() output (true_params + data_list), so
#' callers that need the generating values can reach them.
simulate_fixture <- function(model = "two_stage",
                             n = FIXTURE_N,
                             S = FIXTURE_S,
                             M = 2L, P = 1L, K = 2L,
                             g = 2L, gt = 2L, d = 2L,
                             useSpatField = TRUE,
                             seed = 1L) {
  set.seed(seed)

  # M and K are only meaningful for the replicated designs; the simpler model
  # types collapse to one row per site.
  M_vec <- if (model %in% c("occupancy", "two_stage")) rep(M, n) else rep(1L, n)
  N     <- sum(M_vec)
  K_vec <- if (model == "two_stage") rep(K, N * P) else rep(1L, N * P)

  datasettings <- list(n = n, S = S, g = g, M = M_vec, P = P, K = K_vec,
                       ncov_psi = 2L, ncov_theta = 1L)

  jsdmParams <- list(gt = gt, d = d, ds = 0,
                     sigma_b = 0.5, sigma_bs = 0.5, sigma_ts = 0.5,
                     sigma_h = 1, sigma_s = 0.5, l_s = 0.3,
                     tau = rep(1, S),
                     useSpatField = useSpatField)

  params <- list(p              = matrix(stats::runif(P * S, 0.3, 0.6), P, S),
                 q              = matrix(stats::runif(P * S, 0.01, 0.05), P, S),
                 theta0         = stats::runif(S, 0.02, 0.1),
                 theta_baseline = stats::runif(S, 0.2, 0.5))

  simulateOccJSDMData(datasettings, params, jsdmParams, model = model)
}

#' Fit a simulated fixture with tier-1 MCMC settings.
fit_fixture <- function(sim,
                        spatial = TRUE,
                        d = 2L,
                        traits = TRUE,
                        MCMCparams = FIXTURE_MCMC) {
  dat <- sim$data_list
  if (!traits) dat$traits <- NULL

  listParams <- list(n_factors = d)
  spatCovariates <- NULL
  if (spatial) {
    spatCovariates <- fixture_spat_covariates()
    listParams$n_supportpoints <- FIXTURE_KNOTS
  }

  suppressMessages(suppressWarnings(
    runOccJSDM(dat,
               listParams      = listParams,
               occCovariates   = fixture_occ_covariates(),
               spatCovariates  = spatCovariates,
               MCMCparams      = MCMCparams)
  ))
}

# Cached fits -------------------------------------------------------------
#
# Several tests want the same fitted object and each fit costs ~0.5-1 s. Built
# lazily, so a run that skips those tests pays nothing.

.fixture_cache <- new.env(parent = emptyenv())

cached_fixture <- function(key, build) {
  if (is.null(.fixture_cache[[key]])) {
    .fixture_cache[[key]] <- build()
  }
  .fixture_cache[[key]]
}

#' Two-stage, spatial, multi-primer (P = 2).
#'
#' The workhorse fixture. P > 1 is required by the primer-subsetting
#' regression test (TODO.Rmd Fixed bugs 21), which cannot detect pooling
#' across primers when there is only one.
fixture_twostage <- function() {
  cached_fixture("twostage", function() {
    fit_fixture(simulate_fixture(model = "two_stage", P = 2L), spatial = TRUE)
  })
}

#' Continuous-response fit.
#'
#' Needed because tau is only sampled when model == "continuous"
#' (R/runOccJSDM.R:1130); for every other model type it is fixed at 1 and
#' tau_output is NA by design. See the Fixed bugs 8 test.
fixture_continuous <- function() {
  cached_fixture("continuous", function() {
    fit_fixture(simulate_fixture(model = "continuous", M = 1L, P = 1L, K = 1L),
                spatial = FALSE)
  })
}
