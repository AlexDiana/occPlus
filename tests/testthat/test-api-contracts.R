# Tier 1: argument handling and error paths.
#
# These guard the user-facing contract: that arguments are honoured, and that
# wrong input fails loudly rather than silently doing something else. Several
# of the fixed bugs were of exactly that second kind -- plotFPTPStage2Rates()
# quietly ignoring primerName being the clearest.

test_that("an unknown primerName errors and names the available primers", {
  fit <- fixture_twostage()
  expect_error(plotFPTPStage2Rates(fit, primerName = 99),
               "not found")
  # The message must be actionable: silently plotting the pooled result is how
  # the original bug (Fixed bugs 21) went unnoticed for so long.
  expect_error(plotFPTPStage2Rates(fit, primerName = 99),
               "Available primers")
})

test_that("primerName accepts the stored type and its character form", {
  # primerNames is stored as an *integer* vector, so a natural
  # primerName = "2" would fail an identity match. Both forms must work.
  fit <- fixture_twostage()
  p2 <- fit$infos$primerNames[2]
  expect_equal(plotFPTPStage2Rates(fit, primerName = p2)$data,
               plotFPTPStage2Rates(fit, primerName = as.character(p2))$data)
})

test_that("plotFPTPStage2Rates defaults to the first primer", {
  fit <- fixture_twostage()
  expect_equal(plotFPTPStage2Rates(fit)$data,
               plotFPTPStage2Rates(fit,
                                   primerName = fit$infos$primerNames[1])$data)
})

test_that("idx_species subsets the species actually plotted", {
  fit <- fixture_twostage()
  expect_equal(nrow(plotFPTPStage2Rates(fit, idx_species = 1:2)$data), 2L)
  expect_equal(nrow(plotFPTPStage2Rates(fit)$data), as.integer(FIXTURE_S))
})

test_that("an unrecognised covariate name is rejected", {
  sim <- simulate_fixture(model = "two_stage")
  expect_error(
    suppressMessages(suppressWarnings(
      runOccJSDM(sim$data_list,
                 occCovariates = "no_such_covariate",
                 MCMCparams = FIXTURE_MCMC)
    )),
    "not in data"
  )
})

test_that("n_supportpoints overrides the GP knot default", {
  # The default, max(30, floor(n * 0.2)), crashes below 31 sites and is a flat
  # 30 below 150 (TODO.Rmd group B item 3). The override is what makes small-n
  # spatial fits possible at all, and the simulation study depends on it
  # (PLAN.md 5.4), so it is worth a guard.
  sim <- simulate_fixture(model = "two_stage", n = 20L)
  expect_no_error(
    suppressMessages(suppressWarnings(
      runOccJSDM(sim$data_list,
                 listParams = list(n_factors = 2L, n_supportpoints = 6L),
                 occCovariates  = fixture_occ_covariates(),
                 spatCovariates = fixture_spat_covariates(),
                 MCMCparams = FIXTURE_MCMC)
    ))
  )
})
