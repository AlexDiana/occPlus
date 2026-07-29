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
  # The default used to be max(30, floor(n * 0.2)), which crashed below 31
  # sites and was a flat 30 below 150. Fixed in 42198d9 -- it is now
  # min(floor(n * 0.2), n - 1) (TODO.Rmd Fixed bugs 29) -- so this no longer
  # guards a crash. Kept because the simulation study pins n_supportpoints
  # rather than trusting any default (PLAN.md 5.4), and that override needs to
  # keep working.
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

# --- GAM covariate-effect exports ----------------------------------------
#
# Added in Alex's 29 July pull. Three functions were exported with no test
# reference at all, so nothing caught a signature change or a broken return.
# Smoke level only: these assert shape and that the call runs, not numbers.

test_that("returnCovariateEffect() returns one row per grid point per species", {
  fit <- fixture_twostage()
  eff <- returnCovariateEffect(fit, covName = "X_psi.EnvCov.1",
                               idx_species = 1:2, confidence = 0.95)
  expect_s3_class(eff, "data.frame")
  expect_gt(nrow(eff), 0)
  # a fitted value and an interval, whatever they end up being called
  expect_gte(ncol(eff), 4)
  expect_true(any(grepl("species|sp", names(eff), ignore.case = TRUE)))
})

test_that("plotCovariateEffect() builds a plot for each covariate asked for", {
  fit <- fixture_twostage()
  p <- plotCovariateEffect(fit, covNames = "X_psi.EnvCov.1",
                           idx_species = 1:2, confidence = 0.95)
  expect_type(p, "list")
  expect_length(p, 1)
})

test_that("the GAM effect functions reject a covariate that is not in the fit", {
  fit <- fixture_twostage()
  expect_error(
    returnCovariateEffect(fit, covName = "NotACovariate", idx_species = 1:2)
  )
})
