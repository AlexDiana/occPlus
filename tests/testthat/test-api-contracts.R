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
  # min(floor(n * 0.2), n - 1) (TODO.md Fixed bugs 29) -- so this no longer
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

# --- predictNewSites() argument handling ---------------------------------
#
# X_psi and X_s had no defaults, so the is.null() guards the author wrote
# could never fire: R raised "argument is missing, with no default" first.
# Worse, the guards used `&` rather than `&&`, which evaluates both sides, so
# the missing promises were forced even when the caller had asked for neither
# term. There was no way to call the function without supplying both matrices.
# Fixing that exposed three further defects downstream, each covered below.

test_that("predictNewSites() names the argument it needs", {
  fit <- fixture_twostage()
  expect_error(predictNewSites(fit), "X_psi is required")
})

test_that("predictNewSites() runs with only the term it is asked for", {
  # Previously impossible: useSpatial = FALSE reached C++ that sliced Ks_all
  # and Bs_output unconditionally, aborting with "Cube::slice(): index out of
  # bounds". And useEnvCov = FALSE gave n = 0 rows, because the site count was
  # read from X unconditionally.
  sim <- simulate_fixture(model = "two_stage", useSpatField = TRUE)
  fit <- fit_fixture(sim, MCMCparams = FIXTURE_MCMC)
  X_psi <- sim$data_list$info[1:5, fixture_occ_covariates(), drop = FALSE]
  X_s <- sim$data_list$info[1:5, fixture_spat_covariates(), drop = FALSE]

  full <- suppressMessages(predictNewSites(fit, X_psi = X_psi, X_s = X_s))
  no_spatial <- suppressMessages(
    predictNewSites(fit, X_psi = X_psi, useSpatial = FALSE))
  no_envcov <- suppressMessages(
    predictNewSites(fit, X_s = X_s, useEnvCov = FALSE))

  # (quantiles, new sites, species) in all three cases
  for (p in list(full, no_spatial, no_envcov)) {
    expect_equal(dim(p)[1], 3L)
    expect_equal(dim(p)[2], 5L)
    expect_equal(dim(p)[3], fit$infos$S)
    expect_true(all(p >= 0 & p <= 1, na.rm = TRUE))
    expect_true(all(p[1, , ] <= p[2, , ] & p[2, , ] <= p[3, , ], na.rm = TRUE))
  }

  # Each term must actually contribute, or the switches are cosmetic and the
  # dimension checks above would pass on a function that silently ignored them.
  expect_false(isTRUE(all.equal(full, no_spatial)))
  expect_false(isTRUE(all.equal(full, no_envcov)))
})

test_that("predictNewSites() rejects a request with no site count", {
  # With both terms off nothing says how many new sites to predict for. This
  # used to return a zero-row array rather than complain.
  fit <- fixture_twostage()
  expect_error(
    suppressMessages(predictNewSites(fit, useEnvCov = FALSE, useSpatial = FALSE)),
    "Nothing determines the new sites")
})

test_that("predictNewSites() handles a fit with no spatial field", {
  # Defaults used to hard-stop here ("Cannot use spatial if it was not
  # estimated initially") despite the docs promising the term was ignored when
  # absent. Asking for it explicitly should still be an error.
  sim <- simulate_fixture(model = "two_stage", useSpatField = FALSE)
  fit <- fit_fixture(sim, spatial = FALSE, MCMCparams = FIXTURE_MCMC)
  X_psi <- sim$data_list$info[1:5, fixture_occ_covariates(), drop = FALSE]

  expect_no_error(suppressMessages(predictNewSites(fit, X_psi = X_psi)))
  expect_error(
    suppressMessages(predictNewSites(fit, X_psi = X_psi, useSpatial = TRUE)),
    "not estimated in this fit")
})

test_that("the GAM effect functions reject a fit predating infos$X0_psi", {
  # X0_psi holds the raw occupancy covariates and was added by e60e3ad. A fit
  # made before that lacks it, and the grid construction then does min(NULL)
  # -> Inf and seq(Inf, -Inf, ...) -> "'from' must be a finite number", which
  # says nothing about the cause. That is exactly how the shipped
  # sampleresults dataset fails and how the vignette build breaks, so the
  # guard is aimed at a real object, not a hypothetical one.
  fit <- fixture_twostage()
  expect_false(is.null(fit$infos$X0_psi))   # a current fit must carry it

  stale <- fit
  stale$infos$X0_psi <- NULL
  expect_error(returnCovariateEffect(stale, "X_psi.EnvCov.1", 1:2),
               "raw occupancy covariates")
  expect_error(plotCovariateEffect(stale, "X_psi.EnvCov.1", 1:2),
               "raw occupancy covariates")

  # and the message must name the function the user actually called
  expect_error(returnCovariateEffect(stale, "X_psi.EnvCov.1", 1:2),
               "returnCovariateEffect")
  expect_error(plotCovariateEffect(stale, "X_psi.EnvCov.1", 1:2),
               "plotCovariateEffect")
})

test_that("symbols used by live code are actually imported into the namespace", {
  # Guards TODO.md group D item 4 phase 1. These three are called by live code
  # but were missing from NAMESPACE, which R CMD check reported only as a NOTE
  # about undefined globals.
  #
  # It was not only a NOTE. pivot_longer is reached in the *categorical*
  # covariate branch of returnCovariateEffect_base() (R/jsdmfun.R:429) and
  # plotCovariateEffect_base() (:612); the numeric branch summarises to
  # quantiles and never calls it. Verified against a properly installed copy
  # with tidyr unattached: pivot_longer did not resolve from the namespace at
  # all, so both exported GAM functions would fail with "could not find
  # function" on any categorical occupancy covariate.
  #
  # Asserted on the imports environment rather than by calling the functions,
  # deliberately. A functional test would pass under devtools::test() whether
  # or not the import exists, because load_all() resolves unimported symbols
  # through the global environment. That is exactly how this gap survived a
  # suite that has caught a lot else. Checking the import directly holds under
  # both load_all() and R CMD check.
  imports <- parent.env(asNamespace("occJSDM"))
  for (fn in c("pivot_longer", "setNames", "rnbinom")) {
    expect_true(exists(fn, envir = imports, inherits = FALSE), info = fn)
  }
})
