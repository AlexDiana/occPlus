# Tier 1: every advertised model configuration completes a fit.
#
# This is the tier that would have caught the sample_ls regression (TODO.Rmd
# Fixed bugs 22), which broke every non-spatial fit and survived because
# nothing exercised that path -- the vignette passes spatCovariates, so it
# would not have caught it either.
#
# Assertions are deliberately structural (see helper-fixtures.R): a completed
# fit with correctly shaped output. Numeric recovery belongs in tiers 2-3.

expect_valid_fit <- function(fit, n_species = FIXTURE_S) {
  expect_type(fit, "list")
  expect_true(all(c("results_output", "infos") %in% names(fit)))
  expect_identical(fit$infos$S, n_species)

  # B0_output (species intercepts) is asserted rather than psi_output because
  # it is the only posterior block every model type populates: `binary` and
  # `continuous` fits return jsdm_output and WAIC only -- no psi_output,
  # z_output or theta_output. Verified, not assumed.
  B0 <- fit$results_output$jsdm_output$B0_output
  expect_false(is.null(B0))
  expect_identical(nrow(B0), n_species)
  expect_true(any(is.finite(B0)))

  expect_true(is.finite(fit$results_output$WAIC))
  invisible(fit)
}

test_that("two-stage model fits, with and without a spatial field", {
  sim <- simulate_fixture(model = "two_stage")
  expect_valid_fit(fit_fixture(sim, spatial = TRUE))
  # The non-spatial case is the regression guard: see Fixed bugs 22.
  expect_valid_fit(fit_fixture(sim, spatial = FALSE))
})

test_that("occupancy model (field replicates, no PCR stage) fits", {
  expect_valid_fit(fit_fixture(simulate_fixture(model = "occupancy")))
})

test_that("binary model (pure JSDM, no replicates) fits", {
  sim <- simulate_fixture(model = "binary")
  expect_valid_fit(fit_fixture(sim, spatial = TRUE))
  expect_valid_fit(fit_fixture(sim, spatial = FALSE))
})

test_that("continuous model fits", {
  expect_valid_fit(fit_fixture(
    simulate_fixture(model = "continuous", M = 1L, P = 1L, K = 1L),
    spatial = FALSE))
})

test_that("multiple primers are handled", {
  fit <- fit_fixture(simulate_fixture(model = "two_stage", P = 3L))
  expect_valid_fit(fit)
  # Detection rates are estimated per primer, so p_output carries a primer
  # margin. Guards the P > 1 path behind Fixed bugs 2. unname() because the
  # dim vector carries a name here.
  expect_equal(unname(dim(fit$results_output$p_output)[1]), 3L)
})

test_that("model fits with no latent factors (d = 0)", {
  expect_valid_fit(fit_fixture(simulate_fixture(model = "two_stage", d = 0L),
                               d = 0L))
})

test_that("model fits without species traits", {
  expect_valid_fit(fit_fixture(simulate_fixture(model = "two_stage"),
                               traits = FALSE))
})
