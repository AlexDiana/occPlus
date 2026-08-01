# Tier 1: R/diagnostics.R.
#
# This file closes the coverage gap recorded under the diagnostics item in
# TODO.md: R/diagnostics.R had eight functions and no tests, three of them
# exported and user-facing.
#
# Most assertions here run on hand-built arrays rather than on a fitted model.
# That is deliberate and it is stronger, not weaker, than fitting: the shape
# and NA contracts below are properties of the diagnostics code alone, so
# testing them on synthetic input keeps them free of the RNG
# non-reproducibility described in helper-fixtures.R, and costs no fit time.
#
# Where a fitted object is genuinely needed (the two wrappers), the cached
# tier-1 fixtures are reused.

# Helper: a [dim1, dim2, niter, nchain] array of independent draws.
diag_array <- function(dim1 = 2L, dim2 = 3L, niter = 50L, nchain = 2L,
                       seed = 42L) {
  set.seed(seed)
  array(stats::rnorm(dim1 * dim2 * niter * nchain),
        dim = c(dim1, dim2, niter, nchain))
}

# as4d() -------------------------------------------------------------------
#
# Every other function in the file starts by calling this, so its padding rule
# is load-bearing for all of them.

test_that("as4d() pads 2- and 3-dimensional arrays to 4 dimensions", {
  # [niter, nchain]: a scalar parameter, padded on both index dimensions.
  scalar_par <- array(1:20, dim = c(10L, 2L))
  expect_equal(dim(occJSDM:::as4d(scalar_par)), c(1L, 1L, 10L, 2L))

  # [dim1, niter, nchain]: a species-indexed parameter, padded on dim2 only.
  vector_par <- array(1:60, dim = c(3L, 10L, 2L))
  expect_equal(dim(occJSDM:::as4d(vector_par)), c(3L, 1L, 10L, 2L))

  # Already 4d: returned untouched.
  full <- diag_array()
  expect_identical(occJSDM:::as4d(full), full)
})

test_that("as4d() preserves draws in the order the diagnostics expect", {
  # Padding must not permute the data: element [x, , i] of a 3d array has to
  # land at [x, 1, , i] of the padded one, or every Rhat and ESS below is
  # computed on the wrong series.
  vector_par <- array(stats::rnorm(3 * 10 * 2), dim = c(3L, 10L, 2L))
  padded <- occJSDM:::as4d(vector_par)
  expect_equal(padded[2L, 1L, , 2L], vector_par[2L, , 2L])
})

test_that("as4d() rejects input that is not an array of 2 to 4 dimensions", {
  expect_error(occJSDM:::as4d(1:10), "2, 3, or 4 dimensions")
  expect_error(occJSDM:::as4d(array(1:8, dim = c(2L, 2L, 2L, 1L, 1L))),
               "2, 3, or 4 dimensions")
})

# computeRhat() ------------------------------------------------------------

test_that("computeRhat() returns one Rhat per parameter element", {
  rhat <- occJSDM:::computeRhat(diag_array(dim1 = 2L, dim2 = 3L))
  expect_true(is.matrix(rhat))
  expect_equal(dim(rhat), c(2L, 3L))
})

test_that("computeRhat() is all NA with a single chain", {
  # Documented behaviour: the Gelman-Rubin statistic is undefined for one
  # chain. It must return NA rather than erroring, because
  # returnConvergenceDiagnostics() calls it unconditionally.
  rhat <- occJSDM:::computeRhat(diag_array(nchain = 1L))
  expect_true(all(is.na(rhat)))
})

test_that("computeRhat() returns NA for a parameter that never moved", {
  # coda::gelman.diag() errors on a zero-variance chain. A fixed parameter is
  # normal in this package (tau outside the continuous model, for one), so the
  # guard has to hold or diagnostics fall over on ordinary fits.
  arr <- diag_array(dim1 = 2L, dim2 = 1L)
  arr[1L, 1L, , ] <- 3.7
  rhat <- occJSDM:::computeRhat(arr)
  expect_true(is.na(rhat[1L, 1L]))
  expect_false(is.na(rhat[2L, 1L]))
})

test_that("computeRhat() finds no divergence between chains from one process", {
  # Independent draws from the same distribution: Rhat should sit near 1. The
  # bound is loose on purpose. This asserts the statistic is wired up and
  # oriented correctly, not that it hits a particular value.
  rhat <- occJSDM:::computeRhat(diag_array(niter = 400L, nchain = 3L))
  expect_true(all(is.finite(rhat)))
  expect_true(all(rhat < 1.1))
})

test_that("computeRhat() detects chains that disagree", {
  # The complement of the test above: offsetting one chain far from the other
  # must push Rhat up. Without this, a computeRhat() that returned 1 for
  # everything would pass the whole file.
  arr <- diag_array(dim1 = 1L, dim2 = 1L, niter = 200L, nchain = 2L)
  arr[1L, 1L, , 2L] <- arr[1L, 1L, , 2L] + 50
  expect_gt(occJSDM:::computeRhat(arr)[1L, 1L], 1.1)
})

# computeESSparams() -------------------------------------------------------

test_that("computeESSparams() returns one ESS per parameter element", {
  ess <- occJSDM:::computeESSparams(diag_array(dim1 = 2L, dim2 = 3L))
  expect_true(is.matrix(ess))
  expect_equal(dim(ess), c(2L, 3L))
  expect_true(all(ess > 0))
})

test_that("computeESSparams() pools chains", {
  # ESS is computed across the pooled chains, so doubling the number of chains
  # of independent draws should roughly double it. Asserted as an inequality
  # for the reason given in helper-fixtures.R.
  one <- occJSDM:::computeESSparams(
    diag_array(dim1 = 1L, dim2 = 1L, niter = 200L, nchain = 1L))
  two <- occJSDM:::computeESSparams(
    diag_array(dim1 = 1L, dim2 = 1L, niter = 200L, nchain = 2L))
  expect_gt(two[1L, 1L], one[1L, 1L])
})

test_that("computeESSparams() accepts a padded scalar parameter", {
  ess <- occJSDM:::computeESSparams(array(stats::rnorm(100L), dim = c(50L, 2L)))
  expect_equal(dim(ess), c(1L, 1L))
})

# paramOutputToLong() ------------------------------------------------------

test_that("paramOutputToLong() emits one row per draw per chain per element", {
  long <- occJSDM:::paramOutputToLong(
    diag_array(dim1 = 2L, dim2 = 3L, niter = 10L, nchain = 2L))
  expect_equal(nrow(long), 2L * 3L * 10L * 2L)
  expect_named(long,
               c("param", "label1", "label2", "chain", "iter", "value"))
})

test_that("paramOutputToLong() honours supplied dimension labels", {
  long <- occJSDM:::paramOutputToLong(
    diag_array(dim1 = 2L, dim2 = 2L, niter = 5L),
    param_name = "beta_psi",
    dimnames1 = c("cov_a", "cov_b"),
    dimnames2 = c("sp1", "sp2"))
  expect_setequal(long$label1, c("cov_a", "cov_b"))
  expect_setequal(long$label2, c("sp1", "sp2"))
  expect_true(all(long$param == "beta_psi"))
})

test_that("paramOutputToLong() defaults labels to the index positions", {
  long <- occJSDM:::paramOutputToLong(
    diag_array(dim1 = 2L, dim2 = 1L, niter = 5L))
  expect_setequal(long$label1, c("1", "2"))
})

# summarisePosterior() -----------------------------------------------------

test_that("summarisePosterior() returns the documented columns, one row per element", {
  out <- occJSDM:::summarisePosterior(diag_array(dim1 = 2L, dim2 = 3L))
  expect_named(out, c("param", "idx1", "idx2", "label1", "label2",
                      "mean", "sd", "q2.5", "q97.5", "rhat", "ess"))
  expect_equal(nrow(out), 2L * 3L)
})

test_that("summarisePosterior() summarises the draws it was given", {
  # Deterministic given the input array, so this is safe to assert exactly:
  # the arithmetic is the function's own, not the sampler's.
  arr <- diag_array(dim1 = 2L, dim2 = 2L, niter = 30L, nchain = 2L)
  out <- occJSDM:::summarisePosterior(arr)
  row <- out[out$idx1 == 2L & out$idx2 == 1L, ]
  draws <- as.vector(arr[2L, 1L, , ])
  expect_equal(row$mean, mean(draws))
  expect_equal(row$sd, stats::sd(draws))
  expect_equal(row$q2.5, as.numeric(stats::quantile(draws, .025)))
  expect_equal(row$q97.5, as.numeric(stats::quantile(draws, .975)))
})

test_that("summarisePosterior() credible interval brackets the mean", {
  out <- occJSDM:::summarisePosterior(diag_array(niter = 100L))
  expect_true(all(out$q2.5 <= out$mean))
  expect_true(all(out$mean <= out$q97.5))
})

test_that("summarisePosterior() carries NA rhat through from a single chain", {
  out <- occJSDM:::summarisePosterior(diag_array(nchain = 1L))
  expect_true(all(is.na(out$rhat)))
  expect_true(all(is.finite(out$ess)))
})

# plotTraceplot() ----------------------------------------------------------

test_that("plotTraceplot() returns a ggplot over the supplied draws", {
  p <- plotTraceplot(diag_array(dim1 = 2L, dim2 = 2L, niter = 10L))
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 2L * 2L * 10L * 2L)
})

test_that("plotTraceplot() facets on one dimension when the second is absent", {
  # dim2 == 1 is the common case (a species-indexed parameter). Faceting on
  # both would produce a column of empty strips.
  one_dim <- plotTraceplot(diag_array(dim1 = 3L, dim2 = 1L, niter = 10L))
  two_dim <- plotTraceplot(diag_array(dim1 = 3L, dim2 = 2L, niter = 10L))
  expect_length(one_dim$facet$params$facets, 1L)
  expect_length(two_dim$facet$params$facets, 2L)
})

test_that("plotTraceplot() labels the y axis with the parameter name", {
  p <- plotTraceplot(diag_array(dim1 = 1L, dim2 = 1L, niter = 10L),
                     param_name = "beta0_psi")
  expect_equal(p$labels$y, "beta0_psi")
})

# returnConvergenceDiagnostics() -------------------------------------------

cached_convergence <- function() {
  cached_fixture("convergence_twostage",
                 function() returnConvergenceDiagnostics(fixture_twostage()))
}

test_that("returnConvergenceDiagnostics() returns the documented table", {
  out <- cached_convergence()
  expect_s3_class(out, "tbl_df")
  expect_named(out, c("param", "idx1", "idx2", "label1", "label2",
                      "mean", "sd", "q2.5", "q97.5", "rhat", "ess"))
  expect_gt(nrow(out), 0L)
})

test_that("returnConvergenceDiagnostics() covers the two-stage parameter blocks", {
  # A two-stage eDNA fit has the detection blocks as well as the occupancy
  # ones. If a block silently stopped being assembled, the table would still
  # be well formed, so the contents have to be asserted separately.
  out <- cached_convergence()
  expect_setequal(unique(out$param),
                  c("beta0_psi", "beta_psi", "beta_theta", "p", "q", "theta0"))
})

test_that("returnConvergenceDiagnostics() labels rows with species names", {
  fit <- fixture_twostage()
  out <- cached_convergence()
  beta0 <- out[out$param == "beta0_psi", ]
  expect_equal(nrow(beta0), as.integer(FIXTURE_S))
  expect_setequal(beta0$label1, as.character(fit$infos$speciesNames))
})

test_that("returnConvergenceDiagnostics() skips blocks a fit does not have", {
  # Documented behaviour: components absent from results_output are silently
  # dropped rather than erroring. The continuous fit is the case that matters,
  # since it has no two-stage detection blocks.
  out <- returnConvergenceDiagnostics(fixture_continuous())
  expect_gt(nrow(out), 0L)
  expect_true("beta0_psi" %in% out$param)
  expect_false(any(c("p", "q") %in% out$param))
})

# computeDiagnostics() -----------------------------------------------------
#
# This one reports by printing, so its contract is its console output. Capture
# messages and warnings together: the tier-1 fixture trips the convergence
# thresholds by construction (see the low-ESS test), so a test that captured
# only messages would leak a dozen warnings into the run and bury real ones.

capture_diagnostics <- function(results_output) {
  warns <- character()
  msgs <- withCallingHandlers(
    testthat::capture_messages(computeDiagnostics(results_output)),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(messages = paste(msgs, collapse = ""), warnings = warns)
}

# Cached via the helper-fixtures.R cache: several tests below want the same
# console capture, and it costs a gelman.diag() call per parameter element.
# Measured at 1.0 s for the whole file either way, so this is tidiness rather
# than a budget fix -- the tier-1 total is 23 s and test-recovery.R is 13 s of
# it. Do not assume removing the cache would cost anything material.
cached_diagnostics <- function() {
  cached_fixture("diagnostics_twostage",
                 function() capture_diagnostics(fixture_twostage()$results_output))
}

test_that("computeDiagnostics() returns NULL invisibly", {
  fit <- fixture_twostage()
  expect_invisible(
    suppressWarnings(suppressMessages(computeDiagnostics(fit$results_output))))
  out <- suppressWarnings(suppressMessages(
    computeDiagnostics(fit$results_output)))
  expect_null(out)
})

test_that("computeDiagnostics() reports every sampled block of a two-stage fit", {
  # Both the detection blocks carried at the top level of results_output and
  # the JSDM blocks nested under jsdm_output must appear. They are assembled
  # by two different code paths, so one can go missing without the other.
  out <- cached_diagnostics()
  for (block in c("beta_theta", "p", "q", "theta0",
                  "B0", "B", "L", "sigmah")) {
    expect_match(out$messages, paste0("Parameter Set: ", block, " "),
                 fixed = TRUE, info = block)
  }
  expect_match(out$messages, "ESS")
})

test_that("computeDiagnostics() names blocks differently from returnConvergenceDiagnostics()", {
  # Documenting a real footgun rather than asserting it is good. The two
  # exported diagnostics functions use different vocabularies for the same
  # parameters: computeDiagnostics() strips "_output" off the storage name and
  # prints "B0", while returnConvergenceDiagnostics() relabels it "beta0_psi".
  # Anyone matching output from one against the other needs to know. If these
  # are ever reconciled, this test should fail and be deleted.
  out <- cached_diagnostics()
  expect_match(out$messages, "Parameter Set: B0 ", fixed = TRUE)
  expect_false(grepl("Parameter Set: beta0_psi", out$messages, fixed = TRUE))
  expect_true("beta0_psi" %in% cached_convergence()$param)
})

test_that("computeDiagnostics() skips the latent-state blocks", {
  # z, w, psi, theta and the variance partition are site- or draw-level latent
  # quantities, not parameters to assess for convergence. Printing them would
  # bury the blocks that matter under thousands of rows.
  out <- cached_diagnostics()
  for (skipped in c("z", "w", "psi", "theta", "idx_ls", "varPart")) {
    expect_false(
      grepl(paste0("Parameter Set: ", skipped, " "), out$messages, fixed = TRUE),
      info = skipped)
  }
})

test_that("computeDiagnostics() warns when effective sample size is low", {
  # The tier-1 fixture runs 20 iterations on 2 chains, so ESS is far below the
  # 50 threshold by construction. This asserts the warning path fires at all,
  # which is the half of the function a user actually depends on.
  out <- cached_diagnostics()
  expect_true(any(grepl("Low ESS", out$warnings)))
})

test_that("computeDiagnostics() warns on a block that has not converged", {
  # The other threshold, max R-hat > 1.1. Asserted on the same short fixture
  # for the same reason: 20 iterations do not converge.
  out <- cached_diagnostics()
  expect_true(any(grepl("Convergence issue", out$warnings)))
})
