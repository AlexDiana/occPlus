# Tier 2: per-species occupancy recovery. Skipped on CRAN; runs locally and in CI.
#
# test-recovery.R asks whether the parameters are recovered, pooled over every
# block and every species. It never looks at occupancy itself, and it never
# looks at one species on its own. This file fills both gaps: it prints one row
# per species showing how well psi and z were recovered together with the
# coefficients driving them, and it asserts on quantities derived from that
# table. The table is the thing you read; the assertions are what fails.
#
# Occupancy is scored as point-estimate agreement, not coverage, because there
# are no draws to build an interval from: runOccJSDM() defaults to
# summarisedLatentPresences = TRUE, so psi_output and z_output come back as
# n x S running means with no iteration or chain dimension.
#
# Two metrics per species, deliberately different in kind:
#
#   psi_cor  Pearson correlation across sites between the estimated psi and the
#            true plogis(eta). Null value 0. Invariant to the level of psi, so
#            it asks only whether the model put occupancy in the right places.
#   auc      Probability that a randomly chosen truly-occupied site is ranked
#            above a randomly chosen unoccupied one. Null value 0.5. Its ceiling
#            is well below 1 because z is a coin flip given psi, so even a
#            perfect psi cannot rank every site correctly.
#
# Both are printed. Only psi_cor is asserted on. Measured over 80 species-fits
# they correlate at 0.94, so the AUC column earns its place by being readable
# on its own scale, not by carrying independent information.

# --- Why the assertions are pooled and counted, never per species ---------
#
# Configuration: base scenario at n = 100, S = 10, n_supportpoints = 20, with
# nchain 2, nburn 500, niter 500. Eight independent seed sets gave:
#
#                          min    median   max     floor   margin
#   mean psi_cor            0.498   0.602   0.753    0.30    0.20
#   species with cor >= 0.3 7      8.5     10        5       2
#   mean AUC                0.748   0.779   0.855    --      --
#   driver coverage         0.800   0.867   1.000    0.55    0.25
#   driver cor(est, truth)  0.647   0.741   0.900    0.30    0.35
#
# The single per-species psi_cor over those 80 species-fits ranged from -0.401
# to 0.913. That range is the whole reason there is no "every species clears X"
# assertion here: any floor high enough to be diagnostic fails on ordinary
# seeds, and any floor low enough to pass is indistinguishable from the null.
# Repeating one seed four times moved the worst species from 0.107 to 0.302, so
# this is not a property of particular datasets that a fixed seed would pin
# down -- the sampler is not bit-reproducible (randinvg() draws from R's global
# RNG inside an OpenMP loop), and the minimum of ten noisy correlations is the
# least stable summary available.
#
# What is stable is the mean (0.642 to 0.659 across those same four repeats)
# and, less so, the count above a threshold (8 to 10). So the two assertions
# are a tight floor on the mean, which catches recovery degrading across the
# board, and a loose floor on the count, which catches several species
# collapsing while the average holds. They fail for different reasons.
#
# The smaller configuration used by test-recovery.R (n = 50, S = 5) was
# measured first and rejected: mean psi_cor fell to 0.326 on the worst of eight
# seeds, leaving no room under a 0.30 floor. The larger fit costs ~13 s.

.occrec_cache <- new.env(parent = emptyenv())

occrec_scenario <- function() {
  sc <- simstudy_scenarios()[[1]]   # base
  sc$n <- 100L
  sc$S <- 10L
  sc$n_supportpoints <- 20L
  sc
}

occrec_mcmc <- list(nchain = 2, nburn = 500, niter = 500, nthin = 1)

#' Probability a random occupied site outranks a random unoccupied one.
#'
#' The rank form of the Mann-Whitney U, so ties contribute 0.5 rather than
#' being dropped. NA when a species is entirely occupied or entirely absent,
#' since there is then nothing to discriminate between.
occrec_auc <- function(score, label) {
  pos <- score[label == 1]
  neg <- score[label == 0]
  if (length(pos) == 0L || length(neg) == 0L) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}

#' One row per species: occupancy recovery and the drivers behind it.
#'
#' Cached because the fit costs ~13 s and both tests read it.
occrec_table <- function() {
  if (!is.null(.occrec_cache$tbl)) return(.occrec_cache$tbl)

  scenario <- occrec_scenario()
  truth <- draw_truth(scenario, simstudy_seed_for(scenario, 1L))
  sim <- simstudy_simulate(truth)
  fit <- simstudy_fit(sim, scenario, occrec_mcmc)

  S <- scenario$S
  psi_true <- stats::plogis(sim$true_params$jsdmParams_true$eta)
  z_true <- sim$true_params$z_true
  psi_est <- fit$results_output$psi_output
  z_est <- fit$results_output$z_output

  # The occupancy drivers: the species intercept B0 and its environmental
  # slopes B. Reuse simstudy_rows() rather than reaching into the posterior
  # here, so these inherit the length-mismatch guard in helper-simstudy.R.
  rows <- simstudy_rows(fit, sim, truth, scenario, replicate = 1L)
  drv <- rows[rows$block %in% c("B0", "B"), ]

  # Element to species. Truth is flattened column-major by
  # posterior_draw_matrix(), so B0 (length S) gives one element per species and
  # B (ncov_psi x S) gives ncov_psi consecutive elements per species. Asserted
  # rather than assumed: a wrong species mapping would not fail loudly, it
  # would quietly attribute every species' coefficients to its neighbour.
  ncov <- scenario$ncov_psi
  stopifnot(nrow(drv) == S * (1L + ncov))
  drv$species <- ifelse(drv$block == "B0",
                        drv$element,
                        ceiling(drv$element / ncov))
  stopifnot(all(tabulate(drv$species, nbins = S) == 1L + ncov))

  tbl <- data.frame(
    species  = seq_len(S),
    occ_true = colMeans(z_true),
    occ_est  = colMeans(z_est),
    psi_cor  = vapply(seq_len(S), function(s)
                      stats::cor(psi_est[, s], psi_true[, s]), numeric(1)),
    auc      = vapply(seq_len(S), function(s)
                      occrec_auc(psi_est[, s], z_true[, s]), numeric(1)),
    drv_cov  = vapply(seq_len(S), function(s)
                      mean(drv$covered[drv$species == s]), numeric(1)),
    drv_err  = vapply(seq_len(S), function(s) {
                      d <- drv[drv$species == s, ]
                      mean(abs(d$post_mean - d$truth))
                    }, numeric(1)),
    stringsAsFactors = FALSE
  )

  .occrec_cache$tbl <- tbl
  .occrec_cache$drv <- drv
  tbl
}

occrec_drivers <- function() {
  occrec_table()
  .occrec_cache$drv
}

#' Emit the table so it is visible when the suite runs.
#'
#' message() rather than print(): it goes to stderr and survives testthat's
#' reporters, which is the point of the file.
occrec_show <- function(tbl) {
  txt <- utils::capture.output(print(tbl, digits = 3, row.names = FALSE))
  message("\nper-species occupancy recovery ",
          "(psi_cor null 0, auc null 0.5, drv_cov nominal 0.95)\n",
          paste(txt, collapse = "\n"), "\n")
}

test_that("occupancy is recovered per species (tier 2)", {
  skip_on_cran()

  tbl <- occrec_table()
  occrec_show(tbl)

  expect_equal(nrow(tbl), occrec_scenario()$S)
  # A NA here means a species was entirely occupied or entirely absent, which
  # makes its row uninterpretable rather than merely bad.
  expect_true(all(is.finite(tbl$psi_cor)))
  expect_true(all(is.finite(tbl$auc)))

  # 1. Pooled. The stable summary: measured 0.498-0.753 across eight seeds and
  #    0.642-0.659 across four repeats of one seed. This is what catches
  #    recovery degrading across the board.
  expect_gt(mean(tbl$psi_cor), 0.30)

  # 2. Counted. Measured 7-10 of 10; the floor is half the species. This is
  #    what catches several species collapsing while the mean is held up by the
  #    rest -- the blind spot a pooled check has by construction.
  expect_gte(sum(tbl$psi_cor >= 0.30), ceiling(occrec_scenario()$S / 2))

  # 3. Discrimination is printed, not asserted per species, for the reason in
  #    the header. One pooled check that it is above its null, which a
  #    site-order or species-order transposition would break outright.
  expect_gt(mean(tbl$auc), 0.60)
})

test_that("the coefficients driving occupancy are recovered per species (tier 2)", {
  skip_on_cran()

  tbl <- occrec_table()
  drv <- occrec_drivers()

  # B0 plus ncov_psi slopes per species, so 3 elements here. Three Bernoulli
  # draws cannot support a per-species coverage assertion: the only values
  # available are 0, 1/3, 2/3 and 1, and the observed per-species minimum over
  # eight seeds was 1/3 on a model that is covering at 0.87 overall. Likewise
  # the per-species correlation of estimate against truth is computed on three
  # points and went as low as -0.877. Per species these columns are for
  # reading; the assertions below are pooled.
  expect_equal(nrow(drv), occrec_scenario()$S * (1L + occrec_scenario()$ncov_psi))
  expect_true(all(is.finite(tbl$drv_err)))

  # Coverage of the occupancy-driver elements. Nominal 0.95, measured
  # 0.800-1.000 over eight seeds. The floor is loose because 30 elements from
  # one fit is a small sample, but a mapping error would drive this to near 0.
  expect_gt(mean(drv$covered), 0.55)

  # Point estimates track truth. Measured 0.647-0.900. Complements coverage in
  # the same way as test-recovery.R item 3: an interval can contain the truth
  # while the point estimate says nothing.
  expect_gt(stats::cor(drv$post_mean, drv$truth), 0.30)
})
