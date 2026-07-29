# Tier 2: recovery canary. Skipped on CRAN; runs locally and in CI.
#
# NOT a calibration check -- R = 5 is far too few for that (PLAN.md 9). This
# is a smoke signal that parameter recovery has not grossly broken between
# full tier-3 runs, and in particular that the truth/posterior mapping in
# helper-simstudy.R still lines up. A mapping error is the failure this is
# best placed to catch: it drives coverage and the truth/estimate correlation
# to ~0 rather than merely degrading them.
#
# Thresholds are set from measurement, not taste. Re-derived 27 July under
# ds = 2 (the grid previously used ds = 0, which at the time produced no
# spatial field at all -- since fixed, TODO.md Fixed bugs 30). Eight independent seed sets at this
# configuration gave:
#
#                          min    median   max     floor   margin
#   overall coverage       0.85    0.88    0.90     0.70     0.15
#   lowest block coverage  0.72    0.78    0.80     0.40     0.32
#   cor(post_mean, truth)  0.50    0.61    0.71     0.30     0.20
#
# The floors sit well beneath the observed minima, so a pass is not luck and
# a failure is not noise. The narrowest margin is on the correlation check;
# it was measured over eight seed sets rather than three for that reason.
#
# Switching from ds = 0 to ds = 2 moved these very little -- overall coverage
# fell ~0.02 and the watch-list blocks were unchanged -- so the floors did not
# need revising. Recorded because "unchanged" is itself the useful finding: a
# real spatial field did not degrade recovery of the non-spatial parameters.

# Both tests below read the same 5-replicate run. Computed once: it costs
# ~30 s, and running it per test_that() would double tier 2 for nothing.
.recovery_cache <- new.env(parent = emptyenv())

recovery_rows <- function() {
  if (is.null(.recovery_cache$rows)) {
    scenario <- simstudy_scenarios()[[1]]   # base
    scenario$n <- 50L
    scenario$S <- 5L
    scenario$n_supportpoints <- 10L
    .recovery_cache$rows <- simstudy_scenario(
      scenario, R = 5L, verbose = FALSE,
      MCMCparams = list(nchain = 2, nburn = 300, niter = 300, nthin = 1))
  }
  .recovery_cache$rows
}

test_that("parameter recovery has not grossly broken (tier 2 canary)", {
  skip_on_cran()

  rows <- recovery_rows()
  expect_false(is.null(rows))
  summary_tbl <- simstudy_summarise(rows)

  # 1. Overall coverage. Measured ~0.88; nominal is 0.95. The floor is loose
  #    because R = 5 makes this noisy, and because two blocks are known to sit
  #    below nominal (see the test below).
  expect_gt(mean(rows$covered), 0.70)

  # 2. No block at or near zero. This is the mapping-error signature: if truth
  #    and posterior were misaligned for a block -- a transpose, a wrong
  #    index, a parameter that means something different on each side -- that
  #    block alone would collapse while the others looked fine. Exactly how
  #    sigma_bs was caught.
  expect_gt(min(summary_tbl$coverage), 0.40)

  # 3. Estimates track truth. Complements coverage: an interval can contain
  #    the truth while the point estimate is uninformative.
  expect_gt(stats::cor(rows$post_mean, rows$truth), 0.30)
})

test_that("known-undercovered blocks have not got worse", {
  skip_on_cran()

  # Across eight seed sets under ds = 2: beta_theta 0.72-0.84 (median 0.78)
  # and resid_cor 0.72-0.85 (median 0.80), against a nominal 0.95. Consistent
  # enough that chance is unlikely, but R = 5 is too small to call it -- the
  # effective SE is ~8%. Essentially unchanged from the earlier ds = 0
  # measurement, so the undercoverage is not something the null spatial field
  # was causing.
  #
  # Recorded here so that (a) the canary above is not tuned around them
  # silently, and (b) a further deterioration is caught. The full tier-3 run
  # at R = 100 is what should settle whether this is real; see
  # dev/simstudy/PLAN.md 10.3.
  summary_tbl <- simstudy_summarise(recovery_rows())

  watched <- summary_tbl[summary_tbl$block %in% c("beta_theta", "resid_cor"), ]
  expect_gt(nrow(watched), 0)

  # Floor of 0.55 against a measured 0.72-0.85: room for sampling noise at
  # R = 5, but a real regression would break through.
  expect_gt(min(watched$coverage), 0.55)
})
