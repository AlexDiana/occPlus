# Tier 2: recovery canary. Skipped on CRAN; runs locally and in CI.
#
# NOT a calibration check -- R = 5 is far too few for that (PLAN.md 9). This
# is a smoke signal that parameter recovery has not grossly broken between
# full tier-3 runs, and in particular that the truth/posterior mapping in
# helper-simstudy.R still lines up. A mapping error is the failure this is
# best placed to catch: it drives coverage and the truth/estimate correlation
# to ~0 rather than merely degrading them.
#
# Thresholds are set from measurement, not taste. Three independent seed sets
# at this configuration gave:
#
#   overall coverage      0.89, 0.89, 0.88
#   lowest block coverage 0.78, 0.76, 0.76
#   cor(post_mean, truth) 0.62, 0.56, 0.68
#
# The floors below sit well beneath those, so a pass is not luck and a
# failure is not noise.

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

  # beta_theta and resid_cor sat at 0.78-0.80 and 0.76-0.79 across three
  # independent seed sets -- consistent enough that it is unlikely to be
  # chance, but R = 5 is too small to call it. Recorded here so that (a) the
  # canary above is not tuned around them silently, and (b) a further
  # deterioration is caught. The full tier-3 run at R = 100 is what should
  # settle whether this is real; see dev/simstudy/PLAN.md 10.
  summary_tbl <- simstudy_summarise(recovery_rows())

  watched <- summary_tbl[summary_tbl$block %in% c("beta_theta", "resid_cor"), ]
  expect_gt(nrow(watched), 0)

  # Floor of 0.55 against a measured 0.76-0.80: room for sampling noise at
  # R = 5, but a real regression would break through.
  expect_gt(min(watched$coverage), 0.55)
})
