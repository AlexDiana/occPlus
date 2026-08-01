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

# The floors, in one place, so the printed report and the expectations below
# cannot drift apart. `swept` is the min-to-max range observed over the eight
# seed sets described in the header; it is text because it is a record of a
# measurement, not something recomputed at run time.
OCCREC_SPECIES_CUT <- 0.30   # per-species psi_cor threshold that the count uses
OCCREC_AUC_CUT <- 0.60       # per-species auc threshold, flagged but not counted

OCCREC_FLOORS <- list(
  mean_psi_cor = list(floor = 0.30, swept = "0.498 - 0.753",
    what = "mean per-species correlation of estimated psi with the truth",
    why  = "The stable summary. Varied by 0.017 across four repeats of one seed, so a failure here is not noise. Catches recovery degrading across every species at once."),
  n_species_ok = list(floor = NULL, swept = "7 - 10 of 10",
    what = "species whose psi correlation clears 0.30",
    why  = "Floor is half the species. Catches several species collapsing while the mean is held up by the rest, which is the blind spot a pooled check has by construction."),
  mean_auc = list(floor = 0.60, swept = "0.748 - 0.855",
    what = "mean per-species discrimination of the true occupancy states",
    why  = "Pooled only. Correlates with psi_cor at 0.94 over 80 species-fits, so it adds little independently, but a site-order or species-order transposition would drive it to its 0.5 null outright."),
  drv_cov = list(floor = 0.55, swept = "0.800 - 1.000",
    what = "fraction of B0 and B elements inside their 95% interval",
    why  = "Nominal is 0.95. Loose floor because 30 elements from one fit is a small sample, but a wrong truth-to-posterior mapping drives this to near 0 rather than merely degrading it."),
  drv_cor = list(floor = 0.30, swept = "0.647 - 0.900",
    what = "correlation of B0 and B posterior means with their truths",
    why  = "Complements coverage the same way test-recovery.R item 3 does: an interval can contain the truth while the point estimate says nothing.")
)

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

#' The five asserted quantities, computed once.
#'
#' Both the printed report and the expectations read this, so a floor cannot be
#' explained one way and tested another.
occrec_checks <- function() {
  if (!is.null(.occrec_cache$chk)) return(.occrec_cache$chk)
  tbl <- occrec_table()
  drv <- occrec_drivers()

  .occrec_cache$chk <- list(
    mean_psi_cor = mean(tbl$psi_cor),
    n_species_ok = sum(tbl$psi_cor >= OCCREC_SPECIES_CUT),
    n_species_floor = ceiling(occrec_scenario()$S / 2),
    mean_auc = mean(tbl$auc),
    drv_cov = mean(drv$covered),
    drv_cor = stats::cor(drv$post_mean, drv$truth)
  )
  .occrec_cache$chk
}

#' Emit the table so it is visible when the suite runs.
#'
#' message() rather than print(): it goes to stderr and survives testthat's
#' reporters, which is the point of the file. The legend and the flags are here
#' so that a row can be read without opening this source file.
occrec_show <- function(tbl) {
  sc <- occrec_scenario()

  flag <- ifelse(tbl$psi_cor < OCCREC_SPECIES_CUT,
                 sprintf("<- psi_cor below %.2f", OCCREC_SPECIES_CUT),
                 ifelse(tbl$auc < OCCREC_AUC_CUT,
                        sprintf("<- auc below %.2f", OCCREC_AUC_CUT), ""))
  shown <- cbind(tbl, ` ` = format(flag, justify = "left"))
  txt <- utils::capture.output(print(shown, digits = 3, row.names = FALSE))

  message(sprintf(
"
per-species occupancy recovery -- %s scenario, n = %d sites, S = %d species, %d chains x %d iter

  occ_true  fraction of sites where the species is truly present
  occ_est   fraction estimated present (posterior mean of the latent z)
  psi_cor   correlation of estimated psi with plogis(eta) across sites (null 0)
  auc       P(occupied site ranked above unoccupied site) (null 0.5)
  drv_cov   fraction of this species' B0 and B inside their 95%% interval (nominal 0.95)
  drv_err   mean |posterior mean - truth| for those same %d coefficients

%s
",
    sc$label, sc$n, sc$S, occrec_mcmc$nchain, occrec_mcmc$niter,
    1L + sc$ncov_psi, paste(txt, collapse = "\n")))
}

#' Emit what the assertions check, what they caught, and why they are shaped
#' the way they are.
occrec_show_checks <- function(which) {
  chk <- occrec_checks()
  lines <- vapply(which, function(nm) {
    f <- OCCREC_FLOORS[[nm]]
    floor_val <- if (is.null(f$floor)) chk$n_species_floor else f$floor
    op <- if (is.null(f$floor)) ">=" else ">"
    fmt <- if (is.null(f$floor)) "%s%9.0f  %s %-6.0f  swept %s\n      %s\n      %s"
           else "%s%9.3f  %s %-6.2f  swept %s\n      %s\n      %s"
    sprintf(fmt, formatC(nm, width = -24), chk[[nm]], op, floor_val,
            f$swept, f$what, f$why)
  }, character(1))

  message("\nassertions (floors set below the worst of eight seed sets; ",
          "'swept' is the range seen there)\n\n  ",
          paste(lines, collapse = "\n\n  "), "\n")
}

test_that("occupancy is recovered per species (tier 2)", {
  skip_on_cran()

  tbl <- occrec_table()
  chk <- occrec_checks()
  occrec_show(tbl)
  occrec_show_checks(c("mean_psi_cor", "n_species_ok", "mean_auc"))

  expect_equal(nrow(tbl), occrec_scenario()$S)
  # A NA here means a species was entirely occupied or entirely absent, which
  # makes its row uninterpretable rather than merely bad.
  expect_true(all(is.finite(tbl$psi_cor)))
  expect_true(all(is.finite(tbl$auc)))

  # The three checks printed above. Rationale lives in OCCREC_FLOORS so that
  # the run explains itself; keep it there rather than duplicating it here.
  expect_gt(chk$mean_psi_cor, OCCREC_FLOORS$mean_psi_cor$floor)
  expect_gte(chk$n_species_ok, chk$n_species_floor)
  expect_gt(chk$mean_auc, OCCREC_FLOORS$mean_auc$floor)
})

test_that("the coefficients driving occupancy are recovered per species (tier 2)", {
  skip_on_cran()

  tbl <- occrec_table()
  drv <- occrec_drivers()
  chk <- occrec_checks()
  occrec_show_checks(c("drv_cov", "drv_cor"))

  # B0 plus ncov_psi slopes per species, so 3 elements here. Three Bernoulli
  # draws cannot support a per-species coverage assertion: the only values
  # available are 0, 1/3, 2/3 and 1, and the observed per-species minimum over
  # eight seeds was 1/3 on a model that is covering at 0.87 overall. Likewise
  # the per-species correlation of estimate against truth is computed on three
  # points and went as low as -0.877. Per species these columns are for
  # reading; the assertions below are pooled.
  expect_equal(nrow(drv), occrec_scenario()$S * (1L + occrec_scenario()$ncov_psi))
  expect_true(all(is.finite(tbl$drv_err)))

  # The two checks printed above; see OCCREC_FLOORS for what each one catches.
  expect_gt(chk$drv_cov, OCCREC_FLOORS$drv_cov$floor)
  expect_gt(chk$drv_cor, OCCREC_FLOORS$drv_cor$floor)
})
