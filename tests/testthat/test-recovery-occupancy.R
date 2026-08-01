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

  mean_psi_cor = list(
    floor = 0.30, digits = 3, swept = "0.498 to 0.753",
    what = "For each species, how closely the estimated chance of being present tracks the true chance across the 100 sites, then averaged over the 10 species. A value of 0 would mean the estimates carry no information about where the species really is.",
    why  = "This is the steadiest number in the file. Refitting the same dataset four times moved it by only 0.017, so when it drops, the cause is the model and not chance. It is the check that notices the model getting worse at all ten species at once, which is the kind of failure that no single species would make obvious."),

  n_species_ok = list(
    floor = NULL, digits = 0, swept = "7 to 10 out of 10",
    what = "How many of the 10 species had that same tracking measure reach 0.30 or better. The floor is half the species.",
    why  = "An average can stay healthy while one or two species are recovered badly, because the good species pull it back up. This check looks at the species one at a time so that kind of failure cannot hide. It counts species rather than demanding that every species pass, because a single species is a noisy thing to measure: across the eight test datasets the worst species ranged from -0.401 to 0.913, so any threshold strict enough to be meaningful would fail on ordinary data."),

  mean_auc = list(
    floor = 0.60, digits = 3, swept = "0.748 to 0.855",
    what = "The chance that the model rates a site where the species really is present above a site where it really is absent, averaged over the 10 species. A coin flip would score 0.5.",
    why  = "It measures much the same thing as the tracking check above, and the two agreed closely when measured over 80 species (0.94), so it is not doing independent work. It is here because it is easy to read on its own, and because it would fall to a coin flip outright if sites or species were ever accidentally shuffled out of order."),

  drv_cov = list(
    floor = 0.55, digits = 3, swept = "0.800 to 1.000",
    what = "Each species has 3 coefficients that determine where it is likely to be found: a baseline and two environmental effects. The model reports a range of plausible values for each. This is the share of all 30 such ranges that actually contain the true value.",
    why  = "The ranges are built to contain the truth 95% of the time, so this should sit near 0.95. The floor is much lower because 30 values from one fit is a small sample and will bounce around. It is still worth checking, because if the true values were ever lined up against the wrong coefficients this number would collapse towards 0 rather than merely drift."),

  drv_cor = list(
    floor = 0.30, digits = 3, swept = "0.647 to 0.900",
    what = "How closely the model's best guess at those same 30 coefficients tracks their true values.",
    why  = "The check above asks whether the truth falls inside the model's stated range. This one asks whether the model's actual answer is close. Both are needed: a range wide enough to contain almost anything will pass the first check while telling you nothing, and only this one would notice.")
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

OCCREC_WIDTH <- 78L

#' Wrap a paragraph to the report width at a fixed indent.
occrec_para <- function(txt, indent = 6L) {
  paste(strwrap(txt, width = OCCREC_WIDTH, indent = indent, exdent = indent),
        collapse = "\n")
}

#' A wrapped paragraph with a short label in a hanging left column.
#'
#' Wrapped first and labelled after, because strwrap() collapses runs of
#' spaces and would eat the padding if the label went in ahead of it.
occrec_labelled <- function(label, txt, lab_w = 10L, indent = 2L) {
  body <- strwrap(txt, width = OCCREC_WIDTH - lab_w - indent)
  out <- paste0(strrep(" ", indent), formatC(label, width = -lab_w), body[1])
  # Guard the one-line case: paste0() recycles a zero-length vector to "",
  # which would emit a line of nothing but padding.
  if (length(body) > 1L) {
    out <- c(out, paste0(strrep(" ", lab_w + indent), body[-1]))
  }
  paste(out, collapse = "\n")
}

#' Emit the table so it is visible when the suite runs.
#'
#' message() rather than print(): it goes to stderr and survives testthat's
#' reporters, which is the point of the file. The legend and the flags are here
#' so that a row can be read without opening this source file.
occrec_show <- function(tbl) {
  sc <- occrec_scenario()
  n_drv <- 1L + sc$ncov_psi

  flag <- ifelse(tbl$psi_cor < OCCREC_SPECIES_CUT,
                 sprintf("<- psi_cor below %.2f", OCCREC_SPECIES_CUT),
                 ifelse(tbl$auc < OCCREC_AUC_CUT,
                        sprintf("<- auc below %.2f", OCCREC_AUC_CUT), ""))
  shown <- cbind(tbl, ` ` = format(flag, justify = "left"))
  txt <- sub("\\s+$", "",
             utils::capture.output(print(shown, digits = 3, row.names = FALSE)))

  legend <- c(
    occrec_labelled("occ_true",
      sprintf("Fraction of the %d sites where the species really is present.", sc$n)),
    occrec_labelled("occ_est",
      "Fraction of sites where the model thinks it is present."),
    occrec_labelled("psi_cor",
      "How closely the estimated chance of being present tracks the true chance, across sites. Runs from -1 to 1; 0 means the estimates say nothing about where the species is."),
    occrec_labelled("auc",
      "Chance that the model rates a site where the species really is present above one where it really is absent. A coin flip scores 0.5."),
    occrec_labelled("drv_cov",
      sprintf("This species has %d coefficients setting where it is likely to be found: a baseline and %d environmental effects. The model gives a plausible range for each. This is the share of those %d ranges that contain the true value; it should be near 0.95.",
              n_drv, sc$ncov_psi, n_drv)),
    occrec_labelled("drv_err",
      sprintf("Average distance between the model's best guess and the true value, over those same %d coefficients. Smaller is better; it is on the scale of the coefficients themselves, so there is no fixed good value.",
              n_drv)))

  message(sprintf(
    "\nper-species occupancy recovery\n\n%s\n\n%s\n\n%s\n",
    occrec_para(sprintf(
      "One row per species from a single fit of the '%s' scenario: %d sites, %d species, %d chains of %d iterations after %d burn-in. A row is flagged when it falls below a threshold one of the checks below uses.",
      sc$label, sc$n, sc$S, occrec_mcmc$nchain, occrec_mcmc$niter,
      occrec_mcmc$nburn), indent = 2L),
    paste(legend, collapse = "\n"),
    paste(txt, collapse = "\n")))
}

#' Emit what each assertion measures, what it produced, and why it is shaped
#' the way it is.
occrec_show_checks <- function(which) {
  chk <- occrec_checks()

  blocks <- vapply(which, function(nm) {
    f <- OCCREC_FLOORS[[nm]]
    is_count <- is.null(f$floor)
    floor_val <- if (is_count) chk$n_species_floor else f$floor
    rule <- sprintf(if (is_count) "must be at least %.0f" else "must exceed %.2f",
                    floor_val)
    head <- sprintf("  %s%s   %-20s   swept %s",
                    formatC(nm, width = -16),
                    formatC(chk[[nm]], format = "f", digits = f$digits, width = 7),
                    rule, f$swept)
    paste(head,
          occrec_para(paste("What it measures:", f$what)),
          occrec_para(paste("Why it is here:", f$why)),
          sep = "\n\n")
  }, character(1))

  message(sprintf("\nchecks\n\n%s\n\n%s\n",
    occrec_para(paste(
      "Each check below shows the value this run produced, the floor it has to",
      "clear, and a range labelled 'swept'. Swept means the whole simulate-and-fit",
      "was repeated on eight separately generated datasets before these floors",
      "were chosen, and that range is the smallest and largest value seen across",
      "those eight runs. Every floor sits below the smallest of them, so passing",
      "is not luck and failing is not ordinary random variation. If the model or",
      "the settings above change, sweep again and reset the floors rather than",
      "nudging one until the suite goes green."), indent = 2L),
    paste(blocks, collapse = "\n\n")))
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
