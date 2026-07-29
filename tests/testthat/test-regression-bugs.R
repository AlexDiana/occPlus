# Tier 1: one test per entry in TODO.md's "Fixed bugs" record.
#
# Every one of those fixes was verified by reading code or by a one-off manual
# run; none had a test. These are the tests. Each is named for the entry it
# guards, so a failure points straight at the history.
#
# Structural assertions only -- see helper-fixtures.R for why.

test_that("Fixed bugs 22: a fit with no spatial covariates runs", {
  # Regression introduced by enabling sample_ls() (Fixed bugs 12) and fixed in
  # d6b70b1. update_jSDMcoef() gated the sample_ls() call on ncol(Bs), the
  # *species* count, which is never 0 -- so it ran even with no spatial field
  # and read all-NA grids from precomputeSORmatrices(), dying with "missing
  # value where TRUE/FALSE needed". The guard is now nrow(Bs).
  #
  # This broke the pure-JSDM and classical-occupancy configurations the
  # listserv announcement advertises, and was masked in the vignette because
  # that passes spatCovariates.
  expect_no_error(fit_fixture(simulate_fixture(model = "two_stage"),
                              spatial = FALSE))
})

test_that("Fixed bugs 8: tau is sampled for continuous models", {
  # tau_output was returned entirely NA. Note the fix is only observable for
  # model == "continuous": sample_tau() is called only on that branch
  # (R/jsdmfun.R:1185, R/runOccJSDM.R:1130). For every other model type tau is
  # fixed at 1 and tau_output is NA *by design*, so asserting non-NA on a
  # two-stage fit would be wrong.
  tau_cont <- fixture_continuous()$results_output$jsdm_output$tau_output
  expect_false(all(is.na(tau_cont)))
  expect_true(any(is.finite(tau_cont)))

  # The complementary half: still NA where tau is not a free parameter.
  tau_ts <- fixture_twostage()$results_output$jsdm_output$tau_output
  expect_true(all(is.na(tau_ts)))
})

test_that("Fixed bugs 21: plotFPTPStage2Rates() honours primerName", {
  # p_output/q_output are [primer, species, iteration, chain]; the function
  # summarised them with apply(p_output, 2, ...), pooling over primer, so the
  # plot was identical whatever primerName was passed.
  fit <- fixture_twostage()
  primers <- fit$infos$primerNames
  expect_gt(length(primers), 1)

  by_primer <- lapply(primers, function(p) {
    d <- plotFPTPStage2Rates(fit, primerName = p)$data
    d[order(d$Species), c("Species", "p1", "p2", "q1", "q2")]
  })

  # Different primers must give different intervals. Pooling made these equal.
  expect_false(isTRUE(all.equal(by_primer[[1]], by_primer[[2]])))

  # And each primer's interval must equal the quantiles of that primer's own
  # slice -- "different" alone would not rule out subsetting the wrong margin.
  p_out <- fit$results_output$p_output
  for (i in seq_along(primers)) {
    manual <- t(apply(p_out[i, , , , drop = FALSE], 2,
                      function(x) stats::quantile(x, c(0.025, 0.975))))
    got <- as.matrix(by_primer[[i]][, c("p1", "p2")])
    dimnames(manual) <- dimnames(got) <- NULL
    expect_equal(got, manual)
  }
})

test_that("Fixed bugs 20 and 23: computeSpeciesDetected uses all chains", {
  # Two bugs in one function. Fixed bugs 20: B was referenced but never
  # assigned, so it errored with "object 'B' not found". Fixed bugs 23: the
  # draws used were a prefix of the flattened iter-within-chain array, i.e.
  # chain 1 only whenever per-chain niter >= 500.
  #
  # Tested behaviourally with synthetic posteriors: chain 1 says "never
  # collected" and chain 2 says "always collected". If only chain 1 were used
  # the detection count would be ~0; spanning both chains it must be well
  # above 0.
  S <- 3L; niter <- 40L; nchain <- 2L; P <- 1L

  beta_theta <- array(0, dim = c(1L, S, niter, nchain))
  beta_theta[, , , 1] <- -20   # logistic(-20) ~ 0: never collected
  beta_theta[, , , 2] <-  20   # logistic( 20) ~ 1: always collected

  p_out <- array(0.99, dim = c(P, S, niter, nchain))

  out <- occJSDM:::computeSpeciesDetected(beta_theta, p_out,
                                          M = 2L, K = 2L,
                                          primer = 0, alpha = 0.95)

  expect_false(anyNA(out))                    # guards Fixed bugs 20
  expect_gt(max(out), 0)                      # chain 2 reached: Fixed bugs 23

  # Sanity check the test itself: with *both* chains saying "never collected",
  # the count must be 0. Without this, expect_gt above could pass for the
  # wrong reason.
  beta_never <- array(-20, dim = c(1L, S, niter, nchain))
  out_never <- occJSDM:::computeSpeciesDetected(beta_never, p_out,
                                                M = 2L, K = 2L,
                                                primer = 0, alpha = 0.95)
  expect_equal(max(out_never), 0)
})

test_that("Fixed bugs 9: the WAIC iteration counter agrees with niter", {
  # The WAIC running means divided by the raw MCMC iteration index rather than
  # the number of WAIC accumulations. The fix added a dedicated
  # currentWAICiter plus an assertion, stop("Current wAIC iter wrong"), which
  # trips if the two disagree. Thinning is the case most likely to expose an
  # off-by-one, so exercise nthin > 1.
  expect_no_error(
    fit_fixture(simulate_fixture(model = "two_stage"),
                MCMCparams = list(nchain = 2, nburn = 20, niter = 10, nthin = 2))
  )
})

test_that("Fixed bugs 1 and 11: row order of the input data does not matter", {
  # data$info was re-sorted internally but data$OTU was not, so responses were
  # silently permuted against covariates; and collection covariates were
  # grouped by Sample alone rather than by (Site, Sample).
  #
  # Asserted on the prepared design matrices, not on the posterior: covariate
  # preparation is deterministic, whereas fitted values are not reproducible
  # across platforms while group A item 2 is open (see helper-fixtures.R).
  sim <- simulate_fixture(model = "two_stage", P = 2L)

  shuffled <- sim
  set.seed(99)
  perm <- sample.int(nrow(sim$data_list$info))
  shuffled$data_list$info <- sim$data_list$info[perm, , drop = FALSE]
  shuffled$data_list$OTU  <- sim$data_list$OTU[perm, , drop = FALSE]

  fit_sorted   <- fit_fixture(sim)
  fit_shuffled <- fit_fixture(shuffled)

  expect_equal(fit_shuffled$X_psi,   fit_sorted$X_psi)
  expect_equal(fit_shuffled$X_theta, fit_sorted$X_theta)

  # The internal data_info must be canonical, i.e. independent of input order.
  #
  # ignore_attr = TRUE covers row.names, and that exclusion was checked rather
  # than assumed -- asserting on them made this test fail against Alex's
  # 42198d9 when nothing behavioural had changed. What was verified:
  #
  #   * The (Site, Sample, Primer) sequence, the covariate columns and both
  #     design matrices are identical between sorted and shuffled input.
  #   * 150 of 320 row.names do differ, but only by the ordering of the K
  #     replicates *within* a (Site, Sample, Primer) group -- the same set of
  #     names, permuted inside each group.
  #   * Critically, row.names are ORDER-DEPENDENT BUT PAIRING-PRESERVING.
  #     Re-run with meaningful input names (SAMP001, SAMP002, ...) instead of
  #     R's auto-disambiguated "38", "38.1", both sorted and shuffled input
  #     kept every name attached to its own OTU counts. So a user who looks up
  #     by name -- data_info["SAMP042", ], OTU["SAMP042", ] -- gets the right
  #     row either way. Only the order in which rows appear differs.
  #   * Nothing in the package reads these row.names. The only rownames()
  #     reads in R/ are on data$traits (species lookup); elsewhere they are
  #     assigned, never consulted.
  #
  # So this is not merely "cosmetic": the pairing that a user could rely on is
  # intact, which is the property worth protecting. Do NOT re-tighten this to
  # compare row.names -- it would fail on within-group replicate ordering,
  # which is exchangeable (see the OTU note below) and carries no information.
  #
  # The case that WOULD break is positional joining, e.g.
  # cbind(fit$infos$data_info, my_metadata) assuming the caller's row order.
  # That breaks for sorted input too, since the fit always sorts internally,
  # so it is neither caused by nor detectable from this test.
  expect_equal(fit_shuffled$infos$data_info, fit_sorted$infos$data_info,
               ignore_attr = TRUE)

  # The response must travel with its covariates. Note this is deliberately
  # NOT expect_equal() on infos$OTU row-for-row: (Site, Sample, Primer) is not
  # a unique key -- the K PCR replicates within a group share it -- and the
  # internal sort does not fix their relative order, so shuffling the input
  # permutes replicates *within* a group. Those replicates are exchangeable
  # (conditionally iid given w), so that permutation is meaningless, and
  # asserting exact row equality would fail on a difference that does not
  # exist statistically. Verified: the within-group rows are the same multiset
  # either way.
  #
  # What must hold is that no row crosses a group boundary -- that is what the
  # original bug did.
  grp <- with(fit_sorted$infos$data_info,
              interaction(Site, Sample, Primer, drop = TRUE))
  expect_equal(rowsum(fit_shuffled$infos$OTU, grp),
               rowsum(fit_sorted$infos$OTU,   grp))
})

test_that("Fixed bugs 12: the GP length-scale sampler is reached", {
  # sample_ls() was wrapped in if(F){...}, so it never executed and idx_ls
  # stayed at its hard-coded starting value of 3. The fix was to enable the
  # call, so that is exactly what this asserts.
  #
  # Deliberately NOT asserting that idx_ls_output varies. Two reasons, both
  # measured rather than assumed:
  #
  #   (a) It would be stochastic -- acceptance is a Metropolis step, so a
  #       short chain may not move (PLAN.md 5.2).
  #   (b) More importantly it is currently constant for a *different* reason:
  #       the chain walks to idx_ls = 10, the top of l_s_grid, and stays
  #       there, for every true l_s tried (0.30, 0.15, 0.08) across seeds. So
  #       the sampler runs and accepts moves, but the length-scale is not
  #       recovered -- see the note filed against sigma_s being hard-coded to
  #       1 at the sample_ls() call site (R/jsdmfun.R:1254) while the
  #       simulator generates the field with sigma_s = 0.5. Asserting
  #       "varies" would encode that defect; asserting "is reached" stays
  #       correct once it is fixed.
  called <- FALSE
  suppressMessages(trace("sample_ls", where = asNamespace("occJSDM"),
                         tracer = function() called <<- TRUE, print = FALSE))
  on.exit(try(suppressMessages(
    untrace("sample_ls", where = asNamespace("occJSDM"))), silent = TRUE),
    add = TRUE)

  invisible(fit_fixture(simulate_fixture(model = "two_stage"), spatial = TRUE))
  expect_true(called)
})

# --- Alex's fixes of 28-29 July (Fixed bugs 24-27) ------------------------

test_that("Fixed bugs 24: sigma_h is sampled, not stuck at 1", {
  # sigma_h was initialised to 1 and passed through update_jSDMcoef()
  # unchanged, so sigmah_output was constant. A new sample_sigmah() is now
  # called each iteration. Inert for the fitted model -- U at training sites
  # is drawn under a hard-coded unit-variance prior regardless -- but it is
  # the factor-score SD predictNewSites() uses at NEW sites.
  sh <- fixture_twostage()$results_output$jsdm_output$sigmah_output
  expect_false(is.null(sh))
  expect_gt(length(unique(as.vector(sh))), 1)
  expect_true(all(sh > 0, na.rm = TRUE))
})

test_that("Fixed bugs 25: collection-covariate priors are centred on 0", {
  # b_betatheta was rep(1, ncov_theta) with only element 1 (the intercept)
  # overwritten, so every covariate slope carried a prior centred on +1 and
  # was dragged upward. Asserted on the prior itself rather than on recovery:
  # recovery is stochastic and belongs in tiers 2-3, whereas the prior is a
  # fact about the code and is what actually regressed.
  f <- test_path("..", "..", "R", "runOccJSDM.R")
  skip_if_not(file.exists(f), "package source not reachable (installed-package check)")
  line <- grep("b_betatheta\\s*<-\\s*rep\\(", readLines(f, warn = FALSE), value = TRUE)
  expect_length(line, 1)
  expect_match(line, "rep\\(0,")     # NOT rep(1, ...)
})

test_that("Fixed bugs 26: a fixed seed reproduces a fit", {
  # Fixed bugs 26 replaced randinvg()'s use of R's global RNG inside an OpenMP
  # loop with thread_local engines, closing the data race. But those engines
  # never read R's RNG state, so set.seed() still did not control the sampler
  # (get_rng() was seeded from the literal 12345, mvrnormArmaQuick_TS() from
  # std::random_device). runOccJSDM() now draws a base seed from R via
  # setOccJSDMSeed(); see src/rng.h.
  #
  # This is the one place tier 1 asserts numeric equality. The blanket ban in
  # PLAN.md 5.1 exists *because* the sampler used to be irreproducible, so the
  # test that it no longer is must necessarily be an equality test.
  #
  # Scope: reproducibility holds for a given thread count. Each thread derives
  # its own stream from the base seed, so if this package is ever built with
  # OpenMP actually enabled, changing the thread count changes which stream
  # produces which element. That is inherent to per-thread streams, not a bug,
  # but it does mean "same seed" is only half of the contract -- see the note
  # in src/rng.h. (As built here, SHLIB_OPENMP_CXXFLAGS is empty and the
  # pragmas compile to nothing, so the point is currently moot.)
  sim <- simulate_fixture(model = "two_stage")
  mc <- list(nchain = 2, nburn = 20, niter = 20, nthin = 1)

  set.seed(4242); a <- fit_fixture(sim, MCMCparams = mc)
  set.seed(4242); b <- fit_fixture(sim, MCMCparams = mc)
  expect_equal(a$results_output$jsdm_output$B0_output,
               b$results_output$jsdm_output$B0_output)
  expect_equal(a$results_output$p_output, b$results_output$p_output)

  # A different seed must give a different answer, or the test above would
  # also pass if the sampler had simply stopped being random.
  set.seed(9999); d <- fit_fixture(sim, MCMCparams = mc)
  expect_false(isTRUE(all.equal(a$results_output$jsdm_output$B0_output,
                               d$results_output$jsdm_output$B0_output)))
})

test_that("consecutive fits under one seed are independent", {
  # The base seed is drawn with sample.int() rather than being a constant, so
  # R's RNG advances between fits. Without that, every fit in a session would
  # replay the same stream -- which is what the old literal 12345 did across
  # *processes*, correlating the simulation study's replicates (PLAN.md 8).
  sim <- simulate_fixture(model = "two_stage")
  mc <- list(nchain = 2, nburn = 20, niter = 20, nthin = 1)

  set.seed(777)
  f1 <- fit_fixture(sim, MCMCparams = mc)
  f2 <- fit_fixture(sim, MCMCparams = mc)
  expect_false(isTRUE(all.equal(f1$results_output$jsdm_output$B0_output,
                               f2$results_output$jsdm_output$B0_output)))
})

test_that("Fixed bugs 27: listPriors actually reaches the Stage 2 priors", {
  # a_p/b_p/a_q/b_q were documented as settable but hard-coded -- the same
  # failure as prior_beta_psi (Fixed bugs 14), which is why this is worth a
  # test rather than a read-through. Two documented-but-ignored parameters
  # have already shipped.
  #
  # Beta(1, 50) against the default Beta(5, 1) is a large enough shift that
  # the posterior must move even on a short chain, so this does not depend on
  # numeric reproducibility.
  sim <- simulate_fixture(model = "two_stage", P = 2L)
  mc <- list(nchain = 2, nburn = 40, niter = 40, nthin = 1)

  fit_default <- fit_fixture(sim, MCMCparams = mc)
  fit_low <- suppressMessages(suppressWarnings(
    runOccJSDM(sim$data_list,
               listParams = list(n_factors = 2L, n_supportpoints = FIXTURE_KNOTS),
               occCovariates = fixture_occ_covariates(),
               spatCovariates = fixture_spat_covariates(),
               listPriors = list(a_p = 1, b_p = 50),
               MCMCparams = mc)))

  expect_gt(mean(fit_default$results_output$p_output, na.rm = TRUE),
            mean(fit_low$results_output$p_output, na.rm = TRUE))
})

test_that("Fixed bugs 31: plotCollectionRates() runs and labels the right bars", {
  # plotCollectionRates() errored with "object 'Min' not found" for every
  # input. plotSpeciesRates() had been extracted as a shared helper and never
  # wired up: it read columns Min/Max while its only caller passed the
  # "2.5%"/"97.5%" that quantile() produces, it filtered on a Species column
  # the caller never created, and it referenced speciesNames as a free variable
  # that exists in neither its arguments nor the namespace. Three independent
  # breakages in one call path, so nothing had ever exercised it.
  #
  # Asserting the label-to-value pairing, not just absence of error: the helper
  # also ordered on the full species set while filtering to a subset, which
  # silently mismatches names to bars. That is the defect the other three rate
  # plots still have (group B item 1), so this test is the template for those.
  fit <- fixture_twostage()
  nm <- fit$infos$speciesNames

  expect_no_error(plotCollectionRates(fit))

  idx <- c(3L, 1L)
  d <- plotCollectionRates(fit, idx_species = idx)$data
  expect_equal(nrow(d), length(idx))
  expect_equal(as.character(d$Species), nm[idx])
  expect_true(all(d$Min >= 0 & d$Max <= 1))
  expect_true(all(d$Min <= d$Max))

  # Values must follow the species, not the row position.
  d_all <- plotCollectionRates(fit)$data
  expect_equal(d$Min, d_all$Min[idx])
})

test_that("plotCollectionRates() rejects a model with no collection stage", {
  # It used to fail with a subscript error out of beta_theta_output; a JSDM-only
  # fit has no Stage 1 at all, so say so.
  expect_error(plotCollectionRates(fixture_continuous()),
               "collection-stage")
})

test_that("listPriors$b_betatheta_slope_var actually reaches the sampler", {
  # Added 29 July 2026 to test whether beta_theta's slope prior variance was
  # the cause of coverage worsening with M (PLAN.md 13.7-13.9). It was not --
  # tightening it moved coverage by noise -- but the override itself is a
  # real, permanent capability (same pattern as Fixed bugs 27's a_p/b_p), so
  # it gets the same "does it actually move the posterior" test they got.
  #
  # fit_fixture() never passes collCovariates, so it fits an intercept-only
  # Stage 1 even though simulate_fixture()'s truth has a real slope
  # covariate (beta_theta_true is 2 x S; the fitted beta_theta_output would
  # be 1 x S x niter x nchain without collCovariates). b_betatheta_slope_var
  # only has anything to act on if a slope row exists, so this test builds
  # its own runOccJSDM() calls with collCovariates = "X_theta", matching
  # what the simulation study's simstudy_fit() does.
  sim <- simulate_fixture(model = "two_stage", P = 2L)
  mc <- list(nchain = 2, nburn = 40, niter = 40, nthin = 1)
  args <- list(sim$data_list,
              listParams = list(n_factors = 2L, n_supportpoints = FIXTURE_KNOTS),
              occCovariates = fixture_occ_covariates(),
              collCovariates = "X_theta",
              spatCovariates = fixture_spat_covariates(),
              MCMCparams = mc)

  fit_default <- suppressMessages(suppressWarnings(do.call(runOccJSDM, args)))
  fit_tight <- suppressMessages(suppressWarnings(do.call(
    runOccJSDM, utils::modifyList(args, list(listPriors = list(b_betatheta_slope_var = 0.05))))))

  expect_equal(dim(fit_default$results_output$beta_theta_output)[1], 2L)

  slope_default <- fit_default$results_output$beta_theta_output[2, , , , drop = FALSE]
  slope_tight <- fit_tight$results_output$beta_theta_output[2, , , , drop = FALSE]
  expect_lt(sd(slope_tight), sd(slope_default))
})
