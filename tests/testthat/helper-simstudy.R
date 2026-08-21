# Shared machinery for the simulation study. See dev/simstudy/PLAN.md.
#
# Used by tier 2 (test-recovery.R, skip_on_cran) and tier 3
# (test-coverage-study.R, env-gated), and by dev/simstudy/run_study.R. Nothing
# here runs during an ordinary check(): defining functions is free.
#
# The two seams that let SBC be added later without restructuring (PLAN.md 7):
#
#   draw_truth(scenario, seed)      -- fixed values now; prior draws for SBC
#   statistic_coverage(draws, truth) -- CI indicator now; statistic_rank() for SBC
#
# Everything expensive (grid, simulate/fit loop, aggregation) is shared.

# --- Scenario grid (PLAN.md 6.1) -----------------------------------------

SIMSTUDY_BASE <- list(
  model = "two_stage",
  n = 100L, S = 10L, M = 2L, P = 2L, K = 3L,
  g = 2L, gt = 2L, d = 2L,
  # ds = 2 is now a choice, not a requirement. It used to be required: at
  # ds = 0 the simulator's cross-species spatial covariance collapsed to
  # jitter, giving sd(spatField) = 0.0019 against ~1.0 at ds = 2, so every
  # spatial cell would have been testing against a null field. Fixed in
  # 42198d9 (TODO.md Fixed bugs 30); re-measured 29 July at seed 42, ds = 0
  # now gives 0.598. Left at 2 so the grid keeps a clearly non-degenerate
  # field and stays comparable with the 28 July run.
  ds = 2L,
  ncov_psi = 2L, ncov_theta = 1L,
  useSpatField = TRUE,
  n_supportpoints = 20L,   # pinned, never left to the default (PLAN.md 5.4)
  fit_d = NULL,            # NULL = fit at the simulated d
  fit_traits = TRUE,       # FALSE drops data$traits before fitting
  p_range = c(0.3, 0.6),
  theta_baseline_range = c(0.2, 0.5)
)

#' The 10 scenario cells.
#'
#' One-factor-at-a-time from SIMSTUDY_BASE, not factorial. Cell 1 doubles as
#' the traits x spatial interaction, since both write into the same linear
#' predictor via computePsiCoef().
simstudy_scenarios <- function() {
  mk <- function(label, ...) {
    spec <- utils::modifyList(SIMSTUDY_BASE, list(...))
    spec$label <- label
    spec
  }
  list(
    mk("base"),
    # Traits are dropped at *fit* time, not simulated away: simulateOccJSDMData()
    # errors with g = 0 ("length of 'dimnames' [2] not equal to array extent"),
    # so a genuinely trait-free dataset cannot be generated. Verified.
    mk("spatial_isolated",  fit_traits = FALSE),
    mk("traits_isolated",   useSpatField = FALSE, n_supportpoints = NULL),
    mk("primers_3",         P = 3L),
    # Low information: weak detection AND fewer sites. n_supportpoints is set
    # explicitly; the default no longer crashes below 31 sites (TODO.md Fixed
    # bugs 29) but the study pins it regardless (PLAN.md 5.4).
    mk("low_information",   n = 40L, n_supportpoints = 8L,
                            p_range = c(0.1, 0.3),
                            theta_baseline_range = c(0.05, 0.2)),
    mk("d_underfit",        d = 4L, fit_d = 2L),
    mk("d_overfit",         d = 2L, fit_d = 4L),
    # Trimmed from S = 30 to S = 20 on 27 July. At S = 30 this cell cost
    # 1.75x base -- the most expensive in the grid by a wide margin -- and the
    # study is run on a fanless machine that throttles under sustained load.
    # S = 20 is still double the base S = 10, so the cell still answers "does
    # calibration hold as species count grows"; it just does so at 1.2x base
    # instead of 1.75x. Relabelled so the name does not lie about S.
    mk("species_20",        S = 20L),
    mk("occupancy",         model = "occupancy"),
    mk("binary",            model = "binary")
  ,

    # --- Community-feature contrasts (PLAN.md 18) -----------------------
    #
    # These answer "how well does occJSDM handle a community WITH this
    # feature against one WITHOUT it", which the ten cells above cannot: they
    # each draw independent truths, so any difference between them carries the
    # variance of two separate simulated worlds on top of the effect.
    #
    # All five set `seed_label = "base"`, so they pair against the existing
    # `base` cell without altering it -- `base` has no seed_label, so its key
    # is already the string "base" and every historical result stays
    # comparable.
    #
    # The pairing is exact for any arm holding P and S, because draw_truth()
    # consumes RNG sized only by P*S (p, q) and S (theta0, theta_baseline).
    # All five hold both, so all five draw a bit-identical truth. What differs
    # per arm is how much of that truth remains comparable element-wise:
    # detcov_none drops beta_theta's slope rows, and env_none/env_rich change
    # B's dimension, so those blocks cannot be differenced element-wise even
    # though the detection-side truths match exactly. PLAN.md 18.4 has the
    # table.
    #
    # NOT production cells by default, same as the M ladder: name them.

    # The spatial contrast, and the cell whose name matches what it does.
    # `traits_isolated` above is the arm that actually turns the spatial field
    # off; its name says the opposite. Renaming it would break comparability
    # with runs on disk, so this cell is added rather than that one renamed.
    mk("no_spatial",    useSpatField = FALSE, n_supportpoints = NULL,
                        seed_label = "base"),

    # Environmental covariate count. Never varied anywhere before: ncov_psi is
    # 2 in all 21 existing cells, so nothing in the study has ever spoken to
    # how covariate richness affects recovery.
    #
    # env_none needs the occCov guard in simstudy_fit() -- seq_len(0) yields
    # character(0), which is not what runOccJSDM() wants. Whether the
    # *simulator* tolerates ncov_psi = 0 is unverified; g = 0 already errors
    # outright, so treat the first run of this cell as a test of the cell.
    mk("env_none",      ncov_psi = 0L, seed_label = "base"),
    mk("env_rich",      ncov_psi = 4L, seed_label = "base"),

    # Detection covariates as a production contrast. `nocollcov` below sets the
    # same flag but pairs against itself, which answers the B0-bias question it
    # was built for and not this one.
    mk("detcov_none",   ncov_theta = 0L, seed_label = "base"),

    # Weak detection with n held at 100. This is what makes low_information
    # interpretable: that cell drops n to 40 AND weakens detection AND cuts
    # support points, so its damage is unattributable. Against this arm, the
    # remaining difference is sample size alone.
    mk("weak_detection", p_range = c(0.1, 0.3),
                         theta_baseline_range = c(0.05, 0.2),
                         seed_label = "base"),

    # Does more of the same data help? 3x the sites, everything else held.
    #
    # The prediction this tests, and it is not a simple one. Species-level
    # parameters (B0, B, G, the loadings) pool across sites, so they should
    # sharpen, and rare species should gain most: a species at 10% prevalence
    # has ~10 occupied sites at n = 100 and ~30 at n = 300, and the 10 August
    # run put the weakest recovery in exactly that band (psi_cor 0.366 in
    # (0.1, 0.2] against 0.655 mid-range). But the latent site scores U are
    # n x d, one row per site, each informed only by that site's own
    # observations -- so tripling n gives three times as many U rows each
    # estimated from as much data as before, and the site-specific half of
    # eta does not improve at all. Expect the systematic part to sharpen and
    # psi_cor to lift by less than the extra data would naively suggest.
    #
    # Calibration may well go the OTHER way, and that is the more interesting
    # half. The M ladder already ran this experiment on a different axis and
    # it went backwards: beta_theta coverage falls 0.747 -> 0.579 as M goes
    # 2 -> 20, and q falls 0.945 -> 0.614 as K goes 3 -> 30, both with flat
    # bias. That is intervals shrinking around a bias that is not shrinking.
    # Whether n behaves like M and K is unknown -- they enter the detection
    # stage and n enters the occupancy stage -- but "more data fixes
    # calibration" is an assumption this model has already falsified twice.
    #
    # n_supportpoints stays at 20, deliberately, so the model specification is
    # held fixed and only the data grows. Note the cost: 300 sites over the
    # same domain with the same knot count is a relatively coarser spatial
    # approximation, so a null or negative result on the spatial-adjacent
    # blocks should be read with that in mind rather than as "more sites did
    # not help".
    mk("sites_300", n = 300L, seed_label = "base"),

    # Does fitting MORE latent factors than the data contains help the species
    # that need them? base with fit_d = 4 against its simulated d = 2, paired.
    #
    # `d_overfit` above already does this, but unpaired -- it draws its own
    # truths, so a difference against base carries the variance of two
    # simulated worlds. Measured that way on 10 August the aggregate looked
    # like nothing (0.622 against 0.627), but split by how latent-driven each
    # species is, the differences ran -0.018, -0.023, +0.026, +0.035 from the
    # most environment-driven quarter to the most latent-driven. A clean sign
    # flip in order across four bands is not what noise usually looks like,
    # but each difference is only 1.5 to 2 SE and the comparison is unpaired.
    #
    # This cell makes it exact: identical truths to base, differing only in
    # how many factors are fitted. Read it on the latent-share bands, not the
    # aggregate -- the aggregate is where the effect hides, because the gain
    # for latent-driven species and the loss for environment-driven ones are
    # of similar size and cancel.
    mk("fit_d4", fit_d = 4L, seed_label = "base"),

    # --- M ladder (PLAN.md 13) -----------------------------------------
    #
    # Not production cells: these answer whether group B items 4-6 are code
    # defects or Stage 1 under-identification. Select them explicitly with
    # --scenarios=M2,M5,M10,M20,K30.
    #
    # **They are NOT inert, and a previous version of this comment said they
    # were.** run_study.R takes this entire list whenever --scenarios is
    # empty, so a bare `Rscript run_study.R` runs every cell defined here, not
    # the 10 production ones -- many hours instead of 2 h 45 m. That misreading
    # cost a wasted launch on 2 August. The production grid must be named:
    #
    #   --scenarios=base,binary,d_overfit,d_underfit,low_information,
    #               occupancy,primers_3,spatial_isolated,species_20,
    #               traits_isolated
    #
    # `seed_label` makes all five share a truth at each replicate, so the
    # comparison across arms is paired -- see simstudy_seed_for().
    #
    # K30 is the control: identical row count to M20 (12,000) but spent on
    # PCR replicates instead of field samples. Without it, "everything
    # improved" cannot be told from "more data helps".
    mk("M2",  M = 2L,  K = 3L,  seed_label = "mladder"),
    mk("M5",  M = 5L,  K = 3L,  seed_label = "mladder"),
    mk("M10", M = 10L, K = 3L,  seed_label = "mladder"),
    mk("M20", M = 20L, K = 3L,  seed_label = "mladder"),
    mk("K30", M = 2L,  K = 30L, seed_label = "mladder"),

    # --- Tighter beta_theta prior at high M (PLAN.md 13.8) --------------
    #
    # The M ladder found beta_theta coverage falling *with* M (0.747 at M2
    # to 0.579 at M20) while its bias stayed small and flat -- the
    # signature of an overconfident interval, not insufficient data. The
    # prime suspect is B_betatheta's slope variance, hard-coded to 2
    # until runOccJSDM.R added a listPriors override for this test. These
    # arms repeat M10 and M20 with that variance tightened to 0.5 (SD
    # 0.71 against the default's 1.41).
    #
    # Same seed_label as the M10/M20 arms above: paired not only against
    # each other but against the *original* M10/M20 results already on
    # disk (dev/simstudy/results/simstudy-20260729-212525.rds), since
    # nothing about the truth-generation changes here.
    mk("M10_tightprior", M = 10L, K = 3L, seed_label = "mladder",
                         listPriors = list(b_betatheta_slope_var = 0.5)),
    mk("M20_tightprior", M = 20L, K = 3L, seed_label = "mladder",
                         listPriors = list(b_betatheta_slope_var = 0.5)),

    # --- Prior variance at M = 2 (PLAN.md 14) ---------------------------
    #
    # theta0 went from 0.938-0.959 pre-fix to 0.978-0.985 post-fix at the
    # SAME M = 2, so its overcoverage is not "M = 2 is too little data":
    # it would have overcovered before Alex's fixes too. theta0's own
    # prior is untouched, so the suspected route is coupling through
    # b_betatheta -> beta_theta -> collection probability -> latent w,
    # which sample_theta0() conditions on.
    #
    # The 13.8 arms moved theta0 toward nominal only slightly (M10 0.944
    # to 0.940, M20 0.952 to 0.942), as expected: those are the arms where
    # the data dominates the prior. M = 2 is where the prior should
    # dominate, and is the arm 13.8 did not cover.
    #
    # Two variances rather than one, so a null result cannot be blamed on
    # not having pushed hard enough.
    mk("M2_tight",  M = 2L, K = 3L, seed_label = "mladder",
                    listPriors = list(b_betatheta_slope_var = 0.5)),
    mk("M2_vtight", M = 2L, K = 3L, seed_label = "mladder",
                    listPriors = list(b_betatheta_slope_var = 0.1)),

    # --- No collection covariates (PLAN.md 15) --------------------------
    #
    # Alex's discriminator for the B0 bias, sharpened. B0 bias separates
    # perfectly on whether a cell estimates beta_theta: -0.125 to -1.056 in
    # the nine that do, +0.012 in `binary`, the one that does not. But
    # `binary` differs in many ways at once, and `occupancy` is already
    # one-stage yet still estimates beta_theta (bias -0.151), so neither
    # existing cell isolates the cause.
    #
    # This arm is `base` with ncov_theta = 0: the two-stage machinery, the
    # latent w/z and p/q/theta0 all stay, and only the beta_theta *slopes*
    # go. The intercept row, logit(theta_baseline), necessarily remains.
    #
    # Read it against B0 bias, not coverage: B0 covers at 0.94-0.96 in every
    # cell including `binary`, so coverage cannot distinguish anything here.
    mk("nocollcov", ncov_theta = 0L, seed_label = "nocollcov"),

    # --- Continuous, B0 with no beta_theta anywhere (PLAN.md 16) --------
    #
    # Alex's correction to what `nocollcov` was meant to test. Setting
    # ncov_theta = 0 removes the beta_theta *slopes* but necessarily keeps
    # its intercept row, logit(theta_baseline), and in an occupancy model
    # that intercept and B0 are both intercepts on the same linear
    # predictor chain: psi sets how often a site is occupied, theta sets
    # how often an occupied site yields a positive sample, and with M = 2
    # the data barely separates them. So `nocollcov`'s residual -0.0633
    # bias could be that confounding rather than a defect in B0 itself,
    # which is what he was actually proposing to rule out.
    #
    # model = "continuous" is the arm with no beta_theta at all: z is
    # drawn as Normal(eta, tau) and observed directly, so there is no
    # detection stage, no latent w, and B0 is estimated straight from the
    # Gaussian likelihood. `binary` already removes beta_theta too, but it
    # goes through the Polya-Gamma path; continuous is a different
    # likelihood and a different branch of the sampler, so agreement
    # between the two is worth more than either alone.
    #
    # Same caveat as `binary`, stated rather than glossed: this changes
    # many things at once, not one. It is a second independent reading of
    # "B0 with nothing confounded against it", not a controlled contrast.
    mk("continuous", model = "continuous"),

    # --- q cost-of-identifiability arms (TODO group B item 7, PLAN.md 21)
    #
    # The M ladder found q coverage falling 0.945 -> 0.614 as K rises 3 ->
    # 30. Alex sees no implementation defect in sample_pq_cpp(); the live
    # hypothesis is cost of identifiability: the informative Beta(1, 20)
    # prior (mean 0.0476) holds the posterior a fixed distance from the
    # truth, and sharpening intervals with K turn that fixed offset into
    # falling coverage.
    #
    # Step 0 confirmed both cheap signatures from the EXISTING M-ladder run:
    # q's bias is flat and tiny across K while width contracts ~4x, and the
    # per-species error correlates negatively with distance from the prior
    # mean (slope ~ -0.14, flat in K). What remains is the discriminating
    # test: does degradation scale with how far the truth sits from the
    # prior mean?
    #
    # Four cells: base settings, K = 3 vs K = 30, q drawn near vs far from
    # the prior mean. All four share seed_label "qtest", so p, theta0,
    # theta_baseline and every JSDM parameter are paired bit-identically
    # across cells; only the q truths differ (draw order is unchanged --
    # same P*S uniforms). Prediction under the hypothesis: qnear stays at
    # or near nominal coverage at both K; qfar degrades steeply with K,
    # worse than the observed 0.945 -> 0.614.
    mk("qnear_K3",  K = 3L,                        seed_label = "qtest"),
    mk("qnear_K30", K = 30L,                       seed_label = "qtest"),
    mk("qfar_K3",   K = 3L,  q_range = c(.15, .3), seed_label = "qtest"),
    mk("qfar_K30",  K = 30L, q_range = c(.15, .3), seed_label = "qtest"))
}

# --- Seam 1: how the true parameters are chosen --------------------------

#' Fix the generating parameters for one replicate.
#'
#' v1 draws the nuisance detection parameters from the scenario's stated
#' ranges and takes the structural parameters from the spec; the JSDM
#' coefficients themselves are then realised by simulateOccJSDMData().
#'
#' For SBC this is the function to replace: draw from the *model's own
#' priors* instead (all of which are proper -- see PLAN.md 7 for the list).
#' Nothing downstream needs to change.
draw_truth <- function(scenario, seed) {
  set.seed(seed)
  S <- scenario$S
  P <- scenario$P

  M_vec <- if (scenario$model %in% c("occupancy", "two_stage")) {
    rep(scenario$M, scenario$n)
  } else {
    rep(1L, scenario$n)
  }
  N <- sum(M_vec)
  K_vec <- if (scenario$model == "two_stage") {
    rep(scenario$K, N * P)
  } else {
    rep(1L, N * P)
  }

  list(
    model = scenario$model,
    datasettings = list(n = scenario$n, S = S, g = scenario$g,
                        M = M_vec, P = P, K = K_vec,
                        ncov_psi = scenario$ncov_psi,
                        ncov_theta = scenario$ncov_theta),
    jsdmParams = list(gt = scenario$gt, d = scenario$d, ds = scenario$ds,
                      # Near the InvGamma(10, 1) prior mean of sqrt(1/9); see
                      # the note above simstudy_param_blocks() on why a truth
                      # far from the prior makes this block uninterpretable.
                      sigma_b = 0.35, sigma_bs = 0.5, sigma_ts = 0.5,
                      sigma_h = 1, sigma_s = 0.5, l_s = 0.15,
                      tau = rep(1, S),
                      useSpatField = scenario$useSpatField),
    params = list(
      p = matrix(stats::runif(P * S, scenario$p_range[1],
                              scenario$p_range[2]), P, S),
    # q_range is optional and defaults to the historical c(0.01, 0.05).
    # Used by the q cost-of-identifiability arms (PLAN.md 21): those need a
    # truth drawn well clear of the Beta(1, 20) prior mean of 0.0476, to test
    # whether q's coverage loss with K scales with that distance.
    q = matrix(stats::runif(P * S,
                            if (!is.null(scenario$q_range)) scenario$q_range[1] else 0.01,
                            if (!is.null(scenario$q_range)) scenario$q_range[2] else 0.05),
               P, S),
      theta0 = stats::runif(S, 0.02, 0.1),
      theta_baseline = stats::runif(S, scenario$theta_baseline_range[1],
                                    scenario$theta_baseline_range[2])
    )
  )
}

# --- Simulate and fit ----------------------------------------------------

simstudy_simulate <- function(truth) {
  simulateOccJSDMData(truth$datasettings, truth$params, truth$jsdmParams,
                      model = truth$model)
}

#' Fit one simulated dataset.
#'
#' collCovariates is passed deliberately: without it the fit carries only an
#' intercept for theta, so beta_theta_output would be (1 x S) against a
#' (ncov_theta+1 x S) truth and the comparison would be impossible. Likewise
#' n_lattrait is pinned to the simulated gt, because the default
#' floor(sqrt(min(S, ncov_psi))) will not generally match.
simstudy_fit <- function(sim, scenario,
                         MCMCparams = list(nchain = 2, nburn = 1000,
                                           niter = 1000, nthin = 1)) {
  d_fit <- scenario$fit_d %||% scenario$d

  listParams <- list(n_factors = d_fit, n_lattrait = scenario$gt)
  spatCovariates <- NULL
  if (isTRUE(scenario$useSpatField)) {
    spatCovariates <- c("Xs.1", "Xs.2")
    listParams$n_supportpoints <- scenario$n_supportpoints
  }

  # ncov_psi = 0 gives seq_len(0), hence character(0), which is not the same
  # thing as "no covariates" to runOccJSDM(). Mirrors the collCov branch just
  # below, which already handles the analogous ncov_theta = 0 case. Needed by
  # the env_none cell (PLAN.md 18.6).
  occCov <- if (scenario$ncov_psi > 0)
    paste0("X_psi.EnvCov.", seq_len(scenario$ncov_psi)) else NULL
  # A scenario with ncov_theta = 0 has no X_theta column at all: the simulator
  # produces a 1 x S beta_theta (the intercept row only). Passing
  # collCovariates regardless errors with "Covariate names provided not in
  # data$info". Needed for the ncov_theta = 0 arm in PLAN.md 15.
  collCov <- if (scenario$model %in% c("occupancy", "two_stage") &&
                 scenario$ncov_theta > 0) "X_theta" else NULL

  dat <- sim$data_list
  if (!isTRUE(scenario$fit_traits)) dat$traits <- NULL

  suppressMessages(suppressWarnings(
    runOccJSDM(dat,
               listParams     = listParams,
               listPriors     = scenario$listPriors %||% list(),
               occCovariates  = occCov,
               collCovariates = collCov,
               spatCovariates = spatCovariates,
               MCMCparams     = MCMCparams)
  ))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# --- Seam 2: the per-parameter statistic ---------------------------------

#' Collapse a posterior array to a (parameter element) x (draw) matrix.
#'
#' Posterior blocks always end in (niter, nchain); the leading dimensions
#' index the parameter. Column-major flattening means as.vector(truth) lines
#' up with the rows without any transposing.
posterior_draw_matrix <- function(arr) {
  dn <- dim(arr)
  k <- length(dn)
  stopifnot(k >= 2)
  n_leading <- if (k == 2) 1L else prod(dn[seq_len(k - 2)])
  matrix(as.vector(arr), nrow = n_leading)
}

#' Coverage statistic (v1).
statistic_coverage <- function(draws, truth, level = 0.95) {
  a <- (1 - level) / 2
  qs <- stats::quantile(draws, c(a, 1 - a), names = FALSE, na.rm = TRUE)
  list(truth     = truth,
       lower     = qs[1],
       upper     = qs[2],
       post_mean = mean(draws, na.rm = TRUE),
       covered   = as.integer(truth >= qs[1] && truth <= qs[2]))
}

#' Rank statistic for simulation-based calibration.
#'
#' Not used in v1 and not yet validated. Kept here so the seam is visible:
#' swapping this in requires draw_truth() to sample from the model's priors,
#' and draws to be thinned to near-independence. See PLAN.md 7.
statistic_rank <- function(draws, truth, ...) {
  draws <- draws[!is.na(draws)]
  list(truth = truth, rank = sum(draws < truth), n_draws = length(draws))
}

# --- Which parameters are comparable at all ------------------------------
#
# PLAN.md 5.3. Deliberately absent, each for a reason that is NOT an oversight:
#
#   sigma_h   never sampled (group A item 1) -- would fail by construction
#   l_s       never recovered (group A item 2): sigma_s is hard-coded to 1 at
#             the sample_ls() call site, so the length-scale rails at the top
#             of its grid regardless of truth
#   U, L      rotation/sign invariant -- element-wise coverage is meaningless
#   A, C      same, plus they are a latent-trait decomposition
#   Bs, Gs    NOT comparable: the simulator builds the spatial field directly
#             from sigma_s/l_s as `spatField`, and leaves the SoR-form `Bs`
#             empty -- Bst <- matrix(0, S, ps) means truth is 0 x S
#             *regardless of ds* (re-checked at ds = 2, still 0 x S). The fit
#             meanwhile represents the field as sparse-GP basis coefficients
#             over n_supportpoints knots. Different parameterisations.
#   sigma_bs  NOT comparable, for the same reason one level up. The simulator
#             sets Bst <- matrix(0, S, ps) (R/jsdmfun.R:557), so sigma_bs
#             generates nothing at all on the truth side, while the sampler
#             estimates it from the fitted sparse-GP coefficients. Measured:
#             true 0.5 against a posterior mean of ~1.6, 0/8 coverage. That is
#             a meaningless comparison, not a bug.
#
# sigma_b IS kept, but read it with care: with a tight InvGamma(10, 1) prior
# and only p*S residual coefficients, the posterior sits at the prior mean
# (sqrt(1/9) = 0.333) unless the data are very informative. A fixed truth far
# from that will show near-zero coverage for reasons of prior-data conflict
# rather than sampler error -- measured at true 0.5 -> 0/8. The base scenario
# therefore sets sigma_b near the prior mean. SBC (PLAN.md 7), which draws the
# truth from the prior, is the principled fix.
#
# The residual correlation matrix IS identified and is handled separately by
# simstudy_rescor_rows().

simstudy_param_blocks <- function(fit, sim, truth) {
  ro <- fit$results_output
  jo <- ro$jsdm_output
  jp <- sim$true_params$jsdmParams_true

  list(
    B0         = list(post = jo$B0_output,         truth = jp$B0),
    B          = list(post = jo$B_output,          truth = jp$B),
    G          = list(post = jo$G_output,          truth = jp$G),
    sigma_b    = list(post = jo$sigmab_output,     truth = jp$sigma_b),
    tau        = list(post = jo$tau_output,        truth = jp$tau),
    beta_theta = list(post = ro$beta_theta_output,
                      truth = sim$true_params$beta_theta_true),
    theta0     = list(post = ro$theta0_output,     truth = truth$params$theta0),
    p          = list(post = ro$p_output,          truth = sim$true_params$p_true),
    q          = list(post = ro$q_output,          truth = sim$true_params$q_true)
  )
}

# --- Extraction ----------------------------------------------------------

#' One row per parameter element, for one fitted replicate.
simstudy_rows <- function(fit, sim, truth, scenario, replicate,
                          statistic = statistic_coverage) {
  blocks <- simstudy_param_blocks(fit, sim, truth)
  out <- list()

  for (nm in names(blocks)) {
    b <- blocks[[nm]]

    # Absent by design for this model type (e.g. p/q outside two_stage, tau
    # outside continuous). Skip quietly.
    if (is.null(b$post) || is.null(b$truth)) next
    if (all(is.na(b$post))) next
    if (length(b$truth) == 0L) next

    draws <- posterior_draw_matrix(b$post)
    tv <- as.vector(b$truth)

    # Hard error rather than a silent skip: a length mismatch means the
    # truth/posterior mapping is wrong, and a wrong mapping produces
    # confident nonsense rather than an obvious failure.
    if (nrow(draws) != length(tv)) {
      stop(sprintf(
        "simstudy: '%s' truth has %d elements but posterior has %d (scenario '%s'). Mapping is wrong.",
        nm, length(tv), nrow(draws), scenario$label))
    }

    stats <- lapply(seq_len(nrow(draws)),
                    function(i) statistic(draws[i, ], tv[i]))
    out[[nm]] <- data.frame(
      scenario  = scenario$label,
      replicate = replicate,
      block     = nm,
      element   = seq_len(nrow(draws)),
      do.call(rbind, lapply(stats, as.data.frame)),
      stringsAsFactors = FALSE
    )
  }

  rc <- simstudy_rescor_rows(fit, sim, scenario, replicate)
  if (!is.null(rc)) out[["resid_cor"]] <- rc

  if (length(out) == 0L) return(NULL)
  do.call(rbind, out)
}

#' Residual correlations: identified even though L itself is not.
#'
#' returnResidualCorrelationMatrix() returns a 3 x S x S array of
#' (2.5%, 50%, 97.5%) quantiles, so the interval is available directly.
#' Truth is cov2cor(crossprod(L)). Only the upper triangle is used -- the
#' diagonal is 1 by construction and would inflate coverage.
simstudy_rescor_rows <- function(fit, sim, scenario, replicate) {
  L <- sim$true_params$jsdmParams_true$L
  if (is.null(L) || nrow(L) == 0L) return(NULL)

  rc <- try(returnResidualCorrelationMatrix(fit, confidence = 0.95),
            silent = TRUE)
  if (inherits(rc, "try-error")) return(NULL)

  V <- crossprod(L)
  # A species with an all-zero true loading vector has an undefined residual
  # correlation (0/0), and cov2cor() would warn and return non-finite values.
  # Drop only those species, not the whole matrix -- excluding every pair
  # because one species is degenerate would throw away most of the data.
  ok <- is.finite(diag(V)) & diag(V) > 0
  if (sum(ok) < 2L) return(NULL)
  V <- V[ok, ok, drop = FALSE]
  rc <- rc[, ok, ok, drop = FALSE]
  truth <- stats::cov2cor(V)
  off <- upper.tri(truth)

  data.frame(
    scenario  = scenario$label,
    replicate = replicate,
    block     = "resid_cor",
    element   = seq_len(sum(off)),
    truth     = truth[off],
    lower     = rc[1, , ][off],
    upper     = rc[3, , ][off],
    post_mean = rc[2, , ][off],
    covered   = as.integer(truth[off] >= rc[1, , ][off] &
                           truth[off] <= rc[3, , ][off]),
    stringsAsFactors = FALSE
  )
}

# --- Occupancy recovery --------------------------------------------------
#
# The parameter blocks above ask whether the sampler recovers coefficients.
# They never ask the question a user actually has, which is whether the model
# puts the species in the right places. That is scored here.
#
# These statistics live in helper-simstudy.R rather than in
# test-recovery-occupancy.R, where they were first written, because two
# callers need them and they must not drift apart: the tier-2 test scores one
# fit and asserts on it, and the tier-3 study now scores every replicate so
# the same numbers get an interval instead of resting on a single draw. One
# definition, two consumers.

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

#' One row per species: how well was occupancy itself recovered?
#'
#' Scored as point-estimate agreement, not coverage, because there are no
#' draws to build an interval from: `runOccJSDM()` defaults to
#' `summarisedLatentPresences = TRUE`, so `psi_output` and `z_output` come
#' back as n x S running means with no iteration or chain dimension.
#'
#' Two metrics, deliberately different in kind. `psi_cor` correlates the
#' estimated occupancy probability against the true one across sites; its null
#' value is 0 and it is invariant to the level of psi, so it asks only whether
#' the model put occupancy in the right *places*. `auc` asks the same question
#' more forgivingly, on a scale whose null is 0.5; its ceiling is well below 1
#' because z is a coin flip given psi, so even a perfect psi cannot rank every
#' site correctly.
#'
#' **These rows are deliberately NOT rbind-able onto `simstudy_rows()` output.**
#' That frame carries `truth`/`lower`/`upper`/`covered` for every row and has
#' no NAs anywhere; these statistics have no interval and no coverage, so
#' mixing them would put NAs into columns that `test-recovery.R` reduces with
#' `mean()`, `min()` and `cor()`, turning three live assertions into silent
#' NAs. They travel in their own table for that reason.
simstudy_occupancy_rows <- function(fit, sim, scenario, replicate) {
  S <- scenario$S
  psi_true <- stats::plogis(sim$true_params$jsdmParams_true$eta)
  z_true <- sim$true_params$z_true
  psi_est <- fit$results_output$psi_output
  z_est <- fit$results_output$z_output

  empty <- data.frame(scenario = character(0), replicate = integer(0),
                      species = integer(0), occ_true = numeric(0),
                      occ_est = numeric(0), psi_cor = numeric(0),
                      auc = numeric(0), stringsAsFactors = FALSE)

  # The `binary` and `continuous` arms have no latent occupancy state: presence
  # is observed directly rather than inferred through a detection process, so
  # the fit returns no psi_output or z_output and the simulator no z_true.
  # There is nothing to recover and nothing to score, which is different from
  # being broken -- these cells still contribute their parameter coverage, so
  # return no occupancy rows rather than failing the replicate.
  if (length(psi_est) == 0L || length(z_est) == 0L || is.null(z_true)) {
    return(empty)
  }

  # Guard the shape rather than trusting it. If a future change to
  # summarisedLatentPresences hands back a 4-D array instead of an n x S mean,
  # colMeans() and cor() would still return numbers, just meaningless ones.
  # Distinct from the case above: something is present but the wrong shape,
  # which is a regression worth failing loudly on.
  if (length(dim(psi_est)) != 2L || ncol(psi_est) != S) {
    stop("simstudy_occupancy_rows(): psi_output is not the expected n x S ",
         "matrix (got ", paste(dim(psi_est), collapse = " x "), ", S = ", S,
         "). Has summarisedLatentPresences changed?")
  }

  # What DRIVES each species' occupancy, carried alongside how well it was
  # recovered. Free at this point: the simulator has already computed the
  # partition, and eta is the truth we are scoring against anyway.
  #
  # This exists because "why do species differ in psi_cor?" turned out to have
  # a clear answer that none of the stored columns could express. Measured over
  # 1000 species-fits of `base` on 10 August 2026: psi_cor correlates +0.73
  # with the environmental share and -0.64 with the biotic share, and 0.61 of
  # its variance is explained once sd_eta joins them. Prevalence, the obvious
  # suspect, comes in at -0.01.
  #
  # The mechanism is structural rather than a defect. The environmental part of
  # eta is X %*% B, where X is known and B pools across every site, so it is
  # well identified. The biotic part is U %*% L, and U is one row per site
  # informed only by that site's own observations. A species living mostly in U
  # is being reconstructed from the weakest-identified part of the model.
  #
  # sd_eta earns its place separately: a species whose true occupancy barely
  # varies across sites has no pattern to recover, and psi_cor then correlates
  # noise against noise. That is a property of the statistic, not of the model.
  jp <- sim$true_params$jsdmParams_true
  vp <- jp$varPart
  eta <- jp$eta

  # varPart's first three columns are a genuine decomposition summing to 1.
  # Do NOT substitute var(component)/var(eta): the components are correlated
  # and that ratio exceeded 1 by two orders of magnitude when first tried,
  # because X %*% B is large and cancelled by the rest.
  has_vp <- !is.null(vp) && all(c("Environmental", "Spatial", "Biotic") %in%
                                colnames(vp)) && nrow(vp) == S

  data.frame(
    scenario  = scenario$label,
    replicate = replicate,
    species   = seq_len(S),
    occ_true  = colMeans(z_true),
    occ_est   = colMeans(z_est),
    psi_cor   = vapply(seq_len(S), function(s)
                       stats::cor(psi_est[, s], psi_true[, s]), numeric(1)),
    auc       = vapply(seq_len(S), function(s)
                       occrec_auc(psi_est[, s], z_true[, s]), numeric(1)),
    sd_eta     = if (!is.null(eta)) apply(eta, 2, stats::sd) else NA_real_,
    share_env  = if (has_vp) vp[, "Environmental"] else NA_real_,
    share_spat = if (has_vp) vp[, "Spatial"]       else NA_real_,
    share_lat  = if (has_vp) vp[, "Biotic"]        else NA_real_,
    stringsAsFactors = FALSE
  )
}

# --- Orchestration -------------------------------------------------------

#' Run one replicate end to end.
#'
#' Returns a **list of two data frames**, not one frame. `params` is the
#' per-element coverage table that `simstudy_summarise()` consumes;
#' `occupancy` is the per-species table from `simstudy_occupancy_rows()`.
#' They are kept apart because their schemas genuinely differ -- see the note
#' on that function. Callers that only want the old behaviour take `$params`.
simstudy_replicate <- function(scenario, replicate, MCMCparams = NULL,
                               statistic = statistic_coverage) {
  seed <- simstudy_seed_for(scenario, replicate)
  truth <- draw_truth(scenario, seed)
  sim <- simstudy_simulate(truth)
  fit <- if (is.null(MCMCparams)) simstudy_fit(sim, scenario) else
    simstudy_fit(sim, scenario, MCMCparams)
  list(
    params    = simstudy_rows(fit, sim, truth, scenario, replicate, statistic),
    occupancy = simstudy_occupancy_rows(fit, sim, scenario, replicate)
  )
}

#' Deterministic seed from (scenario label, replicate).
#'
#' Reproducibility matters twice over: a whole run can be repeated, and a
#' single failing replicate can be re-run in isolation without replaying the
#' ones before it.
simstudy_seed <- function(label, replicate) {
  as.integer(sum(as.integer(charToRaw(label))) * 1000L + replicate)
}

#' Seed key for a scenario, honouring a shared-truth override.
#'
#' Normally the seed keys on the scenario label, so every cell gets its own
#' datasets. A scenario may instead set `seed_label`, in which case all
#' scenarios sharing that value draw from the same seed at each replicate.
#'
#' This exists for the M ladder (PLAN.md 13). Those arms differ only in how
#' much data is collected -- M samples per site, or K PCR replicates -- not in
#' the process generating it, so sharing a seed makes the comparison *paired*:
#' the differences between arms stop carrying the variance of independent
#' truths, which is exactly what a ladder is reading.
#'
#' `draw_truth()` supports this because its RNG consumption is independent of
#' M and K: it draws only `p`, `q`, `theta0` and `theta_baseline`, sized by
#' `P * S` and `S`. Anything sized by `N = sum(M_vec)` would desynchronise the
#' stream and silently break the pairing. **If you add a draw to
#' `draw_truth()`, keep it independent of M and K, or this override becomes a
#' lie.** `test-regression-bugs.R` asserts the pairing holds.
simstudy_seed_for <- function(scenario, replicate) {
  key <- if (!is.null(scenario$seed_label)) scenario$seed_label else scenario$label
  simstudy_seed(key, replicate)
}

#' Run R replicates of one scenario.
#'
#' Returns the same two-component list as `simstudy_replicate()`, with each
#' component stacked over replicates.
simstudy_scenario <- function(scenario, R = 100L, MCMCparams = NULL,
                              statistic = statistic_coverage,
                              verbose = TRUE) {
  reps <- vector("list", R)
  for (r in seq_len(R)) {
    if (verbose) message(sprintf("[%s] replicate %d/%d", scenario$label, r, R))
    reps[[r]] <- simstudy_replicate(scenario, r, MCMCparams, statistic)
  }
  simstudy_bind(reps)
}

#' Stack a list of two-component replicate results into one such result.
#'
#' Used by both `simstudy_scenario()` and the parallel runner, so the two
#' cannot disagree about how the pieces are assembled.
simstudy_bind <- function(reps) {
  reps <- Filter(Negate(is.null), reps)
  list(
    params    = do.call(rbind, lapply(reps, `[[`, "params")),
    occupancy = do.call(rbind, lapply(reps, `[[`, "occupancy"))
  )
}

#' Aggregate per-element rows into per-block coverage, bias and RMSE.
simstudy_summarise <- function(rows) {
  if (is.null(rows) || nrow(rows) == 0L) return(NULL)
  key <- paste(rows$scenario, rows$block, sep = "\r")
  sp <- split(seq_len(nrow(rows)), key)

  out <- lapply(sp, function(ix) {
    d <- rows[ix, , drop = FALSE]
    err <- d$post_mean - d$truth
    data.frame(
      scenario  = d$scenario[1],
      block     = d$block[1],
      n_element = length(unique(d$element)),
      R         = length(unique(d$replicate)),
      n_checked = nrow(d),
      coverage  = mean(d$covered),
      bias      = mean(err),
      rmse      = sqrt(mean(err^2)),
      stringsAsFactors = FALSE
    )
  })
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res[order(res$scenario, res$block), ]
}

#' Aggregate per-species occupancy rows to one row per scenario.
#'
#' **Not aggregated by species index, deliberately.** Every replicate draws its
#' own truth, so "species 3" in replicate 1 and "species 3" in replicate 2 are
#' unrelated draws with different prevalences. Averaging down the species
#' column would be averaging over an arbitrary label. What is meaningful is the
#' distribution over all replicate-by-species pairs, which is what this
#' returns: at R = 100 and S = 10 that is 1000 species-fits per cell, against
#' the single fit the tier-2 table reports.
#'
#' `psi_cor_q05`/`q95` are the interval to quote. `share_ok` uses the same 0.30
#' cut as the tier-2 assertion so the two are comparable.
#'
#' `cor_prev` correlates a species' true prevalence with how well it was
#' recovered. **Do not read it as "rare species are the hard ones".** The
#' 10 August R = 100 run measured it at between -0.07 and +0.005 in every cell,
#' while the prevalence bands differ substantially: recovery climbs from the
#' rarest band, peaks in the middle, then eases off again for the commonest
#' species. A linear correlation across a humped relationship cancels to
#' nothing. It is kept because a *large* value would be informative and its
#' absence is itself worth recording; use
#' `simstudy_occupancy_by_prevalence()` for the actual shape.
simstudy_summarise_occupancy <- function(occ, cut = 0.30) {
  if (is.null(occ) || nrow(occ) == 0L) return(NULL)
  sp <- split(seq_len(nrow(occ)), occ$scenario)

  qs <- function(x, p) as.numeric(stats::quantile(x, p, na.rm = TRUE))

  out <- lapply(sp, function(ix) {
    d <- occ[ix, , drop = FALSE]
    data.frame(
      scenario      = d$scenario[1],
      R             = length(unique(d$replicate)),
      n_species_fit = nrow(d),
      occ_true      = mean(d$occ_true, na.rm = TRUE),
      occ_est       = mean(d$occ_est, na.rm = TRUE),
      psi_cor_mean  = mean(d$psi_cor, na.rm = TRUE),
      psi_cor_q05   = qs(d$psi_cor, 0.05),
      psi_cor_q50   = qs(d$psi_cor, 0.50),
      psi_cor_q95   = qs(d$psi_cor, 0.95),
      auc_mean      = mean(d$auc, na.rm = TRUE),
      auc_q05       = qs(d$auc, 0.05),
      auc_q95       = qs(d$auc, 0.95),
      share_ok      = mean(d$psi_cor >= cut, na.rm = TRUE),
      cor_prev      = suppressWarnings(
                        stats::cor(d$occ_true, d$psi_cor,
                                   use = "complete.obs")),
      stringsAsFactors = FALSE
    )
  })
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res[order(res$scenario), ]
}

#' Occupancy recovery banded by how common the species actually is.
#'
#' The headline table reports recovery averaged over species, which hides the
#' thing that most affects a user: a species present at 5% of sites is much
#' harder to place than one present at 60%. Banding by true prevalence turns
#' that from an anecdote about one row into a curve with an interval on it,
#' and it is what the validation article plots.
simstudy_occupancy_by_prevalence <- function(occ,
                                             breaks = c(0, .1, .2, .3, .5, 1)) {
  if (is.null(occ) || nrow(occ) == 0L) return(NULL)
  band <- cut(occ$occ_true, breaks = breaks, include.lowest = TRUE)
  sp <- split(seq_len(nrow(occ)), list(occ$scenario, band), drop = TRUE)

  qs <- function(x, p) as.numeric(stats::quantile(x, p, na.rm = TRUE))

  out <- lapply(sp, function(ix) {
    d <- occ[ix, , drop = FALSE]
    data.frame(
      scenario      = d$scenario[1],
      prevalence    = as.character(band[ix][1]),
      n_species_fit = nrow(d),
      occ_true      = mean(d$occ_true, na.rm = TRUE),
      psi_cor_mean  = mean(d$psi_cor, na.rm = TRUE),
      psi_cor_q05   = qs(d$psi_cor, 0.05),
      psi_cor_q95   = qs(d$psi_cor, 0.95),
      auc_mean      = mean(d$auc, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res[order(res$scenario, res$occ_true), ]
}
