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

    # --- M ladder (PLAN.md 13) -----------------------------------------
    #
    # Not production cells: these answer whether group B items 4-6 are code
    # defects or Stage 1 under-identification. Select them explicitly with
    # --scenarios=M2,M5,M10,M20,K30; they are otherwise inert because the
    # runner takes the full list only when --scenarios is empty, so leaving
    # them here would silently add 5 cells to every full run.
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
    mk("nocollcov", ncov_theta = 0L, seed_label = "nocollcov"))
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
      q = matrix(stats::runif(P * S, 0.01, 0.05), P, S),
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

  occCov <- paste0("X_psi.EnvCov.", seq_len(scenario$ncov_psi))
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

# --- Orchestration -------------------------------------------------------

#' Run one replicate end to end.
simstudy_replicate <- function(scenario, replicate, MCMCparams = NULL,
                               statistic = statistic_coverage) {
  seed <- simstudy_seed_for(scenario, replicate)
  truth <- draw_truth(scenario, seed)
  sim <- simstudy_simulate(truth)
  fit <- if (is.null(MCMCparams)) simstudy_fit(sim, scenario) else
    simstudy_fit(sim, scenario, MCMCparams)
  simstudy_rows(fit, sim, truth, scenario, replicate, statistic)
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
simstudy_scenario <- function(scenario, R = 100L, MCMCparams = NULL,
                              statistic = statistic_coverage,
                              verbose = TRUE) {
  rows <- vector("list", R)
  for (r in seq_len(R)) {
    if (verbose) message(sprintf("[%s] replicate %d/%d", scenario$label, r, R))
    rows[[r]] <- simstudy_replicate(scenario, r, MCMCparams, statistic)
  }
  do.call(rbind, rows)
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
