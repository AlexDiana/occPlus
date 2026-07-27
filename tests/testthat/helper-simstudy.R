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
    # Low information: weak detection AND fewer sites. n_supportpoints must be
    # set explicitly here or the default would crash below 31 sites
    # (TODO.Rmd group B item 3).
    mk("low_information",   n = 40L, n_supportpoints = 8L,
                            p_range = c(0.1, 0.3),
                            theta_baseline_range = c(0.05, 0.2)),
    mk("d_underfit",        d = 4L, fit_d = 2L),
    mk("d_overfit",         d = 2L, fit_d = 4L),
    mk("species_30",        S = 30L),
    mk("occupancy",         model = "occupancy"),
    mk("binary",            model = "binary")
  )
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
    jsdmParams = list(gt = scenario$gt, d = scenario$d, ds = 0,
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
  collCov <- if (scenario$model %in% c("occupancy", "two_stage")) "X_theta" else NULL

  dat <- sim$data_list
  if (!isTRUE(scenario$fit_traits)) dat$traits <- NULL

  suppressMessages(suppressWarnings(
    runOccJSDM(dat,
               listParams     = listParams,
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
#             from sigma_s/l_s (truth is empty at ds = 0), while the fit
#             represents it as sparse-GP basis coefficients over
#             n_supportpoints knots. Different parameterisations, verified.
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
  seed <- simstudy_seed(scenario$label, replicate)
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
