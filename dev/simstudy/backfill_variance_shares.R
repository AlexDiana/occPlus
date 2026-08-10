#!/usr/bin/env Rscript
#
# Backfill sd_eta and the variance shares into an existing validation bundle.
#
#   Rscript dev/simstudy/backfill_variance_shares.R [--scenario=base]
#
# WHY THIS IS A SEPARATE, TEMPORARY SCRIPT
#
# simstudy_occupancy_rows() now records sd_eta, share_env, share_spat and
# share_lat for every species-fit, so any run from 10 August 2026 onward
# carries them for every cell and this script is unnecessary. It exists only to
# fill them into bundles produced BEFORE that change, so the validation article
# can show the decomposition without waiting on a fresh multi-hour run.
#
# Delete it once no pre-change results file is still in use.
#
# It is cheap because it needs no fitting. draw_truth() and
# simstudy_simulate() regenerate the exact dataset for a given
# (scenario, replicate) in about a second, and the variance partition is a
# property of the simulated truth alone -- nothing here depends on what the
# sampler did with it.

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}
scenario_label <- getopt("scenario", "base")
bundle_path <- getopt("bundle", "dev/simstudy/validation-data.rds")

suppressMessages(devtools::load_all(".", quiet = TRUE))
source("tests/testthat/helper-simstudy.R")

# The runner's generator. draw_truth() calls set.seed() without pinning `kind`,
# so it inherits the caller's; the PSOCK workers are left on L'Ecuyer-CMRG by
# clusterSetRNGStream(). Standalone that would be Mersenne-Twister, and the
# same seed would draw a DIFFERENT truth. This line is not optional.
RNGkind("L'Ecuyer-CMRG")

v <- readRDS(bundle_path)
occ <- v$occupancy
# Check the REQUESTED scenario, not the bundle as a whole. Checking globally
# meant that once any one cell was filled, every later cell was skipped.
if ("share_env" %in% names(occ)) {
  mine <- occ$scenario == scenario_label
  if (any(mine) && !any(is.na(occ$share_env[mine]))) {
    message("'", scenario_label, "' already carries variance shares; nothing to do")
    quit(save = "no")
  }
}

scenarios <- simstudy_scenarios()
labels <- vapply(scenarios, function(s) s$label, "")
sc <- scenarios[[which(labels == scenario_label)]]

reps <- sort(unique(occ$replicate[occ$scenario == scenario_label]))
message(sprintf("regenerating truths for %s, %d replicates ...",
                scenario_label, length(reps)))

pred <- do.call(rbind, lapply(reps, function(r) {
  truth <- draw_truth(sc, simstudy_seed_for(sc, r))
  sim <- simstudy_simulate(truth)
  jp <- sim$true_params$jsdmParams_true
  vp <- jp$varPart
  # Assert the decomposition rather than trusting it: a hand-rolled
  # var(component)/var(eta) was tried first and exceeded 1 by two orders of
  # magnitude, because the components are correlated and cancel.
  stopifnot(all(abs(rowSums(vp[, c("Environmental", "Spatial", "Biotic")]) - 1)
                < 1e-6))
  data.frame(scenario = scenario_label, replicate = r,
             species = seq_len(ncol(jp$eta)),
             sd_eta = apply(jp$eta, 2, stats::sd),
             share_env  = vp[, "Environmental"],
             share_spat = vp[, "Spatial"],
             share_lat  = vp[, "Biotic"],
             stringsAsFactors = FALSE)
}))

# Cross-check before merging: if the regenerated truths do not match what the
# study stored, this is a different dataset and the merge would be nonsense.
chk <- merge(occ[occ$scenario == scenario_label, c("replicate","species","occ_true")],
             pred, by = c("replicate","species"))
stopifnot(nrow(chk) == nrow(pred))

for (cl in c("sd_eta", "share_env", "share_spat", "share_lat"))
  if (!cl %in% names(occ)) occ[[cl]] <- NA_real_

idx <- match(paste(occ$scenario, occ$replicate, occ$species),
             paste(pred$scenario, pred$replicate, pred$species))
hit <- !is.na(idx)
for (cl in c("sd_eta", "share_env", "share_spat", "share_lat"))
  occ[[cl]][hit] <- pred[[cl]][idx[hit]]

v$occupancy <- occ
v$variance_shares_backfilled <- scenario_label
saveRDS(v, bundle_path, compress = "xz")

message(sprintf("filled %d of %d occupancy rows in %s",
                sum(hit), nrow(occ), bundle_path))
d <- occ[hit, ]
cat("\ncorrelation with psi_cor:\n")
for (cl in c("share_env", "share_lat", "share_spat", "sd_eta", "occ_true"))
  cat(sprintf("  %-11s %+0.3f\n", cl, stats::cor(d[[cl]], d$psi_cor)))
