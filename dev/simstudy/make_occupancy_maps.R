#!/usr/bin/env Rscript
#
# Per-site occupancy surfaces for one replicate, for mapping.
#
#   Rscript dev/simstudy/make_occupancy_maps.R [--scenario=base] [--replicate=1]
#                                              [--out=dev/simstudy/mapdata.rds]
#
# Then render the report beside it:
#
#   Rscript -e 'rmarkdown::render("dev/simstudy/occupancy_maps.Rmd")'
#
# WHY THIS EXISTS
#
# The study records per-species aggregates -- occ_true, occ_est, psi_cor, auc --
# and discards the per-site vectors behind them, because storing n x S doubles
# per replicate across ten cells is a lot of disk for something only ever looked
# at one replicate at a time. So maps cannot be drawn from
# simstudy-*.rds or from the occupancy CSV; nothing in either has a site in it.
#
# Refitting one replicate costs about 90 seconds and reproduces the stored
# numbers exactly, which is cheaper than carrying the arrays around.
#
# TWO THINGS ARE LOAD-BEARING. Both look like noise and are not.
#
# 1. RCPP_PARALLEL_NUM_THREADS=1 must be set when running this. Above one
#    thread the sampler is not reproducible (TODO.md group B items 8 and 9),
#    so the maps change every run and stop matching the study. The script
#    refuses to start without it rather than producing quietly wrong output.
#
# 2. RNGkind("L'Ecuyer-CMRG"), below. run_study.R parallelises with a PSOCK
#    cluster and calls clusterSetRNGStream(), which leaves each worker on
#    L'Ecuyer-CMRG. draw_truth() then calls set.seed(seed) WITHOUT pinning
#    `kind`, so it inherits whatever the caller had. Standalone that is the
#    default Mersenne-Twister, and the same seed draws a different truth.
#    Verified on 10 August 2026: this one line is the difference between
#    reproducing study replicate 1 and silently getting another dataset.

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}

scenario_label <- getopt("scenario", "base")
replicate_id   <- as.integer(getopt("replicate", "1"))
out_path       <- getopt("out", "dev/simstudy/mapdata.rds")
pkg_root       <- normalizePath(getopt("pkg", "."), mustWork = TRUE)

threads <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", "")
if (!identical(threads, "1")) {
  stop("Set RCPP_PARALLEL_NUM_THREADS=1 before running this.\n",
       "  Above one thread the sampler is not reproducible, so the maps would\n",
       "  not match the study and would differ on every run. See the header.\n",
       "  Rerun as:\n",
       "    RCPP_PARALLEL_NUM_THREADS=1 Rscript ",
       "dev/simstudy/make_occupancy_maps.R", call. = FALSE)
}

suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))
source(file.path(pkg_root, "tests", "testthat", "helper-simstudy.R"))
RcppParallel::setThreadOptions(numThreads = 1)

# See note 2 in the header. Do not remove as redundant.
RNGkind("L'Ecuyer-CMRG")

scenarios <- simstudy_scenarios()
labels <- vapply(scenarios, function(s) s$label, "")
if (!scenario_label %in% labels)
  stop("no scenario called '", scenario_label, "'. Available: ",
       paste(labels, collapse = ", "), call. = FALSE)
sc <- scenarios[[which(labels == scenario_label)]]

# The study's own settings, so the refit matches rather than approximates.
MCMCparams <- list(nchain = 2, nburn = 1000, niter = 1000, nthin = 1)

message(sprintf("fitting %s replicate %d (about 90 s) ...",
                sc$label, replicate_id))
truth <- draw_truth(sc, simstudy_seed_for(sc, replicate_id))
sim   <- simstudy_simulate(truth)
fit   <- simstudy_fit(sim, sc, MCMCparams)

psi_true <- stats::plogis(sim$true_params$jsdmParams_true$eta)   # n x S
z_true   <- sim$true_params$z_true
psi_est  <- fit$results_output$psi_output

if (length(psi_est) == 0L)
  stop("scenario '", sc$label, "' has no latent occupancy state, so there is ",
       "nothing to map. The binary and continuous arms observe presence ",
       "directly.", call. = FALSE)

# Site coordinates. Xs is stored per row of the long info frame (one row per
# sample), so collapse to one row per site or the map gets N points, not n.
info <- sim$data_list$info
site <- info$Site %||% info$site
xy <- unique(data.frame(site = site, x = info$Xs.1, y = info$Xs.2))
xy <- xy[order(xy$site), ]
stopifnot(nrow(xy) == nrow(psi_true))

S <- ncol(psi_true)
per_species <- data.frame(
  species  = seq_len(S),
  occ_true = colMeans(z_true),
  occ_est  = colMeans(psi_est),
  psi_cor  = vapply(seq_len(S), function(s)
                    stats::cor(psi_est[, s], psi_true[, s]), numeric(1)),
  auc      = vapply(seq_len(S), function(s)
                    occrec_auc(psi_est[, s], z_true[, s]), numeric(1)))

long <- do.call(rbind, lapply(seq_len(S), function(s) {
  data.frame(species = s, x = xy$x, y = xy$y,
             psi_true = psi_true[, s], psi_est = psi_est[, s],
             z_true = z_true[, s])
}))

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(per_species = per_species, long = long,
             n = nrow(xy), S = S,
             scenario = sc$label, replicate = replicate_id,
             MCMCparams = MCMCparams, when = Sys.time()),
        out_path)

message(sprintf("wrote %s  (%d sites, %d species)", out_path, nrow(xy), S))
print(per_species, row.names = FALSE, digits = 3)

# Cross-check against the study, when a run is available to check against.
# A mismatch here means the refit is not the replicate it claims to be, which
# is exactly the failure the RNGkind line above exists to prevent.
vd <- file.path(pkg_root, "dev", "simstudy", "validation-data.rds")
if (file.exists(vd)) {
  v <- readRDS(vd)
  ref <- v$occupancy[v$occupancy$scenario == sc$label &
                     v$occupancy$replicate == replicate_id, ]
  if (nrow(ref) == nrow(per_species)) {
    ok <- isTRUE(all.equal(ref$occ_true, per_species$occ_true))
    message(if (ok)
      "  cross-check: simulated truths match the study run exactly."
      else paste("  CROSS-CHECK FAILED: simulated truths differ from the",
                 "study run.\n  This refit is not that replicate. See note 2",
                 "in the header."))
  }
}
