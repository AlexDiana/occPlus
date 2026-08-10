#!/usr/bin/env Rscript
#
# Extract the compact bundle the validation article renders from.
#
#   Rscript dev/simstudy/export_validation_data.R [path/to/simstudy-*.rds]
#
# With no argument, takes the newest simstudy-*.rds in dev/simstudy/results.
#
# WHY THIS EXISTS
#
# The article used to be 882 lines of hand-typed HTML: every number in it was
# transcribed by hand from a run. That made a rerun cost hours of retyping,
# gave no way to prove the prose matched the data it described, and let the
# whole document go stale invisibly when six commits landed under it. The
# article now reads its numbers, and this script is what it reads.
#
# It is a separate step rather than part of run_study.R because the two have
# different lifetimes: the full results are large and regenerable, while this
# bundle is small and must be committed.
#
# WHAT IS AND IS NOT INCLUDED
#
# `rows` is dropped. It is the per-element coverage table, 155k rows at
# R = 100, and everything the article says about it is already in `summary`.
# `occupancy` is kept in full, because it is only ~11k rows and because the
# distribution across replicates is the point: it is what turns the
# per-species table from a single draw into an interval.
#
# dev/ is in .Rbuildignore, so this file never ships in the package; but
# pkgdown builds articles from the source tree, so the path resolves at
# build time. dev/simstudy/results/ is gitignored and this is not: see the
# comment in .gitignore.

#   --extra=file1.rds,file2.rds
#
# folds in single-cell runs that are not part of the production grid, under
# out$extra_cells, keyed by scenario label. They are kept OUT of `summary`,
# `accuracy` and `occupancy` on purpose: those drive the "across all ten
# scenarios" sections, and quietly adding an eleventh cell would corrupt every
# one of them. Comparison cells are a different kind of thing and travel
# separately.

args <- commandArgs(trailingOnly = TRUE)
res_dir <- "dev/simstudy/results"

# Positional argument only: anything starting with -- is a flag, and taking
# args[1] blindly meant --extra= was read as the source filename.
positional <- args[!grepl("^--", args)]

src <- if (length(positional) >= 1L) positional[1] else {
  f <- list.files(res_dir, pattern = "^simstudy-[0-9-]+\\.rds$",
                  full.names = TRUE)
  if (length(f) == 0L)
    stop("no simstudy-*.rds found in ", res_dir,
         "; run dev/simstudy/run_study.R first")
  f[order(file.mtime(f), decreasing = TRUE)][1]
}

message("reading ", src)
r <- readRDS(src)

# Guard against the newest-file default silently picking a single-cell
# comparison run instead of the production grid. Once cells like sites_300 and
# fit_d4 exist, "newest" is no longer a safe proxy for "the study".
if (length(unique(r$rows$scenario)) < 2L && !length(positional)) {
  stop("the newest results file (", basename(src), ") holds only the '",
       unique(r$rows$scenario), "' cell, which is a comparison run rather ",
       "than the production grid.\n  Name the study run explicitly, e.g.\n",
       "    Rscript dev/simstudy/export_validation_data.R ",
       "dev/simstudy/results/simstudy-<stamp>.rds", call. = FALSE)
}

need <- c("summary", "occupancy")
missing <- need[!need %in% names(r)]
if (length(missing) > 0L) {
  stop("results file predates the occupancy instrumentation (missing: ",
       paste(missing, collapse = ", "), "). ",
       "Re-run dev/simstudy/run_study.R rather than exporting this one.")
}

# Recompute the two derived tables rather than copying them across, so an
# older results file with a stale aggregation cannot leak into the article.
source(file.path("tests", "testthat", "helper-simstudy.R"))

# One concrete fit, kept as a worked example.
#
# The article leads with the grid-wide distribution, which is the honest
# summary, but a reader needs one table small enough to actually read. This is
# replicate 1 of `base`: the same shape the article always showed, except now
# generated rather than transcribed, and sitting next to the spread across the
# other 99 replicates so a striking row can be checked against the population
# instead of being narrated as if it were typical.
#
# drv_cov and drv_err come from the per-element rows for that replicate, which
# is the only reason `rows` is read here at all -- it is dropped from the
# bundle immediately after.
example_cell <- "base"
example_rep <- 1L
species_example <- local({
  occ <- r$occupancy[r$occupancy$scenario == example_cell &
                     r$occupancy$replicate == example_rep, , drop = FALSE]
  if (nrow(occ) == 0L) return(NULL)

  drv <- r$rows[r$rows$scenario == example_cell &
                r$rows$replicate == example_rep &
                r$rows$block %in% c("B0", "B"), , drop = FALSE]
  if (nrow(drv) == 0L) return(occ)

  # Element to species, the same mapping test-recovery-occupancy.R asserts:
  # truth is flattened column-major, so B0 gives one element per species and
  # B gives ncov consecutive elements per species.
  S <- nrow(occ)
  ncov <- sum(drv$block == "B") / S
  if (ncov != round(ncov)) return(occ)
  drv$species <- ifelse(drv$block == "B0", drv$element,
                        ceiling(drv$element / ncov))

  occ$drv_cov <- vapply(occ$species, function(s)
    mean(drv$covered[drv$species == s]), numeric(1))
  occ$drv_err <- vapply(occ$species, function(s) {
    d <- drv[drv$species == s, ]
    mean(abs(d$post_mean - d$truth))
  }, numeric(1))
  occ
})

# Per-block accuracy: how close are the estimates, as opposed to how honest the
# intervals are. `summary` carries bias and rmse; the two numbers a reader can
# actually interpret are derived here and need the raw rows, which is the other
# reason this step exists before they are dropped.
#
#   nrmse  RMSE divided by the spread of the true values themselves. Below 1 the
#          model beats guessing the average of the truth; above 1 it does not.
#          Scale-free, so blocks on different scales are comparable.
#   r      Correlation between posterior mean and truth. nrmse can look tolerable
#          for an estimator that has simply shrunk to the mean; r is what
#          notices, so both are reported.
accuracy <- local({
  key <- paste(r$rows$scenario, r$rows$block, sep = "\r")
  parts <- split(seq_len(nrow(r$rows)), key)
  do.call(rbind, lapply(parts, function(ix) {
    d <- r$rows[ix, , drop = FALSE]
    sd_truth <- stats::sd(d$truth)
    err <- d$post_mean - d$truth
    data.frame(
      scenario = d$scenario[1], block = d$block[1],
      n = nrow(d),
      bias = mean(err),
      rmse = sqrt(mean(err^2)),
      sd_truth = sd_truth,
      # A block whose truth is constant has no spread to normalise by; report
      # NA rather than an infinity that would render as garbage.
      nrmse = if (is.finite(sd_truth) && sd_truth > 0)
                sqrt(mean(err^2)) / sd_truth else NA_real_,
      r = suppressWarnings(stats::cor(d$post_mean, d$truth)),
      stringsAsFactors = FALSE)
  }))
})
rownames(accuracy) <- NULL

# The resid_cor degeneracy (PLAN.md 17), measured rather than asserted.
#
# The article claims residual-correlation coverage is arithmetically one minus
# the share of true correlations sitting at exactly +/-1, because the intervals
# span nearly the whole of [-1, 1] and only the boundary cases fall outside.
# That is the least obvious claim on the page and it was hand-typed from a
# one-off investigation. Computing it here means the article can show the
# predicted and measured numbers side by side, and the claim falls over
# visibly if it ever stops being true.
degeneracy <- local({
  d <- r$rows[r$rows$block == "resid_cor", , drop = FALSE]
  if (nrow(d) == 0L) return(NULL)
  parts <- split(seq_len(nrow(d)), d$scenario)
  res <- do.call(rbind, lapply(parts, function(ix) {
    x <- d[ix, , drop = FALSE]
    w <- x$upper - x$lower
    data.frame(
      scenario      = x$scenario[1],
      n             = nrow(x),
      share_at_one  = mean(abs(abs(x$truth) - 1) < 1e-8),
      median_width  = stats::median(w),
      share_wide    = mean(w > 1.9),
      coverage      = mean(x$covered),
      stringsAsFactors = FALSE)
  }))
  # The prediction the article states, next to what was actually measured.
  res$predicted_coverage <- 1 - res$share_at_one
  rownames(res) <- NULL
  res
})

out <- list(
  summary                 = r$summary,
  accuracy                = accuracy,
  degeneracy              = degeneracy,
  occupancy               = r$occupancy,
  species_example         = species_example,
  example_cell            = example_cell,
  example_replicate       = example_rep,
  occupancy_summary       = simstudy_summarise_occupancy(r$occupancy),
  occupancy_by_prevalence = simstudy_occupancy_by_prevalence(r$occupancy),
  R                       = r$R,
  MCMCparams              = r$MCMCparams,
  failed_replicates       = r$failed_replicates,
  elapsed_min             = r$elapsed_min,
  when                    = r$when,
  provenance              = r$provenance,
  source_file             = basename(src)
)

# --- optional comparison cells -------------------------------------------
extra_arg <- {
  hit <- grep("^--extra=", args, value = TRUE)
  if (length(hit) == 0) "" else sub("^--extra=", "", hit[1])
}
if (nzchar(extra_arg)) {
  # Needed only on this path: regenerating a truth calls simulateOccJSDMData(),
  # which lives in the package rather than the helper.
  suppressMessages(devtools::load_all(".", quiet = TRUE))
  RNGkind("L'Ecuyer-CMRG")   # the runner's generator; see make_occupancy_maps.R
  scen_all <- simstudy_scenarios()
  scen_lab <- vapply(scen_all, function(s) s$label, "")

  out$extra_cells <- list()
  for (f in strsplit(extra_arg, ",")[[1]]) {
    e <- readRDS(f)
    lab <- unique(e$rows$scenario)
    stopifnot(length(lab) == 1L)
    oc <- e$occupancy

    # Runs predating the varPart columns need their shares regenerating.
    # Simulation only, no fitting: the partition is a property of the truth.
    if (!"share_lat" %in% names(oc) || all(is.na(oc$share_lat))) {
      message("  regenerating variance shares for ", lab, " ...")
      sc <- scen_all[[which(scen_lab == lab)]]
      sh <- do.call(rbind, lapply(sort(unique(oc$replicate)), function(r) {
        vp <- simstudy_simulate(draw_truth(sc, simstudy_seed_for(sc, r)
              ))$true_params$jsdmParams_true$varPart
        data.frame(replicate = r, species = seq_len(nrow(vp)),
                   sd_eta = NA_real_,
                   share_env = vp[, "Environmental"],
                   share_spat = vp[, "Spatial"],
                   share_lat = vp[, "Biotic"])
      }))
      oc <- merge(oc, sh, by = c("replicate", "species"))
    }

    out$extra_cells[[lab]] <- list(
      occupancy = oc, summary = e$summary,
      R = e$R, elapsed_min = e$elapsed_min,
      provenance = e$provenance, source_file = basename(f))
    message(sprintf("  + %s: %d occupancy rows, %.0f min",
                    lab, nrow(oc), e$elapsed_min))
  }
}

dest <- file.path("dev", "simstudy", "validation-data.rds")
saveRDS(out, dest, compress = "xz")

message(sprintf("wrote %s (%.0f KB)", dest, file.size(dest) / 1024))
message(sprintf("  %d summary rows, %d occupancy rows, R = %s, %d cell(s)",
                nrow(out$summary), nrow(out$occupancy), out$R,
                length(unique(out$summary$scenario))))
if (!is.null(out$provenance))
  message(sprintf("  from %s%s, %.0f min",
                  substr(out$provenance$git_sha %||% "?", 1, 8),
                  if (isTRUE(out$provenance$git_dirty)) " (dirty)" else "",
                  out$elapsed_min))
