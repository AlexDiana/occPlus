#!/usr/bin/env Rscript
#
# Tier-3 runner for the simulation study. See dev/simstudy/PLAN.md.
#
#   Rscript dev/simstudy/run_study.R [--R=100] [--cores=8] [--scenarios=a,b]
#                                    [--out=dev/simstudy/results]
#                                    [--nburn=1000] [--niter=1000]
#
# Replicates are parallelised across *processes*, which is free: each is an
# independent dataset and fit, needing no package changes -- unlike the
# in-package chain parallelisation in TODO.Rmd group D item 1.
#
# PSOCK rather than mclapply(), for the reasons in TODO.Rmd D.1: mclapply() is
# fork-based and errors outright on Windows, and forking a session with a live
# OpenMP pool is a known deadlock source. PSOCK works on all three platforms.
#
# At the defaults (R = 100, 10 scenarios) expect ~23 h serial, ~3 h on 8 cores.

suppressMessages(library(parallel))

# --- arguments -----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}

R_reps   <- as.integer(getopt("R", "100"))
n_cores  <- as.integer(getopt("cores", max(1L, detectCores() - 1L)))
out_dir  <- getopt("out", "dev/simstudy/results")
which_sc <- getopt("scenarios", "")
nburn    <- as.integer(getopt("nburn", "1000"))
niter    <- as.integer(getopt("niter", "1000"))

pkg_root <- normalizePath(getopt("pkg", "."), mustWork = TRUE)

# --- load package + helpers ----------------------------------------------

suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))
source(file.path(pkg_root, "tests", "testthat", "helper-simstudy.R"))

scenarios <- simstudy_scenarios()
if (nzchar(which_sc)) {
  keep <- strsplit(which_sc, ",")[[1]]
  scenarios <- Filter(function(s) s$label %in% keep, scenarios)
  if (length(scenarios) == 0) stop("no scenario matched --scenarios=", which_sc)
}

MCMCparams <- list(nchain = 2, nburn = nburn, niter = niter, nthin = 1)

message(sprintf("simstudy: %d scenario(s), R = %d, %d core(s), MCMC %d+%d x 2",
                length(scenarios), R_reps, n_cores, nburn, niter))

# --- run -----------------------------------------------------------------

# Flatten to (scenario, replicate) jobs so the pool stays busy: scenarios have
# very different costs (binary is far cheaper than species_30), and running
# them one at a time would leave cores idle at the tail of each.
jobs <- do.call(rbind, lapply(seq_along(scenarios), function(i) {
  data.frame(scenario_i = i, replicate = seq_len(R_reps))
}))
jobs <- jobs[sample.int(nrow(jobs)), ]   # mix costs across workers

run_one <- function(k) {
  s <- scenarios[[jobs$scenario_i[k]]]
  r <- jobs$replicate[k]
  out <- try(simstudy_replicate(s, r, MCMCparams = MCMCparams), silent = TRUE)
  if (inherits(out, "try-error")) {
    # One bad replicate must not lose the whole run. Record and continue.
    return(data.frame(scenario = s$label, replicate = r, block = "ERROR",
                      element = NA_integer_, truth = NA_real_,
                      lower = NA_real_, upper = NA_real_,
                      post_mean = NA_real_, covered = NA_integer_,
                      stringsAsFactors = FALSE))
  }
  out
}

# Progress reporting. Without this the run is silent until it finishes, which
# for a multi-hour grid means no way to tell "working" from "wedged" and no
# ETA -- a real shortcoming when the first full run overran its estimate by
# 50% on a throttling machine.
#
# parLapply() is blocking and silent, so jobs are submitted in chunks and a
# line is emitted after each. Chunks are a multiple of the core count and use
# parLapplyLB() internally, so load balancing is preserved within a chunk; the
# only cost is that each chunk waits for its own slowest job. With jobs
# shuffled (above) that tail is small.
progress_line <- function(done, total, t0, failed) {
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  frac <- done / total
  eta <- if (frac > 0) el * (1 - frac) / frac else NA_real_
  message(sprintf(
    "  [%3.0f%%] %d/%d replicates | %.0f min elapsed | ~%.0f min left%s",
    100 * frac, done, total, el, eta,
    if (failed > 0) sprintf(" | %d failed", failed) else ""))
  flush(stderr())
}

n_jobs <- nrow(jobs)
t0 <- Sys.time()

if (n_cores > 1L) {
  cl <- makeCluster(n_cores)                      # PSOCK; portable
  on.exit(stopCluster(cl), add = TRUE)
  # Independent L'Ecuyer streams per worker, on top of the per-replicate seed
  # in simstudy_replicate(). Reproducible and non-overlapping.
  clusterSetRNGStream(cl, 20260727L)
  clusterExport(cl, c("scenarios", "jobs", "MCMCparams", "pkg_root"),
                envir = environment())
  clusterEvalQ(cl, {
    suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))
    source(file.path(pkg_root, "tests", "testthat", "helper-simstudy.R"))
    NULL
  })

  # Aim for ~20 progress lines, but never a chunk smaller than the pool (that
  # would idle cores) nor larger than 4x it (that would report too rarely).
  chunk_size <- max(n_cores, min(n_cores * 4L, ceiling(n_jobs / 20)))
  chunks <- split(seq_len(n_jobs), ceiling(seq_len(n_jobs) / chunk_size))
  res <- vector("list", length(chunks))
  done <- 0L; failed_so_far <- 0L

  for (ci in seq_along(chunks)) {
    res[[ci]] <- parLapplyLB(cl, chunks[[ci]], run_one)
    done <- done + length(chunks[[ci]])
    failed_so_far <- failed_so_far +
      sum(vapply(res[[ci]], function(d) any(d$block == "ERROR"), logical(1)))
    progress_line(done, n_jobs, t0, failed_so_far)
  }
  res <- unlist(res, recursive = FALSE)

} else {
  res <- vector("list", n_jobs)
  failed_so_far <- 0L
  for (k in seq_len(n_jobs)) {
    res[[k]] <- run_one(k)
    if (any(res[[k]]$block == "ERROR")) failed_so_far <- failed_so_far + 1L
    if (k %% 10L == 0L || k == n_jobs) progress_line(k, n_jobs, t0, failed_so_far)
  }
}
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

rows <- do.call(rbind, res)
failed <- sum(rows$block == "ERROR", na.rm = TRUE)
rows <- rows[rows$block != "ERROR", , drop = FALSE]
summary_tbl <- simstudy_summarise(rows)

# --- write ---------------------------------------------------------------

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
saveRDS(list(rows = rows, summary = summary_tbl, R = R_reps,
             MCMCparams = MCMCparams, failed_replicates = failed,
             elapsed_min = elapsed, when = Sys.time(),
             sessionInfo = utils::sessionInfo()),
        file.path(out_dir, sprintf("simstudy-%s.rds", stamp)))
write.csv(summary_tbl, file.path(out_dir, sprintf("simstudy-%s.csv", stamp)),
          row.names = FALSE)

message(sprintf("\ndone in %.1f min; %d replicate(s) failed", elapsed, failed))
if (failed > 0) message("  (failures are dropped from the summary, not silently counted as covered)")
print(summary_tbl, row.names = FALSE)
message("\nwritten to ", normalizePath(out_dir))

# Reminders that outlive this session -- see PLAN.md 5.3 and 10.
message("\nRead with care:")
message("  * l_s and sigma_h are NOT in the table: never sampled / not recoverable")
message("    (TODO.Rmd group A items 1 and 2). Spatial cells say nothing about range.")
message("  * sigma_b is prior-dominated; its coverage reflects prior-data")
message("    consistency more than data-driven recovery.")
message("  * R = ", R_reps, " gives a coverage SE of ",
        sprintf("%.1f%%", 100 * sqrt(0.95 * 0.05 / R_reps)),
        "; do not read differences smaller than that.")
