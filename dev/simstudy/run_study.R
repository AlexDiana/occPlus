#!/usr/bin/env Rscript
#
# Tier-3 runner for the simulation study. See dev/simstudy/PLAN.md.
#
#   Rscript dev/simstudy/run_study.R [--R=100] [--cores=8] [--scenarios=a,b]
#                                    [--out=dev/simstudy/results]
#                                    [--nburn=1000] [--niter=1000] [--resume]
#                                    [--threads=1] [--caffeinate]
#
# --cores is how many replicates run at once, as separate processes.
# --threads is how many threads each individual fit uses, and defaults to 1;
# see the note where it is parsed for why that default is not negotiable while
# TODO.md group B items 8 and 9 are open.
#
# --resume picks up from the checkpoint written by a previous invocation with
# the same --scenarios and --out, skipping replicates already completed.
# --caffeinate (macOS only) stops the machine sleeping mid-run.
#
# Replicates are parallelised across *processes*, which is free: each is an
# independent dataset and fit, needing no package changes -- unlike the
# in-package chain parallelisation in TODO.md group D item 1.
#
# PSOCK rather than mclapply(), for the reasons in TODO.md D.1: mclapply() is
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

# Threads *inside* each fit, as distinct from --cores, which is how many
# replicate processes run at once. Default 1, deliberately.
#
# This used to be unset, which meant each of the --cores worker processes also
# spawned a full TBB pool, and the results carried whatever the default
# happened to be. Two reasons that is the wrong default here. Bit-reproducibility
# is only available at one thread, because TBB work-stealing varies which
# thread draws for which species. And TODO.md group B items 8 and 9 are live
# defects that fire only above one thread -- every thread seeding identically,
# and an unsynchronised race on R's global RNG -- so a multi-threaded study
# would be measuring the statistics of a known-broken sampler.
#
# AGENTS.md has long claimed every study number was measured at one thread.
# Nothing enforced it until now; this line is what makes the claim true.
threads  <- as.integer(getopt("threads", "1"))

resume   <- any(args == "--resume")
keepawake<- any(args == "--caffeinate")

pkg_root <- normalizePath(getopt("pkg", "."), mustWork = TRUE)

# Optionally hold the machine awake for the duration. Scoped to this process
# via -w, so the assertion dies with the run and there is no power setting to
# undo, even if the run crashes.
#
# A convenience, not the safety net -- checkpointing is what makes an
# interrupted run recoverable. This only avoids wasting wall-clock: on
# 28 July a sleeping laptop cost ~80 min of elapsed time with no work done,
# and the sleep was invisible in the progress output at the time.
if (keepawake) {
  if (Sys.info()[["sysname"]] == "Darwin") {
    system2("caffeinate", c("-is", "-w", Sys.getpid()), wait = FALSE)
    message("  caffeinate: holding the machine awake until this run exits")
    message("  (a closed lid on battery still sleeps; caffeinate cannot override that)")
  } else {
    message("  --caffeinate is macOS-only; ignored on ", Sys.info()[["sysname"]], ".")
    message("  Linux equivalent: systemd-inhibit --what=idle Rscript ...")
  }
}

# --- load package + helpers ----------------------------------------------

suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))
source(file.path(pkg_root, "tests", "testthat", "helper-simstudy.R"))

# Pin the in-fit thread count. Set the environment variable as well as calling
# setThreadOptions(), because the variable is what propagates to the PSOCK
# workers' fresh R sessions, where the fits actually run; the call alone would
# only bind this master process, which does no fitting.
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = as.character(threads))
RcppParallel::setThreadOptions(numThreads = threads)

scenarios <- simstudy_scenarios()
if (nzchar(which_sc)) {
  keep <- strsplit(which_sc, ",")[[1]]
  scenarios <- Filter(function(s) s$label %in% keep, scenarios)
  if (length(scenarios) == 0) stop("no scenario matched --scenarios=", which_sc)
}

MCMCparams <- list(nchain = 2, nburn = nburn, niter = niter, nthin = 1)

message(sprintf(
  "simstudy: %d scenario(s), R = %d, %d core(s), %d thread(s)/fit, MCMC %d+%d x 2",
  length(scenarios), R_reps, n_cores, threads, nburn, niter))

# --- run -----------------------------------------------------------------

# Flatten to (scenario, replicate) jobs so the pool stays busy: scenarios have
# very different costs (binary is far cheaper than species_30), and running
# them one at a time would leave cores idle at the tail of each.
jobs <- do.call(rbind, lapply(seq_along(scenarios), function(i) {
  data.frame(scenario_i = i, replicate = seq_len(R_reps))
}))
jobs <- jobs[sample.int(nrow(jobs)), ]   # mix costs across workers

# simstudy_replicate() returns list(params=, occupancy=): the per-element
# coverage rows, and the per-species occupancy-recovery rows. They are carried
# separately all the way to disk because their schemas differ (see the note on
# simstudy_occupancy_rows()); nothing here ever rbinds one onto the other.
run_one <- function(k) {
  s <- scenarios[[jobs$scenario_i[k]]]
  r <- jobs$replicate[k]
  out <- try(simstudy_replicate(s, r, MCMCparams = MCMCparams), silent = TRUE)
  if (inherits(out, "try-error")) {
    # One bad replicate must not lose the whole run. Record and continue.
    # Only the params frame carries the ERROR marker, which is what the
    # failure count and the drop below key on.
    return(list(
      params = data.frame(scenario = s$label, replicate = r, block = "ERROR",
                          element = NA_integer_, truth = NA_real_,
                          lower = NA_real_, upper = NA_real_,
                          post_mean = NA_real_, covered = NA_integer_,
                          stringsAsFactors = FALSE),
      occupancy = NULL))
  }
  out
}

# Checkpointing. Results used to exist only in the master's memory until the
# run finished, so a job killed at hour ten lost ten hours of completed
# replicates. Each chunk is now appended to a checkpoint file as it lands, and
# a run started with --resume skips work already recorded there.
#
# Keyed on (scenario, replicate), which is safe because simstudy_replicate()
# seeds deterministically from exactly that pair -- a resumed replicate is the
# same replicate, not a fresh draw.
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ckpt_path <- file.path(out_dir, sprintf("checkpoint-%s.rds",
                                        if (nzchar(which_sc))
                                          substr(gsub("[^A-Za-z0-9]+", "-", which_sc), 1, 40)
                                        else "all"))

ckpt_load <- function() {
  if (!file.exists(ckpt_path)) return(NULL)
  out <- try(readRDS(ckpt_path), silent = TRUE)
  if (inherits(out, "try-error")) {
    # A checkpoint truncated by a hard kill must not abort the run.
    message("  checkpoint unreadable, ignoring: ", ckpt_path)
    return(NULL)
  }
  out
}

ckpt_append <- function(bundle) {
  prev <- ckpt_load()
  saveRDS(list(params    = rbind(prev$params,    bundle$params),
               occupancy = rbind(prev$occupancy, bundle$occupancy)),
          ckpt_path)
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
# ETA comes from a TRAILING WINDOW, not the cumulative average. A cumulative
# average is badly misleading whenever throughput changes: on 28 July the
# machine slept mid-run and the log kept reporting "~28 min left" while the
# recent rate implied several times that.
#
# cpu/wall is printed for the same reason. Compared against `cores` it says
# how much of the wall time was actually spent computing; a ratio collapsing
# toward zero while elapsed climbs is the unambiguous signature of a suspended
# machine, which is otherwise hard to tell from thermal throttling by eye.
.prog <- new.env(parent = emptyenv())
.prog$hist <- numeric(0)

progress_line <- function(done, total, t0, failed, cpu_min = NA_real_) {
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  frac <- done / total
  .prog$hist <- c(.prog$hist, el)

  h <- .prog$hist
  win <- utils::tail(h, 6L)
  per_report <- if (length(win) >= 2L)
    (win[length(win)] - win[1]) / (length(win) - 1L) else el
  reports_left <- (total - done) / (done / length(h))
  eta <- per_report * reports_left

  cpu_txt <- if (is.finite(cpu_min) && el > 0)
    sprintf(" | cpu/wall %.1fx", cpu_min / el) else ""

  message(sprintf(
    "  [%3.0f%%] %d/%d replicates | %.0f min elapsed | ~%.0f min left%s%s",
    100 * frac, done, total, el, eta, cpu_txt,
    if (failed > 0) sprintf(" | %d failed", failed) else ""))
  flush(stderr())
}

# Total CPU-minutes burned by the worker processes so far.
cluster_cpu_min <- function(cl) {
  out <- try(clusterEvalQ(cl, sum(proc.time()[1:2])), silent = TRUE)
  if (inherits(out, "try-error")) return(NA_real_)
  sum(unlist(out)) / 60
}

# Drop jobs already recorded in the checkpoint.
done_bundle <- if (resume) ckpt_load() else NULL
done_rows <- done_bundle$params
if (!is.null(done_rows)) {
  have <- unique(paste(done_rows$scenario, done_rows$replicate, sep = "\r"))
  want <- paste(vapply(jobs$scenario_i, function(i) scenarios[[i]]$label, ""),
                jobs$replicate, sep = "\r")
  keep <- !(want %in% have)
  message(sprintf("  resuming: %d of %d replicates already in %s",
                  sum(!keep), length(keep), basename(ckpt_path)))
  jobs <- jobs[keep, , drop = FALSE]
  if (nrow(jobs) == 0L) message("  nothing left to run")
}

n_jobs <- nrow(jobs)
t0 <- Sys.time()

if (n_cores > 1L) {
  cl <- makeCluster(n_cores)                      # PSOCK; portable
  on.exit(stopCluster(cl), add = TRUE)
  # Independent L'Ecuyer streams per worker, on top of the per-replicate seed
  # in simstudy_replicate(). Reproducible and non-overlapping.
  clusterSetRNGStream(cl, 20260727L)
  clusterExport(cl, c("scenarios", "jobs", "MCMCparams", "pkg_root", "threads"),
                envir = environment())
  clusterEvalQ(cl, {
    # Pin before load_all(), so the setting is in force for every fit this
    # worker runs. Belt and braces with the exported environment variable:
    # whichever RcppParallel consults, both say the same number.
    Sys.setenv(RCPP_PARALLEL_NUM_THREADS = as.character(threads))
    RcppParallel::setThreadOptions(numThreads = threads)
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
    ckpt_append(simstudy_bind(res[[ci]]))
    done <- done + length(chunks[[ci]])
    failed_so_far <- failed_so_far +
      sum(vapply(res[[ci]], function(d) any(d$params$block == "ERROR"),
                 logical(1)))
    progress_line(done, n_jobs, t0, failed_so_far, cluster_cpu_min(cl))
  }
  res <- unlist(res, recursive = FALSE)

} else {
  res <- vector("list", n_jobs)
  failed_so_far <- 0L
  for (k in seq_len(n_jobs)) {
    res[[k]] <- run_one(k)
    ckpt_append(res[[k]])
    if (any(res[[k]]$params$block == "ERROR")) failed_so_far <- failed_so_far + 1L
    if (k %% 10L == 0L || k == n_jobs) progress_line(k, n_jobs, t0, failed_so_far)
  }
}
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

bundle <- simstudy_bind(res)
rows <- bundle$params
occ  <- bundle$occupancy
# Include anything carried over from an earlier, interrupted run.
if (!is.null(done_bundle)) {
  rows <- rbind(done_bundle$params, rows)
  occ  <- rbind(done_bundle$occupancy, occ)
}
failed <- sum(rows$block == "ERROR", na.rm = TRUE)
rows <- rows[rows$block != "ERROR", , drop = FALSE]
summary_tbl <- simstudy_summarise(rows)
occ_summary <- simstudy_summarise_occupancy(occ)
occ_prev    <- simstudy_occupancy_by_prevalence(occ)

# --- provenance ----------------------------------------------------------
#
# Which code produced these numbers, recorded in the results file itself.
#
# This exists because of a specific failure. The validation article was written
# from the 2 August run and then six commits touched R/ and src/ underneath it,
# including a new parallel sampler; nothing in the results file said what it had
# been run against, so the article went stale invisibly and the drift was only
# found by comparing timestamps against the git log by hand. The article now
# prints this block, so a reader can see the answer instead of assuming it.
git_out <- function(...) {
  v <- suppressWarnings(try(
    system2("git", c("-C", pkg_root, ...), stdout = TRUE, stderr = FALSE),
    silent = TRUE))
  if (inherits(v, "try-error") || length(v) == 0L) NA_character_ else v[1]
}

provenance <- list(
  git_sha    = git_out("rev-parse", "HEAD"),
  git_branch = git_out("rev-parse", "--abbrev-ref", "HEAD"),
  # A dirty tree means the SHA above does not fully describe what ran, which
  # matters more than the SHA itself for reproducing a result.
  git_dirty  = {
    st <- suppressWarnings(try(
      system2("git", c("-C", pkg_root, "status", "--porcelain"),
              stdout = TRUE, stderr = FALSE), silent = TRUE))
    if (inherits(st, "try-error")) NA else length(st) > 0L
  },
  pkg_version = as.character(utils::packageVersion("occJSDM")),
  r_version   = paste(R.version$major, R.version$minor, sep = "."),
  platform    = R.version$platform,
  n_cores     = n_cores,
  # Every fit here is single-threaded within itself; the parallelism is across
  # replicate processes. Recorded because it is the reason these results say
  # nothing about the multi-threaded configuration a user gets by default
  # (TODO.md group B items 8 and 9).
  rcpp_parallel_threads = threads,
  # Read back from a worker rather than trusting the master's own setting:
  # the fits happen there, and a pin that failed to propagate is exactly the
  # kind of thing that would otherwise go unnoticed.
  rcpp_parallel_threads_in_worker = if (n_cores > 1L)
    tryCatch(clusterEvalQ(cl, Sys.getenv("RCPP_PARALLEL_NUM_THREADS"))[[1]],
             error = function(e) NA_character_) else as.character(threads),
  scenarios   = vapply(scenarios, function(s) s$label, ""),
  nburn = nburn, niter = niter, R = R_reps,
  elapsed_min = elapsed,
  started  = t0,
  finished = Sys.time()
)

# --- write ---------------------------------------------------------------

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
saveRDS(list(rows = rows, summary = summary_tbl,
             occupancy = occ, occupancy_summary = occ_summary,
             occupancy_by_prevalence = occ_prev,
             R = R_reps,
             MCMCparams = MCMCparams, failed_replicates = failed,
             elapsed_min = elapsed, when = Sys.time(),
             provenance = provenance,
             sessionInfo = utils::sessionInfo()),
        file.path(out_dir, sprintf("simstudy-%s.rds", stamp)))
write.csv(summary_tbl, file.path(out_dir, sprintf("simstudy-%s.csv", stamp)),
          row.names = FALSE)
write.csv(occ_summary,
          file.path(out_dir, sprintf("simstudy-occupancy-%s.csv", stamp)),
          row.names = FALSE)

message(sprintf("\ndone in %.1f min; %d replicate(s) failed", elapsed, failed))
if (failed > 0) message("  (failures are dropped from the summary, not silently counted as covered)")
print(summary_tbl, row.names = FALSE)
message("\nwritten to ", normalizePath(out_dir))

# Reminders that outlive this session -- see PLAN.md 5.3 and 10.
message("\nRead with care:")
message("  * l_s and sigma_h are NOT in the table: never sampled / not recoverable")
message("    (TODO.md group A items 1 and 2). Spatial cells say nothing about range.")
message("  * sigma_b is prior-dominated; its coverage reflects prior-data")
message("    consistency more than data-driven recovery.")
message("  * R = ", R_reps, " gives a coverage SE of ",
        sprintf("%.1f%%", 100 * sqrt(0.95 * 0.05 / R_reps)),
        "; do not read differences smaller than that.")
