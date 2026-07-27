# Tier 3: the full coverage study. Env-gated -- never runs in an ordinary
# check(), on CRAN or otherwise.
#
#   OCCJSDM_SIMSTUDY=1                enable
#   OCCJSDM_SIMSTUDY_R=100            replicates per scenario (default 100)
#   OCCJSDM_SIMSTUDY_OUT=<dir>        where to write the results table
#
# At the defaults this is 10 scenarios x 100 replicates, roughly 23 h serial.
# For a parallel run use dev/simstudy/run_study.R instead, which spreads
# replicates across processes; this test exists so the study is reachable
# through the normal test entry point and so its thresholds live with the
# other tests rather than in a script.

test_that("coverage is near nominal across the scenario grid", {
  skip_on_cran()
  skip_if_not(nzchar(Sys.getenv("OCCJSDM_SIMSTUDY")),
              "set OCCJSDM_SIMSTUDY=1 to run the coverage study")

  R <- as.integer(Sys.getenv("OCCJSDM_SIMSTUDY_R", "100"))
  outdir <- Sys.getenv("OCCJSDM_SIMSTUDY_OUT", tempdir())

  scenarios <- simstudy_scenarios()
  rows <- do.call(rbind, lapply(scenarios, simstudy_scenario, R = R))
  summary_tbl <- simstudy_summarise(rows)

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  saveRDS(list(rows = rows, summary = summary_tbl, R = R,
               when = Sys.time(), sessionInfo = utils::sessionInfo()),
          file.path(outdir, sprintf("simstudy-%s.rds", stamp)))
  utils::write.csv(summary_tbl,
                   file.path(outdir, sprintf("simstudy-%s.csv", stamp)),
                   row.names = FALSE)
  message("simstudy results written to ", outdir)

  # Blocks known to sit below nominal *before* this run, so a failure here
  # reports something new rather than re-reporting a known issue. Measured at
  # R = 5 across three seed sets: beta_theta 0.78-0.80, resid_cor 0.76-0.79.
  # Whether that is real is precisely what a full run should settle -- so they
  # are held to a lower floor rather than excluded, and the summary table is
  # the artefact to actually read.
  known_low <- c("beta_theta", "resid_cor")

  ordinary <- summary_tbl[!summary_tbl$block %in% known_low, ]
  watched  <- summary_tbl[ summary_tbl$block %in% known_low, ]

  # 0.85 against a nominal 0.95 (PLAN.md 5.2): at R = 100 the SE is ~2.2%, so
  # this is >4 SE below nominal -- loose enough never to flake, tight enough
  # to catch anything of the magnitude in the Fixed bugs list, which drove
  # coverage to 0.5-0.8.
  if (nrow(ordinary) > 0) {
    bad <- ordinary[ordinary$coverage < 0.85, ]
    expect_true(
      nrow(bad) == 0,
      info = paste0("blocks below 0.85 coverage: ",
                    paste(sprintf("%s/%s=%.3f", bad$scenario, bad$block,
                                  bad$coverage), collapse = ", ")))
  }
  if (nrow(watched) > 0) expect_gt(min(watched$coverage), 0.55)

  # The table itself is the deliverable. Snapshotting makes a change show up
  # as a reviewable diff rather than a bare pass/fail.
  expect_snapshot_value(
    summary_tbl[order(summary_tbl$scenario, summary_tbl$block),
                c("scenario", "block", "coverage")],
    style = "json2", tolerance = 0.02)
})
