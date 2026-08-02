# Simulation-based validation of `occJSDM`

Plan for the "extensive testing on simulated datasets" item under *MEE paper / Doug to dos* in `TODO.md`. Drafted 27 July 2026.

[Guide to the test suite](https://claude.ai/code/artifact/ad3d46eb-1fd4-49b5-b795-6b71474ef1d5 "Guide to the test suite")

**Status as of 31 July 2026: the suite is built, and the grid plus four targeted experiments have been run.** The pre-fix grid (27-28 July) found four defects. The post-fix re-run (29 July) is the current table in §12 and added two more findings. Four targeted experiments followed: the M ladder (§13), the tighter-`beta_theta`-prior arms (§13.8), the prior-variance arms (§14), and the `ncov_theta = 0` discriminator (§15, run at R = 50 and confirmed at R = 200). Between them they disproved three hypotheses, narrowed one finding to a specific parameter row, and traced most of another. Sections marked *OPEN* still need a decision or a measurement.

**Reading order if you are new to this:** §12 for the headline results, then §15.4-§15.6, which are the most recent and supersede parts of §12's attributions. §13.7 and §14.7 are the intermediate experiments. §1-§9 are the design and are mostly settled.

**Two findings are open and both are now sharply localised.** `beta_theta` (the *slope overconfidence* item in `TODO.md` group B) undercovers only on the collection-covariate *slopes*; its intercept row covers at 0.968 (§15.5). `B0` (the *bias doubled* item) carries a bias of which 70% is downstream of those slopes and the rest is a real residual, 3.3 SE from zero at R = 200 (§15.6). Both now need a sampler-level look rather than more simulation.

**A note on citing TODO items.** Group B items are cited here by *subject*, never by position. Positions have been renumbered four times and every numeric reference broke each time, in some cases silently repointing at an unrelated item. See "Cross-referencing TODO items" in `AGENTS.md`.

| Stage | State |
|----|----|
| 1\. Fixtures + tier 1 | **done** (`6d9526d`) -- 89 assertions, 7 s; 112 as of 29 July |
| 2\. Spatial `n` floor | **done** (`41fb736`) -- floor is 31 sites |
| 3\. `helper-simstudy.R` | **done** (`6123036`) -- validated at R = 8 |
| 4\. Tier 2 canary | **done** (`f63eeeb`) -- \~30 s, `skip_on_cran()` |
| 5\. Tier 3 + runner | **done** (`f63eeeb`) -- runner smoke-tested end to end |
| 6\. Full R = 100 run | **done, run twice** -- pre-fix 27-28 July (1000 fits across two runs), post-fix 29 July (1000 fits, 285 min, one run). 0 failures throughout. Results in `dev/simstudy/results/`. See §12. |

Building it surfaced four corrections to this plan and two package bugs; both are recorded in place below rather than only in commit messages.

------------------------------------------------------------------------

## Running the suite

Six levels, cheapest first. Times are wall clock on Doug's machine and include \~9 s of R startup and `load_all()` compilation, so the test work itself is correspondingly less.

| Level | Command | Time | Covers |
|----|----|----|----|
| One file, while iterating | `Rscript -e 'devtools::test(filter="regression")'` | \~11 s | one tier-1 file |
| **Everything local** | `Rscript -e 'devtools::test()'` | \~45 s | tiers 1 + 2 |
| Exactly what CRAN runs | `NOT_CRAN=false Rscript -e 'devtools::test()'` | \~16 s | tier 1 only |
| Full package check | `Rscript -e 'devtools::check()'` | minutes | tests + examples + vignettes + `R CMD check` |
| One study cell | `Rscript dev/simstudy/run_study.R --cores=5 --scenarios=base --caffeinate` | \~50 min | tier 3, one cell at R = 100 |
| Full study | `Rscript dev/simstudy/run_study.R --cores=5 --caffeinate` | \~4.75 h | tier 3, all 10 cells at R = 100 |

Study times are **measured**, not projected: the 29 July run did all 10 cells, 1000 fits, in **285 min** on 5 cores. Ignore the 474 min recorded for the 28 July grid -- that wall time was inflated by the laptop sleeping mid-run; its compute rate was the same \~4 fits/min. Throughput is flat across the whole run, so the cost scales linearly in R. An earlier version of this table said `--cores=8` and \~2.4 h, which was a guess and roughly 3x optimistic. On a fanless MacBook Air, 5 cores is the practical ceiling; use `--caffeinate` for anything over an hour, because a run that sleeps mid-way is easy to misread as thermal throttling (it was, once -- check the cpu/wall ratio in the progress line, which makes sleep unmistakable). `--resume` picks a run back up from the last checkpoint.

**`devtools::test()` is the day-to-day command** -- run it after any change to `R/` or `src/`. In RStudio, Cmd+Shift+T. `filter` is a regex on the filename minus the `test-` prefix; the tier-1 files are `smoke`, `regression` and `api`.

**After touching C++**, if you have been switching branches or the object files look stale:

```         
Rscript -e 'pkgbuild::clean_dll(); devtools::test()'
```

**Before pushing**, use `devtools::check()` rather than `devtools::test()`. It is the only way to see the tests as CRAN will run them: `test()` uses `load_all()`, whereas `check()` installs the package properly.

**The study runner never runs from `devtools::test()`.** Tier 3 is gated on `OCCJSDM_SIMSTUDY`; setting it runs the study *serially* through the test entry point, which is \~19 h. Use `run_study.R` for any real run -- it parallelises across replicates. `--R=`, `--cores=`, `--scenarios=`, `--out=`, `--nburn=` and `--niter=` are all overridable.

### Two things to know before changing a test

**A failing test names its bug.** The regression tests are titled `"Fixed bugs 22: ..."`, `"Fixed bugs 8: ..."` and so on, keyed to the *Fixed bugs* list in `TODO.md`. Read that entry before debugging -- it says what broke, when, and why it mattered.

**Do not "fix" a loose assertion by tightening it.** Several are deliberately loose and the comment above each says why. The GP length-scale test asserts that `sample_ls()` is *reached*, not that `idx_ls` varies, precisely so it does not encode the open `sigma_s` bug (§10.0); tightening it would make it pass for the wrong reason now and fail confusingly once that bug is fixed. The shuffle-invariance test asserts group totals rather than row-for-row equality because PCR replicates are exchangeable, and ignores `row.names` because those were verified to be **order-dependent but pairing-preserving** -- a user looking up by name gets the right row either way, so re-tightening it would fail on a difference that carries no information. The tier-2 thresholds are set from measurement across three seed sets (§6.3), not from taste.

------------------------------------------------------------------------

## 1. Purpose

Check that `runOccJSDM()` recovers known generating parameters from data produced by `simulateOccJSDMData()`, at the nominal rate, across the configurations the package advertises.

Two distinct goals, which need different designs -- worth keeping separate in the writing-up:

- **Regression.** Guard the \~23 bugs already fixed (`TODO.md`, *Fixed bugs*). None of them currently has a test; every one was verified by reading code, or by a one-off manual run. This is the near-term value.
- **Calibration evidence.** Support a claim in the MEE paper that inference is correctly calibrated. Needs more replicates and a tighter argument (§9).

The suite as specified below serves the first goal well and the second adequately. Read §5 before quoting any coverage number in the paper.

------------------------------------------------------------------------

## 2. Decisions settled

| Question | Decision |
|----|----|
| Vehicle | **testthat suite**, not a vignette or pkgdown article |
| Does it ship to CRAN? | **Tier 1 yes**, tiers 2–3 no |
| Replicates | **R = 100** per scenario |
| Statistic | **Coverage** for v1, structured so **SBC** can be added later (§7) |
| Scenario grid | **10 cells** (§6.1) |

### Why a testsuite rather than a vignette or article

- The `sample_ls` regression (`TODO.md` Fixed bugs 22) broke *every* non-spatial fit and survived because nothing exercised that path. A three-line test would have caught it the moment it landed. The vignette would not have -- it passes `spatCovariates`.
- It discharges Phase 3 step 11 of the CRAN plan in `AGENTS.md` ("write the testthat suite"). Done as of `6d9526d`: `test-placeholder.R` is gone, replaced by four real test files.
- No package-size or vignette-build cost. Both matter here: `sampleresults.rda` at 62 MB is already a CRAN blocker.

What this gives up is *presentation*. Tests answer "does it pass?", not "here is the calibration evidence". See §9.

------------------------------------------------------------------------

## 3. Tiers

| Tier | Gate | Runtime | Runs on |
|----|----|----|----|
| 1\. Regression + smoke | none | \< 30 s | every `check()`, including CRAN |
| 2\. Recovery canary | `skip_on_cran()` | \~2 min | local + CI |
| 3\. Coverage study | env var `OCCJSDM_SIMSTUDY` | \~4.75 h | manual / nightly |

Tier 3 cannot be an ordinary test: CRAN wants the whole check under \~10 minutes.

------------------------------------------------------------------------

## 4. Layout

```         
tests/testthat.R                      (exists)
tests/testthat/
  helper-fixtures.R       small cached fits shared across tier-1 tests
  helper-simstudy.R       scenario defs, simulate -> fit -> extract, metrics
  test-smoke-configs.R    tier 1: every model configuration runs
  test-regression-bugs.R  tier 1: one test per Fixed bugs entry
  test-api-contracts.R    tier 1: argument handling and error paths
  test-recovery.R         tier 2, skip_on_cran()
  test-coverage-study.R   tier 3, env-gated
dev/simstudy/
  PLAN.md                 this file
  run_study.R             tier-3 runner (parallel across replicates)
  results/                output tables, git-ignored
```

------------------------------------------------------------------------

## 5. Design constraints

These are the things most likely to make the suite wrong or annoying. Each has bitten something already.

### 5.1 Never assert exact numeric equality of fitted values in tier 1

**Tier 1 assertions must be structural**: does it run, are dimensions right, is a block non-`NA`, does a value vary. Numeric recovery belongs in tiers 2–3, with tolerances.

**The original reason expired on 29 July 2026, and the rule survives it.** This rule was written because a fixed seed did not reproduce the sampler: `randinvg()` drew from R's global RNG inside an OpenMP loop, so a tier-1 `expect_equal()` on posterior quantities would have passed locally and flaked on CRAN. That is fixed (TODO.md Fixed bugs 28) -- `set.seed()` now controls the whole fit, verified across separate processes.

The rule stands anyway, for a different and more durable reason: **fitted values are not a contract**. Any legitimate change to a sampler, a prior default, or the order of draws moves every number, and a suite full of golden values turns each such change into a wall of failures that says nothing about correctness. Structural assertions survive those changes; numeric ones do not.

The one sanctioned exception is `test-regression-bugs.R`'s reproducibility test, where equality *is* the property under test. Note its scope: reproducibility holds for a given thread count, since threads derive separate streams (see `src/rng.h`). This is currently moot -- `SHLIB_OPENMP_CXXFLAGS` is empty on the dev machine, so the pragmas compile to no-ops -- but it would bite on a Linux build.

### 5.2 Coverage assertions are stochastic and will flake if set tight

With R = 100 the SE of an estimated 95% coverage is \~2.2%, so a *correct* implementation lands outside ±1.96 SE about 5% of the time. `expect_gt(cov, 0.90)` would fail intermittently for no reason, and a suite that cries wolf gets ignored.

Mitigations, all three:

- fix a seed per replicate, so a tier-3 run is reproducible;
- set the hard failure threshold well below nominal -- fail at **\< 0.85**, which catches every bug in the *Fixed bugs* list by a wide margin (§9);
- record the coverage table with `expect_snapshot_value()` so changes surface as a reviewable diff rather than a boolean.

### 5.3 Identifiability limits what can be checked at all

The latent-factor model is invariant to rotation and sign. `reparamFactorModel()` (`R/jsdmfun.R`) imposes a QR-based constraint, but there is no reason that constraint matches the simulator's parameterisation.

- **Directly comparable** (verified against the code while building `helper-simstudy.R`, 27 July): `B0`, `B`, `G`, `beta_theta`, `theta0`, `p`, `q`, `tau`. Two caveats. `beta_theta` requires passing `collCovariates` to the fit -- without it the posterior is `(1 x S)` against a `(ncov_theta+1 x S)` truth. And `n_lattrait` must be pinned to the simulated `gt`, since the default `floor(sqrt(min(S, ncov_psi)))` will not generally match, giving `A`/`C` mismatched dimensions.
- **`sigma_b`: keep, but read with care.** With a tight `InvGamma(10, 1)` prior and only `p x S` residual coefficients the posterior sits at the prior mean (`sqrt(1/9)` = 0.333), so a fixed truth far from that yields near-zero coverage through prior-data conflict, not sampler error. Measured: true 0.5 gave 0/8, true 0.35 gave 8/8. Scenarios set it near the prior mean; SBC is the principled fix.
- **NOT comparable, contrary to the original draft of this section:** `Bs`, `Gs`, `sigma_bs`. The simulator builds the spatial field directly from `sigma_s`/`l_s` and sets `Bst <- matrix(0, S, ps)` (`R/jsdmfun.R:557`), so on the truth side `Bs` is empty *regardless of `ds`* (re-checked at `ds = 2`) and `sigma_bs` generates nothing; the fit meanwhile represents the field as sparse-GP basis coefficients over `n_supportpoints` knots. Different parameterisations. Measured for `sigma_bs`: true 0.5 against a posterior mean of \~1.6, 0/8 coverage -- a meaningless comparison rather than a bug.
- **Only as identified functions:** residual correlation matrix, `eta`/`psi`, variance partitioning
- **Must NOT be checked element-wise:** `U`, `L`, `A`, `C` -- element-wise coverage would fail for reasons that are not bugs
- **Excluded:** `sigma_h` -- it was not sampled when this was written (since fixed, *Fixed bugs* 24), but it remains excluded because `U` at training sites is drawn under a hard-coded unit-variance prior regardless, so the study cannot see it. **Also `l_s`** -- the `sample_ls()` item in group B: `sigma_s` is hard-coded to 1 at the `sample_ls()` call site, so the length-scale absorbs the amplitude misspecification and rails at the top of `l_s_grid` for every true value tried. Any coverage figure for `l_s` is meaningless until that is fixed. Exclude explicitly and with a comment, so a future reader does not "fix" the test.

### 5.4 The GP knot count must be pinned, not left to the default

`getDefaultSupportPoints(n)` used to be `max(30, floor(n * 0.2))` -- a constant 30 for any dataset below 150 sites, and a crash below 31 (§10.1). Fixed in `42198d9`; it is now `min(floor(n * 0.2), n - 1)` (`R/jsdmfun.R:875`, `TODO.md` Fixed bugs 29). The study pins `n_supportpoints` anyway rather than trusting any default.

Two consequences for the design, beyond the crash. The knot count would be **uncontrolled but not constant** across cells -- identical for `n = 40` and `n = 140`, then suddenly proportional above 150 -- so any difference attributed to `n` would be partly an artefact of the spatial approximation changing underneath. And at small `n` it approaches one knot per site, which is a qualitatively different (denser) approximation than at `n = 100`, precisely in the low-information cell where calibration is most at risk.

Set `n_supportpoints` explicitly in every spatial cell so it is a controlled factor.

### 5.5 CRAN limits

Two cores maximum (`_R_CHECK_LIMIT_CORES_`). Tier 3's replicate-level parallelism must be off or capped when not env-gated on.

------------------------------------------------------------------------

## 6. Tier detail

### 6.1 Tier 3 -- the study

Base = `two_stage`, traits on, spatial on, `d = 2`, **`ds = 2`**, `n = 100`, `S = 10`, `M = 2`, `P = 2`, `K = 3`, `n_supportpoints = 20`.

`ds >= 1` *used to be* required or there was no spatial field at all: at `ds = 0` the simulator's cross-species spatial covariance collapsed to jitter, `sd(spatField)` = 0.0019 rather than \~1.0. The grid used `ds = 0` until 27 July, which made every spatial cell a null-field test -- the exact trap the audit flagged for Fixed bugs 7. Fixed in `42198d9` (`TODO.md` Fixed bugs 30); re-measured 29 July, `ds = 0` now gives 0.598. The grid stays at `ds = 2` by choice, for a clearly non-degenerate field and comparability with the 28 July run.

Every spatial cell sets `n_supportpoints` explicitly (§10.1). Left to the default it is a constant 30 for any `n` below 150, so it would neither scale with `n` nor be a controlled factor -- and below 31 sites it crashes outright.

| \# | Cell | Differs from base by | Buys |
|----|----|----|----|
| 1 | base | -- | Flagship config; *is* the traits×spatial interaction |
| 2 | spatial isolated | traits off | `Bs` recovery unconfounded by trait terms (**not** `l_s` -- see §5.3) |
| 3 | traits isolated | spatial off | `G` (fourth-corner) recovery unconfounded by the GP |
| 4 | `P = 3` | `P = 2 -> 3` | Guards Fixed bugs 2; matches shipped `sampledata` |
| 5 | low information | low `p`/`theta_baseline`, smaller `n` (see §10.1 -- must set `n_supportpoints`) | Calibration where data are thin -- the realistic eDNA regime |
| 6 | `d` under-fit | truth `d = 4`, fit `d = 2` | The common real-world misspecification |
| 7 | `d` over-fit | truth `d = 2`, fit `d = 4` | Interacts with rotation invariance; least-trusted case |
| 8 | larger `S` | `S = 10 -> 20` | Does calibration hold as species count grows. Trimmed from 30 on 27 July: at `S = 30` this was 1.75x base, the grid's most expensive cell, and the study runs on a fanless machine that throttles. `S = 20` is still double the base and costs 1.2x. |
| 9 | `occupancy` | no PCR stage | Calibration evidence for the collapsed form |
| 10 | `binary` | no replicates at all | Calibration evidence for pure JSDM |

Design is **one-factor-at-a-time**, not factorial (\~48 cells would be unaffordable). Cost: OFAT cannot see interactions. Cell 1 covers the one interaction judged most likely to matter, since traits and spatial both write into the same linear predictor via `computePsiCoef()`.

Deliberately *excluded* from tier 3, because tier 1 answers them in seconds for free: "does configuration X run at all". Tier 3 hours should buy calibration evidence, which only tier 3 can provide.

Cells 9–10 are **cheaper** (less data per fit: no PCR level, no replication) but also check **fewer parameters** -- `occupancy` has no `p`/`q`; `binary` has no `beta_theta`, `theta0`, `p`, `q`. Their rows in the results table will be sparser, not directly comparable to the two_stage cells.

**Metrics** per scenario × parameter block, pooled across species: coverage of the 95% credible interval, bias, RMSE.

**Output:** a tibble (`scenario`, `block`, `coverage`, `bias`, `rmse`, `n_par`, `R`) written to `Sys.getenv("OCCJSDM_SIMSTUDY_OUT", tempdir())`, intended destination `dev/simstudy/results/`.

### 6.2 Tier 1 -- ships to CRAN

One shared fixture fit built in `helper-fixtures.R` and reused: `sampledata`, `nburn = 50, niter = 50, nchain = 2` (**measured 5.6 s**). Keeps the tier well under 30 s.

`test-smoke-configs.R` -- every configuration completes: `two_stage` / `occupancy` / `binary`; with and without `spatCovariates`; with and without traits; `d = 0` and `d = 2`; `P = 1` and `P = 3`. *This is the tier that would have caught the `sample_ls` regression.*

`test-regression-bugs.R` -- one test per fixed bug, named so a failure points at its `TODO.md` entry:

| Guards           | Assertion                                           |
|------------------|-----------------------------------------------------|
| Fixed bugs 22    | non-spatial fit runs                                |
| Fixed bugs 8     | `tau_output` not all `NA`                           |
| Fixed bugs 21    | `plotFPTPStage2Rates()` differs across primers      |
| Fixed bugs 23    | `computeSpeciesDetected()` indices span both chains |
| Fixed bugs 9     | WAIC counter guard does not trip                    |
| Fixed bugs 1, 11 | row-shuffle invariance                              |

Row-shuffle invariance must be asserted on the **prepared design matrices** (`X_psi`/`X_theta` row alignment against site/sample indices), which is deterministic -- not on the posterior, which would flake under §5.1.

`test-api-contracts.R` -- unknown `primerName` errors with the available list; unknown covariate name errors informatively; `idx_species` subsetting returns the right rows; `predictNewSites()` honours its documented `NULL` defaults.

### 6.3 Tier 2 -- canary

R = 3–5 replicates at the base cell. Not a calibration check: a smoke signal that recovery has not grossly broken between full tier-3 runs (e.g. correlation between estimated and true `B` above a loose floor).

**RESOLVED: tier 2 fails.** Thresholds were set from measurement rather than taste, which is what makes failing safe. Three independent seed sets gave overall coverage 0.89/0.89/0.88, lowest block 0.78/0.76/0.76, and `cor(post_mean, truth)` 0.62/0.56/0.68; the floors are 0.70, 0.40 and 0.30. A pass is therefore not luck and a failure is not noise.

------------------------------------------------------------------------

## 7. SBC seam

Coverage now, simulation-based calibration later, without restructuring. Two swappable functions in `helper-simstudy.R`:

``` r
draw_truth(scenario)        # fixed values now; prior draws for SBC
summarise_fit(fit, truth)   # CI-contains indicator now; rank statistic for SBC
```

Everything expensive -- grid, simulate/fit loop, parallelism, artifact writing -- is shared.

Adding SBC later requires:

1.  A prior-draw wrapper. All priors in the sampler are **proper** (checked): `B0`, `L` \~ N(0,1) (`sample_BBsL()`, `R/jsdmfun.R`); `B` trait-predicted mean with variance `sigma_b^2`; `sigma_b`, `sigma_bs` \~ sqrt-InvGamma(10,1); `tau` \~ InvGamma(5,5) (`R/runOccJSDM.R`); `p`, `q`, `theta0` \~ Beta; `l_s` \~ Gamma(1,1) over a 10-point grid. So SBC is feasible.
2.  Documenting what the sampler *actually* assumes. Not free: the audit found the *documented* `prior_beta_psi`/`prior_beta_psi_sd` were a complete no-op, with the live path hard-coding something different (Fixed bugs 14).
3.  Thinning to near-independent draws, which raises per-replicate cost.

Note `sigma_h` would show a flagrantly non-uniform rank histogram under SBC. That is the method working, not a new bug.

------------------------------------------------------------------------

## 8. Runtime budget

Measured on `sampledata` (100 sites, 10 species, 2 chains): **0.039 s per iteration**, so `nburn = 1000 + niter = 1000` is \~84 s per fit.

|                   | Serial  | 5 cores      |
|-------------------|---------|--------------|
| One cell, R = 100 | \~2.3 h | \~28 min     |
| 10 cells          | \~23 h  | **\~4.75 h** |

**Superseded by measurement**: the 29 July run did all 10 cells in 285 min on 5 cores, i.e. \~4.75 h. The earlier \~2.5–3 h figure assumed 8 cores, which this machine cannot usefully supply -- the M4 has only 4 performance cores, so workers past the fourth land on efficiency cores worth roughly a third to a half as much (see §5.5). Cells 9–10 are cheaper than the average; cell 8 (`species_20`) is dearer.

Replicates are embarrassingly parallel at the *process* level -- independent datasets and fits -- so this needs no package changes, unlike the in-package chain parallelisation proposed under *MEE paper / Alex to dos*, "Performance of `runOccJSDM()`", option A.

**Do not buy replicates by shortening chains.** Each interval endpoint is a tail quantile; with 500 total draws only \~12 land below the 2.5% bound, so the endpoints are noisy and that noise feeds into coverage as bias, not just variance. Set chain length by an ESS target (\>= 400 on monitored parameters, which `returnConvergenceDiagnostics()` already reports) and let R take what is left.

------------------------------------------------------------------------

## 9. What R = 100 does and does not support

SE of estimated coverage is `sqrt(0.95 * 0.05 / R)`. Cost is linear in R; precision improves only as sqrt(R) -- halving the error bar costs 4x the compute.

> **Resolved 29 July 2026.** This previously carried a caveat that the replicates were not independent: `get_rng()` was seeded from a literal, per *process*, so every PSOCK worker started its Polya-Gamma stream at the same position. That is fixed (TODO.md Fixed bugs 28) -- `runOccJSDM()` now derives its C++ seed from R's RNG, and because `draw_truth()` sets a per-(scenario, replicate) seed before each fit, every replicate gets its own stream. Verified: a repeated replicate reproduces exactly, and two replicates differ on every column. `sqrt(0.95 * 0.05 / R)` is therefore the right SE for the re-run.
>
> §12 now reports the 29 July re-run, which was made *after* this fix, so its error bars are sound. It is the superseded pre-fix table at §12.6 that inherits the old correlation: treat those figures as indicative only, and never quote an SE against them.

| R       | Coverage SE | Undercoverage reliably flagged |
|---------|-------------|--------------------------------|
| 50      | 3.1%        | below \~0.890                  |
| **100** | **2.2%**    | **below \~0.907**              |
| 200     | 1.5%        | below \~0.920                  |
| 400     | 1.1%        | below \~0.929                  |

R = 100 is chosen for the **regression** goal. Every bug in the *Fixed bugs* list would drive coverage to 0.5–0.8, not 0.93, so R = 100 catches them with room to spare.

It is **not** enough to assert "coverage is nominal" in the paper -- that claim needs R = 200–500, because it asserts the *absence* of a small deviation rather than detecting a large one. The runner takes `R` as an argument for exactly this reason; do a final publication-grade run at higher R if the paper makes that claim.

Also: pooling coverage across species within a block buys precision, but those indicators share one simulated dataset and one MCMC run, so effective N sits between R and R x S, nearer R. **Do not quote an SE computed as if there were R x S independent trials.**

------------------------------------------------------------------------

## 10. Open items

0.  *OPEN, and it degrades the two spatial cells.* **The GP length-scale is never recovered** (the `sample_ls()` item in `TODO.md` group B). Needs a derivation, not a code tweak. Until then no cell says anything about spatial range.

1.  ~~Minimum `n` for spatial cells.~~ **RESOLVED 27 July.** The floor is **31 unique locations** with default settings.

2.  ~~Tier 2 failing vs advisory.~~ **RESOLVED: it fails**, on thresholds measured across three seed sets (§6.3).

3.  ~~`beta_theta` and `resid_cor` sit below nominal.~~ **RESOLVED and both traced.** `resid_cor` is the `reparamFactorModel()` item in `TODO.md` group B, confirmed by a paired re-run in which only 104 of 49,978 coverage decisions flipped. `beta_theta` was partly the prior mean (*Fixed bugs* 25); the residue is the slope-overconfidence item in group B, and simulation has now taken it as far as it can go. Ruled out: under-identification, prior width, pseudo-replication (§13, §14). Established: the defect sits on the **slopes** only, since the intercept row covers at 0.968 with the slopes removed (§15.5). What is left is a sampler-level question about the Polya-Gamma update, filed for Alex.

4.  *OPEN.* **Presentation of tier-3 results for the paper.** Still deferred. The summary object feeds either a short pkgdown article or the manuscript directly.

7.  *OPEN.* **The `B0` residual** (the *bias doubled* item in `TODO.md` group B). Removing the collection-covariate slopes removes 70% of the bias, but a remainder survives at 3.3 SE from zero at R = 200 (§15.6). So fixing item 4 will recover most of `B0` but not all of it, and the item cannot be closed as purely downstream. Cause not identified.

5.  *OPEN, and the case keeps strengthening.* **SBC** (§7). Three separate instances now of a fixed truth conflicting with an informative prior: `sigma_b` reading 1.000 because it is prior-dominated, `p` collapsing where the true value sits far into a `Beta(5, 1)` tail, and `theta0` overcovering. All are artefacts of choosing truth independently of the prior, and SBC cannot produce them by construction.

6.  *OPEN.* **A publication-grade high-R run**, if the paper claims calibration. R = 100 detects a large deviation but cannot assert the absence of a small one (§9). This should be one planned run across the whole grid, not per-arm confirmations.

## 11. Build order

1.  `helper-fixtures.R` + tier 1. Ships, pays immediately, blocks nothing.

2.  ~~Resolve open item 1 (spatial `n` floor).~~ **Done** -- floor is 31 sites; every spatial cell now pins `n_supportpoints` (§10.1).

3.  `helper-simstudy.R` with the two seams (§7).

4.  Tier 2.

5.  Tier 3 + `run_study.R`.

6.  **Full R = 100 run; inspect the table.** **Done 27-28 July 2026** -- all 10 cells at R = 100, 1000 fits, 0 failures, \~497 min on 5 cores (`base` alone first, 100 fits/22.9 min; then the other 9 cells, 900 fits/474 min). Results in `dev/simstudy/results/`, tabulated in §12. Found three undercovering blocks, each since traced to code; see §12 and TODO.md group A.

    **A re-run is now owed, and it checks two things at once.** Alex's 28-29 July fixes (`sigma_h`, the collection-covariate prior mean, the Stage 2 prior wiring) all change fitted values, so §12 is stale on its own terms. Separately, the replicates in that run were not independent -- see the note under §8 -- so its error bars were understated. The re-run therefore measures whether the fixes worked *and* produces the first numbers with honest SEs.

    ```         
    Rscript dev/simstudy/run_study.R --cores=5 --caffeinate
    ```

    Two cautions carried over from the first attempt. **Do not edit `run_study.R` while it is running** -- Rscript parses incrementally, and an earlier ten-hour run died on a torn read. Use `--resume` if it is interrupted. And an indicative post-fix probe at R = 8 showed `beta_theta` still undercovering at 0.775 with its bias halved, while `B0`'s bias *grew* fivefold to -0.529; R = 8 settles nothing, but it means the re-run should be read for regressions from Alex's 421-line `jsdmfun.R` rewrite as well as for the intended improvements.

------------------------------------------------------------------------

## 12. Results of the full run (29 July 2026, post-fix)

1000 replicates, 10 scenarios, 0 failures, **285 min** on 5 cores. 155,578 individual interval checks. Commit `8af22cd`. Coverage SE at R = 100 is 2.2 points, so treat anything in 0.93-0.97 as indistinguishable from nominal.

**This is a paired comparison with the 28 July pre-fix run.** `draw_truth()` seeds on (scenario, replicate), so the simulated data and true values are *bit-identical* between the two runs -- verified, `max|truth difference| = 0`. Every difference below is attributable to the code, not to sampling variation between runs. That is a much stronger design than two independent runs would give, and it is worth preserving: **do not change `simstudy_seed()`**, or future runs lose comparability with these.

**Read 12.1-12.3 alongside §14.7 and §15, which have since superseded parts of them.** The `beta_theta` finding in 12.1 has been narrowed to the collection-covariate *slopes* (§15.5), and the `B0` finding in 12.2 has been traced: 70% of the bias is downstream of those slopes, with a real residual confirmed at R = 200 (§15.6). The table itself stands; the attributions under it have moved on.

| Block | base | spat | trait | prim3 | lowinfo | d-under | d-over | S=20 | occ | bin |
|----|----|----|----|----|----|----|----|----|----|----|
| `B0` | 0.948 | 0.943 | 0.942 | 0.954 | **0.862** | 0.958 | 0.950 | 0.951 | 0.953 | 0.962 |
| `B` | 0.947 | **0.883** | 0.916 | 0.942 | **0.826** | 0.936 | 0.939 | 0.933 | 0.947 | 0.952 |
| `G` | 0.948 | -- | 0.945 | *0.975* | **0.840** | 0.953 | 0.955 | 0.925 | *0.970* | 0.938 |
| `beta_theta` | **0.763** | **0.762** | **0.771** | **0.718** | *0.983* | **0.750** | **0.752** | **0.727** | **0.709** | -- |
| `theta0` | *0.983* | *0.983* | *0.978* | *0.985* | **0.602** | *0.980* | *0.981* | *0.982* | *0.985* | -- |
| `p` | 0.902 | 0.911 | 0.900 | 0.925 | **0.109** | 0.913 | 0.919 | 0.908 | -- | -- |
| `q` | 0.945 | 0.947 | 0.949 | 0.937 | 0.911 | 0.946 | 0.943 | 0.950 | -- | -- |
| `resid_cor` | **0.763** | **0.758** | **0.758** | **0.768** | **0.764** | *0.980* | **0.755** | **0.752** | **0.761** | **0.752** |
| `sigma_b` | *1.000* | *1.000* | *1.000* | *1.000* | *1.000* | *1.000* | *1.000* | *1.000* | *1.000* | *1.000* |

**Bold = undercovers. Italic = overcovers. `--` = not estimated in that cell.**

Every figure above was re-derived from `dev/simstudy/results/simstudy-20260729-143756.csv` on 1 August 2026 and matches.

### 12.1 What the fixes changed

Averaged over all cells, pre-fix -> post-fix:

| Block | Coverage | Bias | Verdict |
|----|----|----|----|
| `beta_theta` | 0.719 -> **0.766** | +0.112 -> **+0.038** | improved, not fixed |
| `theta0` | 0.901 -> 0.944 | -0.020 -> -0.004 | moved, see 12.3 |
| `resid_cor` | 0.775 -> 0.777 | +0.006 -> +0.002 | **unchanged** |
| `p` | 0.821 -> 0.827 | +0.087 -> +0.086 | unchanged, as expected |
| `B0` | 0.947 -> 0.943 | -0.135 -> **-0.228** | **bias doubled** |
| `B`, `G`, `q` | \~0.92-0.94 | \~0 | stable |

1.  **`beta_theta` improved but is not fixed.** Every cell gained 0.03-0.05, and two-thirds of the bias is gone. Alex's correction of the collection-covariate prior mean from 1 to 0 was a real cause -- but not the only one, because 0.766 against nominal 0.95 is still far out. Whatever remains is flat across model type, primer count, species count and factor misspecification, exactly as the pre-fix version was. **Since narrowed:** §15.5 shows the intercept row is perfectly calibrated at 0.968 once the slopes are removed, so the defect is in slope estimation specifically, not in the block. This is the slope-overconfidence item in `TODO.md` group B.

2.  **`resid_cor` is untouched, and the paired design makes that conclusive.** Only **104 of 49,978** coverage decisions flipped. Identical data, identical truth, coverage unmoved at 0.777. This is `reparamFactorModel()` (the `reparamFactorModel()` item in `TODO.md` group B), which is unfixed and currently disputed. The `d_underfit` exception persists exactly as before: *over*covering at 0.980, because under-fitting the ordination widens intervals enough to mask the bias.

3.  **`p` is unchanged, and that is the correct outcome.** `low_information` sits at 0.109. The informative `Beta(5, 1)` is load-bearing for identifiability -- `p` and `q` enter `sample_pq_cpp()` symmetrically, so the prior is what selects the mode. This cell measures the *cost of that constraint* when true `p` is 0.1-0.3. It is not a defect to be fixed by flattening the prior.

### 12.2 New finding: `B0` bias roughly doubled

Nine of ten cells moved more negative: base -0.113 -> -0.208, `occupancy` -0.024 -> -0.151, `primers_3` -0.031 -> -0.151, `low_information` -0.931 -> -1.056. Only `binary` moved the other way (-0.002 -> +0.012).

**Coverage does not reveal this** -- it holds at 0.943 because the intervals are wide enough to absorb the shift. It is visible only in the bias column, which is the argument for tracking both.

**Since diagnosed, and this entry's guess was wrong.** It attributed the shift to the 421-line `jsdmfun.R` rewrite in the same pull. §15.4 and §15.6 show otherwise: setting `ncov_theta = 0` removes 70% of the bias, so most of it is downstream of the collection-covariate slopes, with a genuine residual surviving at 3.3 SE. This is the `B0`-bias item in `TODO.md` group B, and it stays open on that residual.

### 12.3 New finding: `theta0` now overcovers

The all-cell average of 0.944 is misleading. Per cell it is **0.978-0.985 in nine cells** (was 0.938-0.959) and **0.602 in `low_information`** (was 0.477). So it moved from mildly under to distinctly over, except where information is thin.

Overcoverage is the safe direction. This entry originally guessed the widened `diag(2)` prior variance had overshot; §14.7 disproved that directly, and the M ladder (§13.7) instead points at Stage 1 under-identification, which relaxes as M rises.

### 12.4 What holds

`q`, `B0`, `B` and `G` are at or near nominal in every cell except `low_information`. The quantities most likely to be reported in an ecology paper remain trustworthy -- with the `B0` bias caveat in 12.2, which affects the point estimate rather than the interval.

`low_information` is still compromised across the board (`theta0` 0.602, `p` 0.109, `B` 0.826, `G` 0.840, `B0` bias -1.06). Users with small or low-detection datasets remain the most exposed.

### 12.5 Caveats that travel with this table

`l_s` and `sigma_h` are absent because neither is recoverable (the `sample_ls()` item in `TODO.md` group B), so **no cell says anything about spatial range**. `sigma_b` reads 1.000 everywhere because it is prior-dominated by construction (§5.3), not because it is well estimated. Differences below 2.2 points are noise.

### 12.6 The superseded pre-fix table (28 July 2026)

Kept because the *delta* is the evidence that the fixes worked, not the level. Run in two parts -- `base` on 27 July (100 fits, 22.9 min), the other 9 cells on 28 July (900 fits, 474 min wall, inflated by the laptop sleeping mid-run). Its replicates were not independent (see §8), so **never quote an SE against these figures**.

| block | base | spat.isol | trait.isol | primers3 | low.info | d.under | d.over | sp20 | occ | binary |
|----|----|----|----|----|----|----|----|----|----|----|
| `B` | 0.943 | 0.876 | 0.913 | 0.931 | **0.854** | 0.937 | 0.947 | 0.935 | 0.943 | 0.948 |
| `B0` | 0.945 | 0.953 | 0.948 | 0.955 | 0.892 | 0.954 | 0.955 | 0.956 | 0.956 | 0.942 |
| `beta_theta` | **0.717** | **0.730** | **0.721** | **0.679** | **0.860** | **0.714** | **0.709** | **0.693** | **0.676** | -- |
| `G` | 0.945 | -- | 0.953 | 0.963 | 0.880 | 0.943 | 0.958 | 0.925 | 0.963 | 0.958 |
| `p` | **0.898** | **0.900** | **0.904** | 0.917 | **0.103** | 0.906 | 0.909 | 0.903 | -- | -- |
| `q` | 0.948 | 0.942 | 0.944 | 0.935 | 0.904 | 0.946 | 0.945 | 0.947 | -- | -- |
| `resid_cor` | **0.763** | **0.758** | **0.757** | **0.768** | **0.764** | *0.980* | **0.739** | **0.752** | **0.761** | **0.746** |
| `sigma_b` | *1.000* | *1.000* | *1.000* | *1.000* | *1.000* | *1.000* | *1.000* | *1.000* | *1.000* | *1.000* |
| `theta0` | 0.941 | 0.947 | 0.941 | 0.954 | **0.477** | 0.959 | 0.938 | 0.949 | 0.959 | -- |

------------------------------------------------------------------------

## 13. The M ladder: are B4-B6 defects, or Stage 1 under-identification?

**Status: planned, not yet run.** Commissioned by Doug, 29 July 2026, as the next step on the `beta_theta`, `theta0` and `B0` findings (now the slope-overconfidence, `theta0`-overcoverage and `B0`-bias items in `TODO.md` group B), each of which is annotated CLAUDE TO RUN SIMULATION STUDY WITH M \> 10 AND SEE IF IT FIXES IT.

### 13.1 The hypothesis, and why one lever could explain three findings

All three open findings sit downstream of the *field-collection* stage, and the grid has only ever been run at `M = 2` samples per site:

- **B4**, `beta_theta` undercovering at 0.766, is a Stage 1 covariate slope.
- **B5**, `theta0` overcovering at 0.978-0.985, is the Stage 1 species-specific baseline.
- **B6**, `B0` bias doubling to -0.228, is the occupancy intercept -- which is inferred *through* the latent collection state `w`, so it inherits whatever Stage 1 gets wrong.

At `M = 2` there are two Bernoulli collection draws per site-species, observed only indirectly through PCR detections, from which the sampler must identify a species baseline *and* covariate slopes *and* the latent `w` that `z` then depends on. That is very little information, and it is a plausible single cause for all three: weak identification inflates the prior's influence, which produces exactly the observed signature -- slopes too confident, baselines too diffuse, and an occupancy intercept biased toward under-inferred occupancy (negative, as measured).

**If true, this is not a bug.** It is a design requirement: the model needs `M` above some threshold before Stage 1 is identified. That is a documentation and paper item rather than a code fix, and an important one, because `M = 2` to `3` is entirely realistic in eDNA practice.

### 13.2 Design

A **ladder**, not a single M \> 10 point. One high-M run answers "fixed or not"; a ladder answers "fixed, and from where" -- and distinguishes a deviation that decays smoothly toward nominal (information-limited) from one that plateaus at a non-nominal level (a real defect that more data cannot reach).

| Arm   | M   | K   | Rows  | Purpose                                          |
|-------|-----|-----|-------|--------------------------------------------------|
| `M2`  | 2   | 3   | 1200  | Reproduces the current base cell; internal check |
| `M5`  | 5   | 3   | 3000  | Is the improvement already underway?             |
| `M10` | 10  | 3   | 6000  | Doug's threshold                                 |
| `M20` | 20  | 3   | 12000 | Comfortably past it                              |
| `K30` | 2   | 30  | 12000 | **Control** -- matched rows, wrong stage         |

**`K30` is what makes the result interpretable.** Raising `M` raises the total observation count, so if everything improves, the naive reading "more data helps" is unfalsifiable. `K30` holds the row count identical to `M20` but spends it on PCR replicates -- Stage 2 -- instead of field samples. If `M20` fixes the three findings and `K30` does not, the mechanism is specifically Stage 1 identification. If both fix them, it is merely sample size, and the conclusion is much weaker.

**The arms share a truth, so the ladder is paired.** `simstudy_seed_for()` honours a `seed_label` field, and all five arms set `seed_label = "mladder"`, so replicate *r* draws the same truth in every arm. The differences between arms then stop carrying the variance of independent truths -- which is precisely what a ladder reads.

**How far the pairing extends, measured rather than assumed:**

| Quantity | Paired across M? |   |
|----|----|----|
| `B0`, `B`, `G`, `A`, `C`, `L`, `sigma_b`, `tau` | **yes** | the whole JSDM layer |
| `z_true` | **yes** | latent occupancy |
| `p_true`, `q_true`, `theta0`, `theta_baseline` | **yes** | Stage 2 rates and Stage 1 baseline |
| `beta_theta` intercept row | **yes** | it is `logit(theta_baseline)` |
| `beta_theta` slope rows | **no** | see below |
| `w_true` | n/a | its dimension is `N x S`, so it cannot match |

`beta_theta`'s slopes are drawn at `R/simulateData.R:128` with `sample(c(-1,1,0), ...)` *after* `simulateOccJSDMData()` has built the `N`-sized `X_theta`, so by then the stream has diverged and there is no `list_params` hook to inject them. The residual is bounded -- slope signs from a three-point set, with the intercept row paired -- and the estimand is unchanged, so the ladder stays unbiased for `beta_theta`, just noisier than for the blocks above. **B6 (`B0`) and B5 (`theta0`) get the full paired benefit; B4 (`beta_theta`) gets a partial one.** The `M20` vs `K30` control comparison is likewise unpaired for `beta_theta`, since those arms differ in `N`.

`test-recovery.R` asserts all of this. It fails if anyone adds a draw to `draw_truth()` sized by `N`, `M` or `K`, which would silently unpair the ladder.

Base cell only: all three findings are flat across the other nine cells, so one cell isolates the mechanism at a tenth of the cost.

**`R = 50`, but only for the first stage, and the reason is asymmetric.**

Detecting that the deviation *persists* is easy: testing 0.766 against nominal uses the null SE, `sqrt(0.95 * 0.05 / 50)` = 3.1%, so the 18-point gap is about 6 SE.

Confirming it has *recovered* is much harder, and R = 50 cannot do it. That claim needs an interval around the observed value, where the SE is `sqrt(p(1-p)/R)` at the observed p:

| R   | SE at p \~ 0.94 | 95% CI around an observed 0.94 |
|-----|-----------------|--------------------------------|
| 50  | 3.4%            | [0.87, 1.00]                   |
| 100 | 2.4%            | [0.89, 0.99]                   |
| 200 | 1.7%            | **[0.91, 0.97]**               |

At R = 50 an observed 0.94 is indistinguishable from a true 0.88, which is still meaningful undercoverage -- so "fixed" could not be told from "much better but still broken". This is §9's point restated: asserting the absence of a small deviation is a different and more expensive claim than detecting a large one.

**So: two stages.** R = 50 across the ladder answers the question Doug actually asked -- *is `M` the lever?* -- because the answer is carried by the **shape across arms**, and a monotone rise over four arms is evidence no single point can supply. Then re-run only the decisive arm at **R = 200** if it lands high, since that is the one number that would appear in the paper. Roughly 1.9 h now and 2.4 h later, the second spent only if the first justifies it.

### 13.3 Cost (measured 29 July, not projected)

Per fit at the study's `nburn = 1000 + niter = 1000`, extrapolated from a 200-iteration timing on this machine:

| M   | s/fit | R = 50 on 5 cores |
|-----|-------|-------------------|
| 2   | 40    | \~8 min           |
| 5   | 64    | \~13 min          |
| 10  | 100   | \~20 min          |
| 20  | 179   | \~36 min          |

Plus `K30` at roughly the `M20` cost. **Total \~1.9 h.** Note the cost scales *sub*-linearly in rows -- M = 20 is 4.5x M = 2, not 10x -- because much of the per-iteration work is over species and sites, not samples.

### 13.4 What each outcome means

Decided in advance, so this is a test rather than a fishing trip.

| Outcome | Reading | Action |
|----|----|----|
| All three approach nominal as M rises; `K30` does not | Stage 1 under-identification. Not defects. | **Re-run the winning arm at R = 200 before closing anything** -- R = 50 cannot separate 0.94 from 0.88. Then close B4-B6 as "not a bug", document the `M` requirement, and state it in the MEE paper: it is a real limitation for `M = 2` eDNA designs. |
| B4/B5 recover, B6's bias persists | Two causes. Stage 1 explains the coverage findings; `B0` is a separate regression. | Close B4/B5; keep B6 open and go after the `jsdmfun.R` rewrite. |
| Nothing improves, even at M = 20 | Genuine code defects; `M` is not the lever. | Keep all three open. Next suspect is the widened `diag(2)` prior variance, which is the common edit behind B4 and B5. |
| `M20` and `K30` improve equally | Sample size, not Stage 1. | Weak conclusion; the finding is that all three are information-sensitive. Re-think the design before spending more. |

### 13.5 Implementation

No machinery change needed: `M` and `K` are already scenario fields, and `mk()` overrides `SIMSTUDY_BASE` by name.

1.  ~~Add the five arms to `simstudy_scenarios()`.~~ **Done**: `M2`, `M5`, `M10`, `M20`, `K30`, each with `seed_label = "mladder"`. They are inert for a full run, because the runner takes the whole list only when `--scenarios` is empty -- so select them explicitly. `simstudy_seed_for()` is the seed override; `simstudy_seed()` is unchanged and still keys on the label for every other cell.
2.  Sanity-run at `R = 4` across all five arms (\~10 min) to confirm the arms fit and the row counts are as expected, before committing the real run.
3.  `Rscript dev/simstudy/run_study.R --scenarios=M2,M5,M10,M20,K30 --R=50 --cores=5 --caffeinate`
4.  Report coverage **and bias** for `beta_theta`, `theta0` and `B0`, plus the remaining blocks to confirm nothing degrades.

**Check first:** the `M2` arm should agree with the existing base cell (`beta_theta` 0.763, `theta0` 0.983, `B0` bias -0.208) within its SE. Its seeds differ -- `simstudy_seed()` keys on the label -- so it will not reproduce them exactly, but a disagreement beyond about 2 SE would mean the ladder is measuring something other than what it claims, and should be resolved before reading the rest.

### 13.6 Caveats

`M` and the number of *sites* are not interchangeable. This ladder holds `n = 100` fixed and varies samples per site, which is the quantity Stage 1 is identified by. It says nothing about whether more sites would help, and should not be quoted as if it did.

Raising `M` also raises the number of latent `w` states being sampled, so mixing may differ across arms. If a high-M arm shows worse mixing rather than better coverage, check `returnConvergenceDiagnostics()` before concluding anything about identification.

------------------------------------------------------------------------

## 13.7 Results (29 July 2026)

**Run:** 250 fits, R = 50, 5 arms, 189 min, 0 failures. Validity check passed: the `M2` arm agrees with the existing base cell on all three targets (largest gap 0.016, against an SE of \~0.09 for the difference) -- the ladder is measuring what it claims.

|                       | M2     | M5     | M10    | M20       | K30 (control) |
|-----------------------|--------|--------|--------|-----------|---------------|
| `beta_theta` coverage | 0.747  | 0.655  | 0.603  | **0.579** | 0.706         |
| `theta0` coverage     | 0.986  | 0.966  | 0.944  | 0.952     | 0.996         |
| `B0` bias             | -0.160 | -0.025 | +0.002 | -0.028    | -0.091        |

**None of the four outcomes in §13.4 fits.** The actual result is a mixture across the three items:

**`theta0` (B5): clean confirmation.** Overcoverage falls from 0.986 toward nominal as `M` rises (0.944 at `M10`), while `K30` -- matched row count, wrong stage -- makes it *worse* (0.996). This is exactly outcome 1's signature: `M` is the lever, matched data volume elsewhere is not. **Close as Stage 1 under-identification, not a defect** -- pending the R = 200 confirmation §13.4 requires before closing anything.

**`B0` (B6): partial confirmation.** Bias collapses from -0.160 to near zero at `M10`. But `K30` also improves it (-0.091), just less than `M10` or `M20` -- so the recovery is not cleanly Stage-1-specific; more data of either kind helps. Coverage stays 0.94-0.96 throughout in every arm, consistent with the original finding that coverage does not reveal this bias.

**`beta_theta` (B4) goes the wrong way, and this is the important result.** Coverage falls *monotonically* with `M` -- 0.747 -\> 0.579, a 17-point drop across four points, far beyond the R = 50 noise floor -- while bias stays small and flat (+0.02 to +0.05) throughout. Shrinking intervals around a bias that is not itself shrinking is the signature of a real defect being *exposed* by more information, not resolved by it. `K30` (0.706) beats every M arm above `M2`.

**This rules out under-identification as the explanation for B4.** `beta_theta` needs a different next step than B5/B6: not "wait for more data" but "find what makes the interval overconfident". The prime suspect is the same `diag(2)` prior-variance widening implicated in B5 (Fixed bugs 25) -- but if that edit was meant to fix things and B4 gets *worse* with more information, it may need to be revisited rather than extended.

**Unplanned finding: `q` degrades hard with `K`.** 0.945 (`M2`) -\> 0.614 (`K30`), worse than any M arm, including `M20` (0.742) which shares `K30`'s row count. Not investigated here -- outside what this ladder was built to answer -- but worth its own entry; see `TODO.md`.

**What this does not settle.** R = 50 detects that `beta_theta` is getting worse with high confidence (a 17-point monotone trend dwarfs the 3.1% SE many times over), but per §13.2's own caveat, R = 50 cannot confirm that `theta0` has *reached* nominal -- an observed 0.944 is not distinguishable from a true 0.90. Do not close B5 on this run alone.

------------------------------------------------------------------------

## 13.8 Follow-up: tighter `beta_theta` prior at high M

**Status: running.** Commissioned by Doug, 29 July 2026, directly off 13.7's finding that `beta_theta` coverage worsens with M -- the signature of an overconfident interval rather than insufficient data, which the write-up pointed at `B_betatheta`'s slope variance.

**That variance was hard-coded**, with no `listPriors` hook, unlike `p`/`q`/`theta0`. Added one: `listPriors$b_betatheta_slope_var` (`R/runOccJSDM.R`), defaulting to 2 -- the existing value -- so nothing changes unless a caller sets it. `ALEX TO REVIEW` before treating a non-default value as a real fix rather than a diagnostic; see the `b_betatheta` variance-decision item in `TODO.md` group B.

**Design:** two arms, `M10_tightprior` and `M20_tightprior`, repeating the worst two points on the ladder with `b_betatheta_slope_var = 0.5` (SD 0.71, against the default's 1.41). Same `seed_label = "mladder"` as the original ladder, so these pair not only against each other but against the *already-collected* `M10`/`M20` results -- no need to re-run the default-prior arms.

**Verified before the real run**: refit one dataset under both priors with a short chain, holding data and seed fixed; the slope row's posterior spread shrank under the tighter prior (0.187 -\> 0.179 SD at this toy configuration), confirming the override reaches the sampler. Tier 1/2 unaffected -- the new default preserves current behaviour exactly.

**Reading it:** if coverage moves *toward* nominal at the tighter variance, the diagnosis holds and `B_betatheta`'s width is implicated. If it does not move, or moves the wrong way, the cause is elsewhere and the prior-variance hypothesis from 13.7 is wrong -- worth stating plainly either way, per the lesson already recorded about not forcing a real result into a box it does not fit.

R = 50, for the same reason as 13.2: adequate to detect a move, not to confirm arrival at nominal.

------------------------------------------------------------------------

## 13.9 Results: tighter prior did not move `beta_theta` (30 July 2026)

**Run:** 100 fits, R = 50, 97 min, 0 failures.

|                       | M10 (var=2) | M10 (var=0.5) | M20 (var=2) | M20 (var=0.5) |
|-----------------------|-------------|---------------|-------------|---------------|
| `beta_theta` coverage | 0.603       | 0.600         | 0.579       | 0.578         |
| `beta_theta` bias     | +0.027      | +0.026        | +0.017      | +0.017        |

**The hypothesis from 13.7 is disproved.** A 4x reduction in prior variance (SD 1.41 -\> 0.71) moved coverage by 0.001-0.003 -- noise against the 3.1% SE. Every other block (`theta0`, `B0`, `B`, `G`, `q`, `resid_cor`) was likewise unmoved, confirming the change was isolated and had no side effects worth reporting. `B_betatheta`'s width is not what makes the interval overconfident.

**A second candidate was checked and is also dead.** Before writing this up, checked whether added M-samples might be pseudo-replicated -- i.e. whether `X_theta` repeats within a site, which would make more M look like independent information without actually being any. It does not: `X_theta <- cbind(1, matrix(rnorm(N * ncov_theta), N, ncov_theta))` (`R/simulateData.R:126`) draws the covariate independently per *sample*, not per site. So the narrowing interval is not a design-matrix artefact either.

**Where this leaves item 4.** Two plausible mechanisms are now ruled out: the prior is not too tight, and the covariate is not pseudo-replicated. The overconfidence has to come from somewhere in the likelihood or the sampler's variance computation for `beta_theta` -- plausibly how the latent collection state `w` (or `z`) is aggregated across a site's M samples, or a genuine numerical issue in the Polya-Gamma update. That is a C++/sampler-level investigation, not another prior-tuning experiment, and is a job for whoever wrote `sample_beta_cpp_TS`/`sample_betatheta_cpp_parallel`, not a further simulation-study pass. Recorded in `TODO.md` rather than pursued further here.

------------------------------------------------------------------------

## 14. `theta0` overcoverage: a better experiment than the one planned

**Status: planned, not yet run.** Supersedes the previous plan for the `theta0` finding, which was "re-run the M10 or M20 arm at R = 200". That plan answers the wrong question, for the reason below.

### 14.1 The question was framed wrong

13.7 read `theta0` as Stage 1 under-identification, because coverage falls toward nominal as M rises (0.986 at M2, 0.944 at M10) while the matched `K30` control makes it worse. That reading has a hole in it: **the pre-fix grid had `theta0` at 0.938-0.959, which is nominal, at the same M = 2 where it now sits at 0.978-0.985.** If M = 2 were simply too little information for `theta0`, it would have overcovered before Alex's fixes too. It did not.

So the question is not "is `theta0` fine at high M". We already have four arms saying it is, and high M is not the configuration anyone runs. The question is **what changed at M = 2**, which is the production setting and the one the shipped defaults describe.

### 14.2 What the existing evidence points at

`theta0`'s own prior is untouched: `a_theta0`/`b_theta0` are still `Beta(1, 20)`, last changed in `1c25529` on 23 July, before the pre-fix run. So the cause is not `theta0`'s prior directly.

The likely route is coupling. `b_betatheta`'s variance was widened to `diag(2)` in Fixed bugs 25. That changes `beta_theta`, which drives the collection probability, which drives the latent `w`; and `sample_theta0(z, w, ...)` conditions on `w`. A more uncertain `w` gives a wider `theta0` posterior, which is overcoverage.

**The tighter-prior run already supports this, weakly.** Both arms moved `theta0` toward nominal:

| arm | var = 2 | var = 0.5 | change |
|-----|---------|-----------|--------|
| M10 | 0.944   | 0.940     | -0.004 |
| M20 | 0.952   | 0.942     | -0.010 |

Small, as expected: those are the arms where data dominates the prior. **M2 is the one arm where the prior should dominate, and it is the one arm not yet run at the tighter setting.**

The `K30` result (0.996, worse than M2's 0.986) also fits: more PCR replicates sharpen `w` without adding Stage 1 samples, and sharpening the wrong stage does not help a term whose uncertainty comes from the collection process.

### 14.3 The experiment

Two new arms at M = 2, varying only `b_betatheta_slope_var`, using the `listPriors` hook added for 13.8:

| Arm         | M   | `b_betatheta_slope_var` | Purpose                            |
|-------------|-----|-------------------------|------------------------------------|
| `M2`        | 2   | 2 (default)             | already collected, 13.7            |
| `M2_tight`  | 2   | 0.5                     | matches the 13.8 setting           |
| `M2_vtight` | 2   | 0.1                     | far enough to force a visible move |

Same `seed_label = "mladder"`, so these pair against the existing `M2` arm. R = 50.

**Cost: about 16 minutes.** M2 is the cheapest arm at roughly 40 s per fit, so 100 fits on 5 cores is trivial next to the 97 and 189 minute runs already done. This is a much better use of an hour than an R = 200 confirmation of something already measured four ways.

### 14.4 What each outcome means

| Outcome | Reading | Action |
|----|----|----|
| `theta0` walks 0.986 toward 0.95 as the variance tightens | Coupling confirmed. `theta0` overcoverage is a downstream symptom of the `b_betatheta` widening, not an independent defect. | Close B5 as part of B4. One cause, two symptoms; whatever fixes B4's overconfidence should be checked against `theta0` at the same time. |
| `theta0` does not move | Coupling is not the route. The cause is elsewhere in the `jsdmfun.R` rewrite, or in `sample_theta0()` itself. | Next test is `theta0`'s own prior, which already has `listPriors$a_theta0`/`b_theta0` hooks, at M = 2. |
| `theta0` moves the wrong way | The widening was compensating for something, and reverting it alone would make `theta0` worse. | Stop and hand to Alex with both results; this would mean B4 and B5 pull in opposite directions. |

**Bonus, at no extra cost:** the same run says whether `beta_theta` responds to the prior at M = 2. 13.9 found no response at M10/M20, but those are the arms where data dominates. A response at M = 2 would qualify 13.9's "the prior is not the cause" conclusion, which was drawn only from high-M arms.

### 14.5 Priority, stated honestly

**B5 is the least urgent of the open findings and this plan should not grow.** Overcoverage is the safe direction: intervals wider than they need to be cost power, not correctness. Nobody publishes a wrong number because of it. Compare B4, where undercoverage at 0.58 means genuinely overconfident intervals, and B6, where `B0`'s bias doubled invisibly to coverage.

The case for spending 16 minutes here is that it is cheap and likely resolves B5 as a **side effect of diagnosing B4**, not that B5 is important on its own. If the first outcome above holds, B5 stops being a separate item entirely.

### 14.6 What happened to the R = 200 plan

Dropped, not deferred. Its purpose was to confirm `theta0` has genuinely reached nominal at high M, since R = 50 cannot distinguish 0.944 from 0.90. But that confirmation is only worth buying if the claim matters, and "`theta0` is well calibrated when you take 10 or 20 samples per site" is not a claim the paper needs or that users can act on at M = 2 to 3. If a publication-grade calibration claim is wanted later, it belongs in a single planned high-R run across the whole grid (see 9), not a one-arm confirmation of the least urgent finding.

------------------------------------------------------------------------

## 14.7 Results: hypothesis dead, but it diagnosed B6 instead (30 July 2026)

**Run:** 100 fits, R = 50, 20.4 min, 0 failures. All three arms share `seed_label = "mladder"`, and truth is bit-identical across them (verified), so every comparison below is paired.

|                       | var = 2 (default) | var = 0.5 | var = 0.1  |
|-----------------------|-------------------|-----------|------------|
| `theta0` coverage     | 0.986             | 0.982     | 0.980      |
| `beta_theta` coverage | 0.747             | 0.707     | 0.653      |
| `B0` bias             | -0.160            | -0.106    | **-0.044** |
| `B0` coverage         | 0.946             | 0.950     | 0.946      |

### `theta0`: the hypothesis is disproved

A **20-fold** reduction in the prior variance moved coverage by 0.006. The paired bias change reaches 2.7 SE at var = 0.1, so it is detectable, but it is far too small to matter: reaching nominal from 0.986 by this route would need an extrapolation the data does not support. **Coupling through `b_betatheta` is not what makes `theta0` overcover.**

That is outcome 2 of 14.4, so the next test is `theta0`'s own `Beta(1, 20)` prior, which already has `listPriors$a_theta0`/`b_theta0` hooks and needs no code change. Given 14.5's priority argument, that is worth doing only if it can ride along with another run.

### `B0`: an unplanned diagnosis, and the most useful thing here

`B0`'s bias responds strongly and monotonically to the same knob: -0.160 at var = 2, -0.106 at 0.5, -0.044 at 0.1. Paired, that is +0.054 (2.1 SE) and +0.116 (4.0 SE).

**The history lines up exactly.** `42198d9` widened `B_betatheta` from `diag(1)` to `diag(2)` (and corrected the mean from 1 to 0), and that is precisely when `B0`'s bias doubled, from -0.135 pre-fix to -0.228 post-fix. Turning the variance back down moves the bias back. B6 was filed as "possible regression, cause not yet identified"; **the widening is now the leading candidate, with a dose-response curve behind it.**

`B0` coverage is 0.946-0.950 across all three arms, unchanged. That is why the original grid never flagged this: the intervals are wide enough to absorb the shift, and only the bias column moves. Further support for tracking both.

### `beta_theta`: 13.9's conclusion survives, its reasoning does not

13.9 concluded the prior's width is not the cause, from M10/M20 where tightening moved coverage by 0.001-0.003. **At M = 2 it moves coverage a lot: 0.747 to 0.653.** The conclusion holds, because tightening makes coverage *worse*, so the prior is not a fix. But the reasoning was wrong: the prior is not inert, it was inert *in the arms tested*, which were the arms where data dominates it. Anyone re-reading 13.9 should read this section with it.

The detail worth keeping: bias falls (+0.0195 to +0.0063) while coverage falls faster. The intervals shrink more than the bias does. That is the overconfidence signature, now shown to be **controllable by the prior without being caused by it**.

### This is a trade-off, not a fix, and it is Alex's call

Tightening `b_betatheta`'s slope variance helps `B0`'s bias and hurts `beta_theta`'s coverage. Those are two open items pulling in opposite directions on one knob. Setting it needs someone who knows what that prior is meant to encode, which is why this goes to Alex rather than being applied. See `TODO.md`.

------------------------------------------------------------------------

## 15. Is `B0`'s bias caused by `beta_theta`? Alex's discriminator, sharpened

**Status: planned, not yet run.** Alex's suggestion on the `B0` item: "we could add a simulation study on the one-stage model only, that would reveal whether there is any issue in `B0` (since there is no `beta_theta`). If there is no issue, we could delete this point and be sure that the issue is only `beta_theta`."

The logic is right and it is the cleanest discriminator anyone has proposed for this. Two adjustments below, one because existing data already answers most of it, and one because "one-stage" does not actually remove `beta_theta`.

### 15.1 Existing data already points hard at his answer

`B0` bias by cell, post-fix run, sorted:

| cell | has `beta_theta`? | `B0` bias |
|---|---|---|
| `low_information` | yes | -1.056 |
| `base` | yes | -0.208 |
| `spatial_isolated` | yes | -0.204 |
| `traits_isolated` | yes | -0.175 |
| `d_overfit` | yes | -0.166 |
| `occupancy` | yes | -0.151 |
| `primers_3` | yes | -0.151 |
| `species_20` | yes | -0.141 |
| `d_underfit` | yes | -0.125 |
| **`binary`** | **no** | **+0.012** |

Perfect separation. Every cell that estimates `beta_theta` has a negative `B0` bias; the single cell that does not is the only one at zero, and the only positive one.

A second, independent line agrees. 14.7 varied `b_betatheta`'s slope prior variance at `M = 2` and `B0`'s bias moved monotonically with it: -0.160 at variance 2, -0.106 at 0.5, -0.044 at 0.1. So `B0`'s bias tracks `beta_theta`'s prior as well as `beta_theta`'s presence.

**Neither is conclusive on its own.** `binary` differs from the others in far more than `beta_theta`: no replicates at all, no latent `w`/`z` layer, no `p`/`q`/`theta0`, and much less data. It is a confounded comparison, which is exactly why a designed test is still worth running.

### 15.2 "One-stage only" is not the sharp test

The grid already has a one-stage cell. `occupancy` drops the PCR stage entirely and **still estimates `beta_theta`**, with `B0` bias -0.151, squarely among the two-stage cells. So removing a stage does not remove the collection-covariate coefficients, and running "one-stage only" would reproduce a result we already have.

What distinguishes `binary` is not that it has one stage. It is that it has **no collection covariates**.

### 15.3 The sharp test: hold everything, set `ncov_theta = 0`

One arm, identical to `base` except that the Stage 1 design matrix has no covariate columns. That keeps the two-stage machinery, the latent `w`/`z`, and `p`/`q`/`theta0` all in play, and removes only the `beta_theta` **slopes**; the intercept row, `logit(theta_baseline)`, necessarily remains.

- **`B0`'s bias vanishes** -> it is caused by the collection-covariate slopes. Combined with 14.7, that makes `beta_theta` the single upstream cause, and the `B0` item closes as downstream of it: fix the `beta_theta` slope defect and `B0` follows. This is the outcome Alex predicts, and the one the existing evidence favours.
- **`B0`'s bias persists** -> it is the two-stage latent structure rather than `beta_theta`, and `binary`'s clean result was down to one of its other differences. The `B0` item stays open and needs a different line of attack.
- **Bias shrinks but does not vanish** -> both contribute, and the split tells you how much of it fixing `beta_theta` would buy.

**Feasibility, verified.** The simulator handles `ncov_theta = 0`: `beta_theta_true` becomes `1 x S` and no `X_theta` column is produced. One harness fix is needed first, because `simstudy_fit()` passes `collCovariates = "X_theta"` unconditionally and the fit errors with `Covariate names provided not in data$info`. Make that conditional on `ncov_theta > 0`.

**Cost:** one arm at the `base` configuration, about 40 s per fit, so roughly 8 minutes at R = 50. Pair it against the existing `base` cell with `seed_label`, noting that pairing will be partial in the same way as the M ladder: `draw_truth()` consumes no RNG that depends on `ncov_theta`, so the drawn truths match, but the simulator's own stream diverges once it builds a differently-shaped `X_theta`.

**Read it against bias, not coverage.** `B0` coverage sits at 0.94-0.96 in every cell including `binary`, so it will not distinguish anything here. This is the finding that only the bias column reveals, which is the whole reason it went unnoticed for two days.

### 15.4 Results (31 July 2026)

**Run:** 50 fits, 7.5 min, 0 failures. `beta_theta` came back with 10 elements per replicate rather than 20, confirming the slopes were genuinely absent and only the intercept row remained.

**`B0` bias, arms not paired, so SEs are across replicates:**

| arm | `beta_theta` slopes | `B0` bias | SE |
|---|---|---|---|
| `base` | present | -0.208 | 0.031 |
| **`nocollcov`** | **absent** | **-0.057** | 0.039 |
| `binary` | no `beta_theta` at all | +0.012 | 0.012 |

Removing the slopes removes **73% of the bias**: a shift of +0.151 at 3.0 SE. But it does not eliminate it. The residual -0.057 sits between `base` and `binary`, and is 1.5 SE from zero, so it is suggestive rather than established.

**This is outcome 3 of 15.3, not outcome 1.** Both mechanisms contribute. Most of `B0`'s bias is downstream of `beta_theta`, which is the useful part of Alex's prediction and means fixing `beta_theta` should recover most of `B0`. But something in the two-stage machinery contributes independently of the collection-covariate slopes, so the `B0` item cannot simply be deleted on this evidence.

Worth stating plainly because the temptation is to round 73% up to "solved": if `beta_theta` were fixed tomorrow and `B0`'s bias fell to -0.057, that is still four times `binary`'s +0.012, in a headline JSDM quantity, and coverage would not reveal it.

### 15.5 The unplanned finding: the defect is in the slopes, not the block

`beta_theta` coverage, same two arms:

| arm | elements per replicate | coverage |
|---|---|---|
| `base` | 20 (intercept + slope, x 10 species) | 0.763 |
| `nocollcov` | 10 (intercept only) | **0.966** |

**With the slopes removed, `beta_theta`'s intercept row is perfectly calibrated.** So the overconfidence tracked in `TODO.md` group B is not a property of `beta_theta` as a parameter block: it is specific to the **slope** estimation.

This narrows that investigation considerably. Four candidate causes are now ruled out:

- Stage 1 under-identification (13.7: more data makes it worse, not better)
- The slope prior's width (14.7: tightening 20-fold at `M = 2` moves coverage the wrong way)
- Pseudo-replication in `X_theta` (13.9: drawn per sample, not per site)
- **The intercept path (here: nominal at 0.966 once the slopes are gone)**

What remains is whatever handles the covariate columns specifically in the Polya-Gamma update. That is also the step Alex's `microbenchmark()` profiling identified as the slowest in the sampler, so the two lines of enquiry converge on the same code.

### 15.6 Confirmation at R = 200 (31 July 2026)

**Run:** 200 fits, 39.9 min, 0 failures.

**Validity check passed.** Replicates 1-50 reuse the R = 50 seeds and reproduce that run **bit-for-bit**: 6,541 rows matched with truth *and* posterior means identical. So this is a genuine superset of 15.4 rather than a separate experiment, and the two can be read together. Worth noting that this check is only possible because `set.seed()` controls the sampler (*Fixed bugs* 28), the property judged not worth caring about at the user-facing level.

**`B0` bias:**

| arm | `beta_theta` slopes | bias | SE | from zero |
|---|---|---|---|---|
| `base` | present | -0.2078 | 0.0307 | 6.8 SE |
| **`nocollcov`** | **absent** | **-0.0633** | 0.0192 | **3.3 SE** |
| `binary` | none at all | +0.0122 | 0.0119 | 1.0 SE |

**The residual is real.** At R = 50 it was 1.5 SE and 15.4 called it suggestive; at R = 200 it is 3.3 SE. Removing the slopes removes 70% of the bias, so fixing the `beta_theta` slope defect will recover most of `B0`, but a genuine remainder survives. `binary` at 1.0 SE from zero is what "no bias" looks like for comparison.

**Consequence for the `B0` item: it cannot be closed.** Alex proposed deleting it if the one-stage test showed no issue. The test, in its sharpened form, shows a smaller issue rather than none. The item stays open with its cause now split: most of it is downstream of `beta_theta`, and a measurable part is not.

**The 15.5 narrowing holds at four times the replicates.** `beta_theta` intercept-only coverage is **0.968** with SE 0.013, against 0.763 with the slopes present. Nominal on any reading. The overconfidence is specific to the covariate columns, and that conclusion no longer rests on 50 fits.

**A cost note for future planning.** This took 39.9 min against a 30 min projection scaled linearly from the R = 50 run. The arithmetic no longer holds: RcppParallel now spawns TBB threads *inside* each PSOCK worker, so 5 workers oversubscribe the 4 performance cores. Load average ran at 15 with cpu/wall at 6.0x. Per-fit throughput therefore falls as worker count rises, and "fits/min x cores" over-predicts. Measure rather than extrapolate until the interaction between `--cores` and the internal threading is characterised.

## 16. `B0` with no `beta_theta` at all: Alex's actual discriminator

**Alex's correction to 15, given 2 August 2026:** *"Actually I wasn't too clear, my suggestion was to run the model with `model = "continuous"` since that part of the sampler would use `B0` only. Using the occupancy model, we still sample the intercept of `beta_theta` so indetermination between `B0` and `beta_theta` still affects the estimate."*

### 16.1 Why `nocollcov` did not settle it

15.3 chose `ncov_theta = 0` deliberately, to change one thing rather than many, and 15.6 read the surviving -0.0633 bias (3.3 SE) as a real residual after the slopes were removed. Alex's point is that the arm cannot support that reading, because setting `ncov_theta = 0` removes the slope columns but **necessarily keeps the intercept row**, `logit(theta_baseline)`.

In an occupancy model that intercept and `B0` sit on the same chain: `psi` governs how often a site is occupied, `theta` how often an occupied site yields a positive sample, and the product is close to what the data sees. At `M = 2` there is very little information separating them. So the residual is equally consistent with `B0`/`beta_theta`-intercept confounding as with a defect in `B0` itself, and 15.6 did not distinguish those.

### 16.2 The arm

`model = "continuous"`: `z ~ Normal(eta, tau)` observed directly. No detection stage, no latent `w`, no `p`/`q`/`theta0`, and **no `beta_theta` of any kind**. `B0` is estimated straight from the Gaussian likelihood. Verified at setup: `beta_theta_output` and `theta0_output` come back `NULL`, and the blocks that map are `B0`, `B`, `G`, `sigma_b`, `tau` and `resid_cor`.

**What it is and is not.** `binary` already removes `beta_theta`, and 15.1 discounted it for changing many things at once; `continuous` has exactly the same weakness and it is not claimed otherwise. Its value is that it is a *different likelihood and a different branch of the sampler* than `binary` -- Gaussian conjugate updates rather than the Polya-Gamma path -- so two independent readings agreeing is worth more than either alone. This is corroboration, not a controlled contrast. The controlled contrast Alex's point calls for does not exist in this model family, because there is no way to keep the two-stage structure while removing the `theta` intercept.

### 16.3 Run configuration, and why it is not the default

Run at R = 200 to match 15.6, with **`RCPP_PARALLEL_NUM_THREADS=1`**.

This is not a performance choice. `41abe69` introduced a data race on R's RNG inside `BBSL_Worker` (TODO.md group B, the `BBSL_Worker` item): `sampleB_SoR()` draws via `arma::randn()`, which `ARMA_RNG_ALT` routes to R's global RNG, and it is now called from every TBB worker thread. Measured: two fits under one `set.seed()` differ by 4.34 on `B0` and 0.63 on `p`. The race can perturb the posterior itself, not merely the draw order, which is the same order of magnitude as the -0.0633 this section is trying to measure.

At one thread the sampler is bit-for-bit reproducible again (verified: `max diff 0` through `simstudy_replicate()`), and replicate-level parallelism is unaffected because that is process-level. Throughput did not suffer -- 6 single-threaded workers ran marginally faster than 4 workers with TBB threads inside them, consistent with the oversubscription noted in 15.6.

**Consequence for comparability.** These numbers are not directly comparable to 15.4/15.6, which were produced at the default thread count and therefore under the race. Nothing in this section should be read as a paired contrast against `nocollcov`; the `base`/`nocollcov`/`binary` figures are quoted for orientation only, and re-establishing them race-free is work that belongs with the group B item, not here.

### 16.4 Results (2 August 2026)

**Run:** 200 fits, 20.2 min on 6 single-threaded workers, 0 failures.

**`B0` bias, continuous arm:**

| arm | `beta_theta` | `B0` bias | SE | from zero |
|---|---|---|---|---|
| `base` | slopes + intercept | -0.2078 | 0.0307 | 6.8 SE |
| `nocollcov` | intercept only | -0.0633 | 0.0192 | 3.3 SE |
| `binary` | none at all | +0.0122 | 0.0119 | 1.0 SE |
| **`continuous`** | **none at all** | **+0.0066** | **0.0056** | **1.2 SE** |

**`B0` is unbiased when `beta_theta` is absent entirely.** Two arms with no `beta_theta`, on two different likelihoods and two different branches of the sampler, both land within 1.2 SE of zero. That is the corroboration 16.2 was after: `binary`'s clean result was not an artefact of its many other differences, because `continuous` differs from `binary` in most of those same respects and agrees anyway.

**This supports Alex's reading of the `nocollcov` residual.** The -0.0633 at 3.3 SE was the only evidence for a `B0` defect surviving removal of the slopes, and 16.1 gives a mechanism for it that is not a defect: `B0` and `beta_theta`'s retained intercept are confounded at `M = 2`. With both no-`beta_theta` arms at zero, the confounding explanation is the more economical one, and the `B0` item reduces to "downstream of the `beta_theta` slope defect".

**It does not close the item on its own, for a reason that is procedural rather than statistical.** The `base`, `nocollcov` and `binary` figures above were all produced at the default thread count, i.e. under the `BBSL_Worker` race (16.3). `continuous` was not. Comparing a race-free measurement against three racy ones is not a clean contrast, and the honest way to close this is to re-run `base` and `nocollcov` with `RCPP_PARALLEL_NUM_THREADS=1` once the race is fixed, and check the -0.2078 and -0.0633 reproduce. Each arm is internally valid -- every bias is measured against its own known truths -- so the *direction* of the conclusion is not in doubt; the paired arithmetic across arms is what wants redoing.

### 16.5 An unlooked-for finding: `B0` undercovers in `continuous`

`B0` coverage in this arm is **0.879**, against a nominal 0.95 and a coverage SE of about 1.5% at R = 200. That is roughly 4.7 SE low, and it is the **lowest `B0` coverage of any cell measured**: the production grid runs 0.892 (`low_information`) to 0.956, with `binary` at 0.942.

The bias is zero, so this is not the estimate being wrong; it is the interval being too narrow. Note the direction of surprise -- `continuous` is the *only* model type that also estimates the residual variance `tau`, and `tau` carries its own bias here (+0.0370, coverage 0.921), so the extra variance component is being under-propagated into `B0`'s interval rather than widening it.

Not investigated further, and deliberately not filed as a group B item yet: this is one arm, at one configuration, on a model type the package's users are least likely to be running, and it was found while looking for something else. It wants confirming at a second configuration before it is called a defect. Recorded here so it is not lost.

The other blocks in this arm, for completeness: `B` covers 0.943 (bias 0.0003), `G` 0.954 (-0.0008), `resid_cor` 0.703 (0.0003, consistent with the known `reparamFactorModel()` defect measuring 0.74-0.77 elsewhere), `sigma_b` 1.000 (prior-dominated by construction).
