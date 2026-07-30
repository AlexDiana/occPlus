# Simulation-based validation of `occJSDM`

Plan for the "extensive testing on simulated datasets" item under *MEE
paper / Doug to dos* in `TODO.md`. Drafted 27 July 2026.

[Guide to the test
suite](https://claude.ai/code/artifact/ad3d46eb-1fd4-49b5-b795-6b71474ef1d5 "Guide to the test suite")

**Status as of 29 July 2026: complete, and the study has now been run
twice.** The pre-fix run (27–28 July) found four defects; three were
traced to specific lines of code and one turned out not to be a defect
at all. The post-fix re-run (29 July) is the current table in §12, and
it added two further findings. Sections marked *OPEN* still need a
decision or a measurement.

| Stage | State |
|----|----|
| 1\. Fixtures + tier 1 | **done** (`6d9526d`) -- 89 assertions, 7 s; 112 as of 29 July |
| 2\. Spatial `n` floor | **done** (`41fb736`) -- floor is 31 sites |
| 3\. `helper-simstudy.R` | **done** (`6123036`) -- validated at R = 8 |
| 4\. Tier 2 canary | **done** (`f63eeeb`) -- \~30 s, `skip_on_cran()` |
| 5\. Tier 3 + runner | **done** (`f63eeeb`) -- runner smoke-tested end to end |
| 6\. Full R = 100 run | **done, run twice** -- pre-fix 27-28 July (1000 fits across two runs), post-fix 29 July (1000 fits, 285 min, one run). 0 failures throughout. Results in `dev/simstudy/results/`. See §12. |

Building it surfaced four corrections to this plan and two package bugs;
both are recorded in place below rather than only in commit messages.

------------------------------------------------------------------------

## Running the suite

Six levels, cheapest first. Times are wall clock on Doug's machine and
include \~9 s of R startup and `load_all()` compilation, so the test
work itself is correspondingly less.

| Level | Command | Time | Covers |
|----|----|----|----|
| One file, while iterating | `Rscript -e 'devtools::test(filter="regression")'` | \~11 s | one tier-1 file |
| **Everything local** | `Rscript -e 'devtools::test()'` | \~45 s | tiers 1 + 2 |
| Exactly what CRAN runs | `NOT_CRAN=false Rscript -e 'devtools::test()'` | \~16 s | tier 1 only |
| Full package check | `Rscript -e 'devtools::check()'` | minutes | tests + examples + vignettes + `R CMD check` |
| One study cell | `Rscript dev/simstudy/run_study.R --cores=5 --scenarios=base --caffeinate` | \~50 min | tier 3, one cell at R = 100 |
| Full study | `Rscript dev/simstudy/run_study.R --cores=5 --caffeinate` | \~4.75 h | tier 3, all 10 cells at R = 100 |

Study times are **measured**, not projected: the 29 July run did all 10
cells, 1000 fits, in **285 min** on 5 cores. Ignore the 474 min recorded
for the 28 July grid -- that wall time was inflated by the laptop
sleeping mid-run; its compute rate was the same \~4 fits/min. Throughput
is flat across the whole run, so the cost scales linearly in R. An
earlier version of this table said `--cores=8` and \~2.4 h, which was a
guess and roughly 3x optimistic. On a fanless MacBook Air, 5 cores is
the practical ceiling; use `--caffeinate` for anything over an hour,
because a run that sleeps mid-way is easy to misread as thermal
throttling (it was, once -- check the cpu/wall ratio in the progress
line, which makes sleep unmistakable). `--resume` picks a run back up
from the last checkpoint.

**`devtools::test()` is the day-to-day command** -- run it after any
change to `R/` or `src/`. In RStudio, Cmd+Shift+T. `filter` is a regex
on the filename minus the `test-` prefix; the tier-1 files are `smoke`,
`regression` and `api`.

**After touching C++**, if you have been switching branches or the
object files look stale:

```         
Rscript -e 'pkgbuild::clean_dll(); devtools::test()'
```

**Before pushing**, use `devtools::check()` rather than
`devtools::test()`. It is the only way to see the tests as CRAN will run
them: `test()` uses `load_all()`, whereas `check()` installs the package
properly.

**The study runner never runs from `devtools::test()`.** Tier 3 is gated
on `OCCJSDM_SIMSTUDY`; setting it runs the study *serially* through the
test entry point, which is \~19 h. Use `run_study.R` for any real run --
it parallelises across replicates. `--R=`, `--cores=`, `--scenarios=`,
`--out=`, `--nburn=` and `--niter=` are all overridable.

### Two things to know before changing a test

**A failing test names its bug.** The regression tests are titled
`"Fixed bugs 22: ..."`, `"Fixed bugs 8: ..."` and so on, keyed to the
*Fixed bugs* list in `TODO.md`. Read that entry before debugging -- it
says what broke, when, and why it mattered.

**Do not "fix" a loose assertion by tightening it.** Several are
deliberately loose and the comment above each says why. The GP
length-scale test asserts that `sample_ls()` is *reached*, not that
`idx_ls` varies, precisely so it does not encode the open `sigma_s` bug
(§10.0); tightening it would make it pass for the wrong reason now and
fail confusingly once that bug is fixed. The shuffle-invariance test
asserts group totals rather than row-for-row equality because PCR
replicates are exchangeable, and ignores `row.names` because those were
verified to be **order-dependent but pairing-preserving** -- a user
looking up by name gets the right row either way, so re-tightening it
would fail on a difference that carries no information. The tier-2
thresholds are set from measurement across three seed sets (§6.3), not
from taste.

------------------------------------------------------------------------

## 1. Purpose

Check that `runOccJSDM()` recovers known generating parameters from data
produced by `simulateOccJSDMData()`, at the nominal rate, across the
configurations the package advertises.

Two distinct goals, which need different designs -- worth keeping
separate in the writing-up:

- **Regression.** Guard the \~23 bugs already fixed (`TODO.md`, *Fixed
  bugs*). None of them currently has a test; every one was verified by
  reading code, or by a one-off manual run. This is the near-term value.
- **Calibration evidence.** Support a claim in the MEE paper that
  inference is correctly calibrated. Needs more replicates and a tighter
  argument (§9).

The suite as specified below serves the first goal well and the second
adequately. Read §5 before quoting any coverage number in the paper.

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

- The `sample_ls` regression (`TODO.md` Fixed bugs 22) broke *every*
  non-spatial fit and survived because nothing exercised that path. A
  three-line test would have caught it the moment it landed. The
  vignette would not have -- it passes `spatCovariates`.
- It discharges Phase 3 step 11 of the CRAN plan in `AGENTS.md` ("write
  the testthat suite"). Done as of `6d9526d`: `test-placeholder.R` is
  gone, replaced by four real test files.
- No package-size or vignette-build cost. Both matter here:
  `sampleresults.rda` at 62 MB is already a CRAN blocker.

What this gives up is *presentation*. Tests answer "does it pass?", not
"here is the calibration evidence". See §9.

------------------------------------------------------------------------

## 3. Tiers

| Tier | Gate | Runtime | Runs on |
|----|----|----|----|
| 1\. Regression + smoke | none | \< 30 s | every `check()`, including CRAN |
| 2\. Recovery canary | `skip_on_cran()` | \~2 min | local + CI |
| 3\. Coverage study | env var `OCCJSDM_SIMSTUDY` | \~4.75 h | manual / nightly |

Tier 3 cannot be an ordinary test: CRAN wants the whole check under \~10
minutes.

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

These are the things most likely to make the suite wrong or annoying.
Each has bitten something already.

### 5.1 Never assert exact numeric equality of fitted values in tier 1

**Tier 1 assertions must be structural**: does it run, are dimensions
right, is a block non-`NA`, does a value vary. Numeric recovery belongs
in tiers 2–3, with tolerances.

**The original reason expired on 29 July 2026, and the rule survives
it.** This rule was written because a fixed seed did not reproduce the
sampler: `randinvg()` drew from R's global RNG inside an OpenMP loop, so
a tier-1 `expect_equal()` on posterior quantities would have passed
locally and flaked on CRAN. That is fixed (TODO.md Fixed bugs 28) --
`set.seed()` now controls the whole fit, verified across separate
processes.

The rule stands anyway, for a different and more durable reason:
**fitted values are not a contract**. Any legitimate change to a
sampler, a prior default, or the order of draws moves every number, and
a suite full of golden values turns each such change into a wall of
failures that says nothing about correctness. Structural assertions
survive those changes; numeric ones do not.

The one sanctioned exception is `test-regression-bugs.R`'s
reproducibility test, where equality *is* the property under test. Note
its scope: reproducibility holds for a given thread count, since threads
derive separate streams (see `src/rng.h`). This is currently moot --
`SHLIB_OPENMP_CXXFLAGS` is empty on the dev machine, so the pragmas
compile to no-ops -- but it would bite on a Linux build.

### 5.2 Coverage assertions are stochastic and will flake if set tight

With R = 100 the SE of an estimated 95% coverage is \~2.2%, so a
*correct* implementation lands outside ±1.96 SE about 5% of the time.
`expect_gt(cov, 0.90)` would fail intermittently for no reason, and a
suite that cries wolf gets ignored.

Mitigations, all three:

- fix a seed per replicate, so a tier-3 run is reproducible;
- set the hard failure threshold well below nominal -- fail at **\<
  0.85**, which catches every bug in the *Fixed bugs* list by a wide
  margin (§9);
- record the coverage table with `expect_snapshot_value()` so changes
  surface as a reviewable diff rather than a boolean.

### 5.3 Identifiability limits what can be checked at all

The latent-factor model is invariant to rotation and sign.
`reparamFactorModel()` (`R/jsdmfun.R`) imposes a QR-based constraint,
but there is no reason that constraint matches the simulator's
parameterisation.

- **Directly comparable** (verified against the code while building
  `helper-simstudy.R`, 27 July): `B0`, `B`, `G`, `beta_theta`, `theta0`,
  `p`, `q`, `tau`. Two caveats. `beta_theta` requires passing
  `collCovariates` to the fit -- without it the posterior is `(1 x S)`
  against a `(ncov_theta+1 x S)` truth. And `n_lattrait` must be pinned
  to the simulated `gt`, since the default
  `floor(sqrt(min(S, ncov_psi)))` will not generally match, giving
  `A`/`C` mismatched dimensions.
- **`sigma_b`: keep, but read with care.** With a tight
  `InvGamma(10, 1)` prior and only `p x S` residual coefficients the
  posterior sits at the prior mean (`sqrt(1/9)` = 0.333), so a fixed
  truth far from that yields near-zero coverage through prior-data
  conflict, not sampler error. Measured: true 0.5 gave 0/8, true 0.35
  gave 8/8. Scenarios set it near the prior mean; SBC is the principled
  fix.
- **NOT comparable, contrary to the original draft of this section:**
  `Bs`, `Gs`, `sigma_bs`. The simulator builds the spatial field
  directly from `sigma_s`/`l_s` and sets `Bst <- matrix(0, S, ps)`
  (`R/jsdmfun.R:557`), so on the truth side `Bs` is empty *regardless of
  `ds`* (re-checked at `ds = 2`) and `sigma_bs` generates nothing; the
  fit meanwhile represents the field as sparse-GP basis coefficients
  over `n_supportpoints` knots. Different parameterisations. Measured
  for `sigma_bs`: true 0.5 against a posterior mean of \~1.6, 0/8
  coverage -- a meaningless comparison rather than a bug.
- **Only as identified functions:** residual correlation matrix,
  `eta`/`psi`, variance partitioning
- **Must NOT be checked element-wise:** `U`, `L`, `A`, `C` --
  element-wise coverage would fail for reasons that are not bugs
- **Excluded:** `sigma_h` -- never sampled (group A item 1), so it would
  fail by construction. **Also `l_s`** -- group A item 2: `sigma_s` is
  hard-coded to 1 at the `sample_ls()` call site, so the length-scale
  absorbs the amplitude misspecification and rails at the top of
  `l_s_grid` for every true value tried. Any coverage figure for `l_s`
  is meaningless until that is fixed. Exclude explicitly and with a
  comment, so a future reader does not "fix" the test.

### 5.4 The GP knot count must be pinned, not left to the default

`getDefaultSupportPoints(n)` used to be `max(30, floor(n * 0.2))` -- a
constant 30 for any dataset below 150 sites, and a crash below 31
(§10.1). Fixed in `42198d9`; it is now `min(floor(n * 0.2), n - 1)`
(`R/jsdmfun.R:875`, `TODO.md` Fixed bugs 29). The study pins
`n_supportpoints` anyway rather than trusting any default.

Two consequences for the design, beyond the crash. The knot count would
be **uncontrolled but not constant** across cells -- identical for
`n = 40` and `n = 140`, then suddenly proportional above 150 -- so any
difference attributed to `n` would be partly an artefact of the spatial
approximation changing underneath. And at small `n` it approaches one
knot per site, which is a qualitatively different (denser) approximation
than at `n = 100`, precisely in the low-information cell where
calibration is most at risk.

Set `n_supportpoints` explicitly in every spatial cell so it is a
controlled factor.

### 5.5 CRAN limits

Two cores maximum (`_R_CHECK_LIMIT_CORES_`). Tier 3's replicate-level
parallelism must be off or capped when not env-gated on.

------------------------------------------------------------------------

## 6. Tier detail

### 6.1 Tier 3 -- the study

Base = `two_stage`, traits on, spatial on, `d = 2`, **`ds = 2`**,
`n = 100`, `S = 10`, `M = 2`, `P = 2`, `K = 3`, `n_supportpoints = 20`.

`ds >= 1` *used to be* required or there was no spatial field at all: at
`ds = 0` the simulator's cross-species spatial covariance collapsed to
jitter, `sd(spatField)` = 0.0019 rather than \~1.0. The grid used
`ds = 0` until 27 July, which made every spatial cell a null-field test
-- the exact trap the audit flagged for Fixed bugs 7. Fixed in `42198d9`
(`TODO.md` Fixed bugs 30); re-measured 29 July, `ds = 0` now gives
0.598. The grid stays at `ds = 2` by choice, for a clearly
non-degenerate field and comparability with the 28 July run.

Every spatial cell sets `n_supportpoints` explicitly (§10.1). Left to
the default it is a constant 30 for any `n` below 150, so it would
neither scale with `n` nor be a controlled factor -- and below 31 sites
it crashes outright.

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

Design is **one-factor-at-a-time**, not factorial (\~48 cells would be
unaffordable). Cost: OFAT cannot see interactions. Cell 1 covers the one
interaction judged most likely to matter, since traits and spatial both
write into the same linear predictor via `computePsiCoef()`.

Deliberately *excluded* from tier 3, because tier 1 answers them in
seconds for free: "does configuration X run at all". Tier 3 hours should
buy calibration evidence, which only tier 3 can provide.

Cells 9–10 are **cheaper** (less data per fit: no PCR level, no
replication) but also check **fewer parameters** -- `occupancy` has no
`p`/`q`; `binary` has no `beta_theta`, `theta0`, `p`, `q`. Their rows in
the results table will be sparser, not directly comparable to the
two_stage cells.

**Metrics** per scenario × parameter block, pooled across species:
coverage of the 95% credible interval, bias, RMSE.

**Output:** a tibble (`scenario`, `block`, `coverage`, `bias`, `rmse`,
`n_par`, `R`) written to
`Sys.getenv("OCCJSDM_SIMSTUDY_OUT", tempdir())`, intended destination
`dev/simstudy/results/`.

### 6.2 Tier 1 -- ships to CRAN

One shared fixture fit built in `helper-fixtures.R` and reused:
`sampledata`, `nburn = 50, niter = 50, nchain = 2` (**measured 5.6 s**).
Keeps the tier well under 30 s.

`test-smoke-configs.R` -- every configuration completes: `two_stage` /
`occupancy` / `binary`; with and without `spatCovariates`; with and
without traits; `d = 0` and `d = 2`; `P = 1` and `P = 3`. *This is the
tier that would have caught the `sample_ls` regression.*

`test-regression-bugs.R` -- one test per fixed bug, named so a failure
points at its `TODO.md` entry:

| Guards           | Assertion                                           |
|------------------|-----------------------------------------------------|
| Fixed bugs 22    | non-spatial fit runs                                |
| Fixed bugs 8     | `tau_output` not all `NA`                           |
| Fixed bugs 21    | `plotFPTPStage2Rates()` differs across primers      |
| Fixed bugs 23    | `computeSpeciesDetected()` indices span both chains |
| Fixed bugs 9     | WAIC counter guard does not trip                    |
| Fixed bugs 1, 11 | row-shuffle invariance                              |

Row-shuffle invariance must be asserted on the **prepared design
matrices** (`X_psi`/`X_theta` row alignment against site/sample
indices), which is deterministic -- not on the posterior, which would
flake under §5.1.

`test-api-contracts.R` -- unknown `primerName` errors with the available
list; unknown covariate name errors informatively; `idx_species`
subsetting returns the right rows; `predictNewSites()` honours its
documented `NULL` defaults.

### 6.3 Tier 2 -- canary

R = 3–5 replicates at the base cell. Not a calibration check: a smoke
signal that recovery has not grossly broken between full tier-3 runs
(e.g. correlation between estimated and true `B` above a loose floor).

**RESOLVED: tier 2 fails.** Thresholds were set from measurement rather
than taste, which is what makes failing safe. Three independent seed
sets gave overall coverage 0.89/0.89/0.88, lowest block 0.78/0.76/0.76,
and `cor(post_mean, truth)` 0.62/0.56/0.68; the floors are 0.70, 0.40
and 0.30. A pass is therefore not luck and a failure is not noise.

------------------------------------------------------------------------

## 7. SBC seam

Coverage now, simulation-based calibration later, without restructuring.
Two swappable functions in `helper-simstudy.R`:

``` r
draw_truth(scenario)        # fixed values now; prior draws for SBC
summarise_fit(fit, truth)   # CI-contains indicator now; rank statistic for SBC
```

Everything expensive -- grid, simulate/fit loop, parallelism, artifact
writing -- is shared.

Adding SBC later requires:

1.  A prior-draw wrapper. All priors in the sampler are **proper**
    (checked): `B0`, `L` \~ N(0,1) (`sample_BBsL()`, `R/jsdmfun.R`); `B`
    trait-predicted mean with variance `sigma_b^2`; `sigma_b`,
    `sigma_bs` \~ sqrt-InvGamma(10,1); `tau` \~ InvGamma(5,5)
    (`R/runOccJSDM.R`); `p`, `q`, `theta0` \~ Beta; `l_s` \~ Gamma(1,1)
    over a 10-point grid. So SBC is feasible.
2.  Documenting what the sampler *actually* assumes. Not free: the audit
    found the *documented* `prior_beta_psi`/`prior_beta_psi_sd` were a
    complete no-op, with the live path hard-coding something different
    (Fixed bugs 14).
3.  Thinning to near-independent draws, which raises per-replicate cost.

Note `sigma_h` would show a flagrantly non-uniform rank histogram under
SBC. That is the method working (group A item 1), not a new bug.

------------------------------------------------------------------------

## 8. Runtime budget

Measured on `sampledata` (100 sites, 10 species, 2 chains): **0.039 s
per iteration**, so `nburn = 1000 + niter = 1000` is \~84 s per fit.

|                   | Serial  | 5 cores      |
|-------------------|---------|--------------|
| One cell, R = 100 | \~2.3 h | \~28 min     |
| 10 cells          | \~23 h  | **\~4.75 h** |

**Superseded by measurement**: the 29 July run did all 10 cells in 285
min on 5 cores, i.e. \~4.75 h. The earlier \~2.5–3 h figure assumed 8
cores, which this machine cannot usefully supply -- the M4 has only 4
performance cores, so workers past the fourth land on efficiency cores
worth roughly a third to a half as much (see §5.5). Cells 9–10 are
cheaper than the average; cell 8 (`species_20`) is dearer.

Replicates are embarrassingly parallel at the *process* level --
independent datasets and fits -- so this needs no package changes, unlike
the in-package chain parallelisation in `TODO.md` group D item 1.

**Do not buy replicates by shortening chains.** Each interval endpoint
is a tail quantile; with 500 total draws only \~12 land below the 2.5%
bound, so the endpoints are noisy and that noise feeds into coverage as
bias, not just variance. Set chain length by an ESS target (\>= 400 on
monitored parameters, which `returnConvergenceDiagnostics()` already
reports) and let R take what is left.

------------------------------------------------------------------------

## 9. What R = 100 does and does not support

SE of estimated coverage is `sqrt(0.95 * 0.05 / R)`. Cost is linear in
R; precision improves only as sqrt(R) -- halving the error bar costs 4x
the compute.

> **Resolved 29 July 2026.** This previously carried a caveat that the
> replicates were not independent: `get_rng()` was seeded from a
> literal, per *process*, so every PSOCK worker started its Polya-Gamma
> stream at the same position. That is fixed (TODO.md Fixed bugs 28) --
> `runOccJSDM()` now derives its C++ seed from R's RNG, and because
> `draw_truth()` sets a per-(scenario, replicate) seed before each fit,
> every replicate gets its own stream. Verified: a repeated replicate
> reproduces exactly, and two replicates differ on every column.
> `sqrt(0.95 * 0.05 / R)` is therefore the right SE for the re-run.
>
> §12 now reports the 29 July re-run, which was made *after* this fix,
> so its error bars are sound. It is the superseded pre-fix table at
> §12.6 that inherits the old correlation: treat those figures as
> indicative only, and never quote an SE against them.

| R       | Coverage SE | Undercoverage reliably flagged |
|---------|-------------|--------------------------------|
| 50      | 3.1%        | below \~0.890                  |
| **100** | **2.2%**    | **below \~0.907**              |
| 200     | 1.5%        | below \~0.920                  |
| 400     | 1.1%        | below \~0.929                  |

R = 100 is chosen for the **regression** goal. Every bug in the *Fixed
bugs* list would drive coverage to 0.5–0.8, not 0.93, so R = 100 catches
them with room to spare.

It is **not** enough to assert "coverage is nominal" in the paper -- that
claim needs R = 200–500, because it asserts the *absence* of a small
deviation rather than detecting a large one. The runner takes `R` as an
argument for exactly this reason; do a final publication-grade run at
higher R if the paper makes that claim.

Also: pooling coverage across species within a block buys precision, but
those indicators share one simulated dataset and one MCMC run, so
effective N sits between R and R x S, nearer R. **Do not quote an SE
computed as if there were R x S independent trials.**

------------------------------------------------------------------------

## 10. Open items

0.  *OPEN, and it degrades cells 1 and 2.* **The GP length-scale is
    never recovered** (`TODO.md` group A item 2, re-diagnosed 27 July).

    `sample_ls()` treats the fitted field `SE` as a GP draw and scores
    its density under `N(0, sigma_s^2 K(l_s))` while holding `SE` fixed.
    But under the SoR approximation `SE = Ks(l_s) %*% Bs` is a
    *deterministic function of `l_s`*, already smoothed at the current
    value, so it scores better and better under smoother covariances.
    The likelihood is monotone increasing in `l_s` and `idx_ls` pins at
    the grid maximum -- measured at range 10-10 for every true `l_s`
    tried (0.074, 0.171, 0.300) with real spatial signal present.

    **This supersedes the earlier entry here, which blamed `sigma_s`
    being hard-coded to 1 and expected a one-line fix.** Profiling the
    grid at the field's actual amplitude leaves the likelihood monotone,
    so that fix is disproven. `sigma_bs` is already sampled and already
    is the SoR amplitude; a separate `sigma_s` would be redundant. The
    real fix is a modelling change needing Alex -- see the TODO entry.

    Consequences here: `l_s` stays out of the checked parameters (§5.3),
    and the spatial cells still test `Bs` and the overall spatial
    contribution but say nothing about range recovery. **Not blocking**
    -- those cells remain worth running -- but no table produced before
    this is fixed may be read as evidence about the GP. Re-run cells 1
    and 2 afterwards.

1.  ~~*OPEN, blocking tier 3.***Minimum `n` for spatial cells.**~~
    **RESOLVED 27 July 2026, and no longer blocking.** Measured: the
    floor is **31 unique locations** with default settings. `n <= 29`
    fails with `cannot take a sample larger than the population`,
    `n = 30` with
    `number of cluster centres must lie between 1 and nrow(x)`,
    `n >= 31` runs.

    Cause: `getDefaultSupportPoints(n) <- max(30, floor(n * 0.2))`
    (`R/jsdmfun.R:480`) feeds `kmeans(X_s, centers = ps)`, and the
    hard-coded 30 wins for every dataset under 150 sites, so the knot
    count never scales down. Filed as a bug and **since fixed**
    (`42198d9`, `TODO.md` Fixed bugs 29) -- the default was also
    statistically wrong well before it crashed (30 knots for 31 sites is
    \~1 knot per site, defeating the point of a sparse GP). Now
    `min(floor(n * 0.2), n - 1)`.

    **Consequence for this plan:** cell 5 must set
    `listParams$n_supportpoints` explicitly rather than relying on the
    default. Verified working at `n = 20` with 6 knots. Pinning the knot
    count is the right call for a study cell anyway -- it stops the
    spatial approximation silently changing with `n` across scenarios,
    which would otherwise confound the low-information comparison. Apply
    the same reasoning to **every** spatial cell (1, 2, 4–10), not just
    cell 5: set `n_supportpoints` per scenario so it is a controlled
    factor rather than a function of `n`.

2.  ~~*OPEN.* Tier 2 failing vs advisory (§6.3).~~ **RESOLVED: it
    fails**, on measured thresholds -- see §6.3.

3.  ~~*OPEN.***`beta_theta` and `resid_cor` sit below nominal.**~~
    **RESOLVED 28 July by the full run, and both traced to code.**
    `beta_theta` undercovers in *every* cell (0.676-0.730) -- the
    flatness across model types, primer counts and factor
    misspecification is the signature of a structural cause, and it is
    `TODO.md` group A item 3 (prior mean 1 instead of 0 on the
    collection slopes). `resid_cor` sits at 0.74-0.77 in nine cells and
    is group A item 4 (`reparamFactorModel()` breaking `t(L) %*% L`).
    See §12.

4.  *Deferred.* Presentation of tier-3 results for the paper. Build the
    suite first; the summary object then feeds either a short pkgdown
    article or the manuscript directly. Nothing here forecloses that.

5.  *Deferred, but the case is now stronger.* **SBC** (§7). The full run
    found two separate instances of a fixed truth conflicting with an
    informative prior -- `sigma_b` reading 1.000 everywhere because it
    is prior-dominated, and `p` collapsing to 0.103 where the true value
    sits far into the tail of a `Beta(5, 1)`. Both are artefacts of
    choosing truth independently of the prior. SBC draws truth *from*
    the prior and so cannot produce them, which would separate "the
    sampler is wrong" from "this truth disagrees with this prior"
    automatically rather than by the hand-analysis those two took.

------------------------------------------------------------------------

## 12. Results of the full run (29 July 2026, post-fix)

1000 replicates, 10 scenarios, 0 failures, **285 min** on 5 cores.
155,578 individual interval checks. Commit `8af22cd`. Coverage SE at R =
100 is 2.2 points, so treat anything in 0.93–0.97 as indistinguishable
from nominal.

**This is a paired comparison with the 28 July pre-fix run.**
`draw_truth()` seeds on (scenario, replicate), so the simulated data and
true values are *bit-identical* between the two runs -- verified,
`max|truth difference| = 0`. Every difference below is attributable to
the code, not to sampling variation between runs. That is a much
stronger design than two independent runs would give, and it is worth
preserving: **do not change `simstudy_seed()`**, or future runs lose
comparability with these.

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

**Bold = undercovers. Italic = overcovers. `--` = not estimated in that
cell.**

### 12.1 What the fixes changed

Averaged over all cells, pre-fix → post-fix:

| Block         | Coverage          | Bias                | Verdict                |
|---------------|-------------------|---------------------|------------------------|
| `beta_theta`  | 0.719 → **0.766** | +0.112 → **+0.038** | improved, not fixed    |
| `theta0`      | 0.901 → 0.944     | −0.020 → −0.004     | moved, see 12.3        |
| `resid_cor`   | 0.775 → 0.777     | +0.006 → +0.002     | **unchanged**          |
| `p`           | 0.821 → 0.827     | +0.087 → +0.086     | unchanged, as expected |
| `B0`          | 0.947 → 0.943     | −0.135 → **−0.228** | **bias doubled**       |
| `B`, `G`, `q` | \~0.92–0.94       | \~0                 | stable                 |

1.  **`beta_theta` improved but is not fixed.** Every cell gained
    0.03–0.05, and two-thirds of the bias is gone. Alex's correction of
    the collection-covariate prior mean from 1 to 0 was a real cause --
    but not the only one, because 0.766 against nominal 0.95 is still
    far out. **This is now the clearest open target in group A.**
    Whatever remains is structural: it is flat across model type, primer
    count, species count and factor misspecification, exactly as the
    pre-fix version was.

2.  **`resid_cor` is untouched, and the paired design makes that
    conclusive.** Only **104 of 49,978** coverage decisions flipped.
    Identical data, identical truth, coverage unmoved at 0.777. This is
    `reparamFactorModel()` (group A item 2), which is unfixed and
    currently disputed. The `d_underfit` exception persists exactly as
    before: *over*covering at 0.980, because under-fitting the
    ordination widens intervals enough to mask the bias.

3.  **`p` is unchanged, and that is the correct outcome.**
    `low_information` sits at 0.109. The informative `Beta(5, 1)` is
    load-bearing for identifiability -- `p` and `q` enter
    `sample_pq_cpp()` symmetrically, so the prior is what selects the
    mode (TODO.md group A item 3). This cell measures the *cost of that
    constraint* when true `p` is 0.1–0.3. It is not a defect to be fixed
    by flattening the prior.

### 12.2 New finding: `B0` bias roughly doubled

Nine of ten cells moved more negative: base −0.113 → −0.208, `occupancy`
−0.024 → −0.151, `primers_3` −0.031 → −0.151, `low_information` −0.931 →
−1.056. Only `binary` moved the other way (−0.002 → +0.012).

**Coverage does not reveal this** -- it holds at 0.943 because the
intervals are wide enough to absorb the shift. It is visible only in the
bias column, which is the argument for tracking both.

Nothing in the four fixes should move species intercepts this way. The
likeliest candidate is the 421-line `jsdmfun.R` rewrite that shipped in
the same pull (`42198d9`). Species intercepts are a headline quantity
for a JSDM, so this needs a cause before the paper uses `B0`. An R = 8
probe on 29 July hinted at it; R = 100 confirms it.

### 12.3 New finding: `theta0` now overcovers

The all-cell average of 0.944 is misleading. Per cell it is
**0.978–0.985 in nine cells** (was 0.938–0.959) and **0.602 in
`low_information`** (was 0.477). So it moved from mildly under to
distinctly over, except where information is thin.

Overcoverage is the safe direction, but the pattern suggests the widened
`diag(2)` prior variance on the collection coefficients overshot. Worth
a look alongside item 1 above, since both trace to the same edit.

### 12.4 What holds

`q`, `B0`, `B` and `G` are at or near nominal in every cell except
`low_information`. The quantities most likely to be reported in an
ecology paper remain trustworthy -- with the `B0` bias caveat in 12.2,
which affects the point estimate rather than the interval.

`low_information` is still compromised across the board (`theta0` 0.602,
`p` 0.109, `B` 0.826, `G` 0.840, `B0` bias −1.06). Users with small or
low-detection datasets remain the most exposed.

### 12.5 Caveats that travel with this table

`l_s` and `sigma_h` are absent because neither is recoverable (group A
item 1), so **no cell says anything about spatial range**. `sigma_b`
reads 1.000 everywhere because it is prior-dominated by construction
(§5.3), not because it is well estimated. Differences below 2.2 points
are noise.

### 12.6 The superseded pre-fix table (28 July 2026)

Kept because the *delta* is the evidence that the fixes worked, not the
level. Run in two parts -- `base` on 27 July (100 fits, 22.9 min), the
other 9 cells on 28 July (900 fits, 474 min wall, inflated by the laptop
sleeping mid-run).

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

## 11. Build order

1.  `helper-fixtures.R` + tier 1. Ships, pays immediately, blocks
    nothing.

2.  ~~Resolve open item 1 (spatial `n` floor).~~ **Done** -- floor is 31
    sites; every spatial cell now pins `n_supportpoints` (§10.1).

3.  `helper-simstudy.R` with the two seams (§7).

4.  Tier 2.

5.  Tier 3 + `run_study.R`.

6.  **Full R = 100 run; inspect the table.** **Done 27-28 July 2026** --
    all 10 cells at R = 100, 1000 fits, 0 failures, \~497 min on 5 cores
    (`base` alone first, 100 fits/22.9 min; then the other 9 cells, 900
    fits/474 min). Results in `dev/simstudy/results/`, tabulated in §12.
    Found three undercovering blocks, each since traced to code; see §12
    and TODO.md group A.

    **A re-run is now owed, and it checks two things at once.** Alex's
    28-29 July fixes (`sigma_h`, the collection-covariate prior mean,
    the Stage 2 prior wiring) all change fitted values, so §12 is stale
    on its own terms. Separately, the replicates in that run were not
    independent -- see the note under §8 -- so its error bars were
    understated. The re-run therefore measures whether the fixes worked
    *and* produces the first numbers with honest SEs.

    ```         
    Rscript dev/simstudy/run_study.R --cores=5 --caffeinate
    ```

    Two cautions carried over from the first attempt. **Do not edit
    `run_study.R` while it is running** -- Rscript parses incrementally,
    and an earlier ten-hour run died on a torn read. Use `--resume` if
    it is interrupted. And an indicative post-fix probe at R = 8 showed
    `beta_theta` still undercovering at 0.775 with its bias halved,
    while `B0`'s bias *grew* fivefold to -0.529; R = 8 settles nothing,
    but it means the re-run should be read for regressions from Alex's
    421-line `jsdmfun.R` rewrite as well as for the intended
    improvements.

------------------------------------------------------------------------

## 13. The M ladder: are B4-B6 defects, or Stage 1 under-identification?

**Status: planned, not yet run.** Commissioned by Doug, 29 July 2026, as
the next step on `TODO.md` group B items 4, 5 and 6, each of which is
annotated CLAUDE TO RUN SIMULATION STUDY WITH M \> 10 AND SEE IF IT
FIXES IT.

### 13.1 The hypothesis, and why one lever could explain three findings

All three open findings sit downstream of the *field-collection* stage,
and the grid has only ever been run at `M = 2` samples per site:

- **B4**, `beta_theta` undercovering at 0.766, is a Stage 1 covariate
  slope.
- **B5**, `theta0` overcovering at 0.978-0.985, is the Stage 1
  species-specific baseline.
- **B6**, `B0` bias doubling to -0.228, is the occupancy intercept --
  which is inferred *through* the latent collection state `w`, so it
  inherits whatever Stage 1 gets wrong.

At `M = 2` there are two Bernoulli collection draws per site-species,
observed only indirectly through PCR detections, from which the sampler
must identify a species baseline *and* covariate slopes *and* the latent
`w` that `z` then depends on. That is very little information, and it is
a plausible single cause for all three: weak identification inflates the
prior's influence, which produces exactly the observed signature --
slopes too confident, baselines too diffuse, and an occupancy intercept
biased toward under-inferred occupancy (negative, as measured).

**If true, this is not a bug.** It is a design requirement: the model
needs `M` above some threshold before Stage 1 is identified. That is a
documentation and paper item rather than a code fix, and an important
one, because `M = 2` to `3` is entirely realistic in eDNA practice.

### 13.2 Design

A **ladder**, not a single M \> 10 point. One high-M run answers "fixed
or not"; a ladder answers "fixed, and from where" -- and distinguishes a
deviation that decays smoothly toward nominal (information-limited) from
one that plateaus at a non-nominal level (a real defect that more data
cannot reach).

| Arm   | M   | K   | Rows  | Purpose                                          |
|-------|-----|-----|-------|--------------------------------------------------|
| `M2`  | 2   | 3   | 1200  | Reproduces the current base cell; internal check |
| `M5`  | 5   | 3   | 3000  | Is the improvement already underway?             |
| `M10` | 10  | 3   | 6000  | Doug's threshold                                 |
| `M20` | 20  | 3   | 12000 | Comfortably past it                              |
| `K30` | 2   | 30  | 12000 | **Control** -- matched rows, wrong stage         |

**`K30` is what makes the result interpretable.** Raising `M` raises the
total observation count, so if everything improves, the naive reading
"more data helps" is unfalsifiable. `K30` holds the row count identical
to `M20` but spends it on PCR replicates -- Stage 2 -- instead of field
samples. If `M20` fixes the three findings and `K30` does not, the
mechanism is specifically Stage 1 identification. If both fix them, it
is merely sample size, and the conclusion is much weaker.

**The arms share a truth, so the ladder is paired.**
`simstudy_seed_for()` honours a `seed_label` field, and all five arms
set `seed_label = "mladder"`, so replicate *r* draws the same truth in
every arm. The differences between arms then stop carrying the variance
of independent truths -- which is precisely what a ladder reads.

**How far the pairing extends, measured rather than assumed:**

| Quantity | Paired across M? |   |
|----|----|----|
| `B0`, `B`, `G`, `A`, `C`, `L`, `sigma_b`, `tau` | **yes** | the whole JSDM layer |
| `z_true` | **yes** | latent occupancy |
| `p_true`, `q_true`, `theta0`, `theta_baseline` | **yes** | Stage 2 rates and Stage 1 baseline |
| `beta_theta` intercept row | **yes** | it is `logit(theta_baseline)` |
| `beta_theta` slope rows | **no** | see below |
| `w_true` | n/a | its dimension is `N x S`, so it cannot match |

`beta_theta`'s slopes are drawn at `R/simulateData.R:128` with
`sample(c(-1,1,0), ...)` *after* `simulateOccJSDMData()` has built the
`N`-sized `X_theta`, so by then the stream has diverged and there is no
`list_params` hook to inject them. The residual is bounded -- slope
signs from a three-point set, with the intercept row paired -- and the
estimand is unchanged, so the ladder stays unbiased for `beta_theta`,
just noisier than for the blocks above. **B6 (`B0`) and B5 (`theta0`)
get the full paired benefit; B4 (`beta_theta`) gets a partial one.** The
`M20` vs `K30` control comparison is likewise unpaired for `beta_theta`,
since those arms differ in `N`.

`test-recovery.R` asserts all of this. It fails if anyone adds a draw to
`draw_truth()` sized by `N`, `M` or `K`, which would silently unpair the
ladder.

Base cell only: all three findings are flat across the other nine cells,
so one cell isolates the mechanism at a tenth of the cost.

**`R = 50`, but only for the first stage, and the reason is
asymmetric.**

Detecting that the deviation *persists* is easy: testing 0.766 against
nominal uses the null SE, `sqrt(0.95 * 0.05 / 50)` = 3.1%, so the
18-point gap is about 6 SE.

Confirming it has *recovered* is much harder, and R = 50 cannot do it.
That claim needs an interval around the observed value, where the SE is
`sqrt(p(1-p)/R)` at the observed p:

| R   | SE at p \~ 0.94 | 95% CI around an observed 0.94 |
|-----|-----------------|--------------------------------|
| 50  | 3.4%            | [0.87, 1.00]                   |
| 100 | 2.4%            | [0.89, 0.99]                   |
| 200 | 1.7%            | **[0.91, 0.97]**               |

At R = 50 an observed 0.94 is indistinguishable from a true 0.88, which
is still meaningful undercoverage -- so "fixed" could not be told from
"much better but still broken". This is §9's point restated: asserting
the absence of a small deviation is a different and more expensive claim
than detecting a large one.

**So: two stages.** R = 50 across the ladder answers the question Doug
actually asked -- *is `M` the lever?* -- because the answer is carried
by the **shape across arms**, and a monotone rise over four arms is
evidence no single point can supply. Then re-run only the decisive arm
at **R = 200** if it lands high, since that is the one number that would
appear in the paper. Roughly 1.9 h now and 2.4 h later, the second spent
only if the first justifies it.

### 13.3 Cost (measured 29 July, not projected)

Per fit at the study's `nburn = 1000 + niter = 1000`, extrapolated from
a 200-iteration timing on this machine:

| M   | s/fit | R = 50 on 5 cores |
|-----|-------|-------------------|
| 2   | 40    | \~8 min           |
| 5   | 64    | \~13 min          |
| 10  | 100   | \~20 min          |
| 20  | 179   | \~36 min          |

Plus `K30` at roughly the `M20` cost. **Total \~1.9 h.** Note the cost
scales *sub*-linearly in rows -- M = 20 is 4.5x M = 2, not 10x --
because much of the per-iteration work is over species and sites, not
samples.

### 13.4 What each outcome means

Decided in advance, so this is a test rather than a fishing trip.

| Outcome | Reading | Action |
|----|----|----|
| All three approach nominal as M rises; `K30` does not | Stage 1 under-identification. Not defects. | **Re-run the winning arm at R = 200 before closing anything** -- R = 50 cannot separate 0.94 from 0.88. Then close B4-B6 as "not a bug", document the `M` requirement, and state it in the MEE paper: it is a real limitation for `M = 2` eDNA designs. |
| B4/B5 recover, B6's bias persists | Two causes. Stage 1 explains the coverage findings; `B0` is a separate regression. | Close B4/B5; keep B6 open and go after the `jsdmfun.R` rewrite. |
| Nothing improves, even at M = 20 | Genuine code defects; `M` is not the lever. | Keep all three open. Next suspect is the widened `diag(2)` prior variance, which is the common edit behind B4 and B5. |
| `M20` and `K30` improve equally | Sample size, not Stage 1. | Weak conclusion; the finding is that all three are information-sensitive. Re-think the design before spending more. |

### 13.5 Implementation

No machinery change needed: `M` and `K` are already scenario fields, and
`mk()` overrides `SIMSTUDY_BASE` by name.

1.  ~~Add the five arms to `simstudy_scenarios()`.~~ **Done**: `M2`,
    `M5`, `M10`, `M20`, `K30`, each with `seed_label = "mladder"`. They
    are inert for a full run, because the runner takes the whole list
    only when `--scenarios` is empty -- so select them explicitly.
    `simstudy_seed_for()` is the seed override; `simstudy_seed()` is
    unchanged and still keys on the label for every other cell.
2.  Sanity-run at `R = 4` across all five arms (\~10 min) to confirm the
    arms fit and the row counts are as expected, before committing the
    real run.
3.  `Rscript dev/simstudy/run_study.R --scenarios=M2,M5,M10,M20,K30 --R=50 --cores=5 --caffeinate`
4.  Report coverage **and bias** for `beta_theta`, `theta0` and `B0`,
    plus the remaining blocks to confirm nothing degrades.

**Check first:** the `M2` arm should agree with the existing base cell
(`beta_theta` 0.763, `theta0` 0.983, `B0` bias -0.208) within its SE.
Its seeds differ -- `simstudy_seed()` keys on the label -- so it will
not reproduce them exactly, but a disagreement beyond about 2 SE would
mean the ladder is measuring something other than what it claims, and
should be resolved before reading the rest.

### 13.6 Caveats

`M` and the number of *sites* are not interchangeable. This ladder holds
`n = 100` fixed and varies samples per site, which is the quantity Stage
1 is identified by. It says nothing about whether more sites would help,
and should not be quoted as if it did.

Raising `M` also raises the number of latent `w` states being sampled,
so mixing may differ across arms. If a high-M arm shows worse mixing
rather than better coverage, check `returnConvergenceDiagnostics()`
before concluding anything about identification.


------------------------------------------------------------------------

## 13.7 Results (29 July 2026)

**Run:** 250 fits, R = 50, 5 arms, 189 min, 0 failures. Validity check
passed: the `M2` arm agrees with the existing base cell on all three
targets (largest gap 0.016, against an SE of \~0.09 for the difference)
-- the ladder is measuring what it claims.

| | M2 | M5 | M10 | M20 | K30 (control) |
|---|---|---|---|---|---|
| `beta_theta` coverage | 0.747 | 0.655 | 0.603 | **0.579** | 0.706 |
| `theta0` coverage     | 0.986 | 0.966 | 0.944 | 0.952 | 0.996 |
| `B0` bias             | -0.160 | -0.025 | +0.002 | -0.028 | -0.091 |

**None of the four outcomes in §13.4 fits.** The actual result is a
mixture across the three items:

**`theta0` (B5): clean confirmation.** Overcoverage falls from 0.986
toward nominal as `M` rises (0.944 at `M10`), while `K30` -- matched row
count, wrong stage -- makes it *worse* (0.996). This is exactly outcome
1's signature: `M` is the lever, matched data volume elsewhere is not.
**Close as Stage 1 under-identification, not a defect** -- pending the
R = 200 confirmation §13.4 requires before closing anything.

**`B0` (B6): partial confirmation.** Bias collapses from -0.160 to near
zero at `M10`. But `K30` also improves it (-0.091), just less than `M10`
or `M20` -- so the recovery is not cleanly Stage-1-specific; more data of
either kind helps. Coverage stays 0.94-0.96 throughout in every arm,
consistent with the original finding that coverage does not reveal this
bias.

**`beta_theta` (B4) goes the wrong way, and this is the important
result.** Coverage falls *monotonically* with `M` -- 0.747 -> 0.579, a
17-point drop across four points, far beyond the R = 50 noise floor --
while bias stays small and flat (+0.02 to +0.05) throughout. Shrinking
intervals around a bias that is not itself shrinking is the signature of
a real defect being *exposed* by more information, not resolved by it.
`K30` (0.706) beats every M arm above `M2`.

**This rules out under-identification as the explanation for B4.**
`beta_theta` needs a different next step than B5/B6: not "wait for more
data" but "find what makes the interval overconfident". The prime
suspect is the same `diag(2)` prior-variance widening implicated in B5
(Fixed bugs 25) -- but if that edit was meant to fix things and B4 gets
*worse* with more information, it may need to be revisited rather than
extended.

**Unplanned finding: `q` degrades hard with `K`.** 0.945 (`M2`) -> 0.614
(`K30`), worse than any M arm, including `M20` (0.742) which shares `K30`'s
row count. Not investigated here -- outside what this ladder was built to
answer -- but worth its own entry; see `TODO.md`.

**What this does not settle.** R = 50 detects that `beta_theta` is
getting worse with high confidence (a 17-point monotone trend dwarfs the
3.1% SE many times over), but per §13.2's own caveat, R = 50 cannot
confirm that `theta0` has *reached* nominal -- an observed 0.944 is not
distinguishable from a true 0.90. Do not close B5 on this run alone.

------------------------------------------------------------------------

## 13.8 Follow-up: tighter `beta_theta` prior at high M

**Status: running.** Commissioned by Doug, 29 July 2026, directly off
13.7's finding that `beta_theta` coverage worsens with M -- the
signature of an overconfident interval rather than insufficient data,
which the write-up pointed at `B_betatheta`'s slope variance.

**That variance was hard-coded**, with no `listPriors` hook, unlike
`p`/`q`/`theta0`. Added one: `listPriors$b_betatheta_slope_var`
(`R/runOccJSDM.R`), defaulting to 2 -- the existing value -- so nothing
changes unless a caller sets it. `ALEX TO REVIEW` before treating a
non-default value as a real fix rather than a diagnostic; see `TODO.md`
group B item 4.

**Design:** two arms, `M10_tightprior` and `M20_tightprior`, repeating
the worst two points on the ladder with `b_betatheta_slope_var = 0.5`
(SD 0.71, against the default's 1.41). Same `seed_label = "mladder"` as
the original ladder, so these pair not only against each other but
against the *already-collected* `M10`/`M20` results -- no need to re-run
the default-prior arms.

**Verified before the real run**: refit one dataset under both priors
with a short chain, holding data and seed fixed; the slope row's
posterior spread shrank under the tighter prior (0.187 -> 0.179 SD at
this toy configuration), confirming the override reaches the sampler.
Tier 1/2 unaffected -- the new default preserves current behaviour
exactly.

**Reading it:** if coverage moves *toward* nominal at the tighter
variance, the diagnosis holds and `B_betatheta`'s width is implicated. If
it does not move, or moves the wrong way, the cause is elsewhere and the
prior-variance hypothesis from 13.7 is wrong -- worth stating plainly
either way, per the lesson already recorded about not forcing a real
result into a box it does not fit.

R = 50, for the same reason as 13.2: adequate to detect a move, not to
confirm arrival at nominal.

------------------------------------------------------------------------

## 13.9 Results: tighter prior did not move `beta_theta` (30 July 2026)

**Run:** 100 fits, R = 50, 97 min, 0 failures.

| | M10 (var=2) | M10 (var=0.5) | M20 (var=2) | M20 (var=0.5) |
|---|---|---|---|---|
| `beta_theta` coverage | 0.603 | 0.600 | 0.579 | 0.578 |
| `beta_theta` bias | +0.027 | +0.026 | +0.017 | +0.017 |

**The hypothesis from 13.7 is disproved.** A 4x reduction in prior
variance (SD 1.41 -> 0.71) moved coverage by 0.001-0.003 -- noise
against the 3.1% SE. Every other block (`theta0`, `B0`, `B`, `G`, `q`,
`resid_cor`) was likewise unmoved, confirming the change was isolated
and had no side effects worth reporting. `B_betatheta`'s width is not
what makes the interval overconfident.

**A second candidate was checked and is also dead.** Before writing this
up, checked whether added M-samples might be pseudo-replicated -- i.e.
whether `X_theta` repeats within a site, which would make more M look
like independent information without actually being any. It does not:
`X_theta <- cbind(1, matrix(rnorm(N * ncov_theta), N, ncov_theta))`
(`R/simulateData.R:126`) draws the covariate independently per *sample*,
not per site. So the narrowing interval is not a design-matrix artefact
either.

**Where this leaves item 4.** Two plausible mechanisms are now ruled
out: the prior is not too tight, and the covariate is not
pseudo-replicated. The overconfidence has to come from somewhere in the
likelihood or the sampler's variance computation for `beta_theta` --
plausibly how the latent collection state `w` (or `z`) is aggregated
across a site's M samples, or a genuine numerical issue in the
Polya-Gamma update. That is a C++/sampler-level investigation, not
another prior-tuning experiment, and is a job for whoever wrote
`sample_beta_cpp_TS`/`sample_betatheta_cpp_parallel`, not a further
simulation-study pass. Recorded in `TODO.md` rather than pursued further
here.
