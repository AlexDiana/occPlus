# Simulation-based validation of `occJSDM`

Plan for the "extensive testing on simulated datasets" item under *MEE paper / Doug to dos* in `TODO.Rmd`. Drafted 27 July 2026.

**Status as of 27 July 2026: stages 1 and 3-5 built and committed; only the
full R = 100 run (stage 6) remains.** Sections marked *OPEN* still need a
decision or a measurement.

| Stage | State |
|----|----|
| 1. Fixtures + tier 1 | **done** (`6d9526d`) -- 89 assertions, 7 s |
| 2. Spatial `n` floor | **done** (`41fb736`) -- floor is 31 sites |
| 3. `helper-simstudy.R` | **done** (`6123036`) -- validated at R = 8 |
| 4. Tier 2 canary | **done** (`f63eeeb`) -- ~30 s, `skip_on_cran()` |
| 5. Tier 3 + runner | **done** (`f63eeeb`) -- runner smoke-tested end to end |
| 6. Full R = 100 run | **not started** -- `Rscript dev/simstudy/run_study.R --cores=8`, ~2.4 h on 8 cores |

Building it surfaced four corrections to this plan and two package bugs; both
are recorded in place below rather than only in commit messages.

------------------------------------------------------------------------

## 1. Purpose

Check that `runOccJSDM()` recovers known generating parameters from data produced by `simulateOccJSDMData()`, at the nominal rate, across the configurations the package advertises.

Two distinct goals, which need different designs — worth keeping separate in the writing-up:

- **Regression.** Guard the \~23 bugs already fixed (`TODO.Rmd`, *Fixed bugs*). None of them currently has a test; every one was verified by reading code, or by a one-off manual run. This is the near-term value.
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

- The `sample_ls` regression (`TODO.Rmd` Fixed bugs 22) broke *every* non-spatial fit and survived because nothing exercised that path. A three-line test would have caught it the moment it landed. The vignette would not have — it passes `spatCovariates`.
- It discharges Phase 3 step 11 of the CRAN plan in `AGENTS.md` ("write the testthat suite"). Done as of `6d9526d`: `test-placeholder.R` is gone, replaced by four real test files.
- No package-size or vignette-build cost. Both matter here: `sampleresults.rda` at 62 MB is already a CRAN blocker.

What this gives up is *presentation*. Tests answer "does it pass?", not "here is the calibration evidence". See §9.

------------------------------------------------------------------------

## 3. Tiers

| Tier | Gate | Runtime | Runs on |
|----|----|----|----|
| 1\. Regression + smoke | none | \< 30 s | every `check()`, including CRAN |
| 2\. Recovery canary | `skip_on_cran()` | \~2 min | local + CI |
| 3\. Coverage study | env var `OCCJSDM_SIMSTUDY` | 2.5–3 h | manual / nightly |

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

`TODO.Rmd` group A item 3 is still open: `randinvg()` (`src/jsdm.cpp:86`) draws from R's global RNG inside `samplePGvariables()`'s OpenMP loop. A fixed seed therefore **does not** reproduce across platforms — it is deterministic on stock macOS clang (no OpenMP) and racy on Linux/Windows.

A tier-1 test asserting `expect_equal()` on posterior quantities would pass locally and fail intermittently on CRAN. That is the worst available failure mode. **Tier 1 assertions must be structural**: does it run, are dimensions right, is a block non-`NA`, does a value vary. Numeric recovery belongs in tiers 2–3, with tolerances.

Revisit once group A item 3 is closed.

### 5.2 Coverage assertions are stochastic and will flake if set tight

With R = 100 the SE of an estimated 95% coverage is \~2.2%, so a *correct* implementation lands outside ±1.96 SE about 5% of the time. `expect_gt(cov, 0.90)` would fail intermittently for no reason, and a suite that cries wolf gets ignored.

Mitigations, all three:

- fix a seed per replicate, so a tier-3 run is reproducible;
- set the hard failure threshold well below nominal — fail at **\< 0.85**, which catches every bug in the *Fixed bugs* list by a wide margin (§9);
- record the coverage table with `expect_snapshot_value()` so changes surface as a reviewable diff rather than a boolean.

### 5.3 Identifiability limits what can be checked at all

The latent-factor model is invariant to rotation and sign. `reparamFactorModel()` (`R/jsdmfun.R`) imposes a QR-based constraint, but there is no reason that constraint matches the simulator's parameterisation.

- **Directly comparable** (verified against the code while building `helper-simstudy.R`, 27 July): `B0`, `B`, `G`, `beta_theta`, `theta0`, `p`, `q`, `tau`. Two caveats. `beta_theta` requires passing `collCovariates` to the fit — without it the posterior is `(1 x S)` against a `(ncov_theta+1 x S)` truth. And `n_lattrait` must be pinned to the simulated `gt`, since the default `floor(sqrt(min(S, ncov_psi)))` will not generally match, giving `A`/`C` mismatched dimensions.
- **`sigma_b`: keep, but read with care.** With a tight `InvGamma(10, 1)` prior and only `p x S` residual coefficients the posterior sits at the prior mean (`sqrt(1/9)` = 0.333), so a fixed truth far from that yields near-zero coverage through prior-data conflict, not sampler error. Measured: true 0.5 gave 0/8, true 0.35 gave 8/8. Scenarios set it near the prior mean; SBC is the principled fix.
- **NOT comparable, contrary to the original draft of this section:** `Bs`, `Gs`, `sigma_bs`. The simulator builds the spatial field directly from `sigma_s`/`l_s` and sets `Bst <- matrix(0, S, ps)` (`R/jsdmfun.R:557`), so on the truth side `Bs` is empty at `ds = 0` and `sigma_bs` generates nothing; the fit meanwhile represents the field as sparse-GP basis coefficients over `n_supportpoints` knots. Different parameterisations. Measured for `sigma_bs`: true 0.5 against a posterior mean of \~1.6, 0/8 coverage — a meaningless comparison rather than a bug.
- **Only as identified functions:** residual correlation matrix, `eta`/`psi`, variance partitioning
- **Must NOT be checked element-wise:** `U`, `L`, `A`, `C` — element-wise coverage would fail for reasons that are not bugs
- **Excluded:** `sigma_h` — never sampled (group A item 1), so it would fail by construction. **Also `l_s`** — group A item 2: `sigma_s` is hard-coded to 1 at the `sample_ls()` call site, so the length-scale absorbs the amplitude misspecification and rails at the top of `l_s_grid` for every true value tried. Any coverage figure for `l_s` is meaningless until that is fixed. Exclude explicitly and with a comment, so a future reader does not "fix" the test.

### 5.4 The GP knot count must be pinned, not left to the default

`getDefaultSupportPoints(n) <- max(30, floor(n * 0.2))` (`R/jsdmfun.R:480`) is a constant 30 for any dataset below 150 sites, and crashes below 31 (§10.1, and `TODO.Rmd` group B item 3).

Two consequences for the design, beyond the crash. The knot count would be **uncontrolled but not constant** across cells — identical for `n = 40` and `n = 140`, then suddenly proportional above 150 — so any difference attributed to `n` would be partly an artefact of the spatial approximation changing underneath. And at small `n` it approaches one knot per site, which is a qualitatively different (denser) approximation than at `n = 100`, precisely in the low-information cell where calibration is most at risk.

Set `n_supportpoints` explicitly in every spatial cell so it is a controlled factor.

### 5.5 CRAN limits

Two cores maximum (`_R_CHECK_LIMIT_CORES_`). Tier 3's replicate-level parallelism must be off or capped when not env-gated on.

------------------------------------------------------------------------

## 6. Tier detail

### 6.1 Tier 3 — the study

Base = `two_stage`, traits on, spatial on, `d = 2`, `n = 100`, `S = 10`, `M = 2`, `P = 2`, `K = 3`, `n_supportpoints = 20`.

Every spatial cell sets `n_supportpoints` explicitly (§10.1). Left to the default it is a constant 30 for any `n` below 150, so it would neither scale with `n` nor be a controlled factor — and below 31 sites it crashes outright.

| \# | Cell | Differs from base by | Buys |
|----|----|----|----|
| 1 | base | — | Flagship config; *is* the traits×spatial interaction |
| 2 | spatial isolated | traits off | `Bs` recovery unconfounded by trait terms (**not** `l_s` — see §5.3) |
| 3 | traits isolated | spatial off | `G` (fourth-corner) recovery unconfounded by the GP |
| 4 | `P = 3` | `P = 2 -> 3` | Guards Fixed bugs 2; matches shipped `sampledata` |
| 5 | low information | low `p`/`theta_baseline`, smaller `n` (see §10.1 — must set `n_supportpoints`) | Calibration where data are thin — the realistic eDNA regime |
| 6 | `d` under-fit | truth `d = 4`, fit `d = 2` | The common real-world misspecification |
| 7 | `d` over-fit | truth `d = 2`, fit `d = 4` | Interacts with rotation invariance; least-trusted case |
| 8 | larger `S` | `S = 10 -> 30` | Does calibration hold as species count grows |
| 9 | `occupancy` | no PCR stage | Calibration evidence for the collapsed form |
| 10 | `binary` | no replicates at all | Calibration evidence for pure JSDM |

Design is **one-factor-at-a-time**, not factorial (\~48 cells would be unaffordable). Cost: OFAT cannot see interactions. Cell 1 covers the one interaction judged most likely to matter, since traits and spatial both write into the same linear predictor via `computePsiCoef()`.

Deliberately *excluded* from tier 3, because tier 1 answers them in seconds for free: "does configuration X run at all". Tier 3 hours should buy calibration evidence, which only tier 3 can provide.

Cells 9–10 are **cheaper** (less data per fit: no PCR level, no replication) but also check **fewer parameters** — `occupancy` has no `p`/`q`; `binary` has no `beta_theta`, `theta0`, `p`, `q`. Their rows in the results table will be sparser, not directly comparable to the two_stage cells.

**Metrics** per scenario × parameter block, pooled across species: coverage of the 95% credible interval, bias, RMSE.

**Output:** a tibble (`scenario`, `block`, `coverage`, `bias`, `rmse`, `n_par`, `R`) written to `Sys.getenv("OCCJSDM_SIMSTUDY_OUT", tempdir())`, intended destination `dev/simstudy/results/`.

### 6.2 Tier 1 — ships to CRAN

One shared fixture fit built in `helper-fixtures.R` and reused: `sampledata`, `nburn = 50, niter = 50, nchain = 2` (**measured 5.6 s**). Keeps the tier well under 30 s.

`test-smoke-configs.R` — every configuration completes: `two_stage` / `occupancy` / `binary`; with and without `spatCovariates`; with and without traits; `d = 0` and `d = 2`; `P = 1` and `P = 3`. *This is the tier that would have caught the `sample_ls` regression.*

`test-regression-bugs.R` — one test per fixed bug, named so a failure points at its `TODO.Rmd` entry:

| Guards           | Assertion                                           |
|------------------|-----------------------------------------------------|
| Fixed bugs 22    | non-spatial fit runs                                |
| Fixed bugs 8     | `tau_output` not all `NA`                           |
| Fixed bugs 21    | `plotFPTPStage2Rates()` differs across primers      |
| Fixed bugs 23    | `computeSpeciesDetected()` indices span both chains |
| Fixed bugs 9     | WAIC counter guard does not trip                    |
| Fixed bugs 1, 11 | row-shuffle invariance                              |

Row-shuffle invariance must be asserted on the **prepared design matrices** (`X_psi`/`X_theta` row alignment against site/sample indices), which is deterministic — not on the posterior, which would flake under §5.1.

`test-api-contracts.R` — unknown `primerName` errors with the available list; unknown covariate name errors informatively; `idx_species` subsetting returns the right rows; `predictNewSites()` honours its documented `NULL` defaults.

### 6.3 Tier 2 — canary

R = 3–5 replicates at the base cell. Not a calibration check: a smoke signal that recovery has not grossly broken between full tier-3 runs (e.g. correlation between estimated and true `B` above a loose floor).

**RESOLVED: tier 2 fails.** Thresholds were set from measurement rather than
taste, which is what makes failing safe. Three independent seed sets gave
overall coverage 0.89/0.89/0.88, lowest block 0.78/0.76/0.76, and
`cor(post_mean, truth)` 0.62/0.56/0.68; the floors are 0.70, 0.40 and 0.30. A
pass is therefore not luck and a failure is not noise.

------------------------------------------------------------------------

## 7. SBC seam

Coverage now, simulation-based calibration later, without restructuring. Two swappable functions in `helper-simstudy.R`:

``` r
draw_truth(scenario)        # fixed values now; prior draws for SBC
summarise_fit(fit, truth)   # CI-contains indicator now; rank statistic for SBC
```

Everything expensive — grid, simulate/fit loop, parallelism, artifact writing — is shared.

Adding SBC later requires:

1.  A prior-draw wrapper. All priors in the sampler are **proper** (checked): `B0`, `L` \~ N(0,1) (`sample_BBsL()`, `R/jsdmfun.R`); `B` trait-predicted mean with variance `sigma_b^2`; `sigma_b`, `sigma_bs` \~ sqrt-InvGamma(10,1); `tau` \~ InvGamma(5,5) (`R/runOccJSDM.R`); `p`, `q`, `theta0` \~ Beta; `l_s` \~ Gamma(1,1) over a 10-point grid. So SBC is feasible.
2.  Documenting what the sampler *actually* assumes. Not free: the audit found the *documented* `prior_beta_psi`/`prior_beta_psi_sd` were a complete no-op, with the live path hard-coding something different (Fixed bugs 14).
3.  Thinning to near-independent draws, which raises per-replicate cost.

Note `sigma_h` would show a flagrantly non-uniform rank histogram under SBC. That is the method working (group A item 1), not a new bug.

------------------------------------------------------------------------

## 8. Runtime budget

Measured on `sampledata` (100 sites, 10 species, 2 chains): **0.039 s per iteration**, so `nburn = 1000 + niter = 1000` is \~84 s per fit.

|                   | Serial  | 8 cores       |
|-------------------|---------|---------------|
| One cell, R = 100 | \~2.3 h | \~17 min      |
| 10 cells          | \~23 h  | **\~2.5–3 h** |

Cells 9–10 are cheaper than the estimate; cell 8 (`S = 30`) is dearer.

Replicates are embarrassingly parallel at the *process* level — independent datasets and fits — so this needs no package changes, unlike the in-package chain parallelisation in `TODO.Rmd` group D item 1.

**Do not buy replicates by shortening chains.** Each interval endpoint is a tail quantile; with 500 total draws only \~12 land below the 2.5% bound, so the endpoints are noisy and that noise feeds into coverage as bias, not just variance. Set chain length by an ESS target (\>= 400 on monitored parameters, which `returnConvergenceDiagnostics()` already reports) and let R take what is left.

------------------------------------------------------------------------

## 9. What R = 100 does and does not support

SE of estimated coverage is `sqrt(0.95 * 0.05 / R)`. Cost is linear in R; precision improves only as sqrt(R) — halving the error bar costs 4x the compute.

| R       | Coverage SE | Undercoverage reliably flagged |
|---------|-------------|--------------------------------|
| 50      | 3.1%        | below \~0.890                  |
| **100** | **2.2%**    | **below \~0.907**              |
| 200     | 1.5%        | below \~0.920                  |
| 400     | 1.1%        | below \~0.929                  |

R = 100 is chosen for the **regression** goal. Every bug in the *Fixed bugs* list would drive coverage to 0.5–0.8, not 0.93, so R = 100 catches them with room to spare.

It is **not** enough to assert "coverage is nominal" in the paper — that claim needs R = 200–500, because it asserts the *absence* of a small deviation rather than detecting a large one. The runner takes `R` as an argument for exactly this reason; do a final publication-grade run at higher R if the paper makes that claim.

Also: pooling coverage across species within a block buys precision, but those indicators share one simulated dataset and one MCMC run, so effective N sits between R and R x S, nearer R. **Do not quote an SE computed as if there were R x S independent trials.**

------------------------------------------------------------------------

## 10. Open items

0.  *OPEN, and it degrades cells 1 and 2.* **`sigma_s` is never sampled** (`TODO.Rmd` group A item 2, found 27 July while writing tier 1). The spatial field's amplitude is hard-coded to 1 at the `sample_ls()` call site, so the length-scale absorbs the misspecification and rails at the top of `l_s_grid` regardless of the true value — measured at true `l_s` = 0.30, 0.15 and 0.08, across seeds.

    Consequences here: `l_s` is dropped from the checked parameters (§5.3), and the spatial cells still test `Bs` and the overall spatial contribution but cannot speak to range recovery. **Not blocking** — the cells remain worth running — but a coverage table produced before this is fixed must not be read as evidence about the GP. Re-run cells 1 and 2 afterwards.

1.  ~~*OPEN, blocking tier 3.***Minimum `n` for spatial cells.**~~ **RESOLVED 27 July 2026, and no longer blocking.** Measured: the floor is **31 unique locations** with default settings. `n <= 29` fails with `cannot take a sample larger than the population`, `n = 30` with `number of cluster centres must lie between 1 and nrow(x)`, `n >= 31` runs.

    Cause: `getDefaultSupportPoints(n) <- max(30, floor(n * 0.2))` (`R/jsdmfun.R:480`) feeds `kmeans(X_s, centers = ps)`, and the hard-coded 30 wins for every dataset under 150 sites, so the knot count never scales down. Filed as a bug — **`TODO.Rmd` group B item 3** — since the default is also statistically wrong well before it crashes (30 knots for 31 sites is \~1 knot per site, defeating the point of a sparse GP).

    **Consequence for this plan:** cell 5 must set `listParams$n_supportpoints` explicitly rather than relying on the default. Verified working at `n = 20` with 6 knots. Pinning the knot count is the right call for a study cell anyway — it stops the spatial approximation silently changing with `n` across scenarios, which would otherwise confound the low-information comparison. Apply the same reasoning to **every** spatial cell (1, 2, 4–10), not just cell 5: set `n_supportpoints` per scenario so it is a controlled factor rather than a function of `n`.

2.  ~~*OPEN.* Tier 2 failing vs advisory (§6.2).~~ **RESOLVED: it fails**, on measured thresholds -- see §6.3.

3.  *OPEN, and the main empirical question stage 6 should settle.*
    **`beta_theta` and `resid_cor` sit below nominal.** Measured at R = 5
    across three independent seed sets: `beta_theta` 0.78 / 0.80 / 0.78, and
    `resid_cor` 0.79 / 0.76 / 0.76. Everything else landed 0.86-1.00.

    That consistency across independent seed sets makes chance unlikely, but
    R = 5 is far too small to call it -- the effective SE is ~8%, so 0.78 is
    only about 2 SE below nominal. At R = 100 the SE is 2.2% and the question
    is decisive either way.

    Both are held to a lower floor in tiers 2 and 3 rather than excluded,
    precisely so that settling this stays a purpose of the run rather than
    something the thresholds quietly paper over. If the undercoverage is real,
    `beta_theta` points at the Stage 1 collection model and `resid_cor` at the
    factor-loading posterior; neither has an obvious suspect yet.

4.  *Deferred.* Presentation of tier-3 results for the paper. Build the suite first; the summary object then feeds either a short pkgdown article or the manuscript directly. Nothing here forecloses that.

5.  *Deferred.* SBC (§7).

------------------------------------------------------------------------

## 11. Build order

1.  `helper-fixtures.R` + tier 1. Ships, pays immediately, blocks nothing.
2.  ~~Resolve open item 1 (spatial `n` floor).~~ **Done** — floor is 31 sites; every spatial cell now pins `n_supportpoints` (§10.1).
3.  `helper-simstudy.R` with the two seams (§7).
4.  Tier 2.
5.  Tier 3 + `run_study.R`.
6.  **Full R = 100 run; inspect the table.** The only stage left.

    ```
    Rscript dev/simstudy/run_study.R --cores=8
    ```

    ~2.4 h on 8 cores (~19 h serial), extrapolated from 10 s/replicate
    measured at reduced size. The runner has been smoke-tested end to end --
    argument parsing, PSOCK cluster, worker setup, both artefacts, failure
    counting -- but never run at full size.

    Two things to weigh first. **Cells 1 and 2 will be re-run** if `sigma_s`
    (open item 0) is fixed afterwards, so if that fix is near-term, waiting
    saves the duplication. And **one cell settles the open empirical
    question** at a tenth of the cost:

    ```
    Rscript dev/simstudy/run_study.R --cores=8 --scenarios=base
    ```

    ~15 min on 8 cores, and enough to say whether the `beta_theta` /
    `resid_cor` undercoverage is real (§10.3).
