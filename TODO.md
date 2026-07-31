---
title: "TODO"
output: html_document
---

```{=html}
<!-- DO NOT hard-wrap this file. Write each paragraph as a single line.

     There is deliberately no `editor_options: markdown: wrap:` block, and
     `occJSDM.Rproj` sets `MarkdownWrap: None`. Doug edits these docs in
     RStudio's Visual mode, which rewrites the file to canonical form on
     every save. While the wrap was set to 72 columns, anything written by
     hand or by an agent got reflowed on the next open, producing a diff
     nobody authored, over and over.

     Measured 30 July 2026: RStudio's canonicalisation cannot be
     reproduced outside RStudio. Feeding RStudio's own output back through
     pandoc changes 234 of 768 lines, across several flag combinations. So
     matching it by hand or by hook is not achievable, and the fix is to
     remove the line-breaking step entirely: with no wrapping there is no
     algorithm left to mismatch.

     Three related rules:
     - NO EM-DASHES. Use a plain double hyphen. Em-dashes were the single
       largest driver of reflow churn under the old setting, because
       pandoc writes them 2 characters wider and every later wrap point
       shifted.
     - NO PIPE TABLES in this file. Use bullet lists instead. RStudio
       recomputes a table's separator row (|----|----|) from the cell
       contents and the available width, and it switches between
       padding to content width and redistributing across the line at a
       width threshold. So any edit to any cell can rewrite the
       separators, and 873f419 is an example. Lists have no
       width-derived formatting and are stable. PLAN.md keeps its
       tables, and accepts the churn, because the data there earns them.
     - Keep inline `code spans` short. Under the old 72-column wrap, a
       span landing on the boundary could be mangled: in 2242b34 that
       turned `* 2` into `- 2` in the logDetKuu note, silently asserting
       something false. Less dangerous now, but still good practice.

     What this file is FOR: a list of to-dos with enough context to act
     on, addressed to Alex. It is not an investigation log. Evidence and
     superseded reasoning go to AGENTS.md; completed work goes to the
     Fixed bugs section below, which is the record. If an item is
     fixed, write the Fixed bugs entry and DELETE the item; do not leave
     a struck-through stub. -->
```

# **v0.1.0-beta** **Public release**

## A. Review Claude's changes (Alex)

Every code change Claude made, newest first. None has had human review beyond Doug asking for it. All are recoverable from git; revert or rework freely. Each says what to check. Full detail for each is in *Fixed bugs*, which is the record; this section is the queue.

1.  **`main` would not load at all, and two fixes were needed to restore it.** 31 July, `src/Makevars`, `src/Makevars.win`, `DESCRIPTION`, `NAMESPACE`, `R/occJSDM-package.R`, `R/runOccJSDM.R`. After `8f9f315` the test suite went from 167 passing to 25 passing, 3 failures and 33 errors, because the built `.so` had an undefined `RcppParallel::tbbParallelFor` and failed to load, taking every function in the package with it.

    **Change 1, the linking.** `RcppParallel` was in `LinkingTo`, which makes the headers visible but does not link the libraries. Added the documented `RcppParallel::RcppParallelLibs()` call to both `Makevars` files, **and** added `RcppParallel` to `Imports` with an `importFrom`. **Both halves are needed**: the Makevars line alone fixes the missing symbol but the load still fails, because on macOS `libtbb` is referenced through `@rpath` and `devtools` copies the `.so` to a temp directory. Importing the package makes its namespace load first and set up the TBB paths. The reasoning is in a comment in `src/Makevars` so the next person does not stop at the first fix and conclude it did not work.

    **To check:** whether you would rather solve this a different way, for instance an explicit `-Wl,-rpath`. The `Imports` addition is a new hard dependency in `DESCRIPTION`, which is your call.

    **Change 2, your debugging scaffold at `R/runOccJSDM.R:404`.** It had live assignments overwriting every argument with values referencing `occ_data_effort`, a dataset not in the package, plus a live `summarisedLatentPresences`. That block was fully commented before `8f9f315`. **Re-commented rather than deleted**, since you evidently use it, with a note saying it must stay commented and what happened when it did not.

    **To check:** nothing, unless you want it gone entirely.

    **Two things left alone because they are your calls.** `verbose` in `computeNewOutputs()` defaults to `T`, so `predictNewSites()` still prints one line per species unless a caller opts out; it is suppressible now but the default is unchanged. And `sample_z_cpp_parallel` is exported but never called: `runOccJSDM` still uses `sample_z_cpp` at `:1101`, so only the `w` sampler is wired in.

    ALEX TO REVIEW

## **B. Inference-affecting bugs (wrong numbers, silently) (Alex)**

1.  **`sample_ls()` scores the wrong density, so the GP length-scale is never recovered.** `R/jsdmfun.R:1054`. Under the SoR approximation the fitted field `SE = Ks(l_s) %*% Bs` is a deterministic function of `Bs` and `l_s`, but `sample_ls()` treats `SE` as a GP draw and scores it under `N(0, sigma_s^2 K(l_s))` with `SE` held fixed. That is not the conditional posterior, and it is self-defeating: `SE` was already smoothed at the current `l_s`, so it scores better under ever-smoother covariances.

    **Measured:** `idx_ls` rails at the top of `l_s_grid` for every true `l_s` tried (0.074, 0.171, 0.300), with real spatial signal present. The profiled log-likelihood rises monotonically from -376 at `l_s = 0.01` to -154 at `l_s = 0.30`.

    **Impact:** biases the spatial term of every `useSpatField = TRUE` fit and makes `l_s` uninterpretable. Widening the grid would only move the rail.

    **Needs a derivation, not a code tweak.** Either recompute `eta` with `Ks(l_s*)` and score the observation likelihood under the proposal, or integrate out `Bs` and use the SoR marginal likelihood. **Do not attempt blind:** a subtly wrong conditional would be worse than the present state, which at least fails visibly. A cheap interim option is to revert to a documented, user-settable fixed `l_s`.

    Four candidate causes were investigated and ruled out (missing amplitude, wrong amplitude passed, `logDetKuu`, weak data). Detail in AGENTS.

    ALEX TO CHECK

2.  **`reparamFactorModel()` breaks residual covariance = `t(L) %*% L`, inflating reported species correlations.** `R/jsdmfun.R:48`. The rotation preserves `U %*% L` (verified to 4e-16) so the linear predictor is untouched, but it moves scale out of `U` into `L`, and `returnResidualCorrelationMatrix()` computes `cov2cor(t(L) %*% L)` from the reparameterised `L`. Measured `Var(U)` afterwards is `diag(0.23, 2.01)`, not the identity.

    **Measured:** correlations move by up to 0.612, consistently toward the extremes. Across the grid `resid_cor` covers at 0.74-0.77 in nine of ten scenarios.

    **Impact:** `returnResidualCorrelationMatrix()` and `plotResidualCorrelationMatrix()` overstate co-occurrence. This is the headline JSDM output.

    **Fix, two options.** (a) Rotate by `Q` alone, dropping the `diag(diag(R))` scaling, so both the identifiability constraint and the covariance identity hold. (b) Keep the scaling and compute the correlation as `t(L) %*% Var(U) %*% L`. (a) is simpler and preserves the output contract.

    **Your note was that this is a non-issue for logistic models since only the correlation is recoverable.** The premise is right but the measured numbers above are already correlations: `cov2cor()` does not absorb the change, because `diag(1/diag(R))` rescales per *factor* while `cov2cor()` normalises per *species*. If it is nonetheless a non-issue, then the simulation study is measuring the wrong thing and that needs identifying before the `resid_cor` results are quoted. Full exchange in AGENTS.

    ALEX TO MAKE A DECISION

3.  **`beta_theta` intervals are overconfident, and it gets worse with more data.** Coverage 0.77 at the production `M = 2`, falling monotonically to 0.58 at `M = 20`, while bias stays small and flat. Shrinking intervals around a bias that is not shrinking is the signature of a real defect being exposed by more information, not fixed by it.

    **Narrowed 31 July: the defect is in the *slopes*, not in `beta_theta` as a block.** Refitting `base` with `ncov_theta = 0`, so only the intercept row remains, gives `beta_theta` coverage of **0.968** (SE 0.013, R = 200), i.e. nominal, against 0.763 with the slopes present (`PLAN.md` 15.5, 15.6). Whatever is wrong is specific to the covariate columns.

    **Four candidate causes now ruled out**, each by measurement: Stage 1 under-identification (more data makes it worse, not better); the slope prior's width (tightening it 20-fold at `M = 2` moves coverage the wrong way); pseudo-replication in `X_theta` (it is drawn per sample, not per site); and the intercept path (nominal once the slopes are gone).

    **So the cause is in whatever handles the covariate columns in the Polya-Gamma update**, in `sample_beta_cpp_TS`/`sample_betatheta_cpp_parallel`. Note this is the same step Alex's `microbenchmark()` profiling identified as the slowest in the sampler, so the calibration problem and the performance bottleneck sit in the same code. This needs someone who knows it; it is not another prior experiment. Evidence in `PLAN.md` 13, 14 and 15.5.

    ALEX TO INVESTIGATE THE SAMPLER

4.  **Decide `b_betatheta`'s slope prior variance. It trades `B0` bias against `beta_theta` coverage.** Measured at `M = 2`, paired on identical truths, varying only that variance:

    - variance 2, your current default: `B0` bias -0.160, `beta_theta` coverage 0.747
    - variance 0.5: `B0` bias -0.106, `beta_theta` coverage 0.707
    - variance 0.1: `B0` bias -0.044, `beta_theta` coverage 0.653

    **This identified the cause of item 6 below.** `42198d9` widened `B_betatheta` from `diag(1)` to `diag(2)`, which is exactly when `B0`'s bias doubled. Turning it back down moves the bias back, monotonically.

    **But it is a trade, not a fix:** tightening helps `B0` and hurts `beta_theta`. Items 4 and 6 pull opposite ways on one knob. Not known: whether an intermediate value beats both endpoints, whether the trade holds at `M > 2`, and whether fixing item 4 at its source would dissolve it entirely.

    ALEX TO DECIDE THE VALUE (or that item 4 supersedes this)

5.  **`B0`'s bias doubled, and coverage does not show it.** Between the pre- and post-fix runs on identical data, nine of ten scenarios moved more negative: overall -0.135 to -0.228. Coverage held at 0.943 throughout, because the intervals are wide enough to absorb the shift, so this is invisible in the headline table and shows only in the bias column.

    Cause is now most likely the `B_betatheta` widening; see item 5, where the decision lives. `B0` is a headline quantity for a JSDM, so this wants settling before the paper reports species intercepts.

ALEX RESPONSE: WE COULD ADD A SIMULATION STUDY ON THE ONE-STAGE MODEL ONLY, THAT WOULD REVEALE WHETHER THERE IS ANY ISSUE IN B0 (since there is no beta_theta). If there is no issue, we could delete this point and be sure that the issue is only beta_theta.

    **Ran 31 July at R = 50 and confirmed at R = 200 (`PLAN.md` 15.4, 15.6). Alex's hypothesis holds in its main claim and fails in its strong form.** The designed arm was `base` with `ncov_theta = 0`, keeping the two-stage machinery, the latent `w`/`z` and `p`/`q`/`theta0`, removing only the `beta_theta` slopes. `B0` bias:

    - `base`, slopes present: -0.2078 (SE 0.0307), 6.8 SE from zero
    - `nocollcov`, slopes removed: **-0.0633** (SE 0.0192), **3.3 SE from zero**
    - `binary`, no `beta_theta` at all: +0.0122 (SE 0.0119), 1.0 SE from zero

    Removing the slopes removes **70%** of the bias, so fixing item 4 will recover most of `B0`. **But the residual is real**: 1.5 SE at R = 50, 3.3 SE at R = 200. `binary` at 1.0 SE is what no bias looks like by comparison.

    **So this item stays open.** Alex proposed deleting it if the one-stage test showed no issue; the test shows a smaller issue, not none. The cause is now split: most is downstream of item 4, a measurable part is not, and `B0` coverage sits at 0.94-0.96 throughout so it will never surface there.

    The R = 200 run also reproduced the R = 50 replicates bit-for-bit, which is how we know the two are the same experiment rather than two similar ones.

    ALEX: WHAT REMAINS AFTER ITEM 4 IS FIXED IS SMALL BUT NOT ZERO

## **C. Crashes, unreachable code paths, and API bugs (Alex)**

1.  **`thinOutput()`: two of the three original defects remain.** The crash is gone and `thin` is honoured, but (a) the 2-D branch thins **by row** (`R/output.R:57`), silently dropping *sites* from the `psi_output`/`w_output`/`theta_output` posterior-mean matrices, whose rows are sites and not iterations; and (b) the scalar `WAIC` falls through to `print("Dimension not recognised")` and becomes `NULL`.

    Also: it is no longer in `NAMESPACE` but `man/thinOutput.Rd` still ships a `\usage` block, so the manual documents a function users cannot call. Either re-export or drop the `.Rd`.

    Note the CRAN plan's option for shrinking `sampleresults.rda` depends on this working, so it is not simply deletable.

    ALEX TO DISABLE IT FOR NOW, BUT LET'S KEEP IT THERE

2.  **Assorted smaller items.**

    (a) `createDataIdx()` is called with `maxP` (`R/runOccJSDM.R:638`) for `model = "occupancy"` too, where `maxP` was never assigned. It survives only because lazy evaluation never forces the promise; pass `NULL` explicitly.

    (b) `d <- get_param(listParams, "n_factors")` defaults to 0, and the cap at `:716` uses `ncol(OTU)`, which is `NULL` for a single-species vector: `if (d > NULL)` errors.

    (c) `reparamFactorModel(A_output, C_output)` (`:1273`) assumes `qr.Q()` is square; when `ncov_psi < gt` it is not, and `diag(diag(R_current), nrow = d)` recycles.

    (d) The spatial-covariate numeric check at `:540` runs before the "names present in `data$info`" check at `:544`, so a mistyped name gives `undefined columns selected` rather than the intended message.

    (e) `computeSpeciesDetected()`'s roxygen documents the removed Beta-approximation signature instead of its actual arguments.

    ALEX WILL REVIEW THE ABOVE

3.  **`computeNewOutputs()` prints to stdout on every call and cannot be silenced.** `src/jsdm.cpp` runs `Rcout << "Computing species ..."` inside the species loop unconditionally, so every `predictNewSites()` call prints one line per species. `suppressMessages()` does not catch it, because `Rcout` is stdout rather than R's condition system.

    It pollutes any script, vignette chunk or app that predicts in a loop, and unconditional console output from a compute function is the kind of thing CRAN reviewers pick up. **Fix:** a `verbose` argument defaulting to `FALSE`, threaded from `predictNewSites()`. Progress reporting is useful on a slow prediction, so making it opt-in beats deleting it.

    ALEX TO DECIDE AND FIX (touches `src/jsdm.cpp`)

ALEX RESPONSE: Added a verbose argument

## **D. Dead and broken internal code (Alex)**

Ten dead functions were moved to `deprecated/` on 30 July (section A item 1). What remains is below. A systematic scan found 38 dead R functions plus 8 unused `RcppExports` wrappers in total, so most of this is still outstanding.

1.  **The remaining dead functions, once the above is settled.** Roughly 24 more in `R/jsdmfun.R` and `R/mcmcfun.R`, plus `computeMinESS()` in `R/diagnostics.R`, plus 8 unused wrappers in `RcppExports.R`.

    The wrappers need different handling: **do not edit `RcppExports.R`**, it is generated. Remove the `// [[Rcpp::export]]` tag in the C++ and re-run `Rcpp::compileAttributes()`.

    Two things the scan flags that must **not** be deleted: `.onLoad()` (zero callers because R itself calls it; removing it would drop the `mc.cores` cap set for CRAN compliance) and `thinOutput()` (uncallable today, but group C item 1 and the CRAN plan depend on it).

    CLAUDE TO CONTINUE AFTER THE ITEM 1 DECISION

2.  **`globalVariables()` for the data-masked column names.** The `R CMD check` undefined-globals NOTE is down from 84 symbols to 65 as dead code has been removed. What will remain is `dplyr`/`ggplot2` NSE references (`x`, `y`, `Species`, `Min`, `2.5%` and so on), which are false positives and want one `utils::globalVariables()` call in `R/occJSDM-package.R`.

    **Do this last.** Every dead function removed shrinks the list, so enumerating it earlier means writing entries for code about to be deleted. There is no `globalVariables()` anywhere yet, so this sets the convention.

    **Done when** `devtools::check()` reports no NOTE under "checking R code for possible problems", not merely a shorter one.

    CLAUDE TO DO AFTER ITEMS 1 AND 2

3.  **`sample_rnb()` cannot run as written** (`R/jsdmfun.R`). Groundwork for the count-data item, not yet called from anywhere, but it has a scoping bug that will bite when wired up: `r_current <- rnb[s]` reads `rnb` inside the `sapply()` whose result is being assigned to `rnb`, so lookup falls through to the namespace and fails. The current size vector needs to come in as an argument.

    Two more to settle while there: `tune_sd = 5` is a random-walk SD on the *log* scale, so proposals land a factor of `exp(+/-10)` away and acceptance will be near zero (0.1 to 1 is the usual starting range); and the prior terms are stubbed to `0` with the intended `dgamma()` commented out, referencing `prior_shape`/`prior_rate`, which are not defined anywhere. The Metropolis step itself looks right: the `log(r_star) - log(r_current)` Jacobian is the correct correction for a log-scale random walk under a flat prior on `r`.

    ALEX's WORK IN PROGRESS FOR THE COUNTS

## **E. Draft of beta version listserv announcement (Doug)**

1.  Listserv announcement (beta release), drafted July 20 2026:

    > Subject: New R package (beta) - occJSDM, a combined occupancy and joint species distribution model
    >
    > Hi all,
    >
    > Announcing **occJSDM**, an R package for combining occupancy and joint species distribution modelling (<https://github.com/AlexDiana/occJSDM>).
    >
    > occJSDM extends the occPlus two-stage eDNA occupancy model of Ji et al. (2025, *Ecology Letters*, <doi:10.1111/ele.70302>) by adding a JSDM layer. Unusually for an occupancy model, false positives are estimated explicitly at both the field and lab stages and separately for each species and each primer.
    >
    > Note this is **beta software**. We are validating it against simulated data. Environmental effects on occupancy, trait-by-environment interactions, and Stage 2 false-positive rates recover well, but we would treat pairwise species correlations, the spatial term, and collection-covariate effects with caution for now.
    >
    > Highlights:
    >
    > - Occupancy modelling: Accounts for both false-negative and false-positive error at two stages (field and lab), per species. Stage 1: estimates species eDNA collection probability in the field, given true eDNA presence at the site, and contamination probability, given true eDNA absence at the site. Stage 2: estimates species eDNA detection probability in the lab (i.e. successful DNA extraction, PCR, and sequencing), given successful eDNA collection in Stage 1, and contamination probability, given eDNA non-collection in Stage 1. In datasets where multiple primers have been used, each species' detection probability is estimated per primer (allowing one to compare each primer's efficiency for each species), while species occupancies are estimated using information across all primers. Both environmental and detection covariates are supported.
    > - JSDM: Integrates the occupancy model with a JSDM: species fit jointly with nonlinear response curves and latent-factor residual correlations. The JSDM optionally supports species traits shaping occupancy responses (trait x env interactions, aka 'fourth-corner analyses') and spatial autocorrelation (GP kernel) across sites. Occupancies predicted at unsampled sites.
    > - occJSDM not only fits a two-stage occupancy model (both field and PCR replicates required), but if given simpler study designs, can collapse to a classical occupancy model (field replicates only) or to a pure JSDM (no replicates).
    > - MCMC fitting with diagnostics, variance partitioning, ordination, and pairwise residual correlation outputs built in.
    > - occJSDM leverages the taxonomic breadth of eDNA datasets by using ordination (each site's position on the latent axes, and each species' loadings on those axes) to predict species occupancies. Thus, each species' predicted occupancy at a site is informed by the estimated occupancies of the other species at that site, thereby using co-occurrence structure. We also allow species to borrow strength from other species sharing similar traits, including inferred traits, in contrast to the classical approach of having rare species borrow strength from abundant species, as is used in multi-species occupancy models.
    >
    > Vignettes included for data simulation and model fitting/interpretation.
    >
    > Feedback, feature requests, and bug reports are very welcome.

# **MEE paper**

## A. Alex to dos

1.  Trait matrix currently not allowing for categorical variables

2.  Design better model selection criterion (the one currently implemented, which is the same as HMSC, tend to overfit).

3.  ability to analyse count data

4.  scenario for source-sink inference (sites where env covariate coeffcients are negative but spatial covariate coefficients are positive), using an explicit source-sink simulation

5.  remove effect of space on environmental covariates. remove the effect of unobserved environmental on observed environmental covariates and space. thus, adding factors (unobserved env covariates) doesn't change the effect of observed env covariates

6.  **site** variance partitioning to complement the **species** variation partitioning (see Leibold et al., Cai et al.). Consider calling it **variation** partitioning.

7.  **Performance of `runOccJSDM()`**

ALEX NOTE: After manually comparing each MCMC step with microbenchmark(), the slowest step is definitely the sample_betatheta_cpp_parallel. There are few things to note, first the parallelisation does not seem to achieve much speed up, and even in its current state, it is not macos compatible since it uses openMP (rather than RcppParallel). Moreover, the slowest step of the sampler seems to be the sample_Omega_cpp, which samples a very large number of PG variables (N x S). We could consider an alternative PG sampler to speed up the computation.

```         
**Profiled by Alex, 31 July 2026, which answers the "nothing here has been profiled" caveat this list used to carry.** Comparing each MCMC step with `microbenchmark()`:

- **`sample_betatheta_cpp_parallel()` is the slowest step**, decisively.
- **The parallelisation achieves little speedup**, and is inert on macOS because it uses OpenMP rather than RcppParallel. This independently corroborates the measurement below: `SHLIB_OPENMP_CXXFLAGS` is empty on Doug's machine, so every `#pragma omp` compiles to a no-op and the "parallel" samplers run serially. Alex has since moved two samplers to RcppParallel in `8f9f315` for this reason.
- **Within that step, `sample_Omega_cpp()` dominates**: it draws `N x S` Polya-Gamma variables per iteration.
- **Alex's suggestion: consider an alternative Polya-Gamma sampler.** That is the lever with the best expected return, and it is a different kind of work from the parallelisation items below, which redistribute the same cost rather than reducing it.

**This reorders the list.** Items A and B parallelise around a step whose cost is dominated by PG sampling; making the PG draw cheaper would benefit every configuration, including single-core and macOS, where the parallelisation currently does nothing. Worth settling the PG question before investing in either.

The items below are ordered by expected speedup per unit of effort as originally written.

**Before profiling, check whether OpenMP is even on.** Measured 29 July 2026 on the macOS development machine: `R CMD config SHLIB_OPENMP_CXXFLAGS` is *empty*, so the `$(SHLIB_OPENMP_CXXFLAGS)` in `src/Makevars` expands to nothing, every `#pragma omp` compiles to a no-op, and `nm -u src/occJSDM.so | grep -c '__kmpc\|_GOMP'` returns 0. The package is running single-threaded here despite `libomp.dylib` appearing in `otool -L` (pulled in transitively by R's own libraries, not by us).

Consequences: the parallel sections are currently dead weight on this machine, the thread-safety bug of Fixed bugs 26 could never have manifested locally, and any timing measured here says nothing about a Linux build where OpenMP *is* active. Worth confirming what Alex's machine and CRAN's check farm do before investing in items 1 and 2 below -- the payoff differs completely between the two cases.

A.  **Parallelise over chains -- but use a PSOCK cluster, not `mclapply()`.** The `for (chain in 1:nchain)` loop (`R/runOccJSDM.R:895`) is serial and embarrassingly parallel; each chain touches only its own `*_output_chain` arrays, so running the chains in separate *processes* gives close to an `nchain`-fold speedup and sidesteps the RNG thread-safety problem in bug A.2 entirely (each process has its own RNG state). Portability constraints, though:

    - **`parallel::mclapply()` is not an option.** It is fork-based, and R ships a Windows stub whose body is `if (cores > 1L) stop("'mc.cores' > 1 is not supported on Windows")` -- a hard error, not a fallback. That would break the package on Alex's machine and on CRAN's Windows check. Even on Linux/macOS, forking a session with a live threaded-BLAS or OpenMP pool (OpenBLAS on most Linux distros, Accelerate on macOS) is a well-known deadlock source -- and doubly so if A.2 is fixed by *keeping* OpenMP.
    - **`parallel::makeCluster()` (PSOCK) + `parLapply()` works on all three platforms**, and `parallel::clusterSetRNGStream()` gives reproducible, independent L'Ecuyer streams per chain, which is strictly better than the current situation. The cost is that workers are fresh R sessions (`clusterEvalQ(cl, library(occJSDM))`, plus serialising the data out and the per-chain arrays back) -- negligible against a multi-minute MCMC, though the return trip is not free given how large `Bs_output`/`U_output` are (see D.7).
    - **Make it opt-in**: a `cores = 1L` argument, serial by default, so existing behaviour is unchanged and CRAN's two-core limit for examples/tests/vignettes is respected.

    Structurally this needs the chain body extracted into a function returning its own `*_output_chain` arrays -- mechanical, since the loop already writes only to per-chain objects. The one piece of genuinely shared state is the WAIC accumulator, which currently streams across chains sequentially via the single `currentWAICiter` counter introduced by the fix in Fixed bugs 9; parallelising forces one accumulator per chain plus a merge (`mean_lik` is a plain mean; `M2` merges with the standard parallel-variance formula), and the `if (numIters != (currentWAICiter - 1)) stop(...)` guard at `R/runOccJSDM.R:1251` will need to be restated in terms of the summed per-chain counts. The `z_output_mean` / `psi_output_mean` / `w_output_mean` / `theta_output_mean` running means merge by simple addition of the per-chain partials.

B.  **Drop the `options(mc.cores = ...)` call in `.onLoad()`.** `R/zzz.R` sets a global option at load time, which changes the user's session state and affects every other package that reads `mc.cores` -- CRAN policy is explicitly against this. It should be replaced by the `cores` argument in D.1. (Also `parallel::detectCores()` can return `NA`, which `min(2L, NA)` happily propagates.)

C.  **`computePsiCoef()` is called three times per iteration.** `R/jsdmfun.R:984`, `:1052` and `:1074`. Each call recomputes `t(computeBtcoef(...))`, `X %*% B` (`n x p` by `p x S`), `KsBproduct()` and `H %*% L`. The second call exists only to refresh `XB`/`SE` before `sample_U_cpp()`, but neither depends on `U`, so they are unchanged from the first call; and the third differs from the second only in the updated `U`. Call it once at the top, reuse `XB`/`SE`, then recompute just `UL <- U %*% L` and `eta <- XB + SE + UL` at the end -- roughly a two-thirds saving on this block.

D.  **Precompute the constant parts of the `c_imk` update.** `R/runOccJSDM.R:1104` recomputes `y > 0` (an `N3 x S` logical) every iteration although `y` never changes; and `w_all <- w[idx_w_k, , drop = FALSE]` (`:1100`) materialises a second `N3 x S` copy. Hoist `y_pos <- (y > 0)` above the chain loop, and consider folding the `w_all` gather into `sample_pq_cpp()` so the copy never reaches R.

E.  **Make the WAIC accumulation optional, and cheaper.** The three `computeModelLoglik*_cpp()` calls evaluate `R::dbinom` over `n*S + N*S + N3*S` elements at every stored iteration. For binary data `dbinom(y, 1, p, log = TRUE)` is just `y*log(p) + (1-y)*log(1-p)`, several times faster than the general call, and `log(p)` / `log(1-p)` can be tabulated once per iteration (`sample_w_cim_cipp()` already does exactly this at `src/functions.cpp:458-469`). Also add a `computeWAIC = TRUE/FALSE` argument for users who are not doing model comparison.

F.  **Vectorise the starting-value loops.** `R/runOccJSDM.R:946-961` and `:974-981` are triple-nested R loops over `S x n x M` to initialise `w` and `z`. Negligible on the vignette dataset, slow for realistic `S`. Both reduce to grouped "any positive" reductions: `w` from `rowsum(y > 0, idx_w_k) > 0`, `z` from `rowsum(w, idx_z_w) > 0`.

G.  **Do not allocate posterior arrays that are never filled.** The `summarisedLatentPresences = FALSE` half of this was resolved by the fix in Fixed bugs 15 (`w_output_chain` / `theta_output_chain` are now written rather than allocated and abandoned), so what remains is the sizing question: `Bs_output` (`ps x S x niter x nchain`) and `U_output` (`n x d x niter x nchain`) are the two largest components of the 62 MB `sampleresults.rda`; a `keep =` argument selecting which blocks to retain, or storing the latent-factor blocks pre-thinned, would address both the runtime allocation and the CRAN size blocker.

H.  **Reduce the repeated `arma::inv()` calls in the samplers.** `sample_beta_cpp()` (`src/functions.cpp:249-253`), `sampleB()` and `sampleBuniv()` (`src/jsdm.cpp:579-583`, `:600-604`) each call `arma::inv(B)` twice plus `arma::inv(arma::trimatl(L))`, executed `S` or `n` times per iteration. In every caller `B` is diagonal, so this is a dense general inverse of a diagonal matrix. `sampleB_SoR()` (`src/jsdm.cpp:817-839`) already shows the right pattern: take the precision directly as an argument, and draw via a triangular solve against a standard normal instead of forming the inverse Cholesky factor explicitly.
```

## B. Doug to dos

1.  **reproduce all Ecoletts results as a test of the package and decide whether to include in repo**

2.  **extensive testing on simulated datasets** -- **suite built and the R = 100 study run; three things remain.** What exists is summarised under *Completed* below; the authoritative specification and the results table are in `dev/simstudy/PLAN.md`.

    (a) **Re-run once group B items 1, 2 and 4 are fixed.** This is the evidence the fixes worked. Without it they rest on the same code-reading that this exercise showed to be unreliable, and the R = 100 table in `PLAN.md` §12 is already partly stale -- it predates Alex's `42198d9`, which corrected the collection-covariate prior, so `beta_theta` should improve markedly. One command, a few hours:

        ```         
        Rscript dev/simstudy/run_study.R --cores=5 --caffeinate
        ```

    (b) **Decide the replicate count for the paper.** R = 100 was chosen to *detect* defects and did so decisively. Asserting *nominal* coverage in print is a claim about the absence of a small deviation and wants R = 200-500 (`PLAN.md` §9). The runner takes `R` as an argument.

    (c) **Decide how the results are presented** (`PLAN.md` open item 4). Currently a private Claude artifact (<https://claude.ai/code/artifact/ad3d46eb-1fd4-49b5-b795-6b71474ef1d5>), updated 29 July 2026 with the post-fix re-run and the before/after comparison; a pkgdown article is the obvious home -- see item 3.

    **One constraint carried from the bug list:** `l_s` is excluded from coverage checks because it is not recoverable while group B item 1 is open, so no cell of the study speaks to spatial range. Two earlier constraints have since lapsed -- `sigma_h` is now sampled (Fixed bugs 24) and the OpenMP RNG race is closed (Fixed bugs 26), so tier 1's "structural assertions only" rule can be revisited once reproducibility is confirmed on a multi-threaded platform.

3.  **Stand up a pkgdown site.** Discussed 28 July; not started.

    pkgdown turns the package into a static website -- function reference, both vignettes, README, changelog -- and GitHub Pages hosts it. `usethis::use_pkgdown_github_pages()` wires up both plus an Action to rebuild on push. Roughly half an hour, plus time shaping the reference index.

    **Three reasons it earns its place:**

    (a) *It closes a CRAN item.* `DESCRIPTION` has no `URL` or `BugReports` field (CRAN plan item 10 in `AGENTS.md`). A documentation site is a better `URL` than the bare repo.
    (b) *It is a much better landing page for the listserv announcement*, which currently points people at a GitHub file listing. A rendered reference and vignettes make a considerably stronger first impression for a beta release.
    (c) *It gives the validation write-up a home.* The plain-language guide to the test suite and the R = 100 results currently lives as a Claude artifact (<https://claude.ai/code/artifact/ad3d46eb-1fd4-49b5-b795-6b71474ef1d5>, private). As a pkgdown *article* it would sit with the package and be citable as supplementary material for the MEE paper.

    **Two practical constraints.**

    The site would live at `alexdiana.github.io/occJSDM`, since the repo is under Alex's account -- **Alex has to enable Pages**, it is not a setting Doug can change.

    And pkgdown builds every `@examples` block, so the **first build will fail on the group B functions that error unconditionally** (`predictNewSites()` among them). That is the same exposure `\donttest{}` creates under `R CMD check`, which is exactly why the CRAN plan sequences item 8 after the group B fixes. So do this after group B, or expect to `@examples`-guard several functions first.

    **Sequencing caveat:** consider whether to publish a site while group B items 2, 4 and 5 are open. The site would document functions whose credible intervals are currently overconfident, without saying so anywhere a reader would see.

# **Future versions**

1.  use spike-in to estimate abundance change (i.e. eDNAPlus)

2.  model selection of environmental, spatial covariates via regularisation/shrinkage, which would be useful with e.g. geospatial foundation model embeddings as env covariates

3.  parallelisation for speedup

# **Fixed bugs**

This is the authoritative record of audit items that have been resolved. Items 1-10 were fixed in commit `b7b6aa2` (26 July 2026); items 11-20 followed on 27 July 2026. Line numbers point at the *post-fix* code.

Items 16 and 18 are marked **partially fixed**: the crash in each is gone, but part of the original defect remains, and the remainder is tracked as a live item in group B above. Nothing in this list has a test attached -- every entry was verified by reading the code, which is what the *MEE paper* testing item is meant to address.

**Commit b7b6aa2:**

1.  ~~**`data$info` is re-sorted but `data$OTU` is not.**~~ **FIXED.** The sort permutation is now applied to `y`/`OTU` along with `data_info`. `R/runOccJSDM.R:476-558`.

2.  ~~**`sample_pq_cpp()`'s catch-all `else` inflates the Stage 2 FP counts.**~~ **FIXED.** The C++ version now counts the four cases independently, gated on `primerIdx == l`, matching the R reference version. `src/functions.cpp:652-668`.

3.  ~~**`sampleBuniv()` drops the residual precision (operator precedence).**~~ **FIXED.** The expression now correctly evaluates to `1.0 / (sigma * sigma)`. `src/jsdm.cpp:579`.

4.  ~~**`Bs <- list_params$B` in `update_jSDMcoef()`.**~~ **FIXED.** Now correctly reads `list_params$Bs`. `R/jsdmfun.R:940`.

5.  ~~**`sigma_bs` is sampled but never used.**~~ **FIXED.** The spatial block of the prior covariance now correctly uses `sigma_bs^2`. `R/jsdmfun.R:711-713`.

6.  ~~**`useBiotic` is ignored in `computeNewOutputs()`.**~~ **FIXED.** Now tests `useBiotic` where it should (in addition to `useEnvCov`). `src/jsdm.cpp:501`.

7.  ~~**`simulateData()` adds the spatial field *after* the responses are drawn.**~~ **FIXED.** The spatial field is now added before drawing responses. `R/jsdmfun.R:429-461`.

8.  ~~**`tau_output` is returned entirely `NA`.**~~ **FIXED.** The `tau` parameter is now saved in the MCMC loop. `R/runOccJSDM.R:934-1236`.

9.  ~~**WAIC running means divide by the raw iteration counter.**~~ **FIXED.** A dedicated `currentWAICiter` counter is initialised to 1 (`R/runOccJSDM.R:831`) and incremented only inside the WAIC block (`:1200`), so the running means now divide by the number of WAIC accumulations rather than by the raw MCMC iteration index. A guard, `if (numIters != (currentWAICiter - 1)) stop("Current wAIC iter wrong")` (`:1251`), asserts the two agree. Confirmed introduced by `b7b6aa2` via `git log -S currentWAICiter`. *(This fix was not recorded when the other `b7b6aa2` fixes were logged; recovered on 27 July while reconciling cross-references.)*

10. ~~**Non-thread-safe RNG inside the OpenMP loops.**~~ **PARTIALLY FIXED.** `sample_beta_nocov_cpp_TS()` (`src/functions.cpp:308`) now calls the thread-safe `sample_beta_cpp_TS()` rather than the non-TS `sample_beta_cpp()`, which closes the race in `sample_betatheta_cpp_parallel()` (`src/functions.cpp:607`) -- the thread-safe path the audit described as "half-finished" is now wired up on that side. **The `samplePGvariables()` path (`src/jsdm.cpp:398`) is still affected** via `randinvg()`'s `R::rnorm` call; that residual is tracked as **group B item 5** above. *(Also not recorded at the time; recovered 27 July.)*

**Post-audit fixes:**

11. ~~**Collection covariates are grouped by `Sample` alone.**~~ **FIXED.** Now groups by `c("Site", "Sample")` to match the site-then-sample ordering of `X_theta`. `R/runOccJSDM.R:682`.

12. ~~**The GP length-scale is never updated.**~~ **FIXED.** Enabled the `sample_ls()` call in `update_jSDMcoef()`. `R/jsdmfun.R:1063`.

13. ~~**`transformCovariatesMatrix()` never applies the stored factor levels.**~~ **FIXED.** Now correctly uses `df[[col]]` instead of `df[,col]` to re-level the column. `R/runOccJSDM.R:38`.

14. ~~**`listPriors$prior_beta_psi` / `prior_beta_psi_sd` have no effect.**~~ **FIXED.** Removed dead references and cleaned up the prior handling. `R/runOccJSDM.R:738-739`.

15. ~~**`summarisedLatentPresences = FALSE` errors on two-stage models.**~~ **FIXED.** Now correctly allocates and writes `w_output_chain` and `theta_output_chain`. `R/runOccJSDM.R:1152`.

16. ~~**`thinOutput()` cannot run at all.**~~ **PARTIALLY FIXED.** The function still exists (`R/output.R:18`) and now runs: `niter` is read from `fitModel$results_output$jsdm_output$B0_output` instead of the long-gone `beta_ord_output`, and the `thin` argument is honoured rather than hard-coding `by = 5`. **Two of the three original defects remain** -- the 2-D branch still thins by row, and the scalar `WAIC` still falls through to `print("Dimension not recognised")`. Tracked as **group B item 2** above.

17. ~~**`computeDiagnostics()` runs on `psi_output`.**~~ **FIXED.** Added `psi_output` to the skip list to avoid meaningless diagnostics. `R/diagnostics.R:403-404`.

18. ~~**`predictNewSites()` does not honour its documented no-op behaviour.**~~ **PARTIALLY FIXED -- this entry previously overstated what landed; corrected 27 July 2026.** The gating on the presence of covariates was done (`R/output.R:1457,1470`). The **`NULL` defaults were not**: `formals(predictNewSites)` shows `X_psi` and `X_s` still have no defaults at all, so `predictNewSites(fit, X_psi = X)` still fails on the missing `X_s` promise even when `useSpatial = FALSE`. Residual tracked as **group B item 3**.

    Third entry in this list found to overstate a fix (with 16 and 20), which is the argument for the test suite: every one of these was recorded from reading a diff rather than from running the code.

19. ~~**`set.seed(1)` inside `computeSpeciesDetected()`.**~~ **FIXED incidentally in b7b6aa2.** The RNG-resetting line was deleted. `R/output.R:2222`.

20. ~~**`computeSpeciesDetected()` crashes: `B` (bootstrap draw count) is referenced but never assigned.**~~ **PARTIALLY FIXED.** The crash is gone: `idxObs <- unique(floor(seq(1, min(500, niter), length.out = 500)))` and `B <- length(idxObs)` (`R/output.R:2221-2222`) both assign `B` and guard the fewer-than-500-draws case that the original hard-coded `B <- 200` did not. **The sampling is still a prefix, not a random draw** across chains, which was the substantive half of the original critique. **Now fully fixed** -- see Fixed bugs 23.

**Fix of 27 July 2026:**

21. ~~**`plotFPTPStage2Rates()` ignores `primerName`.**~~ **FIXED.** `p_output`/`q_output` are now subset to the requested primer *before* the quantiles are computed (`p_output[idx_primer, , , , drop = FALSE]`, `R/output.R:790-791`), instead of `apply(p_output, 2, ...)` pooling over primer, iteration and chain. `primerName` is matched against `fitModel$infos$primerNames` via `match(as.character(...))`, so both `primerName = 2` and `primerName = "2"` work (`primerNames` is stored as an *integer* vector, so the string form would otherwise silently fail); an unrecognised primer now errors with the list of available primers rather than quietly plotting the pooled result. The dead `idx_speciesprimer <- str_match(...)` line is deleted, and the plot title now reports which primer is shown.

    **Verified against `sampleresults`** (`P = 3`): the three primers now give three different plots (previously identical); each primer's intervals match hand-computed `quantile()` values on that primer's slice; no primer reproduces the old pooled result; the default still equals primer 1; `idx_species` subsetting and rendering are unaffected.

    **Note:** this was the *only* use of `stringr` anywhere in `R/`, so the package no longer needs it. Three `@import stringr` tags remain (`R/output.R:757`, `:857`, `:1008` -- the latter two were already spurious) along with `import(stringr)` in `NAMESPACE` and `stringr` in `DESCRIPTION`. Removing all three is a small, safe cleanup, best folded into the group C dead-code pass.

22. ~~**REGRESSION: `runOccJSDM()` crashes on every fit that has no spatial covariates.**~~ **FIXED** by Alex in `d6b70b1` ("Fixed bug on sample l_s"), the same day it was reported. Any call without `spatCovariates` had been dying in the first MCMC iteration with `missing value where TRUE/FALSE needed`, thrown from the acceptance test in `sample_ls()`, because `precomputeSORmatrices()` leaves `logDetKuu_grid`/`Lm1_grid` entirely `NA` when there is no spatial field but `update_jSDMcoef()` called `sample_ls()` anyway.

    **The actual fix was a one-word correction, and a better one than the fix suggested when this was filed.** The guard was `ps <- ncol(Bs)`; it is now `ps <- nrow(Bs)` (`R/jsdmfun.R`, in `update_jSDMcoef()`'s "read state variables" block). `Bs` is `[ps x S]` -- spatial basis dimensions by species -- so `ncol(Bs)` was the *species* count, which is never 0. `if (ps > 0)` was therefore always true regardless of whether a spatial field existed. `nrow(Bs)` is the spatial dimension and is correctly 0 with no spatial covariates, so `sample_ls()` is simply not reached. This addresses the root cause (a transposed dimension in the guard) rather than the workaround originally proposed here, which was to re-gate on `X_centers > 0`.

    **Verified after the fix** against the shipped `sampledata`: both `runOccJSDM(..., spatCovariates = c("Xs.1","Xs.2"))` and the same call *without* `spatCovariates` now complete. Tracing confirms the mechanism -- `precomputeSORmatrices()` still reports `X_centers = 0, grids all NA = TRUE` in the non-spatial case, but `sample_ls()` is no longer called, so the `NA`s are never read.

    **Still open from the original report:** the small-`n` `kmeans()` failure found while reproducing this (`number of cluster centres must lie between 1 and nrow(x)` at `n = 30` *with* spatial covariates) was not addressed by `d6b70b1` and has not been re-tested. It constrains the feasible grid for the simulation-study vignette, so confirm it before designing that.

23. ~~**`computeSpeciesDetected()` uses a prefix of the draws, not a sample across chains.**~~ **FIXED.** `R/output.R:2229` now reads `idxObs <- unique(floor(seq(1, niter, length.out = 500)))` -- the `min(500, niter)` cap is gone, so the 500 indices are spread evenly across the *entire* flattened draw set rather than landing in its first 500 rows. Since `beta_theta_output` is flattened iter-within-chain at `:2222`, spanning the full range means every chain now contributes. This resolves the substantive half of the original audit critique (see Fixed bugs 20, which fixed only the crash).

    Note the result is a *systematic* thin rather than the random draw from `seq_len(niter * nchain)` originally suggested. That is equivalent for this purpose and is standard practice for thinning an MCMC chain, so no further change is needed.

**Fixes of 28-29 July 2026 (Alex, `42198d9` / `e60e3ad`).** Each verified against the source rather than the commit message.

24. ~~**`sigma_h` is never sampled.**~~ **FIXED.** A new `sample_sigmah(U, a_sigmah, b_sigmah)` (`R/jsdmfun.R:1173`) is now called each iteration (`:1694`), so `sigmah_output` varies instead of sitting at its initialised 1. This was inert for fitted models -- `U` at the training sites is drawn under a hard-coded unit-variance prior regardless -- but it is the factor-score SD `predictNewSites()` uses when simulating `U` at *new* sites, so out-of-sample predictions no longer silently assume unit factor variance.

25. ~~**The prior mean for collection-covariate slopes is 1, not 0.**~~ **FIXED.** `b_betatheta <- rep(0, ncov_theta)` (`R/runOccJSDM.R:757`); previously `rep(1, ...)`, with only element 1 (the intercept) overwritten, so every actual covariate slope carried a prior centred on +1. The prior variance was also widened from `diag(1)` to `diag(2)`.

    This was the most clearly evidenced item in the whole audit: at R = 100 the intercept -- the one coefficient whose prior mean *was* overridden -- covered at 0.937, while the slopes covered at 0.496 with bias +0.113, against true slopes averaging +0.01. Across the full grid `beta_theta` undercovered in **every** cell (0.676-0.730).

    **Confirmed by the 29 July re-run**: `beta_theta` improved in every cell, to 0.709-0.771, with two-thirds of the bias removed. It is still well short of nominal, so a second cause remains -- tracked as group B item 4.

26. ~~**Non-thread-safe RNG in the hottest OpenMP loop.**~~ **FIXED.** `randinvg()` (`src/jsdm.cpp:86`) now draws from the `thread_local` `rnorm()` rather than `R::rnorm`; the old line is commented out beside it. That closes the last hole on the `samplePGvariables()` path, whose Polya-Gamma helpers were already converted in `53c38f1`.

    **Two consequences worth acting on.** A fixed seed should now reproduce on Linux and Windows, which (i) removes the reason `sampleresults.rda` had to be refit on macOS or an OpenMP-disabled build, and (ii) unblocks CRAN item 16. It also means the tier-1 constraint in `dev/simstudy/PLAN.md` §5.1 -- structural assertions only, never numeric equality -- can be revisited; it was adopted solely because of this race. Verify reproducibility on a multi-threaded platform before relying on any of that.

27. ~~**Stage 2 hyperparameters documented as settable but never read.**~~ **FIXED.** `a_p`, `b_p`, `a_q` and `b_q` now read from `listPriors` (`R/runOccJSDM.R:751-754`), matching the behaviour `@param listPriors` already claimed, and the documented defaults were corrected to match the code (5/1 and 1/20). A user running a low-detection study can now override the prior without editing the package.

    **Only the wiring is closed.** What the default *should be* remains open and is a design decision, not a defect -- see group B item 3.

28. **`set.seed()` did not control any of the C++ samplers, so `runOccJSDM()` was not reproducible.** Found 29 July 2026 while writing the regression test for Fixed bugs 26; fixed the same day.

    Fixed bugs 26 replaced `randinvg()`'s use of R's global RNG inside an OpenMP loop with `thread_local` engines, which correctly closed the data race. But neither replacement engine ever read R's RNG state: `get_rng()` was seeded from the literal `12345 + omp_get_thread_num()`, and `mvrnormArmaQuick_TS()` from `std::random_device{}()`. Measured before the fix: two fits of the same fixture under the same `set.seed(4242)` differed by 5.09 on `B0_output`.

    Neither biased the posterior -- both are valid streams -- but users could not reproduce a fit, and the literal `12345` was per *process* rather than per run, so every simulation-study worker started its Polya-Gamma stream at the same position and replicates at the same position in different workers shared random numbers.

    **The fix** (`src/rng.h`, new). Both translation units now share one per-thread `mt19937`, seeded via `std::seed_seq` from a base seed that `runOccJSDM()` draws from R with `sample.int()` and installs through the new `setOccJSDMSeed()`. `set.seed()` therefore controls the whole fit, and consecutive fits stay independent because R's RNG advances between them. Three subtleties, each of which cost a debugging round:

    (a) A `thread_local` engine is constructed once per thread and would keep its original seed for the life of the session, so a generation counter triggers re-seeding when R installs a new base seed.
    (b) The generation must *not* be seed material. Mixing it in made the same base seed yield a different stream on the second fit of a session -- reproducible-looking in a fresh process, broken in a suite.
    (c) `std::normal_distribution` caches the second Box-Muller deviate, and that cache survives `rng.seed()`. Without an explicit `dist.reset()`, the first `rnorm()` of a fit could return a value left over from the previous fit, so reproducibility depended on the parity of normal draws earlier in the process. This is why the test passed standalone and failed inside the suite.

    `arma::randn()` elsewhere in the package was never affected: RcppArmadillo sets `ARMA_RNG_ALT` to route Armadillo's RNG to R's, so those draws already honoured `set.seed()`.

    **Scope.** Reproducibility holds for a given thread count. Threads derive separate streams, so if the package is ever built with OpenMP actually enabled (it is not on the macOS dev machine -- see group D), changing the thread count changes which stream produces which element. Inherent to per-thread streams, not a defect.

    **Tests.** `test-regression-bugs.R` covers same-seed equality, different-seed inequality, and independence of consecutive fits. Verified separately that the simulation harness gives identical output for a repeated replicate and different output across replicates.

    **Reviewed by Alex, 31 July 2026: "We don't care about reproducibility."** Taken as a decision not to invest further, not a request to revert, and the fix stays. Worth recording why: reproducibility here is **load-bearing internally even though it is not a user-facing priority**. The simulation study's paired design depends on it, and that pairing is what produced the strongest evidence in the whole study, that only 104 of 49,978 `resid_cor` coverage decisions flipped between the pre- and post-fix runs on identical truths. Remove the R-derived seeding and every future before/after comparison loses that power. The tier-1 test at `test-regression-bugs.R:243` guards it and should stay.

29. **The sparse-GP knot default no longer floors at 30 or crashes below 31 sites.** Filed 27 July 2026 (then group B item 3) after the simulation study hit it; fixed by Alex in `42198d9`, unlogged. Verified 29 July: `getDefaultSupportPoints()` (`R/jsdmfun.R:875`) is now `min(floor(n * 0.2), n - 1)`. The old `max(30, floor(n * 0.2))` fed `kmeans(X_s, centers = ps)` and so was a constant 30 for any dataset below 150 sites -- roughly one knot per site at n = 31, defeating the point of a sparse GP -- and errored outright below 31. The `n - 1` cap is what removes the crash.

30. **`ds = 0` no longer produces a null spatial field.** Filed 27 July 2026 (then group B item 4); fixed by Alex in `42198d9`, unlogged. The simulator's cross-species spatial covariance used to collapse to jitter at `ds = 0` -- measured `sd(spatField)` of 0.0019 against \~1.0 at `ds = 2` -- so any scenario built at `ds = 0` was silently a null-field test. Verified 29 July at seed 42: `sd(spatField)` is 0.598 at `ds = 0`, 0.467 at `ds = 1`, 0.678 at `ds = 2`. The study grid still uses `ds = 2`, now by choice rather than necessity.

    *Both of these were removed from group B by `42198d9` without a Fixed-bugs entry, which is why they are recorded here late. The check that caught it: `TODO.md`'s group B numbering had a gap at 5 and 6.*

31. **`plotCollectionRates()` errored on every input.** Reported by Doug 29 July 2026, fixed the same day. Failed with `object 'Min' not found` for any `fitModel`, with or without `idx_species`.

    `plotSpeciesRates()` (`R/output.R:819`) had been extracted as a shared helper and never wired up to its only caller. Three independent breakages in the same call path, which is why nothing had ever run it successfully:

    (a) the helper read columns `Min`/`Max`, while `plotCollectionRates()` passed the `2.5%`/`97.5%` that `quantile()` names;
    (b) it filtered on a `Species` column the caller never created -- the code that built one is still there, commented out, from before the extraction;
    (c) it referenced `speciesNames` as a free variable, present in neither its arguments nor the package namespace.

    Its roxygen documented a third state again: columns `Species`, `2.5%`, `97.5%`, and parameters `orderSpecies`/`subset` that the signature does not have.

    **Fixed** by giving the helper an explicit contract -- one row per species with `Min`/`Max`, plus `idx_species` and the *unsubset* `speciesNames` -- and doing the subsetting inside it, so the ordering is computed on exactly the rows plotted. Also added a clear error for models with no collection stage, which previously fell through to a subscript failure.

    **Note, since resolved:** this fixed the species-ordering defect for this function only. `plotOccupancyRates()`, `plotFPTPStage2Rates()` and `plotStage1FPRates()` still ordered on the full species set while indexing a filtered one (this was filed as group C at the time, not group B as originally written here). The test added here asserted label-to-value pairing rather than mere absence of error, and served as the template for fixing those three -- see Fixed bugs 32.

32. **`plotOccupancyRates()`, `plotFPTPStage2Rates()` and `plotStage1FPRates()` shared the species-ordering defect `plotCollectionRates()` had (Fixed bugs 31).** Filed as group C at the time; fixed by Claude 29 July 2026. Each computed `order()` on the *filtered* `idx_species` subset and then used the result to index the *unfiltered* `speciesNames`, so for any `idx_species` other than a prefix `1:k` the factor levels named the wrong species and bars silently vanished.

    `plotOccupancyRates()` and `plotStage1FPRates()` now delegate to the `plotSpeciesRates()` helper fixed in Fixed bugs 31, which subsets first and derives labels from the subset. `plotFPTPStage2Rates()` has a two-interval (`p`/`q`) layout that doesn't fit that helper, so it was fixed inline: order and labels are both now derived from the filtered `data_plot`. `plotStage2FPRates()` was already correct and untouched.

    Verified against a live fit with `idx_species = c(3, 1, 10)`: all four functions now plot exactly that subset with matching labels. Full test suite passes (119/119). `R/output.R`.

33. **`returnCovariateEffect()`/`plotCovariateEffect()` had no `idx_species` default, and fixing that exposed two further bugs in the code they call.** Filed as group C at the time; fixed by Claude 29 July 2026.

    Both functions declared `idx_species` with no default, so `returnCovariateEffect(fit, covName)` errored instead of defaulting to all species -- the same gap as `predictNewSites()` (Fixed bugs 34). Gave both a `NULL` default resolving to all species, matching every other return/plot function in `R/output.R`.

    That default immediately reached two pre-existing bugs that a narrow, prefix `idx_species` had been masking:

    (a) `plotCovariateEffect_base()` (`R/jsdmfun.R`) re-applied `apply(B_output, c(1,2), c)` to arrays its only caller, `plotCovariateEffect()`, had already collapsed. The second `apply()` collapsed the species margin again, leaving the array's third dimension sized by `ncov_psi` instead of `S` -- any species index beyond `ncov_psi` (2 in the test fit) errored with `subscript out of bounds`, exactly the region "all species" reaches immediately. Renamed the formals to `B0_output_vec`/ `B_output_vec` to match the already-collapsed inputs and dropped the redundant `apply()` calls.

    (b) `returnCovariateEffect_base()`'s per-species loop labelled each row `speciesNames[i]` (the loop counter) instead of `speciesNames[sp_idx]` -- the same class of defect as Fixed bugs 32 -- mislabelling every species whenever `idx_species` wasn't the prefix `1:k`.

    Verified against a live fit: the default now facets all 10 species correctly; `idx_species = c(3, 1)` plots and labels `OTU_3`/`OTU_1` correctly (previously mislabelled `OTU_1`/`OTU_2`, then errored once the default was exercised). Full test suite passes (119/119). `R/jsdmfun.R`, `R/output.R`.

34. ~~**`predictNewSites()` could not be called without supplying both `X_psi` and `X_s`.**~~ **FIXED 30 July 2026** (Claude; `R/output.R`, `src/jsdm.cpp`). Residual of the original audit's B.4, previously tracked in group C.

    **The filed defect.** `X_psi` and `X_s` had no defaults, so the `is.null()` guards the author had written could never fire: R raised `argument "X_psi" is missing, with no default` first.

    **Worse than filed.** Those guards used `&`, not `&&`. Since `&` evaluates both sides, the missing promises were forced even when the caller had asked for neither term, so `predictNewSites(fit, useEnvCov = FALSE, useSpatial = FALSE)` failed too. There was no way to call the function at all without both matrices.

    **The fix.** Both default to `NULL`. `useEnvCov` and `useSpatial` adopt the tri-state `useBiotic` already used: `NULL` uses the term if the fit estimated it, `TRUE` uses it and errors if it did not, `FALSE` skips it. That is what their roxygen always claimed, and was unreachable because `useSpatial` hard-stopped on a fit with no spatial field instead of ignoring it. Strictly more permissive: every call that worked before still works identically.

    **Three further defects, exposed by making those paths reachable, all pre-existing and all fixed here:**

    (a) `computeNewOutputs()` sliced `Ks_all` and `Bs_output` unconditionally, so `useSpatial = FALSE` aborted with `Cube::slice(): index out of bounds`. Moved inside the `useSpatial` guard.
    (b) It read the new-site count as `X.n_rows` unconditionally, so `useEnvCov = FALSE` silently returned a zero-row result. Now taken from whichever term is active; both being off is an explicit error, since nothing then determines how many sites to predict for.
    (c) A fit with no spatial field returns `Bs_output` with a zero-length first dimension, which collapses under `apply()` and broke the `aperm()` with `'perm' is of wrong length`. Only reshaped when the spatial term is in play.

    **Tested beyond absence of error** (`test-api-contracts.R`): each term measurably changes the prediction when toggled, probabilities stay in `[0,1]`, and quantiles stay ordered. Without the first of those the switches could have been cosmetic and the shape assertions would still have passed.

    **Noticed, not fixed:** `computeNewOutputs()` prints `Computing species i out of S` to stdout via `Rcout` on every call, unconditionally, and it cannot be silenced. Filed separately in group C.

35. **The vignette could not be built, so `R CMD check` never reached code inspection.** **FIXED 30 July 2026** by Doug regenerating `data/sampleresults.rda`.

    Two failures in sequence, each hidden behind the previous one. First `plotCollectionRates()` errored on every input (Fixed bugs 31). With that fixed, the build failed in `plotCovariateEffect()`, apparently for want of a `covNames` default; naming the covariate in the chunk did not fix it either, and it then failed with `'from' must be a finite number`.

    **The real cause was neither function.** Both read `fitModel$infos$X0_psi`, the raw occupancy covariates, to build a prediction grid. `X0_psi` was added by `e60e3ad` ("Added GAMs"), and the shipped `sampleresults.rda` had last been regenerated on 26 July, before that. `min(NULL, na.rm = TRUE)` is `Inf`, and `seq(Inf, -Inf, ...)` throws. `X0_psi` was the only name a current fit's `infos` carried that the shipped object lacked.

    **Confirmed after the refit:** `X0_psi` is present (100 x 2), `list_X_psi_mat` carries `bs_info`/`target_spline_vars`, both vignette chunks run, and `devtools::check()` completes with **0 errors** for the first time. Remaining status: 2 warnings, 3 notes, tracked separately.

    **The lesson worth keeping**: a stale shipped dataset presents as a bug in whatever function touches it first, and moves to the next function each time one is fixed. Two functions were investigated and one was needlessly suspected before the data was. When an example object fails and a fresh fit does not, suspect the object.

# **Completed work**

Finished work, kept for context rather than as tasks. Bug fixes live under *Fixed bugs* above; this is everything else.

1.  ~~**Purge `data/traitdata_caiwang.rdata` from git history.**~~ **DONE.** Purged via `git filter-repo` and force-pushed the rewritten `main` to `origin`. All commit hashes from `11c5449` onward changed as a result, so **Alex needed to re-clone or hard-reset**; that has since happened.

2.  ~~**Port the model-diagnostics functions from the GLGS-eDNA repo.**~~ **DONE.** Alex's note: *"I ALREADY ADDED DIAGNOSTICS SO NOT NEEDED ANYMORE"*. `computeRhat()`, `summarisePosterior()`, `plotTraceplot()` and `returnConvergenceDiagnostics()` are in `R/diagnostics.R`; the first two were later de-exported (`@noRd`) as redundant with `returnConvergenceDiagnostics()` for user-facing purposes. A "Model diagnostics" vignette section was added to `vignettes/occJSDM.Rmd`.

    **Not ported, and still available if wanted:** `compare_to_true()`/`plot_estimated_vs_true()`, which need `true_params` from a simulation. The simulation suite (item 3 below) now covers that ground more systematically, so these are probably redundant rather than outstanding.

3.  ~~**Build the simulation test suite and run the coverage study.**~~ **DONE 28 July 2026**, except the three follow-ups still listed under *MEE paper / Doug to dos 2*. Full specification and results: `dev/simstudy/PLAN.md`. Provenance and lessons: `AGENTS.md`.

    - **Tier 1** (`6d9526d`): 89 structural assertions, \~7 s, **ships to CRAN**. `helper-fixtures.R`, `test-smoke-configs.R`, `test-regression-bugs.R`, `test-api-contracts.R`. `test-placeholder.R` deleted.

    - **Tier 2** (`f63eeeb`): recovery canary, \~30 s, `skip_on_cran()`, thresholds measured over eight seed sets.

    - **Tier 3 + runner** (`f63eeeb`, `6123036`, `e0718f1`, `32b56a2`): env-gated study plus `dev/simstudy/run_study.R`, with checkpointing, `--resume`, `--caffeinate` and progress reporting.

    - **The run itself:** re-run 29 July 2026 after Alex's fixes and the RNG seeding fix. 10 scenarios x 100 replicates, 1000 fits, 0 failures, 285 min on 5 cores, 155,578 interval checks. Results in `dev/simstudy/results/`, tabulated in `PLAN.md` §12. The 27-28 July pre-fix run is retained for comparison (`PLAN.md` §12.6).

    - **It is a paired comparison, and that is worth protecting.** `draw_truth()` seeds on (scenario, replicate), so the simulated data and true values are bit-identical between the pre- and post-fix runs -- verified, `max|truth difference| = 0`. Every difference is therefore attributable to the code rather than to sampling variation, which is what makes statements like "only 104 of 49,978 `resid_cor` decisions flipped" possible. **Do not change `simstudy_seed()`**, or future runs lose comparability with these two.

    It found six defects -- three while the tests were being written, three from the run -- none of which the static audit had caught. Four of the six are now fixed (Fixed bugs 24-27 and group B); the rest are group B items 1, 2 and 4.
