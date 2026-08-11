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

1.  **`thinOutput()` rewritten.** 4 August 2026 (Claude; `R/output.R`, `man/thinOutput.Rd` deleted, `tests/testthat/test-regression-bugs.R`). Was group C item 1. No *Fixed bugs* entry yet, deliberately: it gets one when you have reviewed it.

    **Three defects, two of them worse than the item recorded.** `jsdm_output` is a *list* of 17 parameter arrays, so `length(dim(x))` was 0, it matched none of the five branches, and fell through to `print("Dimension not recognised")` -- which returns its argument, so the whole list was replaced by that character string. The scalar `WAIC` was destroyed identically, becoming a character string rather than the `NULL` the item predicted. Only the third, 2-D matrices thinned by row and so losing sites, was as described.

    **The fix rests on one rule, and that is the thing worth your eye.** Every array is thinned on its **second-to-last axis**, on the assumption that `runOccJSDM()` always stores iterations there with `nchain` last. True of all 27 arrays today, checked individually, so the five branches collapse to none. Iteration arrays are told apart from posterior means by requiring the last two extents to equal `(niter, nchain)` as a pair; matching only the last would misfire whenever `S == nchain`. If you add an output whose axes are ordered differently, it will be silently left unthinned rather than mangled -- safe, but wrong, and nothing will complain.

    **Verified** on real fits with `niter`, `n`, `N`, `S` and `nchain` all deliberately distinct: an earlier ad-hoc check used `niter == n == 40`, which makes a `[sites x species]` mean shape-indistinguishable from an iteration array and hides the whole bug. All 17 `jsdm_output` elements thin on the correct axis with other extents intact; `WAIC` stays numeric; the mean matrices keep every site; retained draws are the right iterations (max difference 0, not merely the right count); `thin = 1` is exactly the identity. Both `summarisedLatentPresences` settings, since `z_output`/`psi_output` are 2-D means under one and 4-D arrays under the other. Tier-1 regression test added; suite 297 passing.

    **Two decisions left open, both yours.** It is still unexported, per your "disable it for now, but let's keep it there" on the original item, so using it means `occJSDM:::thinOutput()`. And it now has no caller anywhere: the CRAN plan step that wanted it (option (a), thinning `sampleresults.rda`) is superseded by refitting small instead, which fixes the correctness problem thinning cannot. Keeping it is defensible, deleting it is defensible; it should not just drift.

    ALEX TO REVIEW

## **B. Inference-affecting bugs (wrong numbers, silently) (Alex)**

1.  **`sample_ls()` scores the wrong density, so the GP length-scale is never recovered.** `R/jsdmfun.R:1054`. Under the SoR approximation the fitted field `SE = Ks(l_s) %*% Bs` is a deterministic function of `Bs` and `l_s`, but `sample_ls()` treats `SE` as a GP draw and scores it under `N(0, sigma_s^2 K(l_s))` with `SE` held fixed. That is not the conditional posterior, and it is self-defeating: `SE` was already smoothed at the current `l_s`, so it scores better under ever-smoother covariances.

    **Measured:** `idx_ls` rails at the top of `l_s_grid` for every true `l_s` tried (0.074, 0.171, 0.300), with real spatial signal present. The profiled log-likelihood rises monotonically from -376 at `l_s = 0.01` to -154 at `l_s = 0.30`.

    **Impact:** biases the spatial term of every `useSpatField = TRUE` fit and makes `l_s` uninterpretable. Widening the grid would only move the rail.

    **Needs a derivation, not a code tweak.** Either recompute `eta` with `Ks(l_s*)` and score the observation likelihood under the proposal, or integrate out `Bs` and use the SoR marginal likelihood. **Do not attempt blind:** a subtly wrong conditional would be worse than the present state, which at least fails visibly. A cheap interim option is to revert to a documented, user-settable fixed `l_s`.

    Four candidate causes were investigated and ruled out (missing amplitude, wrong amplitude passed, `logDetKuu`, weak data). Detail in AGENTS.

    ALEX TO CHECK

2.  **`reparamFactorModel()` breaks residual covariance = `t(L) %*% L`, inflating reported species correlations.** `R/jsdmfun.R:48`. The rotation preserves `U %*% L` (verified to 4e-16) so the linear predictor is untouched, but it moves scale out of `U` into `L`, and `returnResidualCorrelationMatrix()` computes `cov2cor(t(L) %*% L)` from the reparameterised `L`. Measured `Var(U)` afterwards is `diag(0.23, 2.01)`, not the identity.

    **Measured:** correlations move by up to 0.612, consistently toward the extremes.

    **Impact:** `returnResidualCorrelationMatrix()` and `plotResidualCorrelationMatrix()` overstate co-occurrence. This is the headline JSDM output.

    **Fix, two options.** (a) Rotate by `Q` alone, dropping the `diag(diag(R))` scaling, so both the identifiability constraint and the covariance identity hold. (b) Keep the scaling and compute the correlation as `t(L) %*% Var(U) %*% L`. (a) is simpler and preserves the output contract.

    **REOPENED 2 August 2026: the simulation evidence for this item has been withdrawn.** It previously cited `resid_cor` covering at 0.74-0.77 across the grid, and a paired re-run in which only 104 of 49,978 coverage decisions flipped. The 2 August re-run shows that statistic is degenerate: coverage equals one minus the share of true correlations sitting at exactly ±1, the credible intervals span almost the whole of [-1, 1] (median width 1.999 on a scale bounded by 2), and a truth-determined statistic cannot move when the sampler changes, so the paired result never discriminated. `PLAN.md` §17 has the measurements and §17.6 what would make the statistic informative.

    **This does not refute the item.** The code-level argument above is about a per-factor rescaling surviving `cov2cor()` and is untouched. What has gone is the independent evidence for it, so it now rests on reading the transform alone. Your note was that this is a non-issue for logistic models since only the correlation is recoverable; the premise is right, and the reply had been that the measured numbers were already correlations. That reply is no longer available. Full exchange in AGENTS.

    **Confirmed on a second independent run, 10 August 2026, and on every cell rather than two.** Predicted coverage (one minus the share of true correlations at exactly ±1) equals measured coverage to three decimal places in **nine of the ten** production cells, spanning 0.752 to 0.980. The statistic is degenerate everywhere, not just where it was first noticed.

    **The tenth cell is the one worth knowing about.** `d_overfit` predicts 0.760 and measures 0.756, and it is the only cell whose intervals are not effectively the whole range: median width 1.79 against 2.00 in the other nine. So it is the single cell where `resid_cor` coverage carries any signal at all, and even there it is marginal. Anyone acting on `PLAN.md` §17.6 should start from why that cell differs.

    ALEX TO MAKE A DECISION, now on the code argument alone

3.  **`beta_theta` intervals are overconfident, and it gets worse with more data.** Coverage 0.77 at the production `M = 2`, falling monotonically to 0.58 at `M = 20`, while bias stays small and flat. Shrinking intervals around a bias that is not shrinking is the signature of a real defect being exposed by more information, not fixed by it.

    **Replicated on a third data axis, 10 August 2026, and this is the strongest evidence the item has.** Until now everything rested on the M and K ladders, both of which add data to the *detection* stage. The `sites_300` cell triples the number of sites, which feeds the *occupancy* stage instead, and `beta_theta` coverage falls **0.761 to 0.646** (`PLAN.md` §20). Same direction, comparable magnitude, different pathway.

    **It is not a general degradation, which is what makes it diagnostic.** In the same run, accuracy improves sharply everywhere -- `p` nRMSE 0.822 to 0.418, `q` 0.833 to 0.554, `B` 0.511 to 0.384 -- and coverage of `p` and `G` improves slightly. Only `beta_theta` moves materially the wrong way. More information is making every estimate better while making this one parameter's intervals more confidently wrong, which is the signature the item has claimed from the start.

    **Narrowed 31 July: the defect is in the *slopes*, not in `beta_theta` as a block.** Refitting `base` with `ncov_theta = 0`, so only the intercept row remains, gives `beta_theta` coverage of **0.968** (SE 0.013, R = 200), i.e. nominal, against 0.763 with the slopes present (`PLAN.md` 15.5, 15.6). Whatever is wrong is specific to the covariate columns.

    **Four candidate causes now ruled out**, each by measurement: Stage 1 under-identification (more data makes it worse, not better); the slope prior's width (tightening it 20-fold at `M = 2` moves coverage the wrong way); pseudo-replication in `X_theta` (it is drawn per sample, not per site); and the intercept path (nominal once the slopes are gone).

    **So the cause is in whatever handles the covariate columns in the Polya-Gamma update**, in `sample_beta_cpp_TS`/`sample_betatheta_cpp_parallel`. Note this is the same step Alex's `microbenchmark()` profiling identified as the slowest in the sampler, so the calibration problem and the performance bottleneck sit in the same code. This needs someone who knows it; it is not another prior experiment. Evidence in `PLAN.md` 13, 14 and 15.5.

    ALEX TO INVESTIGATE THE SAMPLER

4.  **Decide `b_betatheta`'s slope prior variance. It trades `B0` bias against `beta_theta` coverage.** Measured at `M = 2`, paired on identical truths, varying only that variance:

    - variance 2, your current default: `B0` bias -0.160, `beta_theta` coverage 0.747
    - variance 0.5: `B0` bias -0.106, `beta_theta` coverage 0.707
    - variance 0.1: `B0` bias -0.044, `beta_theta` coverage 0.653

    **This identified the cause of `B0`'s doubled bias, closed by decision as *Fixed bugs* 46.** `42198d9` widened `B_betatheta` from `diag(1)` to `diag(2)`, which is exactly when the bias doubled; turning it back down moves it back, monotonically. Alex's decision closed the bias question rather than this variance trade, which stands on its own regardless.

    **But it is a trade, not a fix:** tightening helps `B0` and hurts `beta_theta` coverage. This item and the `beta_theta` slope item above pull opposite ways on one knob. Not known: whether an intermediate value beats both endpoints, whether the trade holds at `M > 2`, and whether fixing the `beta_theta` slope defect at its source would dissolve it entirely.

    ALEX TO DECIDE THE VALUE (or that fixing the slope defect supersedes this)

5.  **`B0` coverage undercovers by 4.7 SE in the `continuous` arm, and has not been chased.** 0.879 against nominal 0.95 (`PLAN.md` 16.5), bias zero so the interval is too narrow rather than the estimate wrong. Split out of the `B0` bias item below when that one closed, since Alex's decision addressed the bias, not this. One arm, one configuration, found while looking for something else -- wants confirming at a second configuration before it is called a defect.

    **Correction, 10 August 2026: this item claimed 0.879 was "the lowest `B0` coverage of any cell measured". It was not, and was not when written.** `low_information` sat at 0.865 in the same 2 August run the claim was drawn from, and measures 0.871 in the 10 August re-run. Two independent runs, so the low reading there is real rather than noise.

    **That does not discharge the request above.** `low_information` is degraded across the board -- `p` covers at 0.113 in that cell -- so `B0` sagging there is unsurprising and probably a different phenomenon from `continuous`, where every other block was healthy. What is wanted is still a second *clean* configuration. The correction is to the superlative, not to the item.

    CLAUDE OR ALEX TO CONFIRM AT A SECOND CLEAN CONFIGURATION

6.  **`theta0`'s intervals are \~25% wider than they need to be, and that is the price of its bias being fixed.** Coverage 0.978-0.985 post-fix against 0.938-0.959 pre-fix (`PLAN.md` 12.3). The all-cell average of 0.944 hides it, because `low_information` pulls it down at 0.602.

    **Re-read 2 August from the two saved runs, and the framing above was wrong.** Comparing `simstudy-20260728-175534.rds` (pre-fix) with `simstudy-20260729-143756.rds` (post-fix) on identical data, over all cells except `low_information`:

    - coverage 0.938-0.959 -\> 0.978-0.985
    - mean interval width 0.113 -\> 0.143, i.e. **+25%**
    - mean absolute bias **0.0175 -\> 0.0020**, a factor of nine

    **`theta0`'s point estimate went from clearly biased to essentially unbiased.** That was not recorded anywhere, and it inverts the item. Pre-fix coverage near nominal was a *coincidence*, not health: the `Beta(1, 20)` prior mean of 0.0476 sits below the truth mean of 0.06, so estimates were pulled down, and intervals that were too narrow offset that bias almost exactly. Two errors cancelling. The fixes removed the bias and left the width, so what looks like a regression in the coverage column is a genuine improvement in the bias column with an unaddressed remainder.

    **So this is not "`theta0` was fine and broke".** It is "`theta0` was quietly biased, is no longer, and its intervals have not caught up".

    **One of the two "ruled out" causes is only half ruled out.** The M ladder reading (overcoverage falling toward nominal as M rises, 0.986 at `M2` to 0.944 at `M10`, while the matched `K30` control worsens to 0.996) still has the hole it always had: pre-fix, `theta0` was fine at the *same* M = 2. The coupling hypothesis is the one to reopen. *Fixed bugs* 25 changed `b_betatheta`'s prior **mean** (1 to 0) *and* widened its **variance** (`diag(1)` to `diag(2)`). `PLAN.md` 14.7 tested only the variance -- a 20-fold reduction moved coverage by 0.006 -- and that was read as disproving the coupling. **The mean was never tested**, and it is the half that plausibly matters, since it is also the change that would remove a downward bias.

    **The `theta0`-prior arm is the wrong test and should not be run.** `theta0`'s own prior never changed, so it cannot explain a change in behaviour; and the posterior is not prior-dominated in either run -- width is 0.68 of the prior's 95% width pre-fix and 0.86 post-fix, informative in both. Tightening it would narrow the interval and mechanically improve coverage while explaining nothing.

    **Priority: still the lowest of the open findings.** Overcoverage costs power, not correctness, and the parameter is now unbiased, which is the half that matters for a paper.

    IF THIS IS EVER REVISITED, CLAUDE TO RUN A `b_betatheta` PRIOR **MEAN** ARM, NOT A `theta0` PRIOR ARM. It tests the untested half and would account for the bias improvement and the width increase together.

7.  **`q` (Stage 2 false positives) degrades hard as `K` rises.** Found 29 July 2026 as a side effect of the M-ladder run (`PLAN.md` 13.7), which was not built to look for it. Never investigated beyond the measurement.

    Coverage falls from 0.945 at `M2` (K = 3) to **0.614 at `K30`** (K = 30, same total rows as `M20`). `M20` itself, which holds K = 3 and raises M instead, sits at 0.742. So more PCR replicates make `q` *less* well calibrated, and the effect is larger than the M-driven change in any other item here.

    **Hypothesis, untested:** the same cost-of-identifiability pattern as `beta_theta`. More PCR replicates sharpen the posterior, so if the informative `Beta(1, 20)` prior holds `q` a fixed distance from the truth, sharper intervals show it as worse coverage. If that is what this is, it is not a new bug but a known trade extended to K, and it belongs with the prior-choice decision rather than in the sampler.

    **Alex's response rules out an implementation defect, not the item.** He does not see a cause and reads the `p`/`q` sampler as correctly implemented. That leaves the untested hypothesis above as the live lead -- it does not require the sampler to be wrong, only the prior to be informative relative to how much `K` sharpens the posterior.

    ALEX OR CLAUDE TO TEST THE COST-OF-IDENTIFIABILITY HYPOTHESIS

8.  **Every `rng.h`-based parallel sampler draws an identical stream on every thread, not merely a non-reproducible one.** Found 2 August 2026 while reviewing `522b89e`'s new `BetaThetaWorker`, which calls into this same scheme. `get_rng()` seeds each thread from `seed_seq{base_seed, tid}`, and `tid` comes from `omp_get_thread_num()` -- which returns 0 for every `RcppParallel`/TBB worker thread, on any platform, because none of them ever enters an actual `#pragma omp parallel` region (every such pragma in this codebase is commented out). On this machine specifically it is worse again: `_OPENMP` is not even defined when the package is compiled via its own `Makevars`, confirmed independently three ways. So every thread seeds identically, and their random streams are not merely correlated but literally the same sequence, consumed at different offsets. Verified directly through the package's real build: of 3000 draws across 6 threads, only 872 were distinct.

    **This does not reopen the race *Fixed bugs* 41 closed.** There is no concurrent read-modify-write on shared state; that fix still holds. What it undermines is the claim built on top of it, that per-species draws are then independent. They are not, whenever more than one thread actually runs. Everything reported in `PLAN.md` and the validation article is unaffected, because every run there was pinned to one thread already, for the separate reason of bit-reproducibility. What is affected is any fit run by a user at the package's default (multi-core) thread count, on any sampler that calls into `rng.h` from more than one `RcppParallel` thread.

    **Live scope re-checked 4 August: `sampleB_SoR()` is not currently on the call path.** `43f2342` reverted the live `update_jSDMcoef()` call from `sample_BBsL_parallel()` (the `BBSL_Worker`/`sampleB_SoR()` machinery) back to the serial `sample_BBsL_cpp()`, so this defect does not currently bite there -- it would if that call were switched back. What is live and affected today is `sample_betatheta_cpp_parallel()`'s `BetaThetaWorker`, whose leaf samplers route through `rng.h`.

    **This is the same defect *MEE paper* Alex to-do 8 already proposes fixing**, seeding on the species index rather than the thread. That was filed as a reproducibility improvement; this finding makes it a correctness fix, since keying on `s` sidesteps the broken `tid` computation entirely rather than only making non-reproducibility deterministic. Full verification: `AGENTS.md`.

    ALEX TO DECIDE: fix now, or treat `RCPP_PARALLEL_NUM_THREADS=1` as the interim safety net until *MEE paper* item 8 lands

9.  **`sample_beta_nocov_cpp_TS()` was switched back to the non-thread-safe `sample_beta_cpp()`, reintroducing the exact race `Fixed bugs` 10 closed a week ago.** Introduced `43f2342`, 3 August 2026. `sample_beta_cpp()` draws via `mvrnormArmaQuick()`, which calls `arma::randn()` -- routed to R's single, unsynchronised global RNG, exactly the call *Fixed bugs* 10 replaced with the thread-safe `sample_beta_cpp_TS()`/`mvrnormArmaQuick_TS()` when it closed the identical defect in the same function. This function is the leaf `BetaThetaWorker` calls from every TBB thread, and `sample_betatheta_cpp_parallel()` is on the live default call path, not an unused alternative.

    **Verified directly, not inferred from the diff.** Two fits under the same seed at `RCPP_PARALLEL_NUM_THREADS=6`: `beta_theta_output` differed by up to 3.6. At `RCPP_PARALLEL_NUM_THREADS=1` the same two fits were bit-identical. That is the same diagnostic signature *Fixed bugs* 41 used to confirm its race, and it isolates the defect to concurrency rather than to anything else that changed in `43f2342`.

    **Distinct from item 8 above.** Item 8 is every thread drawing the *same* stream, which is wrong but at least deterministic and does not corrupt shared state. This is a genuine unsynchronised concurrent read-modify-write on R's global RNG -- undefined behaviour, not merely a bad but repeatable answer -- and it is the more urgent of the two.

    ALEX: REVERT `sample_beta_nocov_cpp_TS()` TO CALL `sample_beta_cpp_TS()`, NOT `sample_beta_cpp()`

## **C. Crashes, unreachable code paths, and API bugs (Alex)**

No open items. `thinOutput()` has been fixed by Claude and moved to group A, where it awaits Alex's review; it gets a *Fixed bugs* entry once reviewed. The assorted smaller items this section also held were fixed by Alex and closed as *Fixed bugs* 44 and 47.

## **D. Dead and broken internal code (Alex)**

Ten dead functions were moved to `deprecated/` on 30 July (*Fixed bugs* 37). What remains is below.

**Re-scanned 31 July: 40 dead R functions, up from 38.** `mcmcfun.R` 14, `jsdmfun.R` 12, `RcppExports.R` 11, plus `computeMinESS()` in `R/diagnostics.R`, `thinOutput()` in `R/output.R`, and `.onLoad()` in `R/zzz.R`. The set **grew** because `8f9f315` added three more unused `RcppExports` wrappers, so this cleanup is chasing a moving target while the samplers are being rewritten.

**The `RcppExports` count re-measured 2 August: 11 of 35 wrappers have no caller in `R/`.** It went 11 to 13 as `41abe69` and `46d8804` landed, then to 12 when `sample_z_cpp()` was de-exported (*Fixed bugs* 40), then to 11 when `sampleB_SoR_TS()` was (*Fixed bugs* 42). The moving-target point above is therefore not hypothetical: two commits in one night added two more than this whole item has ever removed. The current 11 are `sample_w_cpp`, `sample_w_cim_cipp`, `sample_betatheta_cpp`, `findClosestPoint`, `dist_matrix`, `gpCovMatrix`, `samplePGvariables`, `convert_to_correlation`, `XsBs`, `XtOmegaX_SoR` and `sample_BBsL_cpp`. One is worth noting rather than batch-deleting: `samplePGvariables` went dead only because `46d8804` replaced it with the parallel version, so it is the serial reference for it.

**Four functions are excluded from all of the below, by decision.** `computePredictiveProbs()`, `partition_r2()`, `returnSpatialEffectMean()` and `plotSpatialEffect()` are dead by the same test as the rest, and were previously listed as a question for Alex. He removed that question in `6722e22` without changing the code, in a commit where he did act on other items, which reads as a decision to keep them. Two have independent reasons to stay: `partition_r2()` relates to the live *MEE paper* item on site variance partitioning, and the `returnSpatialEffectMean()`/`plotSpatialEffect()` pair is the only spatial-field plotting anywhere in the package. `computePredictiveProbs()` looks straightforwardly superseded by `predictNewSites()` and could go whenever Alex says so. **If that reading is wrong, say so and they go with the rest.**

**Timing: this is not urgent and is best done after the sampler rewrite lands.** The payoff is a `R CMD check` NOTE, not a WARNING, and CRAN accepts NOTEs with explanation. Meanwhile the `RcppExports` half needs `src/` edits and the `jsdmfun.R` half needs edits to a file being actively rewritten, so doing either now invites merge conflicts for a cosmetic gain. Alex's profiling note points at replacing the Polya-Gamma sampler, which will change this set again.

1.  **Move the remaining dead functions to `deprecated/`.** About 26 across `R/jsdmfun.R` and `R/mcmcfun.R`, plus `computeMinESS()` in `R/diagnostics.R`, plus 11 unused wrappers in `RcppExports.R`. Excludes the four named above.

    The wrappers need different handling: **do not edit `RcppExports.R`**, it is generated. Remove the `// [[Rcpp::export]]` tag in the C++ and re-run `Rcpp::compileAttributes()`.

    One thing the scan flags that must **not** be deleted: `.onLoad()`, which has zero callers because R itself calls it; removing it would drop the `mc.cores` cap set for CRAN compliance. `thinOutput()` also shows up as dead and is a genuine judgement call rather than an oversight: it is correct and tested as of the `thinOutput()` item in group A, but unexported and with no caller, since the CRAN plan step that wanted it has been superseded by refitting small instead. Delete it or keep it deliberately; do not let it fall out with the batch.

    CLAUDE TO DO AFTER THE SAMPLER REWRITE LANDS

2.  **`globalVariables()` for the data-masked column names.** The `R CMD check` undefined-globals NOTE is down from 84 symbols to 65 as dead code has been removed. What will remain is `dplyr`/`ggplot2` NSE references (`x`, `y`, `Species`, `Min`, `2.5%` and so on), which are false positives and want one `utils::globalVariables()` call in `R/occJSDM-package.R`.

    **Do this last.** Every dead function removed shrinks the list, so enumerating it earlier means writing entries for code about to be deleted. There is no `globalVariables()` anywhere yet, so this sets the convention.

    **Done when** `devtools::check()` reports no NOTE under "checking R code for possible problems", not merely a shorter one.

    CLAUDE TO DO AFTER ITEM 1

3.  **`sample_rnb()` cannot run as written** (`R/jsdmfun.R`). Groundwork for the count-data item, not yet called from anywhere, but it has a scoping bug that will bite when wired up: `r_current <- rnb[s]` reads `rnb` inside the `sapply()` whose result is being assigned to `rnb`, so lookup falls through to the namespace and fails. The current size vector needs to come in as an argument.

    Two more to settle while there: `tune_sd = 5` is a random-walk SD on the *log* scale, so proposals land a factor of `exp(+/-10)` away and acceptance will be near zero (0.1 to 1 is the usual starting range); and the prior terms are stubbed to `0` with the intended `dgamma()` commented out, referencing `prior_shape`/`prior_rate`, which are not defined anywhere. The Metropolis step itself looks right: the `log(r_star) - log(r_current)` Jacobian is the correct correction for a log-scale random walk under a flat prior on `r`.

    ALEX's WORK IN PROGRESS FOR THE COUNTS

4.  **Two `list_jsdmParams` entries do not affect the simulated data, and one of them affects nothing anywhere.** Found 2 August 2026 while commenting that list in `vignettes/simulateOccJSDMData.Rmd`. Both are user-facing: `simulateOccJSDMData()` asks callers to supply them, and the vignette does.

    **`sigma_ts` is wholly inert.** Four occurrences in the whole package, every one of them plumbing: documented in the `@param` at `R/simulateData.R:20`, read into a local at `:59`, passed on at `:85`, received in the signature at `R/jsdmfun.R:909`. No function body references it. It is read, passed, received and discarded.

    **`sigma_bs` generates nothing, but is not simply dead.** In the simulator it appears only in the signature and in the returned `trueParams`; the residual spatial term it would scale is set to an exact zero matrix (`Bst <- matrix(0, S, ps)`), so no draw ever uses it. It *is* live on the fitting side, where `sigma_bs^2` sets a prior variance block. So a caller supplies it as a true value, it generates none of the data, and the sampler then estimates a quantity by that name -- which is exactly why the simulation study excludes `sigma_bs` from its coverage checks (`PLAN.md` 5.3, measured true 0.5 against a posterior mean of \~1.6, 0/8 coverage).

    **Why this is worth a decision rather than a deletion.** These are arguments in an exported function's interface, so removing them is a breaking change, and `sigma_bs` at least has a real meaning on the fitting side that a future simulator could honour by drawing `Bst` properly. The options are: drop `sigma_ts` outright, since nothing anywhere reads it; and for `sigma_bs` either make the simulator use it, or keep it and document in `@param` that it is a fitting-side prior rather than a generating parameter.

    Also stale as a result: `R/simulateData.R:20`'s `@param` lists `sigma_ts` as though it were live, and the vignette prose above the code chunk groups both with the real variance components. The vignette's code comments now say what each one actually does; the roxygen does not.

    ALEX TO REVIEW

## **E. Draft of beta version listserv announcement (Doug)**

**No `NEWS.md` for the beta, decided 10 August 2026.** occJSDM gets one when it is released as 0.2.0. That makes this announcement the only place a user is told what is currently broken, since `TODO.md` and `AGENTS.md` are both excluded from the build and from the site, and the validation article is pkgdown-only and unpublished. The "current limitations" block below is therefore load-bearing rather than throat-clearing, and it has to stay in step with group B. Rationale in `AGENTS.md`, CRAN plan item 22.

1.  Listserv announcement (beta release), drafted July 20 2026; limitations block added 10 August 2026:

    > Subject: New R package (beta) - occJSDM, a combined occupancy and joint species distribution model
    >
    > Hi all,
    >
    > Announcing **occJSDM**, an R package for combining occupancy and joint species distribution modelling (<https://github.com/AlexDiana/occJSDM>).
    >
    > occJSDM extends the occPlus two-stage eDNA occupancy model of Ji et al. (2025, *Ecology Letters*, <doi:10.1111/ele.70302>) by adding a JSDM layer. Unusually for an occupancy model, false positives are estimated explicitly at both the field and lab stages and separately for each species and each primer.
    >
    > Note this is still **beta software**. Feedback, feature requests, and bug reports are very welcome.
    >
    > Highlights:
    >
    > - Occupancy modelling: Accounts for both false-negative and false-positive error at two stages (field and lab), per species. Stage 1: estimates species eDNA collection probability in the field, given true eDNA presence at the site, and contamination probability, given true eDNA absence at the site. Stage 2: estimates species eDNA detection probability in the lab (i.e. successful DNA extraction, PCR, and sequencing), given successful eDNA collection in Stage 1, and contamination probability, given eDNA non-collection in Stage 1. In datasets where multiple primers have been used, each species' detection probability is estimated per primer (allowing one to compare each primer's efficiency for each species), while species occupancies are estimated using information across all primers. Both environmental and detection covariates are supported.
    > - JSDM: Integrates the occupancy model with a JSDM: species fit jointly with nonlinear response curves and latent-factor residual correlations. The JSDM optionally supports species traits shaping occupancy responses (trait x env interactions, aka 'fourth-corner analyses') and spatial autocorrelation (GP kernel) across sites. Occupancies predicted at unsampled sites.
    > - occJSDM not only fits a two-stage occupancy model (both field and PCR replicates required), but if given simpler study designs, can collapse to a classical occupancy model (field replicates only) or to a pure JSDM (no replicates).
    > - MCMC fitting with diagnostics, variance partitioning, ordination, and pairwise residual correlation outputs built in.
    > - occJSDM leverages the taxonomic breadth of eDNA datasets by using ordination (each site's position on the latent axes, and each species' loadings on those axes) to predict species occupancies. Thus, each species' predicted occupancy at a site is informed by the estimated occupancies of the other species at that site, thereby using co-occurrence structure. We also allow species to borrow strength from other species sharing similar traits, including inferred traits, in contrast to the classical approach of having rare species borrow strength from abundant species, as is used in multi-species occupancy models.
    >
    > Current limitations, all three being worked on:
    >
    > - **Spatial field.** The Gaussian-process range parameter is not currently recovered, which biases the spatial term of any fit that uses it. We suggest leaving the spatial field off (`useSpatField = FALSE`) in this release.
    > - **Residual species correlations.** These are currently biased toward the extremes, by up to 0.6 in our simulations. Read `returnResidualCorrelationMatrix()` and `plotResidualCorrelationMatrix()` for the sign and structure of co-occurrence rather than for calibrated magnitudes.
    > - **Collection-covariate slopes.** Credible intervals on the Stage 1 collection covariates (`returnCollectionCovariates()`, `plotCollectionCovariates()`) are narrower than they should be: about 77% coverage against a nominal 95% in simulation, and it worsens as replication increases.
    >
    > Vignettes and articles included on data simulation, model fitting/interpretation, and model performance.

# **MEE paper**

## A. Alex to dos

1.  Trait matrix currently not allowing for categorical variables

2.  Design better model selection criterion (the one currently implemented, which is the same as HMSC, tend to overfit).

3.  ability to analyse count data

4.  scenario for source-sink inference (sites where env covariate coeffcients are negative but spatial covariate coefficients are positive), using an explicit source-sink simulation

5.  remove effect of space on environmental covariates. remove the effect of unobserved environmental on observed environmental covariates and space. thus, adding factors (unobserved env covariates) doesn't change the effect of observed env covariates

6.  **site** variance partitioning to complement the **species** variation partitioning (see Leibold et al., Cai et al.). Consider calling it **variation** partitioning.

7.  **Performance of `runOccJSDM()`**

ALEX NOTE: Most of the MCMC steps have now been parallelised, with the only exception of sample_U_cpp. The rest is mostly stuff. It would also worth investigating a faster to compute the variancePartitioning or the WAIC.

**Verified against the code 4 August 2026, item by item, not from the note or the commit messages.** Live (called from the default `runOccJSDM()` path, uncommented): `sample_z_cpp_parallel()`, `sample_w_cim_cipp_parallel()`, `sample_betatheta_cpp_parallel()`, `sample_pq_cpp_parallel()`, `samplePGvariables_parallel()`. **`sample_U_cpp()` is correctly the one exception**, exactly as the note says. **One correction to the note: `sample_BBsL` is not currently parallelised.** `43f2342` (3 August) switched the live call in `update_jSDMcoef()` from `sample_BBsL_parallel()` back to the serial `sample_BBsL_cpp()`, with the parallel version left commented out beside it -- an apparent revert, cause not recorded. That switch is also what reintroduced group B item 9's RNG race, in the function `sample_BBsL_cpp()`'s neighbour calls into.

**Below, each item A-H is marked against the current code, not against what has generically "been parallelised".** None of them describe MCMC-step threading; they are chain-level parallelism (A) and serial R/C++ inefficiencies (B-H), a different axis of work from the RcppParallel conversions above. **None are done.** Two are worth a closer look regardless: C's specific claim no longer matches the code, and H has an unused, half-built attempt at the fix already sitting in `src/functions.cpp`.

**Profiled by Alex, 31 July 2026, which answers the "nothing here has been profiled" caveat this list used to carry.** Comparing each MCMC step with `microbenchmark()`:

- **`sample_betatheta_cpp_parallel()` is the slowest step**, decisively.
- **At the time this was measured, the parallelisation achieved little speedup and was inert on macOS**, because it used OpenMP rather than RcppParallel. **No longer true of the current code**: `sample_betatheta_cpp_parallel()` now runs via `BetaThetaWorker`+`RcppParallel::parallelFor`, which is genuinely parallel on this machine regardless of the OpenMP question below -- confirmed by the \~5x `cpu/wall` speedups measured throughout the simulation study, all of which go through this exact code path. The `#pragma omp` sections this bullet originally described are all commented out now; nothing in the live sampler still depends on OpenMP compiling.
- **Within that step, `sample_Omega_cpp()` dominates**: it draws `N x S` Polya-Gamma variables per iteration.
- **Alex's suggestion: consider an alternative Polya-Gamma sampler.** That is the lever with the best expected return, and it is a different kind of work from the parallelisation items below, which redistribute the same cost rather than reducing it.

**This reorders the list.** Items A and B parallelise around a step whose cost is dominated by PG sampling; making the PG draw cheaper would benefit every configuration, including single-core, where the parallelisation currently does nothing (macOS is no longer such a case, per above). Worth settling the PG question before investing in either.

The items below are ordered by expected speedup per unit of effort as originally written.

**"OpenMP is inert on this machine" is still true, and still worth knowing, but it no longer means "nothing runs in parallel here."** Measured 29 July 2026: `R CMD config SHLIB_OPENMP_CXXFLAGS` is *empty*, so the `$(SHLIB_OPENMP_CXXFLAGS)` in `src/Makevars` expands to nothing and every literal `#pragma omp` compiles to a no-op. But every sampler actually in use now runs through `RcppParallel`/TBB instead, which does not depend on that flag at all and is genuinely multi-threaded here -- see the `BetaThetaWorker` note above. Do not use this fact to conclude the parallelisation items below are moot; check the specific mechanism a given step uses.

Consequences that still hold: the thread-safety bug of Fixed bugs 26 could never have manifested locally via the literal OpenMP pragmas, and any timing measured here for a step that genuinely still uses `#pragma omp` says nothing about a Linux build where that compiles. Worth confirming what Alex's machine and CRAN's check farm do before investing in items A and B below -- the payoff differs completely between the two cases.

A.  **Parallelise over chains -- but use a PSOCK cluster, not `mclapply()`. NOT DONE, verified 4 August 2026.** `R/runOccJSDM.R`'s `for (chain in 1:nchain)` loop is still plain and serial; no `makeCluster()`, `parLapply()`, `mclapply()` or `cores` argument exists anywhere in the file. This is the simplest available speedup and remains fully open. The `for (chain in 1:nchain)` loop is serial and embarrassingly parallel; each chain touches only its own `*_output_chain` arrays, so running the chains in separate *processes* gives close to an `nchain`-fold speedup and sidesteps every RNG thread-safety problem above entirely (each process has its own RNG state). Portability constraints, though:

    - **`parallel::mclapply()` is not an option.** It is fork-based, and R ships a Windows stub whose body is `if (cores > 1L) stop("'mc.cores' > 1 is not supported on Windows")` -- a hard error, not a fallback. That would break the package on Alex's machine and on CRAN's Windows check. Even on Linux/macOS, forking a session with a live threaded-BLAS pool (OpenBLAS on most Linux distros, Accelerate on macOS) is a well-known deadlock source.
    - **`parallel::makeCluster()` (PSOCK) + `parLapply()` works on all three platforms**, and `parallel::clusterSetRNGStream()` gives reproducible, independent L'Ecuyer streams per chain, which is strictly better than the current situation. The cost is that workers are fresh R sessions (`clusterEvalQ(cl, library(occJSDM))`, plus serialising the data out and the per-chain arrays back) -- negligible against a multi-minute MCMC, though the return trip is not free given how large `Bs_output`/`U_output` are (see group D.7 on that size).
    - **Make it opt-in**: a `cores = 1L` argument, serial by default, so existing behaviour is unchanged and CRAN's two-core limit for examples/tests/vignettes is respected.

    Structurally this needs the chain body extracted into a function returning its own `*_output_chain` arrays -- mechanical, since the loop already writes only to per-chain objects. The one piece of genuinely shared state is the WAIC accumulator, which currently streams across chains sequentially via the single `currentWAICiter` counter introduced by the fix in Fixed bugs 9; parallelising forces one accumulator per chain plus a merge (`mean_lik` is a plain mean; `M2` merges with the standard parallel-variance formula), and the `if (numIters != (currentWAICiter - 1)) stop(...)` guard will need to be restated in terms of the summed per-chain counts. The `z_output_mean` / `psi_output_mean` / `w_output_mean` / `theta_output_mean` running means merge by simple addition of the per-chain partials.

B.  **Drop the `options(mc.cores = ...)` call in `.onLoad()`. NOT DONE, verified 4 August 2026.** `R/zzz.R` still sets a global option at load time, which changes the user's session state and affects every other package that reads `mc.cores` -- CRAN policy is explicitly against this. It should be replaced by the `cores` argument in A. (Also `parallel::detectCores()` can return `NA`, which `min(2L, NA)` happily propagates.)

C.  **`computePsiCoef()` is called three times per iteration. Re-checked 4 August 2026: the described redundancy does not match the current code.** The three calls (now `R/jsdmfun.R` around lines 1530, 1604 and 1627) genuinely need different inputs: `B0`/`G`/`A`/`C`/`Bt` are resampled between call 1 and call 2 (via `sample_BBsL_cpp()`, `sample_GC()`, `sample_A()`), and `U` is resampled between call 2 and call 3. Whether that was always true and the original item's premise was wrong, or the sampling order changed since it was written, is not established -- but as the code stands today there is no free two-thirds saving here to take. Leave open pending someone re-deriving whether any of the three inputs really are unchanged from the previous call.

D.  **Precompute the constant parts of the `c_imk` update. NOT DONE, verified 4 August 2026.** `y_pos <- (y > 0)` still sits inside the chain loop rather than above it, and `w_all <- w[idx_w_k, , drop = FALSE]` is still a separate gather. Hoist `y_pos` above the loop, and consider folding the `w_all` gather into `sample_pq_cpp_parallel()` so the copy never reaches R.

E.  **Make the WAIC accumulation optional, and cheaper.** Not verified changed: `computeModelLoglikJSDM_cpp()`/`FirstStage_cpp()`/`SecondStage_cpp()` still appear to call `R::dbinom` generically rather than the closed-form binary case. A fast C++ alternative to `R::dbinom` exists elsewhere in `src/functions.cpp` (used for a different sampler), so the pattern is available in the codebase if not yet applied here. Also add a `computeWAIC = TRUE/FALSE` argument for users who are not doing model comparison.

F.  **Vectorise the starting-value loops. NOT DONE, verified 4 August 2026.** The triple-nested `for` loops initialising `w` and `z` are unchanged. Both reduce to grouped "any positive" reductions: `w` from `rowsum(y > 0, idx_w_k) > 0`, `z` from `rowsum(w, idx_z_w) > 0`.

G.  **Do not allocate posterior arrays that are never filled. NOT DONE, verified 4 August 2026.** No `keep =` argument or equivalent exists. `Bs_output` (`ps x S x niter x nchain`) and `U_output` (`n x d x niter x nchain`) are still the two largest components of the 62 MB `sampleresults.rda`; a `keep =` argument selecting which blocks to retain, or storing the latent-factor blocks pre-thinned, would address both the runtime allocation and the CRAN size blocker.

H.  **Reduce the repeated `arma::inv()` calls in the samplers. NOT DONE, but a matching implementation already exists unused.** `sampleB()` and `sampleBuniv()` (`src/jsdm.cpp`) are unchanged and still call `arma::inv(B)` twice plus `arma::inv(arma::trimatl(L))` per species/site, on a matrix that is diagonal in every caller. **`sample_beta_cpp_TS_opt()`** (`src/functions.cpp`) already implements exactly the fix this item asks for -- takes the precision `invB` directly as an argument and draws via `mvrnorm_from_chol_prec()`, a triangular solve against a standard normal, matching `sampleB_SoR()`'s pattern -- but it has **no callers anywhere**. Either wire it in (and apply the same pattern to `sampleB()`/`sampleBuniv()`), or delete it if it was abandoned for a reason not recorded here.

<!-- -->

8.  **Make the parallel sampler reproducible at any thread count, by keying the draws on species rather than on thread.** This is a design change in your sampler, which is why it is here rather than being applied.

    **Where it stands now, and this got worse on 2 August.** `sampleB_SoR()` draws from `rnorm()` in `src/rng.h`, whose per-thread `mt19937` is meant to derive from `(base_seed, tid)`. `tid` comes from `omp_get_thread_num()`, which returns 0 for every `RcppParallel`/TBB thread regardless of platform, since none of them ever enters a real OpenMP parallel region. So every thread seeds identically, and their streams are not merely uncoordinated but literally the same sequence. Verified directly: 3000 draws across 6 threads produced only 872 distinct values. This is the item at the bottom of group B; it is a correctness gap at the package's default thread count, not only a reproducibility one, and it applies to every sampler that calls into `rng.h` from more than one thread, not only `sampleB_SoR()`.

    **The change.** Build a generator inside each worker's `operator()` seeded from `(base_seed, s)` for species `s`, and thread it into the sampler as an argument instead of having the sampler reach for a thread-local engine keyed on the broken `tid`. Species `s` then gets the same draws whichever thread runs it, which fixes both problems at once: the stream is no longer shared across threads, and it no longer depends on work-stealing assignment. `src/rng.h` already has the pieces -- the base seed and the `seed_seq` construction -- and the generation-counter and `dist.reset()` traps documented at the top of that file apply unchanged.

    **Why it is worth doing rather than living with.** It is no longer only about bit-reproducibility between runs. At the package's default thread count, species handled by different threads can draw from the same random-number stream, at whatever offset each thread happens to have reached -- a correctness problem for any user who has not manually set `RCPP_PARALLEL_NUM_THREADS=1`. The simulation study itself is unaffected, since every run has been pinned to one thread throughout, for the separate reason of exact reproducibility; but that was never documented as a requirement for statistical validity until now, only for comparing runs against each other.

    **Two things it also unblocks.** `test-regression-bugs.R` currently pins its reproducibility test to one thread and carries a skipping test naming this gap; both can go when this lands. And tier 1's "structural assertions only" rule can finally be revisited, which the testing item under *Doug to dos* already flags as waiting on reproducibility being confirmed on a multi-threaded platform.

    ALEX TO DECIDE AND CLAUDE TO IMPLEMENT (touches `src/jsdm.cpp` and `src/rng.h`)

## B. Doug to dos

1.  **reproduce all Ecoletts results as a test of the package and decide whether to include in repo**

2.  **extensive testing on simulated datasets** -- **suite built and the R = 100 study run; three things remain.** What exists is summarised under *Completed* below; the authoritative specification and the results table are in `dev/simstudy/PLAN.md`.

    (a) **Re-run once the `sample_ls()`, `reparamFactorModel()` and `beta_theta` slope items in group B are fixed.** This is the evidence the fixes worked. Without it they rest on the same code-reading that this exercise showed to be unreliable. Still outstanding: none of the three is fixed.

    **A re-run did happen on 10 August 2026, but not this one.** Its purpose was different: six commits had touched `R/` and `src/` since the 2 August study, including `522b89e`'s new parallel sampler, so the published numbers described code that no longer existed. Result: **nothing moved.** Zero of 83 scenario-by-block cells shifted beyond 2 SE, and every per-block mean coverage change was under 0.003 against a measurement SE of 0.022. The old numbers were stale in provenance, not in fact. Scope: every fit runs at one thread, where the new parallel worker reduces to serial, so this says nothing about the multi-threaded path where group B items 8 and 9 bite.

    **That makes the fix-verification re-run cheaper to read, not redundant.** There is now a clean baseline measured on current code, so the next comparison isolates the fixes instead of confounding them with six commits of drift. Name the production grid explicitly -- a bare invocation takes every cell defined in `helper-simstudy.R`, not the ten production ones, which is hours of wasted compute and has happened once already:

    ````         
    ```
    Rscript dev/simstudy/run_study.R --R=100 --cores=5 --caffeinate \
      --scenarios=base,binary,d_overfit,d_underfit,low_information,occupancy,primers_3,spatial_isolated,species_20,traits_isolated
    ```
    ````

    (b) **Decide the replicate count for the paper.** R = 100 was chosen to *detect* defects and did so decisively. Asserting *nominal* coverage in print is a claim about the absence of a small deviation and wants R = 200-500 (`PLAN.md` §9). The runner takes `R` as an argument.

    (c) ~~**Decide how the results are presented**~~ **SETTLED 10 August 2026.** The write-up is the pkgdown article at `vignettes/articles/validation.Rmd`, and it is now *generated*: every table, figure and number renders from `dev/simstudy/validation-data.rds` rather than being typed in. So presentation is no longer a standing decision -- a re-run plus `export_validation_data.R` refreshes the whole document, and the article cannot silently disagree with the data it describes. Not published yet; see item 3.

    **One constraint carried from the bug list:** `l_s` is excluded from coverage checks because it is not recoverable while the `sample_ls()` item in group B is open, so no cell of the study speaks to spatial range. Two earlier constraints have since lapsed -- `sigma_h` is now sampled (Fixed bugs 24) and the OpenMP RNG race is closed (Fixed bugs 26), so tier 1's "structural assertions only" rule can be revisited once reproducibility is confirmed on a multi-threaded platform.

3.  ~~**Stand up a pkgdown site.**~~ **BUILT 2 August 2026** (`b34b36a`), but deliberately **not published**. `_pkgdown.yml`, the validation article at `vignettes/articles/validation.Rmd`, and `URL`/`BugReports` in `DESCRIPTION` are on `main`. That closes CRAN plan item 10.

    `.github/workflows/pkgdown.yaml` carries only a `workflow_dispatch` trigger, so nothing builds or deploys on push. **How to rebuild locally and how to publish for real are both in `AGENTS.md`, "The documentation site".** Short version: `pkgdown::build_site()` writes to a gitignored `docs/`; publishing needs the workflow run by hand *and* Pages repointed from `main` to `gh-pages`, and neither alone is enough.

    **ALEX: the Pages repoint needs admin**, which Doug does not have. Until then `alexdiana.github.io/occJSDM` serves the README via Jekyll rather than the pkgdown site.

    **One thing to settle before publishing, not two. Corrected 10 August 2026.** This item used to say the first build would fail on functions that error unconditionally, `predictNewSites()` among them. That is wrong twice over: `predictNewSites()` was fixed as *Fixed bugs* 34, and pkgdown does not evaluate `\dontrun{}` blocks, which is 22 of the 24 example blocks in `man/`. The two live ones are `str()`, `head()` and one `plotDetectionRates()` call on shipped data. **There is no example-driven build blocker.**

    **What is still open is the judgement call.** Publishing while the `sample_ls()`, `reparamFactorModel()` and `beta_theta` slope items stand means the site documents functions whose output is currently biased or whose intervals are overconfident. The beta's disclosure now lives in the listserv announcement (group E), which the site does not carry, so publishing puts the documentation somewhere the caveats are not. Either say so on the site or accept the gap knowingly.

# **Future versions**

1.  use spike-in to estimate abundance change (i.e. eDNAPlus)

2.  model selection of environmental, spatial covariates via regularisation/shrinkage, which would be useful with e.g. geospatial foundation model embeddings as env covariates

3.  Let `simulateOccJSDMData()` generate nonlinear responses to environmental covariates, so the GAM/spline fitting path can be checked against a known truth as well as compared against a linear fit. The capability already exists and is simply unreachable: `simulateData()` takes a `usingSplines` argument and spline-expands the covariate matrix when it is true, but `R/simulateData.R:86` hard-codes `usingSplines = F` and no element of the three parameter lists reaches it. By analogy with `useSpatField`, the switch belongs in `list_jsdmParams`. Demonstrate it afterwards in `vignettes/simulateOccJSDMData.Rmd`. Two things to know before scoping this: `listParams$splineVars` is live on the fitting side (`R/runOccJSDM.R:683`) but appears nowhere in the roxygen or `man/runOccJSDM.Rd`, so the feature currently ships undiscoverable as well as unvalidated; and "known truth" for a spline is the fitted response *curve*, not the basis coefficients, so this cannot join the coverage study element-wise the way `B0` and `B` do until a statistic is defined for it.

4.  parallelisation for speedup

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

10. ~~**Non-thread-safe RNG inside the OpenMP loops.**~~ **PARTIALLY FIXED, THEN REGRESSED 3 August -- see group B item 9.** `sample_beta_nocov_cpp_TS()` (`src/functions.cpp:308`) was switched to call the thread-safe `sample_beta_cpp_TS()` rather than the non-TS `sample_beta_cpp()`, which closed the race in `sample_betatheta_cpp_parallel()` (`src/functions.cpp:607`) -- the thread-safe path the audit described as "half-finished" was wired up on that side. **The `samplePGvariables()` path (`src/jsdm.cpp:398`) is still affected** via `randinvg()`'s `R::rnorm` call; that residual was tracked as a group B item and is now closed as *Fixed bugs* 28. *(Also not recorded at the time; recovered 27 July.)* **This entry no longer describes the code.** `43f2342` reverted `sample_beta_nocov_cpp_TS()` back to calling `sample_beta_cpp()`, reopening the exact race this entry records having closed. Verified empirically: two same-seed fits differ by up to 3.6 on `beta_theta_output` at 6 threads, bit-identical at 1. Do not read this entry as current status; read group B item 9.

**Post-audit fixes:**

11. ~~**Collection covariates are grouped by `Sample` alone.**~~ **FIXED.** Now groups by `c("Site", "Sample")` to match the site-then-sample ordering of `X_theta`. `R/runOccJSDM.R:682`.

12. ~~**The GP length-scale is never updated.**~~ **FIXED.** Enabled the `sample_ls()` call in `update_jSDMcoef()`. `R/jsdmfun.R:1063`.

13. ~~**`transformCovariatesMatrix()` never applies the stored factor levels.**~~ **FIXED.** Now correctly uses `df[[col]]` instead of `df[,col]` to re-level the column. `R/runOccJSDM.R:38`.

14. ~~**`listPriors$prior_beta_psi` / `prior_beta_psi_sd` have no effect.**~~ **FIXED.** Removed dead references and cleaned up the prior handling. `R/runOccJSDM.R:738-739`.

15. ~~**`summarisedLatentPresences = FALSE` errors on two-stage models.**~~ **FIXED.** Now correctly allocates and writes `w_output_chain` and `theta_output_chain`. `R/runOccJSDM.R:1152`.

16. ~~**`thinOutput()` cannot run at all.**~~ **PARTIALLY FIXED.** The function still exists and now runs: `niter` is read from `fitModel$results_output$jsdm_output$B0_output` instead of the long-gone `beta_ord_output`, and the `thin` argument is honoured rather than hard-coding `by = 5`. Two of the three original defects survived this fix, plus a third nobody had spotted. All were fixed by Claude on 4 August 2026 and are awaiting review as the `thinOutput()` item in group A; read that item, not this one, for current status.

17. ~~**`computeDiagnostics()` runs on `psi_output`.**~~ **FIXED.** Added `psi_output` to the skip list to avoid meaningless diagnostics. `R/diagnostics.R:403-404`.

18. ~~**`predictNewSites()` does not honour its documented no-op behaviour.**~~ **PARTIALLY FIXED -- this entry previously overstated what landed; corrected 27 July 2026.** The gating on the presence of covariates was done (`R/output.R:1457,1470`). The **`NULL` defaults were not**: `formals(predictNewSites)` shows `X_psi` and `X_s` still have no defaults at all, so `predictNewSites(fit, X_psi = X)` still fails on the missing `X_s` promise even when `useSpatial = FALSE`. **That residual is now closed: see *Fixed bugs* 34.**

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

    **Confirmed by the 29 July re-run**: `beta_theta` improved in every cell, to 0.709-0.771, with two-thirds of the bias removed. It is still well short of nominal, so a second cause remains -- tracked as the slope-overconfidence item in group B, and since narrowed to the collection-covariate slopes.

26. ~~**Non-thread-safe RNG in the hottest OpenMP loop.**~~ **FIXED.** `randinvg()` (`src/jsdm.cpp:86`) now draws from the `thread_local` `rnorm()` rather than `R::rnorm`; the old line is commented out beside it. That closes the last hole on the `samplePGvariables()` path, whose Polya-Gamma helpers were already converted in `53c38f1`.

    **Two consequences worth acting on.** A fixed seed should now reproduce on Linux and Windows, which (i) removes the reason `sampleresults.rda` had to be refit on macOS or an OpenMP-disabled build, and (ii) unblocks CRAN item 16. It also means the tier-1 constraint in `dev/simstudy/PLAN.md` §5.1 -- structural assertions only, never numeric equality -- can be revisited; it was adopted solely because of this race. Verify reproducibility on a multi-threaded platform before relying on any of that.

27. ~~**Stage 2 hyperparameters documented as settable but never read.**~~ **FIXED.** `a_p`, `b_p`, `a_q` and `b_q` now read from `listPriors` (`R/runOccJSDM.R:751-754`), matching the behaviour `@param listPriors` already claimed, and the documented defaults were corrected to match the code (5/1 and 1/20). A user running a low-detection study can now override the prior without editing the package.

    **Only the wiring is closed.** What the default *should be* remains open and is a design decision, not a defect. It was tracked as a group B item until Alex removed it in `093f2bb`; if that removal meant the decision is made, the chosen values should be recorded here, and if not the item needs restoring.

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

36. ~~**Both exported GAM functions failed for any user with a categorical occupancy covariate.**~~ **FIXED 30 July 2026** (Claude; `R/occJSDM-package.R`, `NAMESPACE`, commit `032dcff`). Logged 1 August: it was completed but never given a *Fixed bugs* entry, so it was invisible in this file.

    **The defect.** `tidyr::pivot_longer()` was called but never imported, and `tidyr` was in `DESCRIPTION` `Imports:` with nothing imported from it in `NAMESPACE`. It is reached in the categorical-covariate branch of `returnCovariateEffect()` and `plotCovariateEffect()`, so both would fail with "could not find function" from a properly installed namespace. `stats::setNames` and `stats::rnbinom` were missing the same way.

    **Why it survived the test suite.** `devtools::load_all()` resolves unimported symbols through the global environment, so a functional test passes whether or not the import exists. The regression test therefore asserts on the imports environment directly, which holds under both `load_all()` and `R CMD check`.

    **Not cosmetic, despite arriving as an `R CMD check` NOTE.** The "checking dependencies in R code" NOTE is now gone entirely, 3 notes to 2. `dnbinom` and `bs` remain on the undefined-globals list, correctly: they are reached only by dead code, and importing `bs` would add `splines` to `DESCRIPTION` to support code that should not ship.

    **Reviewed by Alex, 2 August 2026: "Agree with the fix."**

37. ~~**Ten dead functions moved to `deprecated/`.**~~ **DONE 30 July 2026** (Claude; commits `7bae018`, `f2e2701`). Logged 1 August so the deprecation has a stable anchor; it was previously cited only by section-A position, which has since been reused.

    `sample_z`, `sample_w`, `sample_cimk`, `sample_betatheta` from `R/jsdmfun.R`, plus six dead samplers and plots from `R/mcmcfun.R` and `R/output.R`. All were unreachable from any exported entry point. Roughly 26 more remain, tracked in group D, deliberately deferred until the sampler rewrite lands so the two do not collide.

    **Reviewed by Alex, 2 August 2026: "Ok to deprecate these functions."**

38. ~~**`main` would not load at all, and two fixes were needed to restore it.**~~ **FIXED 31 July 2026** (Claude; `src/Makevars`, `src/Makevars.win`, `DESCRIPTION`, `NAMESPACE`, `R/occJSDM-package.R`, `R/runOccJSDM.R`). Closed 2 August 2026. After `8f9f315` the test suite went from 167 passing to 25 passing, 3 failures and 33 errors, because the built `.so` had an undefined `RcppParallel::tbbParallelFor` and failed to load, taking every function in the package with it.

    **Change 1, the linking.** `RcppParallel` was in `LinkingTo`, which makes the headers visible but does not link the libraries. Added the documented `RcppParallel::RcppParallelLibs()` call to both `Makevars` files, **and** added `RcppParallel` to `Imports` with an `importFrom`. **Both halves are needed**: the Makevars line alone fixes the missing symbol but the load still fails, because on macOS `libtbb` is referenced through `@rpath` and `devtools` copies the `.so` to a temp directory. Importing the package makes its namespace load first and set up the TBB paths. The reasoning is in a comment in `src/Makevars` so the next person does not stop at the first fix and conclude it did not work.

    **Change 2, the debugging scaffold at `R/runOccJSDM.R:404`.** It had live assignments overwriting every argument with values referencing `occ_data_effort`, a dataset not in the package, plus a live `summarisedLatentPresences`. That block was fully commented before `8f9f315`. **Re-commented rather than deleted**, since Alex evidently uses it, with a note saying it must stay commented and what happened when it did not.

    **Two things left alone at the time, both since resolved by Alex:**

    - `verbose` in `computeNewOutputs()` defaults to `T`, so `predictNewSites()` still prints one line per species unless a caller opts out. It is suppressible now, which was the harder half, but the default behaviour, the test-output noise and the CRAN-reviewer exposure are all unchanged. **Alex: "We can leave verbose on, people won't think of turning it on."**
    - `sample_z_cpp_parallel()` was exported but not yet called at the time; `runOccJSDM` still used `sample_z_cpp`, so only `sample_w_cim_cipp_parallel()` was wired in and half the parallelisation work was unreachable. **Alex: "sample_z_cpp_parallel now used."** Confirmed: it is called at `R/runOccJSDM.R:1113` as of `41abe69`, which leaves `sample_z_cpp()` itself with no callers -- whether *that* should now be deprecated is a new open item in group A.

    The `verbose`-default choice noted above is the same decision closed separately as *Fixed bugs* 39.

39. ~~**`computeNewOutputs()` prints to stdout on every call and cannot be silenced.**~~ **FIXED** (Alex added a `verbose` argument to `predictNewSites()`, threaded through to the `Rcout` call in `src/jsdm.cpp`). Closed 2 August 2026.

    **The original defect.** `src/jsdm.cpp` ran `Rcout << "Computing species ..."` inside the species loop unconditionally, so every `predictNewSites()` call printed one line per species, and `suppressMessages()` did not catch it because `Rcout` is stdout rather than R's condition system.

    **Verified live, not just from the response.** `verbose` reaches the C++ `if(verbose)` gate around the `Rcout` call (`src/jsdm.cpp:482`) via `computeNewOutputs()` (`R/output.R:1683`). Ran both ways on a fitted model: `verbose = FALSE` suppresses all four per-species lines; `suppressMessages()` alone, with `verbose` left at its default, still lets all four through -- expected, not a gap, given the default chosen below.

    **The item's fix spec asked for a default of `FALSE`; Alex chose `T` instead, deliberately.** `predictNewSites()` still prints by default unless a caller opts out, which is not what was originally asked for -- but it is a design decision, not an unfinished fix. Same decision already recorded under *Fixed bugs* 38 (the RcppParallel-linking fix): **"We can leave verbose on, people won't think of turning it on."** **Closed by decision, not by elimination.**

40. ~~**`sample_z_cpp()` was exported but had no callers.**~~ **DE-EXPORTED 2 August 2026** (Claude; `src/functions.cpp`, `R/RcppExports.R`, `src/RcppExports.cpp`).

    `41abe69` wired `sample_z_cpp_parallel()` into `runOccJSDM()`, leaving the serial version reachable only as an unused wrapper. Tag removed and `Rcpp::compileAttributes()` re-run; the C++ body stays in `src/` as the serial reference implementation of the parallel one. Detail and decision provenance in `AGENTS.md`, "Detail behind Fixed bugs 40-42".

41. ~~**`BBSL_Worker` called the non-thread-safe `sampleB_SoR()`, racing on R's RNG from every TBB thread.**~~ **FIXED 2 August 2026** (Claude; `src/jsdm.cpp`). Introduced by `41abe69`; found the same day while running the `continuous` arm for the `B0` item in group B.

    `sampleB_SoR()`'s `arma::randn()` reached R's global RNG from every TBB worker thread, permitting duplicate and torn draws, so the chain could fail to target the intended posterior. It now draws from `rnorm()` in `src/rng.h` instead. The race is closed and results are statistically valid again.

    **Bit-reproducibility above one thread is not restored**, because TBB work-stealing varies which thread draws for which species even at a fixed thread count. Set `RCPP_PARALLEL_NUM_THREADS=1` whenever a run has to be reproducible; the suite is 248 passing single-threaded, and two tests in `test-regression-bugs.R` are pinned to one thread for this reason, with a companion skip naming the gap. Closing it is *MEE paper* Alex to-do 8. Measurements and full mechanism in `AGENTS.md`, "Detail behind Fixed bugs 40-42".

42. ~~**`sampleB_SoR_TS()` was exported, had no callers, and its name invited the exact error it was written to prevent.**~~ **DEPRECATED 2 August 2026** (Claude; `src/jsdm.cpp`, `deprecated/jsdm-sampleB_SoR_TS.cpp`, `R/RcppExports.R`, `src/RcppExports.cpp`). Closes the decision left open at the end of *Fixed bugs* 41.

    Written as the thread-safe variant of `sampleB_SoR()` for the race above, and left purposeless when that race was closed another way. Worse than merely dead: it seeds from OS entropy, so anything built on it would ignore `set.seed()` entirely, despite the `_TS` name reading as the endorsed choice. De-exported and moved to `deprecated/`, not deleted. Group D's dead-wrapper count drops to 11 of 35, re-measured rather than decremented; the suite is unchanged at 248 passing. Detail in `AGENTS.md`, "Detail behind Fixed bugs 40-42".

43. ~~**The `runOccJSDM()` debugging scaffold went live for the fourth time.**~~ **FIXED 2 August 2026** (Claude; `R/runOccJSDM.R`, `tests/testthat/test-regression-bugs.R`). First occurrence `8f9f315` (*Fixed bugs* 38); second `46d8804`, fixed in `e5d0105`; third some time before `11981a1`, which Alex fixed himself; fourth `522b89e`, minutes later, apparently as a side effect of editing nearby lines rather than a deliberate change.

    Re-commented as before. This time also added a static regression test, since four occurrences of the identical defect is not a coincidence to keep fixing by hand: it reads `R/runOccJSDM.R` and asserts the live `data = occ_data_effort` line is not present, catching the defect before any fit is attempted rather than after. Verified: tier 1 passes clean, 144 of 144.

44. ~~**Two of the "assorted smaller items" in group C.**~~ **FIXED** (Alex; `R/runOccJSDM.R`, `R/output.R`). Closed 2 August 2026, Alex marked both in `f4a59e4`.

    (a) `createDataIdx()`'s `maxP` for `model = "occupancy"`. Verified moot rather than fixed as originally described: `maxP` is used nowhere outside `if (model == "two_stage")` blocks, so the described crash path is not reachable as the code now stands. Confirmed live via the tier-1 occupancy smoke test, which passes.

    (b) `computeSpeciesDetected()`'s roxygen no longer documents any signature at all -- title, description and `@noRd` only -- so the stale Beta-approximation `@param` block is simply gone.

45. ~~**`listPriors$b_betatheta_slope_var`, a new prior hook.**~~ **APPROVED 2 August 2026** (Alex: "Happy with the new fix"). Closes group A item 1.

    Exposes `B_betatheta`'s previously hard-coded slope variance, default 2 unchanged. Built as a diagnostic for the `beta_theta` slope item, which is still open; approving the hook is not a decision on the value, which remains the `b_betatheta` variance item in group B. Alex's `522b89e` also removed the explanatory comment at the call site that had marked it as provisional, consistent with treating it as a permanent feature now.

46. ~~**`B0`'s bias doubled between the pre- and post-fix runs.**~~ **CLOSED BY DECISION 2 August 2026** (Alex: "we can close the point on B0 and assume that the bias is only due to the confounding with beta_theta"). Not closed by elimination -- the quantitative link was never established, only that `B0` is unbiased in both arms tested with `beta_theta` absent entirely (`binary` +0.0122 at 1.0 SE, `continuous` +0.0066 at 1.2 SE, `PLAN.md` 16.4), against -0.0633 at 3.3 SE with only the intercept present. Alex accepted that as sufficient. Full investigation, including the two hypotheses tested and ruled out along the way, in `AGENTS.md`.

    The coverage sub-finding this item also carried (`continuous` at 0.879 against nominal, `PLAN.md` 16.5) was not addressed by this decision and is carried forward as its own item in group B.

47. ~~**Three assorted smaller items (C2b, C2c, C2d), all verified fixed (Alex).**~~

    **(b) `d > NULL` errors on single-species input.** `get_param(listParams, "n_factors")` returns 0 by default; the old cap `if (d > ncol(OTU))` errored when `OTU` was a single-species vector because `ncol()` returns `NULL`. Fixed by branching first on `ncol(OTU) > 1`; the `else` branch sets `d <- 0` with a message. `R/runOccJSDM.R:759-770`.

    **(c) `reparamFactorModel(A_output, C_output)` fails when `gt > ncov_psi`.** `C_output` is `[gt x ncov_psi]`; when `gt > ncov_psi`, `qr.Q()` returns a non-square `[gt x ncov_psi]` matrix and `Q %*% diag(diag(R), nrow = gt)` fails on conformability. Fixed by `gt_default <- floor(sqrt(min(S, ncov_psi)))`, which ensures the default `gt <= ncov_psi`. A caller who manually sets `listParams$n_lattrait > ncov_psi` could still trigger it, but that is a usage error. `R/runOccJSDM.R:773`.

    **(d) Spatial-covariate `is.numeric` check ran before the names-present check.** A mistyped `spatCovariates` name gave `undefined columns selected` rather than the intended "Covariate names provided not in data\$info". Fixed by reordering: the `%in% colnames(data$info)` check (`:570`) now precedes `sapply(data_info[,spatCovariates], is.numeric)` (`:574`). `R/runOccJSDM.R:570-576`.

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

    It found six defects -- three while the tests were being written, three from the run -- none of which the static audit had caught. Four of the six are now fixed (Fixed bugs 24-27 and group B); the rest are the `sample_ls()`, `reparamFactorModel()` and `beta_theta` slope items in group B.
