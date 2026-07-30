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

     Two related rules:
     - NO EM-DASHES. Use a plain double hyphen. Em-dashes were the single
       largest driver of reflow churn under the old setting, because
       pandoc writes them 2 characters wider and every later wrap point
       shifted.
     - Keep inline `code spans` short. Under the old 72-column wrap, a
       span landing on the boundary could be mangled: in 2242b34 that
       turned `* 2` into `- 2` in the logDetKuu note, silently asserting
       something false. Less dangerous now, but still good practice. -->
```

# **v0.1.0-beta** **Public release**

## A. Review Claude fixes (Alex)

1.  **Review the RNG seeding fix, and answer two questions it raises.** **Written by Claude** on 29 July 2026: `src/rng.h` (new), `src/functions.cpp`, `src/jsdm.cpp`, `R/runOccJSDM.R`. Flagging it rather than leaving it to surface in a diff, since it changes files you own and you have not seen the design. Revert or rework it freely. It is tested (see below) but it has had no human review beyond Doug agreeing it should be fixed before the simulation study is re-run.

    **What it fixes.** Pointing `randinvg()` at the `thread_local` engine correctly closed the OpenMP data race, but neither engine it left behind read R's RNG state (`get_rng()` was seeded from the literal `12345 + omp_get_thread_num()`, `mvrnormArmaQuick_TS()` from `std::random_device{}()`), so `set.seed()` had no effect on the sampler. Two fits of the same data under one seed differed by 5.09 on `B0_output`. Not a bias -- both are valid streams -- but users could not reproduce a fit, and the literal seed was per *process*, so simulation-study replicates at the same index in different workers shared random numbers. Full detail, including three non-obvious traps for anyone reworking it, is in **Fixed bugs 28**; the third (`std::normal_distribution` caching the second Box-Muller deviate across `rng.seed()`) is the one most likely to bite.

    **Question 1: is thread-count invariance worth having?** Reproducibility currently holds *for a given thread count*. Threads derive separate streams, so with OpenMP genuinely active, changing the thread count changes which stream produces which element. Making the fit invariant needs a different scheme -- a counter-based RNG keyed on the element index rather than the thread id. Doable, but a bigger change than the bug warranted, so it was left as your call.

    **Question 2: what does `R CMD config SHLIB_OPENMP_CXXFLAGS` return on your machine?** On Doug's machine it is *empty*, so every `#pragma omp` compiles to a no-op and the built `.so` contains no OpenMP runtime calls at all -- the package runs single-threaded there. If yours differs, your experience of both the thread-safety bug and of performance differs from his, and the answer decides whether the group D parallelism items are worth anything. See group

    ALEX TO REVIEW THE FIX

2.  **Review the `plotCollectionRates()` fix, and decide how far to take it.** **Written by Claude** (AI agent, in a session with Doug) on 29 July 2026, not by Doug. Touches `R/output.R` only. Reverted or reworked freely; it has had no human review beyond Doug reporting the error.

    **What was broken.** The function errored on *every* input with `object 'Min' not found`. `plotSpeciesRates()` had been extracted as a shared helper and never wired up to its only caller: it read columns `Min`/`Max` while the caller passed `quantile()`'s `2.5%`/`97.5%`, filtered on a `Species` column the caller never created, and referenced `speciesNames` as a free variable that exists in neither its arguments nor the namespace. Full detail in **Fixed bugs 31**.

    **The part that matters most: the vignette calls this.** `vignettes/occJSDM.Rmd:515` calls `plotCollectionRates(fitmodel, idx_species = 1:10)` in a plain evaluated chunk -- no `eval = FALSE`. So the vignette cannot have been rebuilt since this broke, and `R CMD check` would have failed on it. That makes it a CRAN blocker, not a cosmetic plotting bug. **Worth rebuilding the vignette to confirm it now passes**, and worth asking what else has gone stale in it for the same reason.

    **Three decisions.**

    (a) **Is the helper's new contract the one you want?** `plotSpeciesRates(data_plot, idx_species, speciesNames)` now takes one row per species with `Min`/`Max` columns plus the *unsubset* names, and does the subsetting itself, so ordering is computed on exactly the rows plotted. If the helper was meant to serve the other three rate plots, this is the contract they would adopt -- so it should be a shape you are happy with.

    (b) **Do you want the other three fixed the same way?** `plotOccupancyRates()`, `plotFPTPStage2Rates()` and `plotStage1FPRates()` still order on the full species set while indexing a filtered one, which silently mismatches labels to bars (group B item 1). The test added for `plotCollectionRates()` asserts label-to-value pairing against an independent recomputation and is the template.

    (c) **Delete the commented-out block?** The original inline implementation is still sitting in `plotCollectionRates()`, below the return. It predates the extraction and is now doubly dead.

        ALEX TO REVIEW THE FIX

3.  **Investigate why `beta_theta`'s posterior gets overconfident as M rises, and review the new `listPriors` hook added to test it.** Two experiments this week, both by Claude, both negative -- the cause is not found, and it needs someone who knows the sampler.

    **The finding.** `beta_theta` coverage falls *monotonically* with M (0.747 at M2 down to 0.579 at M20) while bias stays small and flat throughout -- a shrinking interval around a bias that is not itself shrinking. That is the signature of overconfidence, not of insufficient data: more information makes it worse, not better. Full numbers in group B item 4 and `PLAN.md` 13.7.

    **Two candidate causes, both ruled out.** (a) `B_betatheta`'s slope prior variance, hard-coded to `diag(2)` -- tightened it 4x (SD 1.41 to 0.71) and coverage moved by 0.001-0.003, noise against the 3.1% SE. (b) Pseudo-replication in `X_theta` -- checked whether the M samples at a site share a covariate value, which would make added M look like independent information without being any. It does not: `X_theta` is drawn independently per sample, not per site (`R/simulateData.R:126`). Detail in `PLAN.md` 13.8-13.9.

    **So the cause is in the likelihood or the sampler's variance computation for `beta_theta`.** Worth a look: how the latent `w`/`z` state is aggregated across a site's M samples in `sample_beta_cpp_TS()` / `sample_betatheta_cpp_parallel()` (`src/functions.cpp`), or a numerical issue in the Polya-Gamma update they call. This is a sampler-code question -- the two obvious prior and data-generation explanations are both closed off.

    **What to review alongside it.** Testing (a) needed a real `listPriors` hook, since `B_betatheta`'s variance had none, unlike `p`/`q`/`theta0`. Added `listPriors$b_betatheta_slope_var` to `R/runOccJSDM.R`, default unchanged at 2, so nothing changes unless a caller sets it. Tested in `test-regression-bugs.R` and it reaches the sampler correctly, but has had no review from you.

    **Qualified 30 July (`PLAN.md` 14.7).** Candidate (a) was ruled out on evidence from M10/M20 only, which are the arms where data dominates the prior. At M = 2 the prior is *not* inert: tightening it moves coverage 0.747 to 0.653. It remains ruled out as a **fix**, since that is the wrong direction, but the same knob turns out to control `B0`'s bias, which is now a decision in its own right. See item 5 below before changing anything here.

    ALEX TO INVESTIGATE THE SAMPLER

4.  **Review the `predictNewSites()` fix, and confirm two API decisions in it.** **Written by Claude** on 30 July 2026, not by Doug. Touches `R/output.R` and **`src/jsdm.cpp`**. Revert or rework freely; it has had no human review beyond Doug asking for the fix. Full account in Fixed bugs 34.

    **What was broken.** `X_psi` and `X_s` had no defaults, so the `is.null()` guards you had written could never fire. Underneath that, the guards used `&` rather than `&&`, which evaluates both sides, so the missing promises were forced even when the caller had asked for neither term. `predictNewSites(fit, useEnvCov = FALSE, useSpatial = FALSE)` failed too. There was no way to call the function without supplying both matrices.

    Fixing the guards made three code paths reachable for the first time, and all three were broken. Two of those fixes are in your C++, which is the main reason this needs your eyes:

    (a) `computeNewOutputs()` sliced `Ks_all` and `Bs_output` before the `useSpatial` check, so `useSpatial = FALSE` aborted with `Cube::slice(): index out of bounds`. Moved inside the guard.
    (b) It read the new-site count as `X.n_rows` unconditionally, so `useEnvCov = FALSE` returned a silently empty result. Now taken from whichever term is active.
    (c) R-side: a fit with no spatial field returns `Bs_output` with a zero-length first dimension, which collapses under `apply()` and broke the `aperm()`. Only reshaped when the spatial term is in play.

    **Decision 1: is the tri-state the API you want?** `useEnvCov` and `useSpatial` now default to `NULL` and follow the pattern `useBiotic` already used: `NULL` uses the term if the fit estimated it, `TRUE` uses it and errors if it did not, `FALSE` skips it. This matches what their roxygen always promised, which was unreachable because `useSpatial` hard-stopped instead of ignoring. It is strictly more permissive, so no call that worked before changes behaviour. But it does change two documented defaults from `TRUE` to `NULL`.

    **Decision 2: should "no covariates and no spatial" be an error, or should the function take `n_new`?** With both terms off, nothing in the arguments says how many new sites to predict for; the answer would be the intercept plus the biotic term, identical at every site. It currently errors. The alternative is an explicit `n_new` argument, which would make that a legitimate call. I chose the error because it is the conservative reading, but it is a design question, not a bug.

    ALEX TO REVIEW THE FIX

5.  **Decide `b_betatheta`'s slope prior variance. It trades `B0`'s bias against `beta_theta`'s coverage, and both are open findings.** Measured 30 July 2026 (`PLAN.md` 14.7). This is a modelling decision, not a bug: nothing here is broken in a way that has a right answer without knowing what the prior is meant to encode. Nothing has been changed; the default is still `diag(2)` as you set it.

    **What the knob does.** All three arms below are M = 2, R = 50, paired on identical truths, varying only `b_betatheta_slope_var`:

    | variance                 | `B0` bias  | `beta_theta` coverage | `theta0` coverage |
    |--------------------|-----------------|------------------|-----------------|
    | 2 (your current default) | -0.160     | 0.747                 | 0.986             |
    | 0.5                      | -0.106     | 0.707                 | 0.982             |
    | 0.1                      | **-0.044** | **0.653**             | 0.980             |

    Paired, `B0`'s bias improvement is 2.1 SE at 0.5 and 4.0 SE at 0.1. `B0` coverage is 0.946-0.950 throughout, unchanged.

    **Why this matters for `B0` (group B item 6).** That item was filed as a possible regression with no known cause. It now has one: `42198d9` widened `B_betatheta` from `diag(1)` to `diag(2)` while correcting the mean from 1 to 0, and that is exactly when `B0`'s bias doubled, from -0.135 to -0.228. Turning the variance back down moves the bias back, monotonically. The `jsdmfun.R` rewrite in the same pull is no longer the leading suspect.

    **Why it is not simply "turn it back down".** The same change degrades `beta_theta`'s coverage, from 0.747 to 0.653 at variance 0.1. Group B items 4 and 6 pull in opposite directions on one parameter. Tightening buys a better point estimate for species intercepts at the cost of more overconfident collection-covariate intervals.

    **What is not known.** Whether some intermediate value is better than both endpoints, since only three points were measured; whether the trade looks the same at M \> 2, since this was run only at the production setting; and whether the correct move is to keep `diag(2)` and fix `beta_theta`'s overconfidence at its source instead, which is item 3 above. If the sampler-level cause of item 3 turns out to be the real problem, this trade may dissolve.

    **The hook exists**: `listPriors$b_betatheta_slope_var`, added for this test, defaults to 2 so current behaviour is unchanged. Reverting it is a one-line change if you would rather it stayed hard-coded.

    ALEX TO DECIDE THE VALUE (or to decide that item 3 supersedes this)

6.  **Review two C++ changes made to clear the install WARNING, one of which was wider than asked.** **Written by Claude** on 30 July 2026. Touches `src/functions.cpp` and `src/jsdm.cpp`. Revert freely.

    **Change 1, requested: `&` to `&&` on boolean operands, 8 sites.** `src/functions.cpp:670,672,674,676` (in `sample_pq_cpp`'s counting loop) and `src/jsdm.cpp:203,219,235,252` (the four `isPointInBand*` helpers). `clang` reported these as `-Wbitwise-instead-of-logical`.

    **This was safe, and here is the argument.** In C++ `==` and the relational operators bind tighter than `&`, so `a == b & c == d` already parsed as `(a == b) & (c == d)`, which is what was intended: the results were never wrong. Every operand is a plain comparison or an Armadillo element read, with no side effects, so short-circuiting cannot skip anything that mattered. The change is therefore semantics-preserving and only stops evaluating operands whose answer is already determined. 167 tests pass unchanged.

    **Worth contrasting with the R-side version of the same habit** (Fixed bugs 34): there `&` *was* a real bug, because the unevaluated operand was a missing-argument promise, so forcing it turned an intended error message into `argument "X_psi" is missing`. Same idiom, different consequence, because R has lazy evaluation and C++ does not.

    **Change 2, NOT requested, please check you are happy with it.** Removed `TRUNC` and `TRUNC_RECIP` from the top of both files (4 declarations). They were declared and referenced nowhere, and drew `-Wunused-const-variable`, which was the *only* remaining package-attributable cause of the install WARNING once change 1 landed. I took them out rather than leave the CRAN blocker half-fixed, and left a comment in each file recording the values (0.64 and 1/0.64, the truncation point from Windle 2013) and saying to restore them if the Polya-Gamma sampler is ever changed to use the truncated form. **If you would rather keep them, say so and I will restore them with a `// [[maybe_unused]]` or equivalent instead.**

    **What this achieved, measured:** the package now contributes **zero** compiler warnings. `R CMD check` still reports an install WARNING, but it is no longer ours: the sole remaining "significant warning" is from **R's own header**, `R_ext/Boolean.h:62`, where a `#pragma` names a warning group this clang version does not recognise. That is an R-build and toolchain artifact, is environment-specific, and is not fixable from this package. Do not chase it.

    ALEX TO REVIEW, ESPECIALLY CHANGE 2

## **B. Inference-affecting bugs (wrong numbers, silently) (Alex)**

1.  **`sample_ls()` evaluates the wrong density, so the GP length-scale is never recovered and rails at the top of its grid.** Found 27 July 2026 while writing the test suite.

    **This entry previously blamed `sigma_s` being hard-coded to 1 and proposed sampling it. That diagnosis was wrong and the proposed fix is disproven** -- see "What was ruled out" below. The real defect is deeper and is a modelling change, not a parameter addition.

    **The defect.** Under the SoR (subset-of-regressors) approximation the fitted spatial field is `SE = Ks(l_s) %*% Bs` (`computePsiCoef()`), i.e. a *deterministic function of `Bs` and `l_s`*. But `sample_ls()` (`R/jsdmfun.R:1054`) treats `SE` as a **GP draw** and scores its density under `N(0, sigma_s^2 K(l_s))` while holding `SE` fixed. That is not the conditional posterior of `l_s`, and it is self-defeating: `SE` was already smoothed at the current `l_s`, so it scores progressively better under ever-smoother covariances. The result is a likelihood that is **monotone increasing in `l_s`**, with nothing to penalise the upper bound.

    **Measured.** `idx_ls` sits at 10 -- the maximum of `l_s_grid = seq(0.01, 0.3, length.out = 10)` -- with range 10-10 across the whole post-burn-in run, for **every** true `l_s` tried (0.074, 0.171, 0.300) and with genuine spatial signal present (`ds = 2`). Profiling the grid confirms the cause: the log-likelihood rises monotonically from -376 at `l_s = 0.01` to -154 at `l_s = 0.30`, so the true value is never competitive.

    **What was ruled out** (recorded so it is not re-investigated):

    (a) *Amplitude is not missing from the model.* `sigma_bs` **is** sampled (`R/jsdmfun.R:1211`) and **is** the SoR amplitude, since `Bs` carries prior variance `sigma_bs^2` and the kernel is built at unit amplitude (`K2(..., 1, l_s)`). A separate `sigma_s` would be redundant and would create an identifiability problem.
    (b) *Passing the correct amplitude does not help.* Profiling the grid at the fitted field's actual amplitude (`sd(SE) = 0.85`) instead of 1 leaves the likelihood monotone -- it still rises to the grid maximum. So the one-line fix this entry used to propose would not work.
    (c) *`logDetKuu` is correct.* Note the factor is `* 2`, not `- 2`: `sum(log(rcppeigen_get_diag(K_xx))) * 2` matches the true log-determinant at all ten grid points, as it should for a Cholesky factor. `rcppeigen_get_diag()` returns the **Cholesky** diagonal despite the name.
    (d) *Not an artefact of weak data.* Identical behaviour with `ds = 2` (real spatial signal) as with `ds = 0`.

    **Impact.** Biases the spatial term of every fit using `useSpatField = TRUE`, and makes `l_s` uninterpretable. Note widening `l_s_grid` alone would simply move the rail.

    **Fix (needs Alex -- a derivation, not a code tweak).** Evaluate `l_s` where it actually acts, on the data: either recompute `eta` with `Ks(l_s*)` and score the observation likelihood under the proposal, or integrate out `Bs` and use the standard SoR marginal likelihood. The first is the smaller change but costs an `eta` recomputation per proposal; the second is cleaner and harder to derive. Roughly half a day if the intended conditional is already clear, longer if it needs settling. **Do not attempt blind** -- a subtly wrong conditional would be worse than the present state, which at least fails visibly.

    **Cheap interim alternative (\~1 h).** Since `sample_ls()` does not work, and enabling it also caused the non-spatial crash (Fixed bugs 22 / group B item 1), reverting to a *documented, user-settable, fixed* `l_s` is arguably more honest than the current state: it removes a moving part that is not moving correctly and exposes the assumption instead of burying it in a call-site literal.

    **Not visible before `b7b6aa2`.** While `sample_ls()` sat inside `if(F){...}` the length-scale was openly fixed, which was the original audit's complaint. Enabling it exposed that the update itself was never correct.

    **Test coverage:** `tests/testthat/test-regression-bugs.R` asserts that `sample_ls()` is *reached* rather than that `idx_ls` varies, precisely so the test does not encode this defect. Revisit that test once this is fixed.

    ALEX TO CHECK

2.  **`reparamFactorModel()` breaks the identity that residual covariance = `t(L) %*% L`, so reported species correlations are inflated.** Found 28 July 2026 by the coverage study; confirmed numerically. **Disputed -- see Alex's note at the end of this item.**

    `reparamFactorModel()` (`R/jsdmfun.R:48`) rewrites the stored `U` and `L` each iteration as `L_new = diag(1/diag(R)) %*% R` and `U_new = U %*% Q %*% diag(diag(R))`, where `L = QR`. That **preserves `U %*% L`** -- verified to 4e-16 -- so the linear predictor is untouched and the rotation is legitimate for identifiability.

    But `returnResidualCorrelationMatrix()` computes the correlation as `cov2cor(t(L) %*% L)` from the *stored, reparameterised* `L`. That identity holds only while `U` is standardised, and the reparameterisation moves scale out of `U` into `L`: measured `Var(U)` afterwards is `diag(0.23, 2.01)`, not the identity.

    **Measured**, `d = 2`, `S = 6`. These are `cov2cor()` output -- already correlations, not covariances:

    | pair | original | after reparam |
    |------|----------|---------------|
    | 1    | +0.713   | **+0.926**    |
    | 2    | -0.811   | **-0.958**    |
    | 3    | -0.168   | **-0.780**    |
    | 4    | +0.443   | **+0.875**    |
    | 5    | +0.585   | **+0.285**    |

    Maximum change **0.612**, consistently toward the extremes. Inflated point estimates with correctly-sized intervals is exactly what produces undercoverage, and across the full grid `resid_cor` covers at 0.74-0.77 in nine of ten scenarios (`PLAN.md` §12). The exception is diagnostic: `d_underfit` *over*covers at 0.980, so under-fitting the ordination widens the intervals enough to mask the bias -- a model with too few factors looks better calibrated than it is.

    **Impact.** `returnResidualCorrelationMatrix()` and `plotResidualCorrelationMatrix()` overstate species co-occurrence. This is the headline JSDM output and the basis of the ordination interpretation.

    **Fix, two options.** (a) Make the reparameterisation purely orthogonal -- rotate by `Q` alone, dropping the `diag(diag(R))` scaling. Then `t(L_new) %*% L_new = t(L) %*% L` and `Var(U_new) = I`, so both the identifiability constraint and the covariance identity hold. (b) Keep the scaling and compute the correlation as `t(L) %*% Var(U) %*% L`. (a) is simpler and preserves the existing output contract.

    ------------------------------------------------------------------------

    **Alex's note:** *"THIS IS NOT AN ISSUE FOR LOGISTIC MODELS (as the covariance matrix is not recoverable anyway, only the correlation one)."*

    **Response -- the premise is right but does not resolve it.** In a logistic latent-variable model the scale is indeed unidentified and only the correlation is meaningful. But the table above is already correlations, and they still move by 0.612. `cov2cor()` does not absorb the change because `diag(1/diag(R))` rescales **per factor** while `cov2cor()` normalises **per species** -- a per-factor rescaling is not a global scale change, so it survives the normalisation.

    **If it is nonetheless a non-issue, then the simulation study is measuring the wrong thing and that needs identifying**: it compares `cov2cor(crossprod(L_true))` from the simulator's own `L` against the same function of the stored `L`. Either that comparison is invalid, or the correlations really are inflated. Worth settling before the `resid_cor` row of `PLAN.md` §12 is quoted anywhere.

    **Confirmed by the 29 July re-run, in a paired design.** `draw_truth()` seeds on (scenario, replicate), so the simulated data and true values are *bit-identical* between the pre- and post-fix runs -- verified, `max|truth difference| = 0`. Across that pair only **104 of 49,978** `resid_cor` coverage decisions flipped, and coverage was unmoved at 0.777 (0.752-0.768 in nine of ten cells). Same data, same truth, same answer: the undercoverage is not sampling noise and does not depend on the RNG or on any of the four fixes. `d_underfit` still *over*covers at 0.980, which is diagnostic -- under-fitting the ordination widens the intervals enough to hide the bias, so the defect's visible severity depends on the fitted factor dimension. See `PLAN.md` 12.1.

    ALEX TO MAKE A DECISION

3.  **Choose the informative priors for `p`, `q` and `theta0`.** Alex's note in this file: *"NEED TO CHOOSE THE INFORMATIVE PRIORS"*. **A design task, not a bug** -- an earlier version of this entry called the prior a defect, which was wrong.

    **An informative prior is required here, not optional.** `p` and `q` enter the sampler perfectly symmetrically --

    ```         
    p(l,s) = rbeta(a_p + [w=1 & detected], b_p + [w=1 & not detected])
    q(l,s) = rbeta(a_q + [w=0 & detected], b_q + [w=0 & not detected])
    ```

    -- and **no `p > q` constraint is enforced anywhere in the code** (checked). The augmented likelihood is therefore invariant under swapping *(collected, `p`)* with *(not collected, `q`)*: the label-switching multimodality of false-positive occupancy models. The only thing selecting the correct mode is the prior, `Beta(5, 1)` pushing `p` high against `Beta(1, 20)` pushing `q` low; `theta0` gets the same treatment at Stage 1. Flatten them and the model is unidentified.

    **What the study measured, correctly interpreted.** In `low_information` (true `p` in 0.1-0.3 against a prior mean of 0.833) `p` covers at **0.103** with bias **+0.49**, against 0.90-0.92 elsewhere. That is not a defect firing -- it is **the cost of the identifiability constraint, quantified**, in the regime where constraint and data disagree most.

    **The decision.** How informative should the default be?

    - Too weak: the chain can flip to the mirror mode, which presents as catastrophic nonsense rather than a warning.
    - Too strong: low-detection studies -- the eDNA norm, and what `data/sampledata.rda` was deliberately regenerated to represent -- are biased upward by up to 0.5.

    Worth informing with the Ji et al. (2025) posteriors, which are real detection rates from the system the model was built for.

    **Already settled, no decision needed** (Alex, `42198d9`): all four hyperparameters now read from `listPriors`, so a user can override them, and the documented defaults match the code.

    **Whatever is chosen, document the tradeoff** under `@param listPriors`. Someone running a low-detection study needs to know to set `a_p`/`b_p`, and nothing currently tells them.

    ALEX RESPONSE: WE'LL ACCEPT THIS PRIOR AND THE BIAS

4.  **`beta_theta` still undercovers at 0.766 after the prior-mean fix, so a second cause remains.** Measured by the 29 July R = 100 re-run (`PLAN.md` 12.1).

    Alex's correction of the collection-covariate prior mean from 1 to 0 (Fixed bugs 25) was a real cause and did real work: every cell gained 0.03-0.05 coverage and two-thirds of the bias went (+0.112 -\> +0.038). But 0.766 against a nominal 0.95 is still far outside the 2.2-point SE, and the residue has the same signature as before -- **flat across model type, primer count, species count and factor misspecification** (0.709-0.771 in nine cells), which says structural rather than conditional.

    `low_information` is the one cell that now *over*covers, at 0.983, having been 0.860 pre-fix. A prior that no longer pulls, combined with the widened `diag(2)` variance, appears to overshoot when the data are thin -- see item 5.

    **This is the clearest open target in group A**: it is the largest remaining deviation that is definitely a defect (unlike `p`, which is the identifiability constraint working as intended).

    **M-ladder result, 29 July 2026 (`PLAN.md` 13.7): the opposite of the hypothesis, and this rules out under-identification.** Doug's directive was to re-run at M \> 10 and see if it fixes this. Coverage instead falls *monotonically* with M -- 0.747 at M2, 0.655 at M5, 0.603 at M10, **0.579 at M20** -- while bias stays small and flat throughout (+0.02 to +0.05). The matched control (`K30`: same row count as `M20`, spent on PCR replicates instead) covers at 0.706, beating every M arm above M2.

    Shrinking intervals around a bias that is not itself shrinking is the signature of a real defect being *exposed* by more information, not resolved by it. More field samples make this worse, so the next step is not more data but finding what makes the interval overconfident.

    **Tighter-prior result, 30 July 2026: the prior is not the cause, but the first test was run in the wrong place.** Tightening `B_betatheta`'s slope variance at M10/M20 moved coverage by 0.001-0.003 (`PLAN.md` 13.9). **At M = 2 it moves it a lot: 0.747 to 0.653** (`PLAN.md` 14.7). The conclusion survives, since tightening makes coverage *worse* and so is not a fix, but the prior is not inert as 13.9 implied: it was inert only in the arms tested, which were the ones where data dominates it. Bias falls while coverage falls faster, so the intervals shrink more than the bias does. A second candidate, pseudo-replication in `X_theta`, is also ruled out. **Full detail and the review item are above, under "Review Claude fixes" item 3**: this needs sampler-level investigation.

    ALEX TO REVIEW A3 ABOVE (THIS IS REDUNDANT WITH A3).

5.  **`theta0` now overcovers at 0.978-0.985, having been near nominal.** Measured by the same re-run (`PLAN.md` 12.3).

    Pre-fix it sat at 0.938-0.959, i.e. fine. Post-fix it is 0.978-0.985 in nine cells. The all-cell average of 0.944 hides this, because `low_information` pulls it down (0.602, up from 0.477 but still the worst cell in the table).

    Overcoverage is the safe direction and this is not urgent, but the pattern points at the same edit as item 4 -- the widened `diag(2)` prior variance on the collection coefficients looks to have overshot. Worth examining the two together, since one change plausibly produced both.

    **M-ladder result, 29 July 2026 (`PLAN.md` 13.7).** Overcoverage falls toward nominal as M rises (0.986 at M2, 0.944 at M10), while the matched control (`K30`: same row count, spent on PCR replicates instead) makes it worse at 0.996.

    **That was read as Stage 1 under-identification. The reading has a hole in it, found 30 July.** Pre-fix, `theta0` sat at 0.938-0.959 at the *same* M = 2 where it now sits at 0.978-0.985. If M = 2 were simply too little information for `theta0`, it would have overcovered before Alex's fixes too. It did not. So the question is not whether `theta0` is fine at high M, which four arms already say it is, but **what changed at M = 2** -- the setting users actually run.

    **Likely route: coupling, not `theta0`'s own prior.** `a_theta0`/`b_theta0` are untouched, still `Beta(1, 20)`, last changed 23 July before the pre-fix run. But `b_betatheta`'s variance was widened to `diag(2)` in Fixed bugs 25, which changes `beta_theta`, which drives the collection probability, which drives the latent `w`; and `sample_theta0(z, w, ...)` conditions on `w`. The tighter-prior run supports this weakly: `theta0` moved toward nominal in both arms tested (M10 0.944 to 0.940, M20 0.952 to 0.942), small because those are the arms where data dominates the prior. M2, where the prior should dominate, is the one arm not yet run at the tighter setting.

    **Plan: `PLAN.md` 14.** Two arms at M = 2 varying only `b_betatheta_slope_var` (0.5 and 0.1 against the default 2), R = 50, about 16 minutes. If `theta0` walks toward 0.95 as the variance tightens, this closes as a downstream symptom of item 4 rather than a separate defect, and the two should be fixed together.

    **The previous plan, "re-run the M10 or M20 arm at R = 200", is dropped.** It would confirm that `theta0` reaches nominal at high M, but that is not a claim the paper needs or that users can act on at M = 2 to 3. See `PLAN.md` 14.6.

    **Priority: lowest of the open findings, and this should not grow.** Overcoverage is the safe direction; it costs power, not correctness. The case for the 16 minutes is that it likely resolves this as a side effect of diagnosing item 4, not that it matters on its own.

    **Result, 30 July 2026 (`PLAN.md` 14.7): the coupling hypothesis is disproved.** A 20-fold reduction in `b_betatheta`'s slope variance moved `theta0` coverage by 0.006 (0.986, 0.982, 0.980 at variance 2, 0.5, 0.1). Paired, the bias change is detectable at 2.7 SE but far too small to matter. `b_betatheta` is not the route.

    **Next test, if it is worth doing at all:** `theta0`'s own `Beta(1, 20)` prior, which already has `listPriors$a_theta0`/`b_theta0` hooks and needs no code change. Given this is the least urgent open finding and overcoverage is the safe direction, it is worth running only alongside another study, not on its own.

    CLAUDE TO PIGGYBACK A theta0-PRIOR ARM ON THE NEXT RUN

6.  **`B0` bias roughly doubled, and coverage does not show it.** Measured by the same re-run (`PLAN.md` 12.2). **Possible regression, cause not yet identified.**

    Nine of ten cells moved more negative between the pre- and post-fix runs on identical data: base -0.113 -\> -0.208, `occupancy` -0.024 -\> -0.151, `primers_3` -0.031 -\> -0.151, `low_information` -0.931 -\> -1.056. Only `binary` moved the other way. Overall -0.135 -\> -0.228.

    **Coverage stays at 0.943**, because the intervals are wide enough to absorb the shift -- so this is invisible in the headline table and shows only in the bias column. That is the argument for tracking both.

    Nothing in the four fixes should move species intercepts this way. The likeliest candidate is the 421-line `jsdmfun.R` rewrite that shipped in the same pull (`42198d9`). An R = 8 probe on 29 July hinted at it and R = 100 confirmed it, so it is not noise. **`B0` is a headline quantity for a JSDM -- this needs a cause before the MEE paper reports species intercepts.**

    **M-ladder result, 29 July 2026 (`PLAN.md` 13.7): partial confirmation, mechanism ambiguous.** Bias collapses from -0.160 at M2 to near zero at M10 (+0.002), which looks like item 5's clean pattern. But the matched control (`K30`) also improves it (-0.091), just less than M10/M20 do -- so the recovery is not cleanly Stage-1-specific; more data of either kind helps somewhat. Coverage stays 0.94-0.96 in every arm, consistent with the original finding that coverage does not reveal this bias.

    Still needs the `jsdmfun.R` rewrite investigated as a candidate cause (see below), since the ladder does not rule it out -- it only shows that *some* of the effect is an M/data-volume story.

    **Cause identified, 30 July 2026 (`PLAN.md` 14.7), found while testing something else.** `B0`'s bias responds strongly and monotonically to `b_betatheta`'s slope variance: -0.160 at the current default of 2, -0.106 at 0.5, -0.044 at 0.1. Paired against the same truths, that is +0.054 (2.1 SE) and +0.116 (4.0 SE).

    **The history matches exactly.** `42198d9` widened `B_betatheta` from `diag(1)` to `diag(2)` while correcting the mean from 1 to 0, and that is precisely when `B0`'s bias doubled from -0.135 to -0.228. Turning the variance back down moves the bias back. So the widening, not the `jsdmfun.R` rewrite, is the leading candidate.

    `B0` coverage is 0.946-0.950 across all three variances, unchanged, which is why the original grid never flagged this.

    **But it is a trade-off, not a fix.** Tightening that variance helps `B0`'s bias and *hurts* `beta_theta`'s coverage (0.747 to 0.653 at variance 0.1, see item 4). Two open items pull in opposite directions on one knob. Setting it needs someone who knows what the prior is meant to encode.

    ALEX TO DECIDE: filed as "Review Claude fixes" item 5, with the dose-response table and the trade-off against item 4.

7.  **`q` (Stage 2 false positives) degrades hard as `K` rises.** Found 29 July 2026, as a side effect of the M-ladder run (`PLAN.md` 13.7) -- not something that run was built to look for.

    Coverage: 0.945 at `M2` (K = 3) down to **0.614 at `K30`** (K = 30, same total rows as `M20`). `M20` itself, which keeps K = 3 and raises M instead, sits at 0.742 -- worse than `M2` but far better than `K30`. So more PCR replicates make `q` less well calibrated, not more, and the effect is bigger than the M-driven change in any of items 4-6.

    Not investigated beyond this measurement. Worth checking against the same label-switching mechanism noted in item 3 above -- more PCR replicates sharpen the posterior, and if the informative prior is pulling `p`/`q` away from the true values by a fixed amount, sharper intervals would show it as *worse* coverage, exactly as seen here for `beta_theta` in item 4. If so this is not a new bug but the same cost-of-identifiability story, extended to K.

    ALEX TO DIAGNOSE THE CAUSE AND DECIDE WHAT TO DO ABOUT THIS

## **C. Crashes, unreachable code paths, and API bugs (Alex)**

1.  **`thinOutput()`: two of the three original defects remain** (residual of the original audit's B.2; see Fixed bugs 16). The crash is gone and `thin` is now honoured, but

    (a) the 2-D branch still thins **by row** (`x[idx_thinned,,drop=F]`, `R/output.R:57`), which silently drops *sites* from the `psi_output` / `w_output` / `theta_output` posterior-mean matrices under the default `summarisedLatentPresences = TRUE` -- those matrices are `n x S`, so their rows are sites, not iterations; and (b) the scalar `WAIC` still falls through to `print("Dimension not recognised")` (`:62`) and becomes `NULL`. Also worth settling: `thinOutput()` is no longer in `NAMESPACE` (de-exported), but `man/thinOutput.Rd` still ships a `\usage` block for it, so the manual documents a function users cannot call. Either re-export it or drop the `.Rd`.

        ALEX TO DISABLE IT FOR NOW, BUT LET'S KEEP IT THERE

2.  **Assorted smaller items.**

    (a) `createDataIdx()` is called with `maxP` (`R/runOccJSDM.R:638`) for `model = "occupancy"` too, where `maxP` was never assigned -- it survives only because R's lazy evaluation never forces the promise; pass `NULL` explicitly.

    (b) `d <- get_param(listParams, "n_factors")` defaults to 0, and the cap at `:716` uses `ncol(OTU)`, which is `NULL` for a single-species vector `OTU` -- `if (d > NULL)` errors.

    (c) `reparamFactorModel(A_output, C_output)` (`:1273`) assumes `qr.Q()` is square; when `ncov_psi < gt` it is not, and `diag(diag(R_current), nrow = d)` recycles.

    (d) The spatial-covariate numeric check at `:540` runs *before* the "names present in `data$info`" check at `:544`, so a mistyped spatial covariate name gives `undefined columns selected` rather than the intended message.

    (e) `simulateOccJSDMData()`'s `@details` still says "For two-stage eDNA data specifically, use `simulateOccJSDMData()`" -- a self-reference left over from the merge of `simulateOccJSDMDataGeneral()`.

    (f) `computeSpeciesDetected()`'s roxygen still documents the removed Beta-approximation signature (`ab_p`, `K`, `primer`, `alpha`) instead of its actual arguments.

        ALEX WILL REVIEW THE ABOVE

3.  **`plotSpeciesResponseCurve()` is exported but has an internal helper's signature.** It takes `species_name, target_cov, beta_mcmc_j, list_matrix, raw_df, n_points` -- raw MCMC pieces and a raw data frame, with no `fitModel` argument -- while every other plotting function in `R/output.R` takes the fitted object. A user holding a `runOccJSDM()` result has no documented way to assemble these arguments. Either give it a `fitModel` signature like its neighbours, or drop it from `NAMESPACE` and let `plotCovariateEffect()` call it internally. Until this is settled it is untestable, which is why the new smoke tests cover the other two GAM exports and not this one.

    WAIT FOR ALEX

4.  **`computeNewOutputs()` prints to stdout on every call, and it cannot be silenced.** Found 30 July 2026 while testing the `predictNewSites()` fix (Fixed bugs 34).

    `src/jsdm.cpp:476` runs `Rcout << "Computing species " << j + 1 << " out of " << S << std::endl;` inside the species loop, unconditionally. Every `predictNewSites()` call therefore prints one line per species, with no way to turn it off: `suppressMessages()` does not catch it, because `Rcout` is stdout rather than R's condition system, and `suppressWarnings()` does not either. `capture.output()` works but forces the caller to discard everything.

    **Why it matters beyond tidiness.** It pollutes test output, which is how it was noticed: the new `predictNewSites()` tests interleave `Computing species 1 out of 4` with testthat's progress dots. It would do the same inside any user script, vignette chunk, or Shiny app that predicts in a loop. Unconditional console output from a compute function is also the kind of thing CRAN reviewers pick up.

    **Fix.** Add a `verbose` argument, defaulting to `FALSE`, threaded from `predictNewSites()` through to `computeNewOutputs()`, and guard the `Rcout`. Progress reporting is genuinely useful for a slow prediction over many species, so removing it outright would be a loss; making it opt-in keeps it.

    **Related, same anti-pattern, already tracked:** `thinOutput()` uses `print("Dimension not recognised")` at `R/output.R:43` and `:63` where it should warn or error. That is part of item 1 above, not a separate fix, but worth doing in the same pass since it is the same class of problem.

    ALEX TO DECIDE AND FIX (touches `src/jsdm.cpp`)

## **D. Dead and broken internal code (Alex)**

None of this is reachable from an exported function, but it will draw `R CMD check` "no visible binding" notes and is a trap for anyone reading the source. Suggest moving the genuinely dead functions to `deprecated/` and wiring up the two missing imports.

### Plan for items 1-3, written 30 July 2026

**The catalogue in items 1-3 is incomplete.** It names 14 functions. A systematic scan (every function defined in `R/`, counting call sites in `R/`, `tests/` and `vignettes/`, excluding roxygen comment lines, and checking `NAMESPACE`) finds **38 dead R functions plus 8 unused `RcppExports.R` wrappers**. Deprecating only the 14 named would leave two thirds of the problem and most of the `R CMD check` NOTE behind.

Counts by file: `mcmcfun.R` 17, `jsdmfun.R` 17, `output.R` 3, `diagnostics.R` 1, `RcppExports.R` 8, `zzz.R` 1.

**Three different treatments, and conflating them is the trap:**

**(a) Move to `deprecated/`.** The genuinely dead R functions. `deprecated/` is already tracked, already excluded by `.Rbuildignore`, and already holds `functions_old.cpp` and the old vignette, so the precedent and the mechanism both exist. Suggested layout: `deprecated/R/mcmcfun-dead.R`, `deprecated/R/jsdmfun-dead.R`, `deprecated/R/output-dead.R`, each with a header saying what it was and when it was moved.

**(b) Remove the `// [[Rcpp::export]]` tag in C++, then re-run `Rcpp::compileAttributes()`.** The 8 `RcppExports.R` wrappers (`XsBs`, `XtOmegaX_SoR`, `convert_to_correlation`, `dist_matrix`, `findClosestPoint`, `gpCovMatrix`, `sample_betatheta_cpp`, `sample_w_cpp`). **Do not edit `RcppExports.R`**: it is generated, and hand edits are silently overwritten by the next `document()`. If the C++ function itself is also unused, it can go to `deprecated/functions_old.cpp` alongside the existing dead C++.

**(c) Leave alone. Two false positives the scan flags that must not be deleted:**

- **`.onLoad()`** (`R/zzz.R`). Zero callers because **R itself calls it** on namespace load. Deleting it would silently drop the `mc.cores` cap set for CRAN compliance.
- **`thinOutput()`** (`R/output.R`). Genuinely uncallable today (not exported, no callers), but the CRAN plan's option (a) for shrinking `sampleresults.rda` depends on it, and it is group C item 1. **Fix it or decide against it before deprecating it**; do not sweep it up here.

**Doug's directives, and one ambiguity to settle.** Item 1 says "deprecate, but keep for reference", which is exactly what `deprecated/` is for. Item 3 says "deprecate both". **Item 2 says "deprecate the sample functions not used in the MCMC", which is narrower than the item's own text**: item 2 also lists `computePredictiveProbs()`, `partition_r2()`, `returnSpatialEffectMean()` and `plotSpatialEffect()`, which are not `sample_*` functions. All four are dead by the same test. Confirm whether they go too, or stay pending a decision.

**Verified before proposing any of this:** no dead name is referenced as a string anywhere in `R/`, `src/` or `tests/`, so nothing is reachable by `do.call()`, `get()` or `match.fun()`. That is the failure mode a caller-count scan would otherwise miss, and it is why the scan alone was not treated as sufficient.

**Execution order:**

1.  Move group (a) file by file, running `devtools::test()` after each file rather than at the end, so a break is attributable to one move.
2.  Handle group (b), then `Rcpp::compileAttributes()` and `devtools::document()`.
3.  Re-run `devtools::check()` and record how far the undefined-globals NOTE falls. **Measured so far: 84 symbols down to 65**, a 23% reduction, from items 1, 3 and the four `sample_*` of item 2 alone. Removed: `Freq`, `Frequency`, `M`, `N3`, `OTU`, `Occupancy`, `Type`, `Var1`, `Var2`, `data_info`, `p`, `primerIdx`, `sumM`, `w`, `z`, plus `U`, `gt`, `gts`, `computeBscoef` earlier. The remaining `jsdmfun.R` functions and the 8 `RcppExports` wrappers are still to go.
4.  Only then start item 4 phase 2 (`globalVariables()`), which is blocked on this precisely so the surviving list is enumerated once.

**What "done" looks like:** 167 tests still passing, the NOTE reduced to data-masked column names only, and no entry in it naming a function that no longer exists.

**Risk, stated plainly:** this deletes a lot of Alex's code from the package build. It is recoverable from `deprecated/` and from git, and none of it is reachable, but it is his call whether "not currently reachable" means "not wanted". Worth his sign-off before step 1 rather than after step 4.

1.  ~~`R/mcmcfun.R` dead samplers~~ **DONE 30 July 2026** (Claude). `sample_z()`, `sample_w()`, `sample_cimk()` and `loglik_sigma1()` moved to `deprecated/R/mcmcfun-dead-samplers.R`, per "DEPRECATE, BUT KEEP FOR REFERENCE". The three samplers referenced undeclared `M`, `sumM`, `n`, `z`, `w`, `y` and are superseded by the `_cpp` versions; `loglik_sigma1()` was an unfinished stub whose entire body was `p[primerIdx[]]`.

2.  `R/jsdmfun.R`. **The four `sample_*` functions are DONE, 30 July 2026** (Claude): `sample_BCsL()`, `sample_U()`, `sample_Br()` and `sample_BL_fixed()` moved to `deprecated/R/jsdmfun-dead-samplers.R`, per Doug's scoping to sampler code only. All four were unreachable (no callers in `R/`, `tests/` or `vignettes/`, not in `NAMESPACE`, no string references so nothing could reach them via `do.call()`/`get()`/`match.fun()`) and all four referenced undefined globals, so none could have run as written. Superseded by the `_cpp` implementations. Measured effect: `R CMD check`'s undefined-globals NOTE fell from 84 symbols to 80, losing exactly `U`, `gt`, `gts` and `computeBscoef`. 167 tests unchanged.

    **Still open in this item, deliberately.** **Four more functions** listed here are dead by the same test but are **not** sampler code, so Doug's directive did not cover them: `computePredictiveProbs()` (references `fitModel`, `psi_output`, `X_ord`, `beta_ord_output`, none of which exist), `partition_r2()` (calls `pseudo_R2()` with one argument; it takes two), and the pair `returnSpatialEffectMean()` / `plotSpatialEffect()` (reference a global `Xs_centers` and a nonexistent `returnSpatialEffect()`).

    Two of them look like unfinished intent rather than abandoned code, which is why they were not swept up: `partition_r2()` relates to the live *MEE paper* item on site variance partitioning, and the spatial pair is the only spatial-field plotting anywhere in the package. `computePredictiveProbs()` looks straightforwardly superseded by `predictNewSites()`.

    ALEX/DOUG TO DECIDE ON THE REMAINING FOUR

3.  ~~`R/output.R` dead plots~~ **DONE 30 July 2026** (Claude). `plotReadIntensity()` and `plotOccupancyStates()` moved to `deprecated/R/output-dead-plots.R`, per "DEPRECATE BOTH". The first read `results_output$mu1_output` / `mu0_output` / `sigma1_output` / `sigma0_output` and `infos$maxexplogy1`, none of which `runOccJSDM()` produces any more; the second referenced undefined `data_info` / `OTU`. Their roxygen blocks moved with them, so nothing re-attached to the following function, and `man/plotReadIntensity.Rd` was removed by `document()` as a consequence. `NAMESPACE` is unchanged: neither was exported.

4.  **Missing imports, and the wider `R CMD check` undefined-globals NOTE.** Plan written 30 July 2026 from a full `R CMD check` run, which corrected this entry's premise.

    **This entry named the wrong functions.** It cited `plotCovariateTrend()` as the live `tidyr` user and `createSplinesObjects()` as the live `splines` user. Both are **dead**: zero call sites, not exported. Verified by mapping every symbol in the check NOTE to its enclosing function and counting callers.

    **What is actually live and missing an import:**

    | symbol | needs | live caller |
    |------------------------|------------------------|------------------------|
    | `pivot_longer` | `tidyr` | `returnCovariateEffect_base`, `plotCovariateEffect_base` (Alex's GAM functions, both exported, both used in the vignette) |
    | `setNames` | `stats` | `create_covariates_matrix`, `transform_new_covariates` (the latter is on `predictNewSites()`'s path) |
    | `rnbinom` | `stats` | `simulateData`, reached from `simulateOccJSDMData()` |

    **Dead-code-only, so do not import these; delete the code instead** (items 1-3 and 5 above): `bs` from `splines` in `createSplinesObjects()`, `dnbinom` from `stats` in `sample_rnb()`, and the third `pivot_longer` in `plotCovariateTrend()`. `splines` is not in `DESCRIPTION` at all, so importing it would add a dependency purely to support code that should not ship.

    **The NOTE is three different problems wearing one hat.** Roughly 90 symbols are listed, and treating them uniformly is the trap:

    (a) **Genuine missing imports.** The three rows above. `R CMD check` names `dnbinom`, `rnbinom` and `setNames` itself in its "Consider adding" line; it cannot name `pivot_longer` because it does not know which package that is meant to come from.
    (b) **Data-masked column names, which are false positives.** `x`, `y`, `species`, `Min`, `Max`, `2.5%`, `97.5%`, `Output`, `value`, `upper`, `lower`, `iter`, `chain` and so on are `dplyr`/`ggplot2` NSE references, not undefined objects. These are silenced with `utils::globalVariables()`, or avoided with the `.data$` pronoun. Silencing is cosmetic; the NOTE is the only symptom.
    (c) **Genuinely undefined objects in dead code, which are real bugs.** `computeBscoef`, `computePsiOutput`, `transformCoefficients`, `returnSpatialEffect`, `sample_beta_cpp`, `sample_beta_nocov_cpp` are functions that do not exist. `Ks_new`, `Xs_centers`, `X_ord`, `Tr`, `data_info`, `gt`, `gts` and the rest are undefined variables. All sit inside the functions already catalogued in items 1-3. These do not want imports; they want deleting.

    **Sequencing matters, and it is the main thing to get right.** Do items 1-3 (deprecate the dead code) **first**. Category (c) then disappears entirely, and category (b) shrinks to whatever the surviving functions use. Doing item 4 first means writing `globalVariables()` entries for code that is about to be deleted, and then deleting them again.

    **Upgraded 30 July: this is not only a NOTE. There is a latent runtime failure behind it.** Verified against a properly installed copy, with `tidyr` not attached: `pivot_longer` does not resolve from the package namespace at all (absent from both the namespace and its imports environment). It is reached only in the **categorical covariate** branch of `returnCovariateEffect_base()` (`R/jsdmfun.R:429`) and `plotCovariateEffect_base()` (`:612`); the numeric branch summarises to mean and quantiles and never calls it. So both exported GAM functions will fail with `could not find function "pivot_longer"` on any categorical occupancy covariate, for a real installed user.

    Nobody has hit it because every test and every vignette chunk uses numeric covariates. `create_covariates_matrix()` carries `cat_levels` and `is_numeric`, so categorical covariates are otherwise supported: this is a genuine gap, not an unreachable branch.

    **Phase 1: DONE 30 July 2026.** Added `setNames` and `rnbinom` to the `stats` list and a new `@importFrom tidyr pivot_longer`, both in `R/occJSDM-package.R`. **Verified by re-running `devtools::check()`: the "checking dependencies in R code" NOTE is gone entirely** (3 notes down to 2), and `pivot_longer`, `setNames` and `rnbinom` have left the undefined-globals list. `dnbinom` and `bs` remain in it, correctly, since they were deliberately not imported. Steps as executed:

    1.  Add to the consolidated tag in `R/occJSDM-package.R`, which is where the package's `@importFrom stats ...` line already lives: `setNames` and `rnbinom` to the existing `stats` list, and a new `@importFrom tidyr pivot_longer`.
    2.  **Do not import `dnbinom` or `bs`.** `dnbinom` is used only by `sample_rnb()` (item 5, not wired up) and `bs` only by `createSplinesObjects()` (dead). `bs` would additionally mean adding `splines` to `DESCRIPTION`, a new declared dependency existing purely to support code that should not ship. Both should go when that code goes.
    3.  Added a test, but asserting on the **imports environment** rather than calling the functions. A functional test of the categorical branch would pass under `devtools::test()` whether or not the import existed, because `load_all()` resolves unimported symbols through the global environment: that is precisely how this gap survived the suite. Checking `exists(fn, envir = parent.env(asNamespace("occJSDM")), inherits = FALSE)` holds under both `load_all()` and `R CMD check`. 167 tests passing.
    4.  Re-run `devtools::check()` and confirm the `Namespace in Imports field not imported from: tidyr` NOTE is gone, and that `pivot_longer`, `setNames` and `rnbinom` have left the undefined-globals list.

    **Phase 2: `globalVariables()` for the data-masked names. Do this after items 1-3, not before.**

    Measured 30 July: **38 functions generate undefined-global complaints, and 19 of them are dead** (zero callers, not exported): `computePredictiveProbs`, `loglik_sigma1`, `plotCovariateTrend`, `plotOccupancyStates`, `plotSpatialEffect`, `returnSpatialEffectMean`, and the twelve unused `sample_*` functions. Removing the dead code halves the complaint sources at a stroke, and takes the whole of category (c) with it.

    Doing Phase 2 first means enumerating NSE names for functions that are about to be deleted, then deleting the entries again. Wait.

    When it is time: one `utils::globalVariables()` call in `R/occJSDM-package.R`, next to the imports, with a comment saying what it is and why. Not scattered across files. There is currently no `globalVariables()` anywhere in the package, so this establishes the convention.

    **Success criterion for the whole item:** `devtools::check()` reports no NOTE in either "checking dependencies in R code" or "checking R code for possible problems". Shrinking the list is not finishing it.

    PHASE 1 DONE. CLAUDE TO DO PHASE 2 AFTER ITEMS 1-3

5.  **`sample_rnb()` cannot run as written** (new in `0abb104`, `R/jsdmfun.R:581-614`). Groundwork for the count-data item under *MEE paper*, not yet called from anywhere, but it has a scoping bug that will bite the moment it is wired up: `r_current <- rnb[s]` (`:590`) reads `rnb` inside the `sapply()` at `:588` whose result is *being assigned to* `rnb`, so at that point `rnb` does not exist in the function frame and lookup falls through to the namespace and fails with `object 'rnb' not found`. The current size vector needs to come in as an argument, e.g. `sample_rnb(z, eta, rnb, tune_sd = ...)`. Two more things to settle while there: `tune_sd = 5` is a random-walk SD on the *log* scale, so proposals land a factor of `exp(+/-10)` away and acceptance will be near zero (something in the 0.1-1 range is the usual starting point); and the prior terms are stubbed to `0` with the intended `dgamma()` commented out, referencing `prior_shape`/`prior_rate`, which are not defined anywhere. The Metropolis step itself looks right -- the `log(r_star) - log(r_current)` Jacobian is the correct correction for a log-scale random walk under a flat prior on `r`.

    ALEX's WORK IN PROGRESS FOR THE COUNTS

## **E. Draft of beta version listserv announcement (Doug)**

1.  Listserv announcement (beta release), drafted July 20 2026:

    > Subject: New R package (beta) — occJSDM, a combined occupancy and joint species distribution model
    >
    > Hi all,
    >
    > Announcing **occJSDM**, an R package for combining occupancy and joint species distribution modelling (<https://github.com/AlexDiana/occJSDM>). The package is \`eDNA-aware' and thus can be used on metabarcoding data.
    >
    > Note this is **beta software**, under active development.
    >
    > occJSDM extends the occPlus two-stage eDNA occupancy model of Ji et al. (2025, *Ecology Letters*, <doi:10.1111/ele.70302>) by adding a JSDM layer.
    >
    > Highlights:
    >
    > - Occupancy modelling: Accounts for both false-negative and false-positive error at two stages (field and lab), per species. Stage 1: estimates species eDNA collection probability in the field, given true eDNA presence at the site, and contamination probability, given true eDNA absence at the site. Stage 2: estimates species eDNA detection probability in the lab (i.e. successful DNA extraction, PCR, and sequencing), given successful eDNA collection in Stage 1, and contamination probability, given eDNA non-collection in Stage 1. In datasets where multiple primers have been used, each species' detection probability is estimated per primer (allowing one to compare each primer's efficiency for each species), while species occupancies are estimated using information across all primers. Both environmental and detection covariates are supported.
    > - JSDM: Integrates the occupancy model with a full-featured JSDM: species fit jointly, with latent-factor residual correlations. The JSDM optionally supports species traits shaping occupancy responses (trait X env interaction, aka \`fourth-corner analyses') and spatial autocorrelation (GP kernel) across sites.
    > - occJSDM not only fits a full-featured two-stage occupancy model (both field and PCR replicates required), but if given simpler study designs, can collapse to a classical occupancy model (field replicates only) or to a pure JSDM (no replicates).
    > - MCMC fitting with diagnostics, variance partitioning, ordination, and pairwise residual correlation outputs built in.
    > - occJSDM leverages the taxonomic breadth of eDNA datasets by using ordination (each site's position on the latent axes, and each species' loadings on those axes) to predict species occupancies. Thus, each species' predicted occupancy at a site is informed by the estimated occupancies of the other species at that site, thereby using co-occurrence structure. We also allow species to borrow strength from other species sharing similar traits, in contrast to the classical approach of having rare species borrow strength from abundant species, as is used in multi-species occupancy models.
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

    Ordered by expected speedup per unit of effort. Nothing here has been profiled -- worth running `profvis::profvis()` on a vignette- sized fit first to confirm where the time actually goes.

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

## B. Doug to dos

1.  **reproduce all Ecoletts results as a test of the package and decide whether to include in repo**

2.  **extensive testing on simulated datasets** -- **suite built and the R = 100 study run; three things remain.** What exists is summarised under *Completed* below; the authoritative specification and the results table are in `dev/simstudy/PLAN.md`.

    (a) **Re-run once group A items 1-3 are fixed.** This is the evidence the fixes worked. Without it they rest on the same code-reading that this exercise showed to be unreliable, and the R = 100 table in `PLAN.md` §12 is already partly stale -- it predates Alex's `42198d9`, which corrected the collection-covariate prior, so `beta_theta` should improve markedly. One command, a few hours:

        ```         
        Rscript dev/simstudy/run_study.R --cores=5 --caffeinate
        ```

    (b) **Decide the replicate count for the paper.** R = 100 was chosen to *detect* defects and did so decisively. Asserting *nominal* coverage in print is a claim about the absence of a small deviation and wants R = 200-500 (`PLAN.md` §9). The runner takes `R` as an argument.

    (c) **Decide how the results are presented** (`PLAN.md` open item 4). Currently a private Claude artifact (<https://claude.ai/code/artifact/ad3d46eb-1fd4-49b5-b795-6b71474ef1d5>), updated 29 July 2026 with the post-fix re-run and the before/after comparison; a pkgdown article is the obvious home -- see item 3.

    **One constraint carried from the bug list:** `l_s` is excluded from coverage checks because it is not recoverable while group A item 1 is open, so no cell of the study speaks to spatial range. Two earlier constraints have since lapsed -- `sigma_h` is now sampled (Fixed bugs 24) and the OpenMP RNG race is closed (Fixed bugs 26), so tier 1's "structural assertions only" rule can be revisited once reproducibility is confirmed on a multi-threaded platform.

3.  **Stand up a pkgdown site.** Discussed 28 July; not started.

    pkgdown turns the package into a static website -- function reference, both vignettes, README, changelog -- and GitHub Pages hosts it. `usethis::use_pkgdown_github_pages()` wires up both plus an Action to rebuild on push. Roughly half an hour, plus time shaping the reference index.

    **Three reasons it earns its place:**

    (a) *It closes a CRAN item.* `DESCRIPTION` has no `URL` or `BugReports` field (CRAN plan item 10 in `AGENTS.md`). A documentation site is a better `URL` than the bare repo.
    (b) *It is a much better landing page for the listserv announcement*, which currently points people at a GitHub file listing. A rendered reference and vignettes make a considerably stronger first impression for a beta release.
    (c) *It gives the validation write-up a home.* The plain-language guide to the test suite and the R = 100 results currently lives as a Claude artifact (<https://claude.ai/code/artifact/ad3d46eb-1fd4-49b5-b795-6b71474ef1d5>, private). As a pkgdown *article* it would sit with the package and be citable as supplementary material for the MEE paper.

    **Two practical constraints.**

    The site would live at `alexdiana.github.io/occJSDM`, since the repo is under Alex's account -- **Alex has to enable Pages**, it is not a setting Doug can change.

    And pkgdown builds every `@examples` block, so the **first build will fail on the group B functions that error unconditionally** (`predictNewSites()` among them). That is the same exposure `\donttest{}` creates under `R CMD check`, which is exactly why the CRAN plan sequences item 8 after the group B fixes. So do this after group B, or expect to `@examples`-guard several functions first.

    **Sequencing caveat:** consider whether to publish a site while group A items 3, 4 and 6 are open. The site would document functions whose credible intervals are currently overconfident, without saying so anywhere a reader would see.

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

10. ~~**Non-thread-safe RNG inside the OpenMP loops.**~~ **PARTIALLY FIXED.** `sample_beta_nocov_cpp_TS()` (`src/functions.cpp:308`) now calls the thread-safe `sample_beta_cpp_TS()` rather than the non-TS `sample_beta_cpp()`, which closes the race in `sample_betatheta_cpp_parallel()` (`src/functions.cpp:607`) -- the thread-safe path the audit described as "half-finished" is now wired up on that side. **The `samplePGvariables()` path (`src/jsdm.cpp:398`) is still affected** via `randinvg()`'s `R::rnorm` call; that residual is tracked as **group A item 5** above. *(Also not recorded at the time; recovered 27 July.)*

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

    **Confirmed by the 29 July re-run**: `beta_theta` improved in every cell, to 0.709-0.771, with two-thirds of the bias removed. It is still well short of nominal, so a second cause remains -- tracked as group A item 4.

26. ~~**Non-thread-safe RNG in the hottest OpenMP loop.**~~ **FIXED.** `randinvg()` (`src/jsdm.cpp:86`) now draws from the `thread_local` `rnorm()` rather than `R::rnorm`; the old line is commented out beside it. That closes the last hole on the `samplePGvariables()` path, whose Polya-Gamma helpers were already converted in `53c38f1`.

    **Two consequences worth acting on.** A fixed seed should now reproduce on Linux and Windows, which (i) removes the reason `sampleresults.rda` had to be refit on macOS or an OpenMP-disabled build, and (ii) unblocks CRAN item 16. It also means the tier-1 constraint in `dev/simstudy/PLAN.md` §5.1 -- structural assertions only, never numeric equality -- can be revisited; it was adopted solely because of this race. Verify reproducibility on a multi-threaded platform before relying on any of that.

27. ~~**Stage 2 hyperparameters documented as settable but never read.**~~ **FIXED.** `a_p`, `b_p`, `a_q` and `b_q` now read from `listPriors` (`R/runOccJSDM.R:751-754`), matching the behaviour `@param listPriors` already claimed, and the documented defaults were corrected to match the code (5/1 and 1/20). A user running a low-detection study can now override the prior without editing the package.

    **Only the wiring is closed.** What the default *should be* remains open and is a design decision, not a defect -- see group A item 3.

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

    It found six defects -- three while the tests were being written, three from the run -- none of which the static audit had caught. Four of the six are now fixed (Fixed bugs 24-27 and group B); the rest are group A items 1-3.
