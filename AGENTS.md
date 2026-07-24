# occJSDM

R package for fitting occupancy models to eDNA metabarcoding reads data (joint species distribution model / JSDM framework), accounting for a two-stage detection process (field collection + PCR amplification/lab detection), species traits, and optional spatial autocorrelation.

## Package layout

- `R/runOccJSDM.R` -- main model-fitting entry point (`runOccJSDM()`), plus data-prep helpers (`process_covariates()`, `createDataIdx()`, `get_param()`).
- `R/jsdmfun.R` -- core JSDM machinery, including `simulateData()` (the lower-level data simulator) and coefficient/variance-partitioning helpers (`computeBtcoef()`, `computePsiCoef()`, `computeVariancePartitioning()`).
- `R/simulateData.R` -- exports `simulateOccJSDMData()`, which now takes an explicit `model` argument (`"binary"`, `"occupancy"`, `"continuous"`, or `"two_stage"`) and directly returns `data_list` already in the `runOccJSDM()`-ready `info`/`OTU`/`traits` shape (as of Alex's `a9700ab` refactor, below). This merged the former `simulateOccJSDMDataGeneral()` into `simulateOccJSDMData()` and **removed** the separate `toRunOccJSDMFormat()` conversion step entirely -- neither `simulateOccJSDMDataGeneral()` nor `toRunOccJSDMFormat()` exist in the source anymore (dropped from `NAMESPACE` too), even though `tests/testthat/test-toRunOccJSDMFormat.R` and a `\link{toRunOccJSDMFormat}` in `R/data.R`'s roxygen docs still reference the removed function. **This is a live regression** -- the old two-step simulate-then-convert workflow described in `vignettes/simulateOccJSDMData.Rmd` may now be stale; check/update that vignette and the orphaned test before relying on either.
- `R/mcmcfun.R`, `R/output.R`, `R/diagnostics.R` -- MCMC sampling and post-processing/output/plotting functions used after `runOccJSDM()` (e.g. `returnOccupancyCovariates()`, `plotVariancePartitioning()`, `returnLatentPresences()`, `plotResidualCorrelationMatrix()`, `returnOrdinationScores()`/`plotOrdinationScores()`/`returnFactorLoadings()`/`plotFactorLoadings()`; see NAMESPACE for the full exported list).
- `src/jsdm.cpp`, `src/functions.cpp` (via `RcppExports.R`) -- Rcpp/Armadillo backend, including the spatial kernel `K2()` (squared-exponential GP kernel over site coordinates) and `sample_w_cpp()` (samples latent collection state `w` from continuous read intensities; currently unused in `runOccJSDM()` since `threshold < 1` is unsupported).
- `vignettes/occJSDM.Rmd` -- walkthrough of fitting a model with `runOccJSDM()` and using the output/plotting functions.
- `vignettes/simulateOccJSDMData.Rmd` -- walkthrough of `simulateOccJSDMData()`: parameter lists, simulating with/without spatial autocorrelation via `useSpatField`, checking variance partitioning, visualizing occupancy and the linear predictor spatially, and how trait covariates (`Tr`) shape species-specific occupancy coefficients via the `G` matrix. **Likely stale**: as of Alex's `a9700ab` refactor, `simulateOccJSDMData()` takes a `model` argument and returns `runOccJSDM()`-ready data directly, so any vignette content describing a separate `toRunOccJSDMFormat()` conversion step needs updating.
- `TODO.Rmd` -- structured outstanding feature list, organized as **v0.1.0-beta Public release** (Alex to dos / Doug to dos), **MEE paper** (Doug to dos / Alex to dos), and **Future versions**. See "Current work status" below for the current item list.
- `analysis/analysis.R` -- ad hoc analysis script (not part of package build).
- `tests/testthat/` -- unit tests (testthat edition 3, set up via `usethis::use_testthat(3)`); `test-toRunOccJSDMFormat.R` covers `toRunOccJSDMFormat()` (dimensions, column names, intercept handling, mismatched-dimension errors), but that function no longer exists in `R/` as of Alex's `a9700ab` refactor -- this test file is now orphaned and will fail/error until reconciled with the new `simulateOccJSDMData(model = ...)` API.
- `data/` -- `sampledata.rda` / `sampledata_orig.rda` (example dataset used in `occJSDM.Rmd`), `sampleresults.rda` (precomputed fit used by the vignette, `nchain=2, nburn=5000, niter=5000, nthin=1`; **refit and overwritten this session**, see "Current work status" below), `data_out_env_trait.rdata`, and `CaiWang_data/` (a data directory the user is actively reorganizing; `traitdata_caiwang.rdata` was recently removed in favor of this).
- `CITATION.cff` -- citation metadata; primary citation is the Ji et al. (2025) *Ecology Letters* methods paper, secondary is the software repo (`https://github.com/AlexDiana/occJSDM`). No CRAN/Zenodo listing exists, so this positioning follows standard practice for a research tool wrapping a published method.
- `README.md` -- minimal readme with install instructions, vignette pointers, and a "How to cite" section (since GitHub's citation button is easy to miss).
- `DESCRIPTION` -- `Authors@R` lists Alex Diana (`Alex.diana92@yahoo.it`, maintainer/`"cre"`) and Douglas W. Yu (`dougwyu@mac.com`, `"aut"`).

## Data structures

**Model input (`runOccJSDM()`)**: a list with `info` (data.frame, one row per PCR replicate, with `Site`, `Sample`, `Primer`, occupancy covariates `X_psi.*`, and collection covariates `X_theta.*`) and `OTU` (matrix of read counts, rows matching `info`, columns = species).

The `threshold` argument to `runOccJSDM()` controls how `OTU` is interpreted: `threshold >= 1` (default `1`) truncates reads to binary presence/absence (`OTU >= threshold` -\> detection) and samples `w` via `sample_w_cim_cipp()`. `threshold < 1` (including `threshold = 0`) is currently unsupported and raises an error (`stop("Threshold has to be greater than 0")`). Note that while `simulateOccJSDMData()` simulates continuous read counts for two-stage models via a log-normal mixture (`mu1`, `sigma1`, `mu0`, `sigma0`), `runOccJSDM()` truncates these to binary observations at fit time.

**Simulated data (`simulateOccJSDMData()`)**, as of Alex's `a9700ab` refactor, now takes an explicit `model` argument (one of `"binary"`, `"occupancy"`, `"continuous"`, `"two_stage"`) and returns `list(true_params, data_list)` where **`data_list` is already in `runOccJSDM()`-ready form**: `list(info, OTU, traits)` (no separate conversion step needed anymore -- see "Package layout" above for what changed and what's now stale as a result). `true_params` still carries the true occupancy/detection states and coefficients (`jsdmParams_true`, `beta_theta_true`, `z_true`, and for `"two_stage"` also `w_true`/`p_true`/`q_true`), and read counts for `"two_stage"` are generated the same way as before: `y = round(exp(logy1) - 1)` where `logy1` is drawn per-replicate from `Normal(mu1, sigma1)` (true detections), `Normal(mu0, sigma0)` (false-positive/contamination), or `0` (no reads) -- `mu1`/`sigma1`/`mu0`/`sigma0` still come from `list_params` (defaulting to `5`, `1`, `1.5`, `1`). **This section needs re-verification against the current source** (`R/simulateData.R`) next time simulated data is used, since the exact argument list (`list_datasettings`, `list_params`, `list_jsdmParams`) may have shifted along with the `model` argument addition -- confirm before relying on details above.

## Known issues in existing code

- **`model`-used-before-assignment bug in `runOccJSDM()` is fixed** (superseding the earlier note this replaced): `model <- inferDataModel(data)` is now assigned at line \~439, before its first use. Confirmed by successfully running `runOccJSDM()` live this session on both a working two-stage-derived occupancy dataset and, separately, a JSDM-only dataset that got well past this point before hitting the *different* `M`-not-found bug below. `TODO.Rmd` item 1.2 already marked fixed.
- **`runOccJSDM()object 'M' not found` bug for JSDM-only data is fixed** (Alex, commit `c1f9cec` "Fix M"): the "samples per site" `else` branch (`model %in% c("binary","continuous")`, no `Site` column) now sets `M <- NULL` alongside the existing `siteNames <- 1:n` (`R/runOccJSDM.R` \~565); `K` was already `NULL` in this case via the separate "pcr per marker" block's own `else` branch. The final `infos` list references both unconditionally but tolerates `NULL`. Same commit also added an internal `computeSpeciesDetected_M()` helper (`R/output.R`, `@noRd`) alongside the pre-existing `computeSpeciesDetected()`, plus a `set.seed(1)` inside `computeSpeciesDetected()` for reproducibility -- purpose/caller of the new function not yet investigated. TODO.Rmd's former item 1.7 for this bug has been removed accordingly (items renumbered; the remaining open simulate-side bug is now item 1.8, formerly 1.8 too but content unchanged).
- **`simulateOccJSDMData(model = "occupancy")` errors with `arguments imply differing number of rows: 0, N`.** Its `data_info` construction (`R/simulateData.R`) binds `Site = idx_z_k`/`X_psi = X_psi[idx_z_k,]`/`Xs = Xs[idx_z_k,]`, but `idx_z_k` is only populated by `createDataIdx()` when `twostage = TRUE` -- for `model == "occupancy"` (`twostage = FALSE`), `idx_z_k` is `NULL`. Confirmed live. This is now `TODO.Rmd` item 1.2 (renumbered as of `1c25529`, which pruned several other resolved Alex-to-do items -- see "Current work status"; the only other remaining Alex-to-do item, 1.1, is "Allow for different primers per sample", unrelated to this bug). Fix: use `idx_z_w` (populated regardless of `twostage`) instead of `idx_z_k`. **Workaround**: simulate via `model = "two_stage"` with trivial `P = 1`, `K = rep(1, N)`, which routes through the working `idx_z_k`-populated two-stage path and, once fed to `runOccJSDM()`, gets classified as classical occupancy (`Site` repeats, `Sample` doesn't) since `K=P=1` collapses detection to one Bernoulli draw per sample. **Still open** as of this session (July 24 2026) -- not addressed by `1c25529` or any other recent commit.
- **`computeSpeciesDetected_M()` no longer exists** (superseding the earlier note below about its purpose being uninvestigated): as of `1c25529`, `plotCumulativeSpeciesDetections()` was reworked to support both `M` and `K` on the x-axis via a single new `computeSpeciesDetected(beta_theta_output, p_output, M, K, primer, alpha)` that bootstraps Stage 1 (`w`, from `beta_theta_output`) and Stage 2 (detection, from `p_output`) draws directly per species/sample/replicate (`B = 200`), replacing the old Beta-approximation approach entirely. The function now requires `fitModel$infos$model == "two_stage"` (errors otherwise). See "Current work status" for full detail.
- `predictOccupancyProbs()` was renamed `predictNewSites()` (plus a new helper `createSpatialPredMatrix()`) in Alex's `a9700ab` refactor. **Status of the previously-flagged `X_ord` bug is now unclear** -- see the more specific note below (this bullet superseded).
- `computeAverageCollectionProbs()` and `computeConditionalSamplePresenceProbs()` were confirmed working (via live testing against fitted model objects) as long as the fitted model actually populated `results_output$theta_output` / `results_output$w_output` (always true for a fresh `runOccJSDM()` fit with default `summarisedLatentPresences = TRUE`).
- `computeMinESS()` (`R/diagnostics.R`) -- the `ESS_beta0psi`-not-populated bug (introduced in `af5ebe3`, TODO.Rmd item 1.5) **has been fixed** by Alex in `a9700ab`: `ESS_beta0psi` is now computed from `beta0_psi_output` directly and included in the final `min()`, and the function also guards against `NULL`/zero-row edge cases (`beta_theta_output` absent, `beta_psi_output` with 1 row, `L_output` with 0 rows). TODO.Rmd item removed accordingly.
- `plotFPTPStage2Rates()` (`R/output.R:748-824`) **ignores its own `primerName` argument**: `p_output`/`q_output` are 4-D arrays `[primer, species, niter, nchain]`, but the function does `apply(p_output, 2, quantile, ...)`, which pools over primer/iteration/chain regardless of `primerName` (which is only ever used to set a default, never to index/subset the array). There's also dead code -- `idx_speciesprimer <- stringr::str_match(...)` extracts species/primer indices from rownames but the result is never used. Net effect: the plotted CI is always pooled across all primers no matter what `primerName` is passed. Confirmed live (changing `primerName` produces an identical plot). Not yet in `TODO.Rmd` -- worth adding. Fix: subset `p_output[primer_idx, , , ]`/`q_output[primer_idx, , , ]` (matching `primerName` to `fitModel$infos$primerNames`) before computing quantiles.
- Trait-reading fragility in `runOccJSDM()`: it checks `data$traits` via partial name-matching against `data$traitsMatrix` (relying on `$`'s partial matching, since `traitsMatrix` is the only element of `data` starting with `"traits"`), not the unused `traitsMatrix` function argument. Still flagged as TODO.Rmd item (Alex to dos).
- Ordination (TODO.Rmd former item 1.1) has been addressed: `returnOrdination()`/`plotOrdinationScores()` were renamed/cleaned up as `returnOrdinationScores()`/`plotOrdinationScores()`, and new `returnFactorLoadings()`/`plotFactorLoadings()` (plus internal `returnFactorLoadings_jsdm()`/`plotFactorLoadings_jsdm()` in `R/jsdmfun.R`) now cover species loading scores, which were previously missing. Removed from TODO.Rmd accordingly.
- Counts model (`data_type == "counts"`, auto-detected from `OTU` values) is unsupported downstream: `stop("Counts model not supported yet")`. No explicit user-facing `count=` argument exists (see TODO.Rmd item under MEE paper / Alex to dos, "ability to analyse count data").
- MCMC diagnostics (Rhat/ESS/trace plots) ported from `~/src/Cowork/GLGS-eDNA` into `R/diagnostics.R` -- see "MCMC diagnostics" section below. This is now **committed and documented** (vignette section + `NAMESPACE`/`man/` regenerated), not stashed/dropped as an earlier version of this memory file claimed -- see "Current work status" below.
- `returnLatentPresences()`'s roxygen `@note` (`R/output.R` \~line 1671-1675) still describes an old bug (referencing an undefined `varPart_output`), but the current function body (verified by reading it directly) does **not** actually reference `varPart_output` anywhere and runs successfully in the vignette (`returnLatentPresences(fitmodel, idx_species = 1)`, confirmed via live rendering). The bug appears already fixed; only the stale `@note` comment needs deleting.
- `predictNewSites()`'s current roxygen docs (`R/output.R:1401-1422`) no longer contain an `@note` about the `X_ord` bug previously flagged in this memory file -- unclear whether the underlying bug is actually fixed or just undocumented; **needs live testing to confirm** before trusting either way.

## MCMC diagnostics

`R/diagnostics.R` has, in addition to the pre-existing unexported `computeESSparams()`/`computeMinESS()`:

- `as4d()` (internal) -- pads 2D (`[niter, nchain]`, scalar params like `theta0`/`sigmab`) or 3D (`[dim1, niter, nchain]`, e.g. species-only `B0_output`) posterior arrays up to the common `[dim1, dim2, niter, nchain]` shape, so every downstream function can treat any parameter uniformly regardless of how many index dimensions it has.
- `computeRhat(param_output)` -- per-element Gelman-Rubin Rhat via `coda::gelman.diag()`; returns `NA` for elements with a constant (zero-variance) chain or fewer than 2 chains, since `gelman.diag()` errors in those cases. **`@noRd` (unexported)** -- de-exported as of this session because it's redundant with `returnConvergenceDiagnostics()` as a user-facing function; remains an internal helper.
- `summarisePosterior(param_output, param_name, dimnames1, dimnames2)` -- tidy tibble (`param`, `idx1`, `idx2`, `label1`, `label2`, `mean`, `sd`, `q2.5`, `q97.5`, `rhat`, `ess`) for every element of a posterior array, pooling draws across chains for the summary stats. **`@noRd` (unexported)**, de-exported for the same reason as `computeRhat()`.
- `paramOutputToLong()` (internal) -- long-format tibble of every draw, feeding `plotTraceplot()`.
- `plotTraceplot(param_output, param_name, dimnames1, dimnames2)` (exported) -- faceted per-chain ggplot trace plots. Expects a 4-D array (`[covariate/primer, species, niter, nchain]`); to plot a single covariate/primer, slice it out with `drop = FALSE` first (see the vignette's `extractCovariateSlice()` helper below).
- `returnConvergenceDiagnostics(fitmodel)` (exported) -- assembles one tidy table across `beta0_psi` (`jsdm_output$B0_output`), `beta_psi` (`jsdm_output$B_output`), `beta_theta`, `p`, `q`, `theta0`, labelled with species/covariate/primer names from `fitmodel$infos$speciesNames`/`colnames(fitmodel$X_psi)`/`colnames(fitmodel$X_theta)`/`fitmodel$infos$primerNames`. Components absent from a given fit (e.g. `beta_theta`/`p`/`q` for a JSDM-only model) are silently skipped. Only exported user-facing function for tabular diagnostics (calls `computeRhat()`/`summarisePosterior()` internally).

**Only `plotTraceplot()` and `returnConvergenceDiagnostics()` are exported/user-facing** (confirmed in `NAMESPACE`); `computeRhat()`/`summarisePosterior()`/`computeESSparams()`/`computeMinESS()`/`as4d()`/`paramOutputToLong()` are all internal.

Not ported: `compare_to_true()`/`plot_estimated_vs_true()` from GLGS-eDNA, since those need `true_params` from a simulation rather than just a fitted model -- left as a possible future item. No testthat test exists yet for `R/diagnostics.R`.

`DESCRIPTION` has `coda`, `purrr`, `tibble` in `Imports` (this also fixed a pre-existing gap: `coda::` was already called in the old `computeESSparams()`/`computeMinESS()` but `coda` was never listed as a dependency).

`vignettes/occJSDM.Rmd`'s **"Model diagnostics" section** (added and committed this session, `882cffe`) demonstrates `returnConvergenceDiagnostics()` (with a table explaining each coefficient block and its `idx1`/`idx2`/`label1`/`label2` columns: `beta0_psi` species intercepts, `beta_psi` occupancy covariate slopes, `beta_theta` Stage 1/field-collection covariate effects including an auto-added `(Intercept)` column, `p`/`q` Stage 2/lab-PCR true-detection vs. false-positive rates per primer, `theta0` per-species Stage 1 false-positive collection rate) and a chunk flagging parameters exceeding `Rhat > 1.01` or `ess < 400`. A separate **"Traceplots for a single covariate" subsection** (moved by Doug in a follow-up commit, `e817bd1`, from right after "Model diagnostics" to the end of the vignette, after the "Latent presence table" section and before "Prediction at new sites") provides a reusable `extractCovariateSlice(param_output, covNames, covName)` helper (matches a covariate/primer name, slices with `drop = FALSE` to keep the array 4-D) plus a lookup table mapping each covariate/parameter type to its column-names source and coefficient array in `results_output`:

| Covariate/parameter | Column names | Coefficient array |
|----|----|----|
| `X_psi.*` (occupancy covariates) | `colnames(fitmodel$X_psi)` | `results_output$jsdm_output$B_output` (`beta_psi`) |
| species intercept | -- (species-only) | `results_output$jsdm_output$B0_output` (`beta0_psi`) |
| `X_theta.*` (collection covariates) | `colnames(fitmodel$X_theta)` | `results_output$beta_theta_output` (`beta_theta`) |
| primer-level rates | `fitmodel$infos$primerNames` (as character) | `results_output$p_output`, `results_output$q_output` |
| species-only false-positive rate | -- (species-only) | `results_output$theta0_output` (`theta0`) |

## Model inference logic (`inferDataModel()`, `R/runOccJSDM.R:110-155`)

`runOccJSDM()` classifies the fitted model type based on **row-level duplication of `Site`/`Sample` in `data$info`**, not directly on M/K/P: - No repeated values in `Site` → JSDM-only (M=K=P=1). - `Site` repeats, `Sample` does not → classical occupancy model (M\>1, K=P=1). - `Sample` repeats (due to K\>1 and/or P\>1) → two-stage eDNA model, including the case M=1 with K\>1 or P\>1 (a single site sampled with multiple PCR replicates/primers).

## Output functions of note

- `returnVariancePartitioning(fitModel)` -- per-species table of (Env, Spatial, Biotic, StDev) variance fractions; feeds `plotVariancePartitioning()`.
- `returnResidualCorrelationMatrix(fitModel, confidence = .95)` -- returns a 3 × S × S array (quantile × species × species) of posterior credible intervals for pairwise residual correlations. The 50% slice gives the median correlation matrix; comparing bounds' signs gives a significance flag.
- `plotCumulativeSpeciesDetections(fitModel, K, primer = 0, alpha = .95)` -- credible interval for cumulative species detected as a function of PCR replicates (K), via bootstrapped per-species Beta-distributed detection probabilities. No M-based (sample/visit-level) equivalent exists yet (TODO.Rmd, MEE paper / Alex to dos).
- Ordination is now more complete (as of Alex's `a9700ab` refactor): `returnOrdinationScores()`/`plotOrdinationScores()` (renamed from `returnOrdination()`) cover site factor scores, and new `returnFactorLoadings()`/`plotFactorLoadings()` cover species loading scores -- both exported with roxygen docs. The former TODO.Rmd item calling for this has been removed accordingly.

## Git and build artifacts

- **`src/*.o`/`src/*.so` tracking is fixed**: confirmed this session (`git ls-files 'src/*.o' 'src/*.so'` returns nothing) -- no longer tracked, despite the earlier note in this file claiming they were. Superseded.
- **`src/occJSDM.dll` untracking is fixed** (this session, `1b4b8ad`): `.gitignore` already had `src/*.dll` (so new dll files won't get re-added), but the file itself was still tracked from before that rule existed (`git check-ignore` doesn't flag already-tracked files even when a matching ignore rule exists). Ran `git rm --cached src/occJSDM.dll` and committed -- the file remains on disk locally (needed for Doug's Windows build) but git no longer tracks it. Still exists in git *history* (introduced by Alex's `a9700ab` commit, \~8.6 MB compiled Windows binary) -- purging history would need the same reclone coordination with Alex as the `traitdata_caiwang.rdata` item below.
- `.gitignore` also changed in `a9700ab`: added a blanket `/deprecated/` rule (replacing the narrower `deprecated/analysis/`).

## ggtern and plotting

- `plotVariancePartitioning()` (via `plotVarPart()` in `R/jsdmfun.R`, line 1466) produces a warning: "Ignoring unknown labels: L, T, R" because it calls `labs(L = "Environment", T = "Biotic", R = "Spatial")`. The `ggplot2::labs()` function doesn't recognize ternary-axis parameters `L`, `T`, `R`. **Fix**: Remove the `labs()` call entirely — the axis labels are correctly set by the aesthetic names (`Env`, `Biotic`, `Spatial`) in the `aes()` mapping and don't require additional specification. The theme elements `tern.axis.title.T`, `tern.axis.title.L`, `tern.axis.title.R` already control styling.

## `man/*.Rd` tracking

`man/*.Rd` files are currently tracked and committed in git (not gitignored), despite an earlier commit (`dd158a9`, prior to this session) whose message claimed intent to "stop generating/tracking man/\*.Rd" -- that intent was never followed through in `.gitignore`. Regenerating docs via `devtools::document()` after adding/changing roxygen tags will produce `man/*.Rd` diffs that should be committed (or the gitignore situation resolved) deliberately.

## Current work status

- **Later the same day (July 24 2026), a second AI-driven pass**: (1) resaved `data/sampleresults.rda` by assigning the live session's `fitmodel` object to `sampleresults` and running `usethis::use_data(sampleresults, overwrite = TRUE)`, at Doug's explicit request -- this **overrides the earlier caution in this file** ("do not regenerate `sampleresults.rda` without first confirming with Doug") since the request came directly from him this time. **Not yet confirmed**: whether the `fitmodel` used here was actually fit against the current (post-`876eb15`, lower-probability) `sampledata`, which would resolve the sampledata/sampleresults mismatch flagged below -- worth checking before trusting vignette output that depends on it. (2) Committed the previously-flagged uncommitted `vignettes/occJSDM.Rmd` diff as `763e738` ("Update occJSDM.Rmd: M+K cumulative detections prose, refreshed WAIC values") -- confirmed via `git diff` immediately before committing that it contained: the `devtools::load_all()` setup-chunk comment; the M+K cumulative-detections prose/chunk updates (removing Alex's "For Doug:" TODO note, adding explanatory sentences for each of the three `plotCumulativeSpeciesDetections()` calls); a reworded Rhat/ess convergence caveat (now "a few species" generically rather than naming `OTU_8`/`OTU_9`/`OTU_10` specifically, and noting "for all but two estimates" Rhat/ess are fine); removal of a stale `traitsMatrix = data$traits` argument from the `fitmodel2 <- runOccJSDM(...)` call (consistent with the trait-reading-fragility item elsewhere in this file -- traits are read from `data$traits` directly, not the `traitsMatrix` argument); and refreshed `extractWAIC()` values in the model-selection section (`fitmodel`: 31967.49 -\> 29920.78; `fitmodel2`: 31912.9 -\> 29896.9), consistent with `sampleresults`/`fitmodel` having been refit against the current lower-probability `sampledata` from `876eb15`. **This means the `sampledata`/`sampleresults` mismatch flagged in the July 24 entry below is very likely now resolved** (the WAIC values changing is consistent with a refit against new data), though this hasn't been independently re-verified end-to-end (e.g. by re-rendering the full vignette). No `TODO.Rmd` changes this pass despite it being open/visible in the editor.
- **This session (July 24 2026): reconciled memory with a batch of commits made directly by Alex/Doug since the last AI session; no AI-driven code changes this turn.** Working tree had one uncommitted change at session start: `vignettes/occJSDM.Rmd` (a `devtools::load_all()` comment added to the setup chunk, and the "Cumulative species detections" section's prose/chunks updated to reflect `plotCumulativeSpeciesDetections()` now supporting both `M` and `K` -- see below -- with Alex's old "For Doug:" TODO comment removed). Commits absorbed (oldest to newest):
  - `1c25529` "Added diagnostics and plot function" (Alex) -- the largest substantive change this batch:
    - **New `computeDiagnostics(results_output)`** (`R/diagnostics.R`, exported) replaces the old `computeMinESS()`-based single-number check at the end of `runOccJSDM()`: it now prints per-parameter-block Rhat/ESS diagnostics to the console (covering `beta0_psi`, `beta_psi`, `beta_theta`, `p`, `q`, `theta0`, skipping latent-state components `z_output`/`w_output`/`theta_output`/`idx_ls_output`/`varPart_output` and any all-constant block), warning if any block's max Rhat \> 1.1 or min ESS \< 50. Called automatically at the end of `runOccJSDM()` (`computeDiagnostics(results_output)`, replacing the old `minESS <- computeMinESS(results_output); if (minESS < 50) print(...)`). This resolves the former `TODO.Rmd` item "Report MCMC diagnostics by blocks of parameters."
    - **`plotReadIntensity()` de-exported** (`@export` tag removed) -- now internal/`@noRd` in effect (still has full roxygen otherwise; not re-verified whether `@noRd` was added explicitly or just `@export` dropped).
    - **`plotCumulativeSpeciesDetections()` reworked to support both `M` (environmental samples) and `K` (PCR replicates)**, not just `K` -- new signature `plotCumulativeSpeciesDetections(fitModel, M, K, primer = 0, alpha = .95, byK = TRUE)`. `byK = TRUE` (default) facets by `M` with `K` on the x-axis; `byK = FALSE` facets by `K` with `M` on the x-axis. Requires `fitModel$infos$model == "two_stage"` (errors otherwise). Internally rewritten: the old Beta-approximation-based `computeSpeciesDetected()`/`computeSpeciesDetected_M()` pair (which fit a Beta distribution to `p_output` moments) was replaced by a single new `computeSpeciesDetected()` with signature `(beta_theta_output, p_output, M, K, primer, alpha)` that directly bootstraps (`B = 200` reps, `set.seed(1)`) collection (`w`, via `beta_theta_output`/Stage 1) and detection (via `p_output`/Stage 2) draws per species/sample/PCR-replicate, rather than approximating via a Beta-Binomial mixture. **This makes the earlier memory note about `computeSpeciesDetected_M()`'s "purpose/caller ... not yet investigated" stale** -- it no longer exists; superseded by this unified function.
    - `process_covariates()` gained input validation: errors on `Inf`/`NaN` values in covariates, and on any covariate column with zero variance (constant values) with an informative message naming the offending columns.
    - `inferDataModel()`/`runOccJSDM()` switched from `print()` to `message()` throughout for status output (more idiomatic for a package -- `message()` respects `suppressMessages()`/goes to stderr, `print()` doesn't).
    - `runOccJSDM()`'s threshold validation tightened: now explicitly checks `is.numeric(threshold)` and `length(threshold) == 1` before the `threshold <= 0` check (previously just `threshold >= 1` vs. the `else` stop-with-message branch); same behavior for valid input, clearer errors for malformed input.
    - `runOccJSDM()` gained new validation: errors on duplicated species names (in `OTU` colnames or in `data$traits` rownames), and errors if `P` (primers per sample) differs across samples ("Cannot have different number of primers for each sample").
    - **`MCMCparams` list elements now have defaults** via `get_param()`: `nchain` defaults to 2, `nburn`/`niter` to 1000, `nthin` to 1, instead of requiring the caller to supply a complete `MCMCparams` list (previously `nchain <- MCMCparams$nchain` etc. with no fallback, which would silently propagate `NULL` if omitted). **This may affect `helper-sim.R`'s planned `minimal_mcmc` object in the not-yet-implemented test suite** -- worth confirming defaults still get overridden correctly when only some elements are supplied.
    - `sample_betatheta_cpp` calls replaced with `sample_betatheta_cpp_parallel` (new Rcpp export, `R/RcppExports.R`/`src/RcppExports.cpp`/`src/functions.cpp`) -- presumably a parallelized version; not yet benchmarked/verified for behavioral equivalence.
    - `TODO.Rmd` heavily pruned in this commit: several previously-open Alex-to-do items (trait-reading fragility fix confirmation, `.dll` gitignore, `predictNewSites()` for non-surveyed sites, `Xs` recycling bug, `computeMinESS()` copy-paste bug) were struck through as done/resolved and then the strikethrough text removed entirely (condensed to just the two remaining open items: "Allow for different primers per sample" and the `simulateOccJSDMData(model = "occupancy")` `idx_z_k` bug, renumbered 1.1/1.2). In the MEE-paper section, "Switch to integrated WAIC" was promoted higher in Alex's list and "`plotCumulativeSpeciesDetections()` for different values of M" was removed as an open item (now done, per the `output.R` change above).
    - Also touched `vignettes/occJSDM.Rmd` (18 lines) -- superseded by Doug's follow-up commits below for the final vignette wording.
  - `4e786e8` "Update TODO.Rmd" (Doug, minor tweak before the above landed/was reconciled).
  - `84a3f20` "Fix multi-line @importFrom roxygen tag (must be single line)" (Doug) -- `R/occJSDM-package.R`'s `@importFrom` block had been split across multiple lines at some point (roxygen requires each `@importFrom pkg fn1 fn2 ...` tag on one line); collapsed back to a single line.
  - `40f9596` "Fix stale computeSpeciesDetected_MK call to match actual function name" (Doug) -- `plotCumulativeSpeciesDetections()` in `1c25529` called a function named `computeSpeciesDetected_MK()`, which doesn't exist (the actual helper, per the diff above, is named `computeSpeciesDetected()`); fixed the call site and regenerated `man/plotCumulativeSpeciesDetections.Rd`.
  - `6ca8dc0` "Ignore scratch data/sampledata_highp_highthetabaseline.rda" (Doug) -- added to `.gitignore`; an ad hoc/scratch simulated dataset variant (higher `p`/`theta_baseline`), not tracked.
  - `901c50c` "Fix copy-pasted roxygen docs for computeDiagnostics()" (Doug) -- `computeDiagnostics()`'s roxygen `@param`/`@return` had been copy-pasted from a different function (describing `param_output`/`param_name`/`dimnames1`/`dimnames2` args and a ggplot return, none of which match its actual single-argument `results_output`, side-effect/console-print signature). Rewrote the docs to match reality and added explicit `invisible(NULL)` at the end of the function body to match the corrected `@return NULL, invisibly` doc. Regenerated `man/computeDiagnostics.Rd`/`NAMESPACE`.
  - `45dada8` "Regenerate occJSDM-package.Rd to include Alex Diana as author" (Doug) -- `DESCRIPTION`'s `Authors@R` already listed Alex Diana, but `man/occJSDM-package.Rd` (auto-generated authors section) hadn't been regenerated since; ran `devtools::document()` to sync it.
  - `876eb15` "new sampledata with lower collection and detection probabilities" (Doug) -- regenerated `data/sampledata.rda` (26,723 -\> 22,593 bytes) using new `list_params` values in `vignettes/simulateOccJSDMData.Rmd`: `p` (Stage 2 true-positive PCR detection rate) lowered from `runif(P*S, 0.6, 0.95)` to `runif(P*S, 0.3, 0.6)`, and `theta_baseline` (Stage 1 baseline field-collection probability given present) lowered from `runif(S, 0.3, 0.7)` to `runif(S, 0.15, 0.4)`. Rationale (per new vignette comments): with the old higher rates, a single PCR replicate/field sample already detected most present species, so cumulative-detection plots (`plotCumulativeSpeciesDetections()`) barely changed with more replicates `K`/samples `M` -- the new lower rates make detection genuinely improve with more effort, which is the point of that plot. **`data/sampleresults.rda` was NOT regenerated in this commit** -- it was fit against the *old* higher-probability `sampledata.rda`, so `vignettes/occJSDM.Rmd` (which loads `sampleresults` as a precomputed fit) may now be showing a fit against stale/mismatched simulated data. **Needs checking**: whether `occJSDM.Rmd` still references `sampledata` directly anywhere in a way that would break, or whether it only uses `sampleresults` (in which case it's merely inconsistent-but-functional, not broken).
- **Not yet addressed this session**: the uncommitted `vignettes/occJSDM.Rmd` diff (still open in the editor, not yet committed as of this memory update) -- just the `devtools::load_all()` comment and the "Cumulative species detections" section's M+K prose/chunk update, consistent with the `1c25529`/`40f9596` code changes above. The `data/sampleresults.rda`-vs-`sampledata.rda` mismatch flagged above is unresolved. Item 4 of the CRAN submission plan (`sampleresults.rda` is 62MB) is presumably still open and now further complicated by this potential mismatch -- **do not regenerate `sampleresults.rda` without first confirming with Doug**, since he may already be tracking this via a still-open editor session.
- **Previous session (July 22 2026): reconciled memory with commits made directly by Alex/Doug since the last AI session; no AI-driven code changes that turn.** Working tree was already clean and `main`/`origin/main` identical at session start. Commits absorbed: `c1f9cec` "Fix M" (Alex) -- see "Known issues" above, this is the fix for the previously-open `object 'M' not found` bug; `96a226d` "Update AGENTS.md" and `13bb50d` "Update TODO.Rmd" (Doug) -- `TODO.Rmd`'s "Alex to dos" sub-list was renumbered (old item 1.7, the `M`-bug, removed/resolved; old item 1.8, the `simulateOccJSDMData(model = "occupancy")` row-mismatch bug, is now item 1.8 with unchanged content) and gained a new first sub-item, "Report MCMC diagnostics by blocks of parameters" (not otherwise elaborated in `TODO.Rmd` or here yet). **Net effect: the only remaining open simulate-side bug is the `simulateOccJSDMData(model = "occupancy")` row-mismatch** (`idx_z_k` is `NULL` when `twostage = FALSE`) -- any references elsewhere in this file to "two known bugs" or "`TODO.Rmd` items 1.7/1.8" (both plural) are now stale and should be read as referring only to the remaining occupancy-model bug (former item 1.8). `vignettes/occJSDM.Rmd` shows as modified/open in the editor this session (edited directly by Doug, per this file's usual pattern) -- diff not yet inspected.
- **This session (July 21 2026): triaged `devtools::check()`'s 1 warning / 3 notes, one item fixed, two deliberately deferred.** Ran a fresh `devtools::check()` and kept the result in the live session as `res` (an `rcmdcheck` object -- `res$warnings`/`res$notes` are useful for pulling exact diagnostic text without re-running the full check). Findings and actions:
  - **LICENSE NOTE**: see CRAN submission plan item 11 below -- deleted then restored the file, ending in "keep + accept the cosmetic NOTE."
  - **`tidyr` unused-import NOTE**: unchanged/still open this session -- discussed the fix (`tidyr::pivot_longer()` qualification inside `plotCovariateTrend()`) but did not implement it.
  - **"no visible binding" NOTE, base-R-function half**: **fixed and committed** (`ef77c4c`, "Import base stats functions to resolve no-visible-binding NOTE") -- added `@importFrom stats binomial cov cov2cor dbinom dgamma dnorm dpois glm kmeans logLik median predict quantile rbeta rbinom reorder rgamma rnorm rpois runif sd var` to `R/occJSDM-package.R`'s package-doc roxygen block (which already had one `@importFrom` tag, for `Rcpp::evalCpp`), then ran `devtools::document()` and confirmed all 22 `importFrom(stats, ...)` lines landed correctly in `NAMESPACE`.
  - **"no visible binding" NOTE, NSE-column half** (e.g. `Species`/`x`/`y`/`Site`/`` `2.5%` `` flagged across many `plot*()` functions in `R/jsdmfun.R`/`R/output.R`): **the user explicitly decided to leave this unfixed** -- discussed but did not implement either standard remedy (`utils::globalVariables()` vs. rewriting call sites with `.data$colname`). Treat as an accepted, permanent NOTE going forward, not a to-do.
  - **C++ bitwise-operator warning**: pinpointed the exact 7 warning sites (see CRAN submission plan item 9) but the fix itself was not started -- a `bash`/`sed` attempt to view the relevant `functions.cpp` lines was declined by the user mid-session.
  - Net effect on `check()` note/warning count not yet re-measured after `ef77c4c` -- **rerun `devtools::check()`** to get a fresh baseline before further triage.
- **Most recent commits** (as of this session, July 20 2026, newest first, all now pushed to `origin/main` -- confirmed `git rev-parse main origin/main` identical): `f594ea4` "Add placeholder test after removing orphaned test-toRunOccJSDMFormat.R" (added `tests/testthat/test-placeholder.R` with a trivial passing test, since `testthat::test_check()` errors with "No test files found" if `tests/testthat/` is completely empty); `00233a1` "Remove orphaned test-toRunOccJSDMFormat.R" (deleted the test file entirely via `git rm` -- it tested the removed `toRunOccJSDMFormat()` function and errored before ever reaching it; comprehensive tests against the current `simulateOccJSDMData(model = ...)` API are a separate future task); `1a898ad` "Fix vignettes to use library(occJSDM) instead of devtools::load_all()" (both vignette setup chunks called `devtools::load_all()`, which has the side effect of attaching Imports to the search path and silently allowing unqualified tidyverse/ggplot2 calls; `R CMD check` rebuilds vignettes in an isolated environment without `devtools` installed, so this failed with "there is no package called devtools" -- fixed by switching to `library(occJSDM)` and qualifying two bare `count()` calls as `dplyr::count()` and two bare `ylim()` calls as `ggplot2::ylim()`; also repaired incidental corruption where `library(occJSDM)` had been merged into `traitsSummary`); `2b2426d` "Fix DESCRIPTION dependency declarations" (added `FastGP` to Imports -- used in `R/output.R`/`R/jsdmfun.R` for GP kernel math but never declared, masked locally since it was already installed; added `patchwork` to Suggests -- used in the vignette for combining plots; bumped `Depends` to `R (>= 4.1.0)` since the code already uses the native pipe `|>` and `\(...)` lambda shorthand); `45afe09` "Remove redundant @import tidyverse (undeclared dependency)" (`R/simulateData.R` had `@import tidyverse` but `tidyverse` itself was never in DESCRIPTION -- DESCRIPTION already declares the specific tidyverse packages actually used, so the tag was just removed and `NAMESPACE` regenerated). Earlier history (`44c0b03` and everything before it, described below) is unchanged from previous sessions' notes.
- **This session's work (July 20 2026)**: two separate efforts. (1) **Git history purge**: used `git filter-repo --path data/traitdata_caiwang.rdata --invert-paths --force` to permanently remove `data/traitdata_caiwang.rdata` from all commits in `main`'s lineage from `11c5449` onward, then force-pushed to `origin` (`AlexDiana/occJSDM`) -- confirmed the file is gone from the working tree, local history, and the remote. **This resolves the previously-open TODO.Rmd item** ("purge traitdata_caiwang.rdata from git history") described in the older TODO.Rmd-structure bullet below, which is now stale on that point. **Action needed**: Alex needs to re-clone or `git fetch origin && git reset --hard origin/main && git clean -fd` his local copy, since his history now diverges from the rewritten remote -- **this message has not yet been sent to him**; still an open task for the user. Note the parallel `src/occJSDM.dll`-in-history item (see "Git and build artifacts" above) was *not* addressed this session and remains open, unpurged. (2) **Fixed `devtools::check()` from erroring out entirely**: after killing \~60 stuck R processes and 143 leftover temp build dirs left over from repeated retries (which had made `check()` merely *look* permanently hung, not actually broken -- `tools::buildVignettes()` and a full `R CMD INSTALL` both completed in \~30s in isolation), found and fixed 5 real packaging bugs (the 5 commits above: `45afe09`, `2b2426d`, `1a898ad`, `00233a1`, `f594ea4`). **Current `devtools::check()` status: 0 errors, 2 warnings, 5 notes.** Warnings: benign C++ compiler warnings, and undocumented Rd `@param`s for `predictNewSites` (missing `useEnvCov`/`useSpatial`/`useBiotic`) and `runOccJSDM` (missing `traitsMatrix`) -- **not yet fixed, low priority**. Notes are cosmetic (file timestamps, LICENSE/DESCRIPTION mismatch, non-standard top-level files, CITATION.cff placement, unused declared imports, no-visible-global-variable) and not yet triaged individually.
- **Rd documentation warnings fixed** (commit `8d36c26`, "Document missing Rd params for predictNewSites and runOccJSDM"): added `@param` docs for `predictNewSites`'s `useEnvCov`/`useSpatial`/`useBiotic` and `runOccJSDM`'s `traitsMatrix` (noting it's currently unused directly -- traits are actually read from `data$traits`, per the trait-reading-fragility item elsewhere in this file). Regenerated `.Rd` files via `devtools::document()`; confirmed via `tools::checkDocFiles()` and a full `devtools::check()` rerun that both "checking for missing documentation entries" and the general documentation warning are gone. **`devtools::check()` improved from 0 errors/2 warnings/5 notes to 0 errors/1 warning/5 notes** (the remaining warning is the pre-existing benign C++ compiler warnings, unrelated to docs).
- **"Unused declared imports" note partially triaged** (commit `2b2089d`, "Update DESCRIPTION", committed directly by Doug in parallel with an identical AI-drafted edit): removed `cli` and `tidyselect` from `DESCRIPTION`'s `Imports` -- `cli` was used nowhere in `R/`, and the only `tidyselect`-style calls (`contains()`, `starts_with()` inside `dplyr::across()`) are actually re-exported by `dplyr` (`import(dplyr)` already covers them), so `tidyselect` wasn't a real direct dependency either. **`tidyr` was deliberately left in `Imports`**, even though it still triggers the note (`Namespace in Imports field not imported from: 'tidyr'`) -- its only usage is inside `plotCovariateTrend()` (`R/jsdmfun.R:1495`), an unexported function that is never called anywhere in the package and is explicitly **in-progress** (per Doug) rather than dead code, so it was left untouched at Doug's request rather than deleted or wired up with a proper `@importFrom tidyr pivot_longer`. `plotCovariateTrend()` computes a GAM/spline-style marginal effect of one occupancy covariate on the *linear predictor* (not probability) scale via `createSplinesMatrixSingleCov()`/`list_ns`, with no posterior-uncertainty band and all species overlaid in one plot -- **not a duplicate of** the exported `plotOccupancyGradient()`/`returnOccupancyGradient()` (which operate on the probability scale with full posterior credible intervals, faceted per species, and pull directly from a fitted model's posterior draws). The spline helpers it depends on (`createSplinesObjects()`/`createSplinesMatrix()`/`createSplinesMatrixSingleCov()`, `R/jsdmfun.R` \~line 240-272) trace back to an internal `usingSplines` argument of the lower-level `simulateData()`, but `runOccJSDM()` itself has no corresponding spline/GAM argument today -- so this looks like scaffolding for a not-yet-exposed spline/GAM occupancy-covariate feature. **Still open**: decide the eventual fate of `plotCovariateTrend()`/spline support (finish it, or delete later) and, until then, the `tidyr`-unused-import note will persist as expected.
- **Outstanding from this session, not yet done**: (a) notify Alex about the git-history rewrite and recovery steps (drafted but unsent); (b) triage the remaining 5 `check()` notes (file timestamps, `LICENSE`/DESCRIPTION mismatch, non-standard top-level files, CITATION.cff placement, the now-narrowed `tidyr`-unused-import note, no-visible-global-variable/function issues e.g. `plotCovariateTrend`'s `pivot_longer`/`x`/`Value`/`Species` bindings); (c) write a comprehensive test suite (the placeholder test only keeps `testthat` scaffolding alive, it doesn't test anything real); (d) the two known bugs below (`object 'M' not found` for JSDM-only models; `simulateOccJSDMData(model = "occupancy")` row-mismatch) are unrelated to this session's work and were not investigated further -- still open per `TODO.Rmd` items 1.7/1.8; (e) `plotCovariateTrend()`'s in-progress spline/GAM feature (see above) has no corresponding `runOccJSDM()` argument yet -- unclear if/when it'll be finished.
- Earlier: `44c0b03` "Deduplicate sites in predictNewSites example; move section earlier in vignette" (this session -- committed a direct edit Doug made outside this session's tool calls: deduplicates `data$info` by `Site` for the `predictNewSites()` example, `1782 -> 100` unique sites, via `site_rows <- !duplicated(data$info$Site)`, and moved the whole "Prediction at new sites" section earlier in the vignette, right after the model-diagnostics/traceplots material, with a note that the function "is still under active development"); `e817bd1` "Update occJSDM.Rmd" (Doug, moved the "Traceplots for a single covariate" subsection from immediately after "Model diagnostics" to the end of the vignette, after "Latent presence table..." and before "Prediction at new sites" -- done directly by Doug outside this session's tool calls); `882cffe` "Add Model diagnostics vignette section; de-export computeRhat/summarisePosterior" (see "MCMC diagnostics" above for full detail); `c4e24c0` "Refit and update sampleresults.rda; fix Xs doc mismatch"; `2a943c7` "Regenerate sampledata.rda with fixed Xs coordinates"; `bf0c3ee` "Mark TODO 1.1 (trait-reading fragility in runOccJSDM()) as fixed"; `cb3cc94` "Fix Xs recycling bug and modernize simulateOccJSDMData vignette"; `0bb5efd` "Regenerate sampledata to new info/OTU/traits shape". Earlier history (`a9700ab` Alex's refactor and everything before it) is unchanged from previous sessions' notes below.
- **This session's work** (analysis + bug-hunting; `44c0b03` committed, `TODO.Rmd` updated with two new items): (1) Investigated why `predictNewSites()` applied to the original 100 sites doesn't exactly reproduce `computePredictiveOccupancyProbs()`'s fitted values (differences up to \~0.26 in probability) -- found no clear species-level pattern tying the discrepancy to factor-loading or spatial-coefficient strength (weak/inconsistent correlations, n=10 species). (2) To disentangle "estimation noise" from "`predictNewSites()` structurally can't reconstruct latent effects," simulated a dataset with a *known* true occupancy probability (`true_psi <- plogis(true_eta)`, from `simulateOccJSDMData()`'s internal linear predictor) and compared both estimators against it: MAE-vs-truth was nearly identical (fitted 0.202 vs. `predictNewSites` median 0.206), with correlated errors (r=0.94) and similar shrinkage-toward-0.5 bias for both -- concluding the fitted/`predictNewSites` mismatch reflects shared estimation uncertainty in a modest-sized dataset (n=50 sites), not a `predictNewSites`-specific flaw. (3) While building that simulation, discovered two new bugs (see "Known issues" above and `TODO.Rmd` items 1.7/1.8): `runOccJSDM()` errors with `object 'M' not found` for JSDM-only (`model="binary"/"continuous"`) data, and `simulateOccJSDMData(model = "occupancy")` errors with a row-mismatch due to a `NULL` `idx_z_k`. Worked around both by simulating via `model = "two_stage"` with trivial `P=1`, `K=rep(1,N)` instead.
- **Working tree**: clean as of the end of this session (`44c0b03` committed). Re-check with `git status` at the start of any new session, since the vignette is open in the editor and edited directly by Doug between AI sessions.
- Two doc-only staleness items found while working in this area this session (see "Known issues" above for detail): `returnLatentPresences()`'s `@note` roxygen comment describes an old bug that the current function body doesn't actually have; `predictNewSites()`'s docs no longer mention the previously-flagged `X_ord` bug (unclear if actually fixed).
- **`data/sampleresults.rda` refit** (from a previous session, `c4e24c0`): re-ran the `vignettes/occJSDM.Rmd` workflow (`data <- sampledata`, `set.seed(3947)`, drop 3 samples down to 297 samples/1,782 rows, fit with the documented `MCMCparams`), saved over `sampleresults` via `usethis::use_data(sampleresults, overwrite = TRUE)`. Also fixed a doc-only bug in `R/data.R`'s `@format` roxygen for `sampleresults` (documented `X_s`, actual element is `Xs`).
- **TODO.Rmd structure** (current as of `1c25529`/this session, heavily condensed from the earlier structure described below -- see file for full detail): organized as **v0.1.0-beta Public release** (1. Alex to dos, now just 2 items: "Allow for different primers per sample" [new], and the `simulateOccJSDMData(model = "occupancy")` `idx_z_k` bug [renumbered 1.2, unchanged content, see "Known issues"]; 2. Doug to dos, both items now marked **done**: purging `traitdata_caiwang.rdata` from git history via `git filter-repo`, and the "model diagnostics functions from GLGS-eDNA repo" item, whose done-note now also covers the vignette section and `computeRhat()`/`summarisePosterior()` de-export, and explicitly notes a testthat test for `R/diagnostics.R` is still not done), **MEE paper** (1. Doug to dos: reproduce Ecology Letters results as a package test; 2. Alex to dos, reordered/updated by `1c25529`: Overleaf math vignette, **"Switch to integrated WAIC" [new, promoted to position 2]**, GAMs for JSDM, count-data support, source-sink inference scenario, remove space's effect on env covariates, site-level "variation" partitioning -- **the former "M-based `plotCumulativeSpeciesDetections()`" item was removed since it's now done**, see "Current work status"), a new **Outreach** section (a drafted listserv announcement for the beta release, dated July 20 2026, not yet sent), and **Future versions** (spike-in abundance-change estimation, model selection via regularisation/shrinkage e.g. for geospatial foundation model embeddings, and a new item: parallelisation for speedup -- plausibly related to the new `sample_betatheta_cpp_parallel` Rcpp export from `1c25529`). `TODO.Rmd` still does **not** mention the `runOccJSDM()` `model`-used-before-assignment bug -- treat as fixed per the "Known issues" section above unless re-broken.
- **`interim_todo.Rmd`** (an untracked scratch file, never in git) proposed where to add several not-yet-exposed functions to the vignette (`returnOccupancyRates`, `returnTraitsCoeff`, `plotFPTPStage2Rates`, `computeAverageCollectionProbs`, `computeConditionalSamplePresenceProbs`, `thinOutput`); its diagnostics-related row is now resolved (the "Model diagnostics" placeholder it referenced has been filled in). **The file no longer exists on disk** (removed by Doug during this session, after this memory file's last read of it) -- if any of its other suggestions (the functions listed above) are still wanted in the vignette, they'll need to be re-derived or asked about again, since that content isn't preserved anywhere else.

## Notes

- When editing files that are also open in the RStudio editor, prefer re-reading the file immediately before and after edits -- a prior session saw a file get corrupted/duplicated after an `edit` call, likely from an editor/disk desync; rewriting the whole file with `write` fixed it.
- This repo's git history shows a pattern of committing by concern (e.g. testthat setup, generated docs, vignette content fixes kept separate from unrelated whitespace-only diffs) rather than one large commit -- worth continuing when making multiple unrelated changes.
- TODO.Rmd tracks near-term and future development priorities; keep synchronized with actual work being done.
- `doc/` (vignette build output, already in `.Rbuildignore` via `^doc$`) was deleted from the working tree this session -- it's regenerated automatically by `devtools::build_vignettes()`/`R CMD build`, so deleting it is safe and has no effect on `git` (it was never tracked).
- `Meta/vignette.rds` and `dev/` are two other directories clarified in a recent session: `Meta/` is a local, gitignored (`.gitignore:4`, `/Meta/`) build artifact created whenever vignettes are built, holding vignette metadata for `vignette()` lookups -- not something to commit or worry about. `dev/` is a tracked (not gitignored) scratch directory of ad hoc dev scripts (`dev/example.R`, `dev/runmodel.R`) outside `R/`, so not part of the built package; `dev/example.R` calls several plotting-function names (`plotDetectionCovariates()`, `plotOrdinationCovariates()`, `plotFPDetectionRates()`, etc.) that don't match current exports and looks stale/pre-refactor, similar to `analysis/analysis.R` -- `dev/runmodel.R` (37KB) has not yet been read/triaged.

## CRAN submission plan

Goal: get occJSDM accepted on CRAN. Investigated `devtools::check()` output plus DESCRIPTION/data/docs directly (July 20 2026 session) and found the items below. Split into blocking issues (CRAN will very likely reject without these) and should-fix issues (not hard blockers but expected/likely to draw reviewer pushback).

### Blocking

1.  ~~**`Description:` field is still template placeholder text**~~ **FIXED** (`d0f1650`, refined `083df56`): rewritten to describe the actual method (occupancy/JSDM for eDNA metabarcoding read data, two-stage detection process, traits, spatial autocorrelation via a Gaussian process, MCMC fitting), citing Ji et al. (2025) `<doi:10.1111/ele.70302>`. Doug hand-edited wording afterward (fixed "reads data" -\> "read data", tightened phrasing) -- final text lives in `DESCRIPTION` lines 17-24.
2.  ~~**`Title:` not in title case**~~ **FIXED** (`083df56`): Doug chose his own title over the AI-drafted one -- current `Title:` is `Simultaneous Fitting of an eDNA-Aware Occupancy and Joint Species Distribution Model` (title-cased at Doug's request; note it's long, \~90 chars, but CRAN doesn't hard-reject on length). The `occJSDM-package.Rd` internal title was separately set to `occJSDM: Occupancy and Joint Species Distribution Models for eDNA Metabarcoding Data` (`8b99219`) -- these two titles (DESCRIPTION's `Title:` vs. the package-doc `\title`) are intentionally worded differently and don't need to match.
3.  ~~**`Remotes: kassambara/ggcorrplot`**~~ **FIXED** (`083df56`): confirmed via web search that `ggcorrplot` is on CRAN (current release 0.1.4.1) and that the only usage in the codebase (`ggcorrplot::ggcorrplot()` in `R/jsdmfun.R`, called with `method`/`type`/`lab`/`lab_size`/`colors`/`title`) is standard API present in the CRAN release -- confirmed live by running `plotResidualCorrelationMatrix(fitmodel)` successfully against the installed `0.1.4.999` GitHub dev build. Removed the `Remotes:` field entirely.
4.  **`data/sampleresults.rda` is 62 MB -- still OPEN, not yet fixed.** Profiled live this session: `object.size(sampleresults)` = 65.1 MB total, of which `results_output` = 64.7 MB and within that `jsdm_output` = 57 MB. The two largest single arrays are `Bs_output` (22.9 MB, dims `[30, 10, 5000, 2]`) and `U_output` (15.3 MB, dims `[100, 2, 5000, 2]`) -- i.e. full raw MCMC draws, 10,000 total (`niter=5000 x nchain=2`, `nthin=1`), stored as double-precision 4-D arrays with no thinning applied before saving. Discussed 5 options with Doug: (a) post-hoc thin the existing arrays along the `niter` axis (e.g. keep 1-in-25 draws, no refit needed, \~62MB -\> \~2.5MB) -- **recommended**, fast, doesn't change the vignette's fitted numbers structurally; (b) refit with much smaller `MCMCparams` (e.g. `niter=200`) -- same size win but requires re-running MCMC (minutes) and re-rendering vignette prose; (c) audit and drop `jsdm_output` sub-components unused by vignette/output functions (`Gs_output`/`Cs_output`/`U_output`/`As_output`/`tau_output`/`C_output`/`sigmab_output`/`sigmabs_output`/`sigmah_output`/`idx_ls_output`, \~25-30MB combined) -- combine with (a); (d) don't ship a precomputed fit at all, fit a tiny model live at vignette build time instead -- cleanest long-term but touches vignette code and the dataset's now-documented `\value`; (e) compression alone (already `LazyDataCompression: xz`) won't meaningfully shrink near-random MCMC draws. **Recommended combo**: (a) + (c). **Not yet implemented** -- next step is to check which `jsdm_output` components the vignette/output functions actually touch, then write a thin+prune script and regenerate `sampleresults.rda`.
    - **Separately discussed**: an external-download alternative for cases where a package genuinely needs a large file. CRAN disallows downloading anything during install/build/`check` (must work offline), but allows a package to ship a user-facing function (e.g. `get_full_sampleresults()`) that downloads+caches a file at runtime via `tools::R_user_dir()`, guarded in examples/vignettes (e.g. `\donttest{}`, `curl::has_internet()` checks) so CRAN's automated checks never depend on network access. Host candidates: GitHub Releases (simple) or Zenodo/OSF (gets a DOI, fits this project's existing paper-citation/reproducibility framing). **Proposed two-track approach** (not yet actioned): ship the thinned/pruned `sampleresults.rda` (above) for CRAN-checkable examples/vignettes, and optionally host the full untouched 62MB fit externally via an optional download helper for exact paper-reproduction purposes -- these are independent and don't block each other.
5.  ~~**`rstan` is a dead, extremely heavy dependency**~~ **FIXED** (`d0f1650`): removed `rstan (>= 2.26.0)` from `Imports:` and deleted the `rstan::rstan_options(auto_write = TRUE)` call from `R/zzz.R`'s `.onLoad()` (which now only sets `options(mc.cores = ...)`). Confirmed via `grep -rn "rstan"` across `R/`, `DESCRIPTION`, `NAMESPACE` that no references remain.
6.  ~~**Missing `\value` tags on 3 Rd files**~~ **FIXED** (`6a7313f`): added `@return` roxygen tags to `sampledata`/`sampleresults` (in `R/data.R`, short restatements of their existing `@format` descriptions) and to `occJSDM-package` (in `R/occJSDM-package.R`, alongside replacing the placeholder "Package documentation." text with a real description mirroring the DESCRIPTION rewrite, then retitled per Doug's wording in `8b99219`). Regenerated via `devtools::document()`; confirmed via `grep -L "\\value" man/*.Rd` that no Rd files remain missing `\value`.

### Should fix before submitting

7.  ~~**Duplicate `ggplot2` line in `Imports:`**~~ **FIXED** (`083df56`): removed the unversioned duplicate, kept `ggplot2 (>= 3.4.0)`.
8.  **17 Rd files with no `\examples`** -- not a hard CRAN rule but routinely requested by reviewers on first submission; given MCMC runtime, wrap slow examples in `\donttest{}` rather than `\dontrun{}`. **Not yet started.**
9.  **Remaining C++ compiler warning** (bitwise-vs-logical operators in `src/functions.cpp`/`src/jsdm.cpp`) -- CRAN's periodic checks use stricter compiler flags/sanitizers than a local `check()`; fix now rather than risk a later "please fix or be archived" email. Pinpointed live (July 21 2026 session) via `res$warnings` on a stored `rcmdcheck` object from `devtools::check()`: exactly 7 "use of bitwise '&' with boolean operands" warnings, at `functions.cpp:569`, `functions.cpp:571`, `functions.cpp:573`, `jsdm.cpp:215`, `jsdm.cpp:231`, `jsdm.cpp:247`, `jsdm.cpp:264` (an 8th warning in the same block, `-Wfixed-enum-extension` from an R header, is unrelated/benign and not fixable from this package). **Not yet fixed** -- a first attempt to view `functions.cpp:560-580` via `bash`/`sed` was declined by the user; needs the `read` tool instead, then change `&` to `&&` where operands are scalar booleans (need to check each site isn't actually intended as an Armadillo elementwise `&` on vectors, which would need different handling).
10. **No `URL`/`BugReports` fields in DESCRIPTION** -- expected practice; add pointing at the GitHub repo (`https://github.com/AlexDiana/occJSDM`). **Not yet started.**
11. ~~**LICENSE file mismatch**~~ **DECIDED/ACCEPTED (not "fixed")** (July 21 2026 session): briefly deleted the top-level `LICENSE` file (`19a1650`/`19a1756`) since `License: Apache License (>= 2)` is self-contained and doesn't reference it, which cleared the `check()` NOTE ("File LICENSE is not mentioned in the DESCRIPTION file"). **Reverted** (`f274f41`) after realizing GitHub's license-detector badge needs the physical `LICENSE` file present at the repo root regardless of the `DESCRIPTION` `License:` field -- restored the original Apache 2.0 full-text file. Net decision: **keep `LICENSE`, deliberately accept the resulting cosmetic `check()` NOTE** rather than add `+ file LICENSE` (confirmed via earlier testing that this triggers a worse "restrictions not permitted" NOTE for this license type). This item is now closed by decision, not by elimination of the NOTE.
12. **No `cran-comments.md`** -- customary for first submissions; should explain the compiler warning (if any remains after item 9) and test environments checked (win-builder / R-hub / mac-builder, given compiled code via Rcpp/RcppArmadillo). **Not yet started.**
13. **Run `R CMD check --as-cran`**, not just `devtools::check()`, as the final local gate -- `--as-cran` adds a few extra checks (e.g. installed size) that plain `check()` skips. **Not yet run since the fixes above landed** -- should rerun to confirm current state before further work.
14. **Multi-platform pre-check** -- submit a test build to win-builder or R-hub before the real CRAN submission, since compiled-code packages (via `LinkingTo: Rcpp, RcppArmadillo`) fail there more often than pure-R ones. **Not yet started.**
15. ~~**README install instructions omitted `build_vignettes = TRUE`**~~ **FIXED** (`3a42058`): a user reported `vignette("occJSDM", package = "occJSDM")`/`vignette("simulateOccJSDMData", package = "occJSDM")` failing after a fresh install. Root cause: `remotes::install_github()` defaults to `build_vignettes = FALSE` (confirmed by inspecting `args(remotes::install_github)` and `remotes:::normalize_build_opts()`, which adds `--no-build-vignettes` unless `build_vignettes = TRUE` is passed explicitly), so the README's plain `remotes::install_github("AlexDiana/occJSDM")` never built vignettes for GitHub installs. Verified the fix live: built the package with `pkgbuild::build(vignettes = TRUE)`, installed the tarball into a clean library via `R CMD INSTALL`, and confirmed `vignette(package = "occJSDM")` then lists both `occJSDM` and `simulateOccJSDMData`. Updated `README.md`'s install snippet to `remotes::install_github("AlexDiana/occJSDM", build_vignettes = TRUE)` with an explanatory note. **This is purely a GitHub-install gotcha, not a CRAN-submission blocker**: `install.packages()` from CRAN (and `R CMD check`/`R CMD build` during the submission process itself) always builds vignettes as a matter of course, since CRAN's build pipeline runs `R CMD build` (which builds vignettes by default, no `remotes`-style opt-in flag involved) before publishing the source tarball -- so once/if `occJSDM` is on CRAN, ordinary `install.packages("occJSDM")` users won't hit this issue regardless of the README wording. The README fix only matters for the current pre-CRAN GitHub-install workflow, but is good practice to keep regardless (e.g. for anyone installing from a fork or a non-CRAN branch later).

### Status summary

Of the 6 blocking items, 5 are fixed (1, 2, 3, 5, 6) and committed (`d0f1650`, `083df56`, `6a7313f`, `8b99219`) -- confirmed `git rev-parse main origin/main` identical, all pushed to `origin/main`. **Item 4 (`sampleresults.rda` size) is the sole remaining blocker** and needs the thin+prune work described above before it's resolved. Of the should-fix items, 7 and 15 are fixed; 11 is closed by deliberate decision (keep `LICENSE`, accept the cosmetic NOTE); 9 is pinpointed (exact 7 warning lines identified) but not yet fixed; 8, 10, 12-14 are all still open. Separately (outside this numbered list, from `devtools::check()`'s general note/warning triage, not CRAN-blocking): the `stats`-import portion of the "no visible binding" NOTE is fixed (`ef77c4c`); the NSE-column portion is deliberately left unfixed by the user's decision; the `tidyr`-unused-import NOTE remains open. Recommended next steps, in order: (i) finish item 4 (the actual CRAN blocker), (ii) fix item 9's 7 pinpointed C++ bitwise-operator warnings, (iii) rerun `devtools::check()`/`R CMD check --as-cran` to get a fresh baseline, (iv) work through remaining should-fix items 8, 10, 12-14.

## Testing Strategy Implementation (July 22, 2026)

**Rationale**: The package currently has only a placeholder test file (`tests/testthat/test-placeholder.R`). A comprehensive test suite is needed to (1) validate core functionality, (2) prevent regression of fixed bugs, (3) ensure CRAN compliance, and (4) support ongoing development.

### Test Architecture

**Organizational structure** (hierarchical, matching package components):

```         
tests/testthat/
├── test-simulation.R    # Data generation (simulateOccJSDMData)
├── test-model-fit.R     # Core modeling (runOccJSDM)
├── test-outputs.R       # Post-processing functions (all ~34 exported functions)
├── test-diagnostics.R   # MCMC checks (returnConvergenceDiagnostics, plotTraceplot)
├── test-integration.R   # Full workflow tests
├── helper-sim.R         # Shared helper: minimal valid simulate() arguments
└── fixtures/            # Pre-fitted model objects (.rds), generated once and
                         # committed; regenerate with tests/testthat/make-fixtures.R
```

**Key design decisions**: - Use minimal datasets (`n=20` sites, `S=5` species, `nchain=1`) to keep total test runtime under 3 minutes - Cache pre-fitted models as `.rds` files in `fixtures/` to avoid repeated MCMC runs; load with `readRDS(testthat::test_path("fixtures/fit_twostage.rds"))`; regenerate with a dedicated `make-fixtures.R` script (not run by `testthat`, only run manually when the model interface changes) - Test all `simulateOccJSDMData()` model types (`"binary"`, `"occupancy"`, `"continuous"`, `"two_stage"`) - Regression tests for known fixed bugs: `M` variable (JSDM-only, fixed `c1f9cec`) - The `model = "occupancy"` `idx_z_k` bug (TODO item 1.8) is **still open** -- its test must use `skip()` until fixed, not `expect_equal()` (a passing assertion for a broken function is worse than no test)

### Minimal Working Simulate Call

All four arguments to `simulateOccJSDMData()` are **required** with no defaults. All `list_datasettings` keys are **lowercase**. A minimal valid call for tests (use this as `helper-sim.R`):

``` r
minimal_sim_args <- function(model = "two_stage") {
  n <- 20; S <- 5; M <- rep(2, n); P <- 1; K <- rep(2, n * P * max(M))
  list(
    list_datasettings = list(
      n = n, S = S, g = 1, M = M, P = P, K = K,
      ncov_psi = 1, ncov_theta = 1
    ),
    list_params = list(
      p       = matrix(0.7,  P, S),
      q       = matrix(0.05, P, S),
      theta0  = rep(0.05, S),
      theta_baseline = 0.4
    ),
    list_jsdmParams = list(
      gt = 1, d = 2, ds = 0, tau = 1,
      sigma_b = 1, sigma_bs = 0, sigma_ts = 0,
      sigma_h = 0, sigma_s = 0, l_s = NULL,
      useSpatField = FALSE
    ),
    model = model
  )
}

minimal_mcmc <- list(nchain = 1, nburn = 250, niter = 500, nthin = 1)
```

### Core Test Cases

**Data Simulation** (`test-simulation.R`):

``` r
test_that("simulateOccJSDMData returns correct top-level structure", {
  for (model in c("binary", "continuous", "two_stage")) {
    args <- minimal_sim_args(model)
    sim  <- do.call(simulateOccJSDMData, args)
    expect_named(sim, c("true_params", "data_list"))
    expect_true(all(c("info", "OTU", "traits") %in% names(sim$data_list)))
    expect_equal(nrow(sim$data_list$info), nrow(sim$data_list$OTU))
  }
})

test_that("simulateOccJSDMData(model='occupancy') bug (TODO 1.8) is still present", {
  # Remove skip() and flip to expect_no_error() once idx_z_k fix is committed
  skip("TODO item 1.8: idx_z_k is NULL for model='occupancy', fix pending")
  args <- minimal_sim_args("occupancy")
  expect_no_error(do.call(simulateOccJSDMData, args))
})

test_that("binary/continuous models produce correct OTU dimensions", {
  for (model in c("binary", "continuous")) {
    args <- minimal_sim_args(model)
    sim  <- do.call(simulateOccJSDMData, args)
    n    <- args$list_datasettings$n
    S    <- args$list_datasettings$S
    expect_equal(dim(sim$data_list$OTU), c(n, S))
  }
})
```

**Model Fitting** (`test-model-fit.R`):

``` r
test_that("runOccJSDM returns expected output structure for two_stage data", {
  args <- minimal_sim_args("two_stage")
  sim  <- do.call(simulateOccJSDMData, args)
  fit  <- runOccJSDM(sim$data_list, MCMCparams = minimal_mcmc)
  expect_type(fit, "list")
  expect_true(all(c("results_output", "infos", "X_psi") %in% names(fit)))
})

test_that("runOccJSDM JSDM-only (binary) no longer errors on M variable", {
  # Regression test for object 'M' not found bug, fixed in commit c1f9cec
  args <- minimal_sim_args("binary")
  sim  <- do.call(simulateOccJSDMData, args)
  expect_no_error(runOccJSDM(sim$data_list, MCMCparams = minimal_mcmc))
})

test_that("runOccJSDM rejects threshold < 1", {
  # runOccJSDM requires threshold >= 1; threshold < 1 raises an error
  args <- minimal_sim_args("two_stage")
  sim  <- do.call(simulateOccJSDMData, args)
  expect_error(
    runOccJSDM(sim$data_list, threshold = 0, MCMCparams = minimal_mcmc),
    "Threshold has to be greater than 0"
  )
})

test_that("runOccJSDM errors informatively on malformed input", {
  expect_error(runOccJSDM(list()), "data_info or OTU missing")
  args <- minimal_sim_args("two_stage")
  sim  <- do.call(simulateOccJSDMData, args)
  bad  <- sim$data_list
  bad$OTU <- bad$OTU[-1, ]   # row mismatch
  expect_error(runOccJSDM(bad), "cannot have different number of rows")
})
```

**Output Functions** (`test-outputs.R`):

Tests in this file load a pre-fitted fixture to avoid re-running MCMC. The fixture covers all output function code paths (two-stage model, traits, no spatial field). The \~34 exported functions in `output.R` are grouped below by what they need from the fixture:

| Group | Functions | Key assertions |
|----|----|----|
| Occupancy covariates | `returnOccupancyCovariates`, `plotOccupancyCovariates`, `returnOccupancyGradient`, `plotOccupancyGradient` | Returns data.frame/ggplot; nrow == n_species |
| Detection/collection | `returnCollectionCovariates`, `plotCollectionCovariates`, `returnOccupancyRates`, `plotOccupancyRates`, `plotCollectionRates` | Returns data.frame; nrow matches species count |
| Stage 2 rates | `plotFPTPStage2Rates`, `plotDetectionRates`, `plotStage1FPRates`, `plotStage2FPRates` | Returns ggplot; `primerName` arg actually subsets (regression for known bug) |
| Traits | `returnTraitsCoeff`, `plotTraitsCoefficients` | Returns data.frame; only run when traits present |
| Residual correlation | `returnResidualCorrelationMatrix`, `plotResidualCorrelationMatrix` | Returns 3×S×S array; median slice is symmetric |
| Ordination | `returnOrdinationScores`, `plotOrdinationScores`, `returnFactorLoadings`, `plotFactorLoadings` | Returns data.frame/ggplot; nrow == n_sites / n_species |
| Variance partitioning | `returnVariancePartitioning`, `plotVariancePartitioning` | Row fractions sum to \~1 per species |
| Prediction | `predictNewSites`, `computePredictiveOccupancyProbs` | Returns matrix with nrow == n_sites; values in [0,1] |
| Latent state | `returnLatentPresences`, `plotLatentPresences`, `computeAverageCollectionProbs`, `computeConditionalSamplePresenceProbs` | Dimensions match n_species / n_sites |
| Diagnostics | `returnConvergenceDiagnostics`, `plotTraceplot` | Columns rhat/ess present; no NA rows for well-behaved chains |
| Utilities | `thinOutput`, `extractWAIC` | Thinned object smaller; WAIC is a scalar |
| Cumulative detections | `plotCumulativeSpeciesDetections` | Returns ggplot |

``` r
# Load fixture once for all output tests
fit_ts <- readRDS(testthat::test_path("fixtures/fit_twostage.rds"))

test_that("returnVariancePartitioning fractions sum to approximately 1", {
  vp <- returnVariancePartitioning(fit_ts)
  frac_cols <- c("Env", "Spatial", "Biotic")
  row_sums <- rowSums(vp[, frac_cols])
  expect_true(all(abs(row_sums - 1) < 0.01))
})

test_that("returnResidualCorrelationMatrix returns 3 x S x S array", {
  S   <- length(fit_ts$infos$speciesNames)
  rcm <- returnResidualCorrelationMatrix(fit_ts)
  expect_equal(dim(rcm), c(3, S, S))
  # Median slice should be symmetric
  expect_equal(rcm[2,,], t(rcm[2,,]))
})

test_that("returnConvergenceDiagnostics has expected columns", {
  diag <- returnConvergenceDiagnostics(fit_ts)
  expect_true(all(c("param", "rhat", "ess") %in% colnames(diag)))
})
```

**Integration Test** (`test-integration.R`):

``` r
test_that("simulate -> fit -> output pipeline completes without error", {
  set.seed(42)
  args <- minimal_sim_args("two_stage")
  sim  <- do.call(simulateOccJSDMData, args)
  fit  <- runOccJSDM(sim$data_list, MCMCparams = minimal_mcmc)

  expect_no_error(returnVariancePartitioning(fit))
  expect_no_error(returnConvergenceDiagnostics(fit))
  expect_no_error(returnResidualCorrelationMatrix(fit))
  expect_no_error(returnLatentPresences(fit, idx_species = 1))
})
```

### CRAN Compliance

- 17 Rd files currently have no `\examples` -- add `\donttest{}` wrappers using the `minimal_sim_args()` helper to keep example runtime acceptable.
- `devtools::check()` must pass with 0 errors and 0 new warnings before any PR merge.
- Individual test: \< 30 seconds; full suite (excluding fixture generation): \< 3 minutes.

### CI/CD Integration

**GitHub Actions** (`.github/workflows/R-CMD-check.yaml`):

``` yaml
name: R-CMD-check
on: [push, pull_request]
jobs:
  check:
    runs-on: ${{ matrix.os }}
    strategy:
      matrix:
        os: [ubuntu-latest, macos-latest, windows-latest]
        r-version: ['4.5.0']
    steps:
      - uses: actions/checkout@v4
      - uses: r-lib/actions/setup-r@v2
        with:
          r-version: ${{ matrix.r-version }}
      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          extra-packages: any::rcmdcheck, any::covr
      - uses: r-lib/actions/check-r-package@v2
        with:
          args: 'c("--no-manual", "--no-build-vignettes")'
      - name: Test coverage (ubuntu only)
        if: matrix.os == 'ubuntu-latest'
        run: Rscript -e 'covr::codecov()'
```

Note: Rcpp/RcppArmadillo compilation is handled automatically by `r-lib/actions/setup-r-dependencies@v2` on all three platforms; no manual `apt-get` step is needed with the current action versions.

### Current Status

**Not yet started**: - All test files listed above - `helper-sim.R` (minimal arg constructor) - `make-fixtures.R` (fixture generation script) - `fixtures/` directory and `.rds` files - CI/CD pipeline configuration

**Predecessor**: `tests/testthat/test-placeholder.R` should remain until at least `test-simulation.R` and `test-model-fit.R` are in place, to keep `testthat::test_check()` from erroring with "No test files found".

**Next steps** (recommended order): 1. Write `helper-sim.R` and verify `minimal_sim_args()` actually runs end-to-end in the console 2. Write and run `make-fixtures.R` to generate `fixtures/fit_twostage.rds` 3. Implement `test-simulation.R` 4. Implement `test-model-fit.R` (including threshold rejection check) 5. Implement `test-outputs.R` using the fixture 6. Implement `test-diagnostics.R` and `test-integration.R` 7. Configure GitHub Actions CI/CD 8. Remove `test-placeholder.R`
