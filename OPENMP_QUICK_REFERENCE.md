---

editor_options: 
  markdown: 
    wrap: sentence
---

# OpenMP Quick Reference for occJSDM

## Current State

| Component | Location | Status | Issue |
|---------------------|-------------------|-----------------|-----------------|
| **Parallel region 1** | `src/jsdm.cpp:398` | Active | RNG not thread-safe |
| **Parallel region 2** | `src/functions.cpp:607` | Active | RNG not thread-safe |
| **Thread capping** | `src/jsdm.cpp:649` | Commented out | Uncapped (CRAN violation) |
| **Global state** | `R/zzz.R:3` | Active | Changes user's session (CRAN violation) |
| **Build config** | `src/Makevars` | Unconditional OpenMP | No guards |

## The Three CRAN Violations

### 1. Uncapped Threads (CRAN Item #16)

- **Problem**: No limit on thread count; can use all 16+ cores on CRAN infrastructure
- **Why it's bad**: CRAN policy: maximum 2 cores during checks
- **Symptom**: Code passes locally (2 cores on Mac), fails on CRAN (16 cores on Linux)
- **Current attempt to fix**: `omp_set_num_threads()` is commented out, never runs

### 2. RNG Thread-Safety (Bug A.6)

- **Problem**: `rpg()` and `arma::randn()` call non-thread-safe R RNG
- **Why it's bad**: Race condition on global RNG state in parallel threads
- **Why it works locally**: macOS clang silently ignores `#pragma omp` without `-fopenmp`
- **Why it fails on CRAN**: Linux/Windows actually parallelize, exposing the race condition
- **Symptom**: Intermittent failures, platform-dependent, hard to reproduce

### 3. Global State Pollution (Item #17a + #21)

- **Problem**: `.onLoad()` sets `options(mc.cores = ...)`, polluting user's session
- **Why it's bad**: CRAN policy: packages shall not modify global state
- **Secondary problem**: `parallel::detectCores()` is undeclared (undeclared import)
- **Symptom**: Every other package's `mc.cores` read gets occJSDM's value; silent behavior change

## Two Ways to Fix

### Option A: Fix OpenMP (Harder, Less Reliable)

```         
1. Add SystemRequirements: OpenMP in DESCRIPTION
2. Uncomment omp_set_num_threads(2) and add guards (#ifdef _OPENMP)
3. Make RNG thread-safe using thread-local seeds
4. Remove .onLoad() call
5. Result: More complex code, still risky on some platforms, slower
```

### Option B: Remove OpenMP, Parallelize Over Chains (Recommended ✅)

```         
1. Extract chain loop into run_mcmc_chain(chain_id, ...)
2. Use parallel::makeCluster(cores, type = "PSOCK") + parLapply()
3. Implement per-chain output merging (WAIC, running means)
4. Remove R/zzz.R entirely
5. Remove #pragma omp from C++ code
6. Remove OpenMP from src/Makevars
7. Result: Cleaner, safer, faster, portable, thread-safe by design
```

## Why Option B is Better

| Criterion | Option A | Option B |
|--------------------------|-----------------------|-----------------------|
| **CRAN compliance** | ✅ If done perfectly | ✅ By construction |
| **RNG safety** | ⚠️ Needs careful code | ✅ Automatic (separate processes) |
| **Platform portability** | ⚠️ OpenMP quirks per compiler | ✅ PSOCK universal |
| **Deadlock risk** | ⚠️ High (fork + BLAS + OpenMP) | ✅ None (separate processes) |
| **Code complexity** | Medium (guards, thread-local RNG) | Medium-Low (just reorganize loop) |
| **Performance** | Serial `samplePGvariables()`, no speedup | Chain-level parallelism, N×speedup |
| **Backward compatibility** | ✅ Full | ✅ Full (opt-in `cores` arg) |
| **Risk if wrong** | ⚠️ Intermittent race conditions | ✅ Deterministic (easy to debug) |

## Implementation Overview (Option B)

### Step 1: Extract Loop

``` r
# New file: R/parallel_mcmc.R
run_mcmc_chain <- function(chain_id, data_list, X_psi, X_theta, ..., seed) {
  # Copy the entire 300-line chain loop body here
  # Change references to *_output -> *_output_chain
  # Return: list(B_output_chain, U_output_chain, ..., waic_parts_chain)
}
```

### Step 2: Implement Merging

``` r
# Merge WAIC across chains using parallel variance formula
# Merge running means using simple summation
# Rebuild full *_output arrays from per-chain arrays
```

### Step 3: Use Parallelism

``` r
if (cores == 1L) {
  # Serial path: existing loop over chains (unchanged)
  for (chain in 1:nchain) { ... }
} else {
  # Parallel path: new code
  cl <- makeCluster(cores, type = "PSOCK")
  clusterSetRNGStream(cl, seed = ...)
  chain_outputs <- parLapply(cl, 1:nchain, run_mcmc_chain, ...)
  stopCluster(cl)
}
```

### Step 4: Remove OpenMP

``` bash
# src/Makevars: remove SHLIB_OPENMP_CXXFLAGS
# src/jsdm.cpp:398: delete #pragma omp parallel for
# src/functions.cpp:607: delete #pragma omp parallel for
# R/zzz.R: delete entire file
```

## Code Locations

### Parallel Regions (To Remove)

**`src/jsdm.cpp:390–410`** — `samplePGvariables()`

``` cpp
#pragma omp parallel for collapse(2)  // ← DELETE THIS
for(int i = 0; i < n1; i++){
  for(int j = 0; j < n2; j++){
    Omega_mat(i,j) = rpg(1, Xbeta(i,j));
  }
}
// Delete commented-out alternative loop too (lines 404–407)
```

**`src/functions.cpp:600–630`** — `sample_betatheta_cpp_parallel()`

``` cpp
#pragma omp parallel for  // ← DELETE THIS
for (int s = 0; s < S; ++s) {
  // loop body
}
```

### Build Config (To Simplify)

**`src/Makevars:1–2`**

``` makefile
# BEFORE:
PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS)
PKG_LIBS = $(SHLIB_OPENMP_CXXFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

# AFTER:
PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
```

### Global State (To Remove)

**`R/zzz.R`**

``` r
# BEFORE:
.onLoad <- function(libname, pkgname) {
  options(mc.cores = min(2L, parallel::detectCores()))
}

# AFTER:
# Delete this file entirely (or leave empty)
```

### Chain Loop (To Extract)

**`R/runOccJSDM.R:895–~1200`** (approximate; exact line depends on version)

``` r
# Current code:
for (chain in 1:nchain) {
  # ~300 lines of MCMC logic
  # Writes to B_output_chain, U_output_chain, etc.
}

# Refactor to:
# 1. Move entire loop body to run_mcmc_chain() function
# 2. Return list of *_output_chain objects
# 3. Replace loop with parLapply() call (or serial loop if cores==1)
# 4. Merge outputs
```

## Testing Checklist

- [ ] **Build without error**: `devtools::load_all()` after removing OpenMP
- [ ] **Serial path works**: `runOccJSDM(..., cores = 1L)` produces identical output to current code
- [ ] **Parallel path works**: `runOccJSDM(..., cores = 2L)` on 2-core machine completes without error
- [ ] **Reproducibility**: Same seed → identical output (serial or parallel)
- [ ] **WAIC is correct**: Manual MCMC on tiny dataset, compare merged WAIC to expected value
- [ ] **Running means correct**: Compare `z_output_mean` to value computed from raw draws
- [ ] **No compiler warnings**: `devtools::check()` should not report OpenMP warnings
- [ ] **Vignette renders**: `rmarkdown::render("vignettes/occJSDM.Rmd")` produces expected output
- [ ] **No CRAN violations**: `R CMD check --as-cran` passes (resolves items #16, A.6, #17a, #21)

## References in TODO.Rmd

- **Item A.6** (line \~200): RNG thread-safety bug
- **Item D.1** (line \~260): Parallelize over chains using PSOCK + parLapply
- **Item D.2** (line \~303): Remove .onLoad() global state
- **Item #16** (CRAN plan, line \~556): OpenMP uncapped threads
- **Item #17a** (CRAN plan, line \~555): Global state pollution
- **Item #21** (CRAN plan, line \~560): Undeclared `parallel` import

## Common Pitfalls

| Pitfall | What Goes Wrong | How to Avoid |
|------------------|------------------------------|-------------------------|
| **Forgetting to return values from `run_mcmc_chain()`** | Chain outputs are lost | Always return a named list with all `*_output_chain` arrays |
| **Not merging WAIC correctly** | WAIC is wrong or missing | Use parallel variance formula (mean of chains' means, var of chain means + mean of chain vars) |
| **Worker doesn't have package** | `clusterEvalQ()` step omitted | Call `clusterEvalQ(cl, library(occJSDM))` after creating cluster |
| **Worker doesn't have data** | Parallel code errors with "object not found" | Use `parLapply(..., SIMPLIFY = FALSE)` with all needed args passed explicitly |
| **RNG seed not set** | Parallel runs are not reproducible | Use `clusterSetRNGStream(cl, seed = ...)` |
| **Cluster not cleaned up** | File handles leak, R session hangs | Use `on.exit(stopCluster(cl))` or `tryCatch(..., finally = stopCluster(cl))` |
| **Forgetting to remove `#pragma omp`** | OpenMP still compiled, warnings remain | Search `src/` for all occurrences of "pragma omp" and "parallel for" |

## Commits Needed

```         
1. [parallel-mcmc] Extract chain loop and implement parLapply path
2. [parallel-merge] Implement WAIC and running mean merging
3. [remove-openmp] Remove OpenMP from src/Makevars and C++ code
4. [cleanup-zzz] Delete R/zzz.R (global state pollution)
5. [docs] Update roxygen for runOccJSDM cores argument
6. [tests] Add unit tests for serial/parallel equivalence
7. [vignette] Update occJSDM.Rmd to show cores argument usage
```

## FAQ

**Q: Will the code be slower without OpenMP?** No. Chain-level parallelism provides near-linear speedup. On a 2-core machine with 2 chains, you get \~2× speedup. Within-iteration OpenMP provided minimal benefit anyway (e.g., 10% per iteration) and is error-prone.

**Q: What if someone only wants one chain?** Serial by default (`cores = 1L`), so no performance change. If they want parallelism, they must pass `cores > 1L` explicitly.

**Q: Will this break CRAN examples?** No. CRAN never passes `cores > 1L`, so examples use serial path (current behavior). Output is identical.

**Q: How long will the refactor take?** 7–11 hours total. Chain loop extraction is the longest part (\~3 hours). Merging logic and parallelization are straightforward once the loop is extracted.

**Q: Can I test this before committing?** Yes. After refactoring: - Set `cores = 1L`, run vignette, compare to current code → should be identical - Set `cores = 2L` on your machine, run vignette, compare → should be identical (same seed) - If identical, the refactor is correct; commit and move on

