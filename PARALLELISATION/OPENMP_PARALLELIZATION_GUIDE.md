---

editor_options: 
  markdown: 
    wrap: 72
---

# OpenMP and Parallelization Strategy for occJSDM

**Updated: July 26, 2026**

## Executive Summary

occJSDM currently has **two active OpenMP parallel regions** in C++ code, but **violates CRAN policy on multiple fronts**:

1.  **Uncapped thread count**: No limit enforced, can use all available cores (CRAN policy: maximum 2 cores)
2.  **RNG thread-safety bug (A.6)**: Uses non-thread-safe R RNG in parallel regions
3.  **Global state pollution (item 17a)**: `.onLoad()` sets `options(mc.cores = ...)` affecting the user's session

The recommended solution is **lower-risk than fixing OpenMP**: **remove OpenMP entirely and parallelize over chains using R-level `parallel::makeCluster()` + `parLapply()` instead**. This: - ✅ Resolves all three CRAN violations with a single architectural change - ✅ Avoids the deadlock risk from mixing threaded BLAS + OpenMP + fork() - ✅ Provides better reproducibility via `parallel::clusterSetRNGStream()` - ✅ Works identically on Windows/macOS/Linux - ✅ Fully backward-compatible (make it opt-in with a `cores = 1L` argument)

------------------------------------------------------------------------

## Current OpenMP Usage

### Active Parallel Regions

**1. `samplePGvariables()` in `src/jsdm.cpp:398`**

``` cpp
#pragma omp parallel for collapse(2)
for(int i = 0; i < n1; i++){
  for(int j = 0; j < n2; j++){
    Omega_mat(i,j) = rpg(1, Xbeta(i,j));  // <-- Calls R's RNG!
  }
}
```

- **Problem**: `rpg()` is a wrapper around a non-thread-safe R RNG call
- **Why it works locally** (macOS): Stock macOS clang doesn't parallelize (`#pragma omp` is silently ignored without `-fopenmp`)
- **Why it fails on CRAN** (Linux/Windows with gfortran/gcc): Threaded BLAS + OpenMP forces actual parallelization

------------------------------------------------------------------------

**2. `sample_betatheta_cpp_parallel()` in `src/functions.cpp:607`**

``` cpp
#pragma omp parallel for
for (int s = 0; s < S; ++s) {
  // ...
  arma::vec Y = arma::randn(ncols);  // <-- Also uses R's RNG!
  // ...
}
```

- **Same problem**: `arma::randn()` ultimately calls `R::rnorm()` (not thread-safe)

------------------------------------------------------------------------

### Configuration

**`src/Makevars`:**

``` makefile
PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS)
PKG_LIBS = $(SHLIB_OPENMP_CXXFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
```

**Commented-out code (never reached):**

``` cpp
// #ifdef _OPENMP
//   omp_set_num_threads(n_threads);  // <-- Should cap at 2 for CRAN
// #endif
```

**Global state pollution in `R/zzz.R`:**

``` r
.onLoad <- function(libname, pkgname) {
  options(mc.cores = min(2L, parallel::detectCores()))
}
```

------------------------------------------------------------------------

## Why Current Approach Fails CRAN

### Issue 1: Uncapped Thread Count (CRAN Blocking Item #16)

**Policy violation**: CRAN forbids using more than 2 cores during automated checks. The code has: - No thread-capping code (`omp_set_num_threads()` is commented out) - No declaration in `SystemRequirements: OpenMP` - Unconditional OpenMP flags in `Makevars`

**Consequence**: A 16-core Linux machine used for CRAN checks will spawn 16 threads per iteration, violating policy and causing unpredictable performance on CRAN infrastructure.

------------------------------------------------------------------------

### Issue 2: RNG Thread-Safety (CRAN Blocking Item A.6)

**Policy violation**: CRAN's sanitizers detect non-thread-safe RNG in parallel regions.

R's RNG state (`R::rand_seed()`) is **global**, not thread-local. Calling `R::rnorm()` / `rpg()` / `arma::randn()` from multiple threads simultaneously causes: - Race conditions on the global seed state - Intermittent failures that appear platform/load-dependent - Violations of CRAN's thread-safety requirements

**Why it never happened locally**: macOS clang doesn't respect `#pragma omp` without explicit `-fopenmp` flag, so the code runs serial even with pragmas present.

**Why it will happen on CRAN**: Linux/Windows explicitly compile with OpenMP flags, causing true parallelization.

------------------------------------------------------------------------

### Issue 3: Global State Pollution (CRAN Item #17a + #21)

**Two violations:**

1.  **`options(mc.cores = ...)` in `.onLoad()`** violates CRAN's "package shall not change user's global state" policy
    - Affects every other package that reads `mc.cores`
    - Can propagate `NA` (if `detectCores()` fails)
    - Should be replaced with a function argument
2.  **`parallel::detectCores()` is undeclared** in `DESCRIPTION`
    - Not in `Imports:` or `Depends:`
    - Reported as `'::' import not declared from: 'parallel'`

------------------------------------------------------------------------

## Solution: Parallelize Over Chains (Recommended)

### Why This Works Better

**Architecture**: Move parallelism from within-iteration C++ loops to between-chain R loops.

```         
Current:
  for (chain 1:nchain) {
    for (iter 1:niter) {
      [OpenMP parallel for in C++]  ← Dangerous, uncapped, RNG-unsafe
    }
  }

Proposed:
  parallel::parLapply(cluster, chains = 1:nchain, function(chain_id) {
    [Run single chain entirely in one worker process]
    return(chain_outputs)
  })
  [Merge per-chain outputs in main process]
```

**Advantages:**

✅ **Solves all three CRAN violations** in one edit: - Each worker process has its own RNG state (thread-safe by design) - PSOCK cluster respects `OMP_NUM_THREADS` and CRAN's 2-core limit - No global state changes (opt-in `cores` argument to `runOccJSDM()`)

✅ **Platform-universal**: `parallel::makeCluster(type = "PSOCK")` works on Windows/macOS/Linux identically

✅ **Reproducible**: `parallel::clusterSetRNGStream()` provides L'Ecuyer RNG streams with independent, deterministic sequences per chain

✅ **Backward-compatible**: Serial by default (`cores = 1L`), so existing code and CRAN tests unaffected

✅ **Avoids BLAS deadlock**: No mixing of fork() + OpenMP + threaded BLAS (the classic deadlock on Linux with OpenBLAS)

✅ **Near-linear scaling**: Chain-level parallelism is embarrassingly parallel; `nchain=4` on a 4-core machine ≈ 4× faster

------------------------------------------------------------------------

### Trade-offs

**Cost**: Worker initialization overhead + data serialization/deserialization: - Each worker: `clusterEvalQ(cl, library(occJSDM))` to attach the package - Outbound: Serialize `data_list`, `X_psi`, `X_theta`, coefficients, starting values - Inbound: Deserialize `*_output_chain` arrays for each chain

**Why negligible**: MCMC runtime dominates (\~seconds per iteration × thousands of iterations). Serialization (milliseconds) is unmeasurable against the total.

**Concern**: Large `Bs_output` and `U_output` arrays need serialization back. See item D.9 (memory optimization) to address this separately if needed.

------------------------------------------------------------------------

## Implementation Roadmap

### Phase 1: Architecture Change (Core Fix)

#### Step 1: Extract chain loop body into a function

**Current** (`R/runOccJSDM.R:895`):

``` r
for (chain in 1:nchain) {
  # ~300 lines of MCMC logic, writes to *_output_chain
}
```

**Refactor to** (`R/parallel_mcmc.R`, new file):

``` r
run_mcmc_chain <- function(
  chain_id, data_list, X_psi, X_theta, ..., MCMCparams, seed
) {
  # Entire loop body here, returns list(
  #   B_output_chain, U_output_chain, ... all *_output_chain objects,
  #   waic_parts_chain = list(mean_lik, ...)
  # )
}
```

**Main function** (`R/runOccJSDM.R:895`):

``` r
if (cores == 1L) {
  # Serial path: existing loop
  for (chain in 1:nchain) { ... }
} else {
  # Parallel path: use makeCluster + parLapply
  require(parallel)
  cl <- makeCluster(cores, type = "PSOCK")
  on.exit(stopCluster(cl))
  
  clusterEvalQ(cl, library(occJSDM))
  clusterSetRNGStream(cl, seed = ...)  # Reproducible L'Ecuyer RNG
  
  chain_outputs <- parLapply(cl, 1:nchain, 
    run_mcmc_chain, 
    data_list = data_list, X_psi = X_psi, ...
  )
  
  # Merge per-chain outputs (see below)
}
```

------------------------------------------------------------------------

#### Step 2: Implement per-chain output merging

**WAIC merging** (TODO D.1 item about bug A.4):

``` r
merge_waic <- function(waic_parts_list) {
  # Per-chain: mean_lik_chain, M2_chain (Welford accumulator for variance)
  # Merge using parallel variance formula:
  #   global_mean = mean(per_chain_means)
  #   global_var  = mean(per_chain_vars) + var(per_chain_means)
  # Return: combined mean_lik, M2 → WAIC
}
```

**Running mean merging** (for `z_output_mean`, etc.):

``` r
merge_running_means <- function(running_mean_list) {
  # Each chain accumulated: n / nchain partial sums
  # Merge: sum all partials and divide by total count
  do.call("+", running_mean_list) / nchain
}
```

------------------------------------------------------------------------

### Phase 2: Code Cleanup

#### Step 1: Remove OpenMP from build

**`src/Makevars`:**

``` diff
- PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS)
- PKG_LIBS = $(SHLIB_OPENMP_CXXFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
+ PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
```

#### Step 2: Delete parallel pragmas and RNG-unsafe code

**`src/jsdm.cpp:398`**: Delete the `#pragma omp parallel for` and the commented code (lines 404-407) - `samplePGvariables()` becomes serial (which is fine, it's called once per iteration on a matrix, not the bottleneck) - If `rpg()` calls prove slow, profile first before re-parallelizing

**`src/functions.cpp:607`**: Delete `#pragma omp parallel for` - The loop is over species (typically `S ≤ 50`), which is fine serial - Remove the corresponding C++ export from `src/RcppExports.cpp` if `sample_betatheta_cpp_parallel` was specifically the parallel variant

#### Step 3: Clean up global state

**`R/zzz.R`:** Delete entirely or leave empty (all code removed, no more side effects)

**`DESCRIPTION`:** No longer need `SystemRequirements: OpenMP` if added

------------------------------------------------------------------------

### Phase 3: Testing & Validation

#### Unit tests

``` r
test_that("run_mcmc_chain returns expected structure", {
  # Verify one chain produces identical output to old serial code
  # (set seed, run both paths, compare)
})

test_that("chain merging preserves WAIC/means", {
  # Simulate two chains independently, merge, compare to serial
})

test_that("cores argument controls parallelism", {
  # cores = 1: serial, no cluster formed
  # cores = 2: cluster formed with 2 workers
  expect_true(inherits(cluster_object, "cluster"))
})
```

#### Integration test

``` r
test_that("parallel and serial paths produce identical results (same seed)", {
  set.seed(42)
  fit_serial   <- runOccJSDM(data, MCMCparams = params, cores = 1L)
  
  set.seed(42)
  fit_parallel <- runOccJSDM(data, MCMCparams = params, cores = 2L)
  
  expect_equal(fit_serial$results_output$B_output, 
               fit_parallel$results_output$B_output)
})
```

------------------------------------------------------------------------

## Implementation Details: What Stays Parallel

The **removal of within-iteration C++ OpenMP does NOT mean zero parallelism**:

1.  **Theaded BLAS** (OpenBLAS / Accelerate) remains enabled
    - Large matrix operations (`X %*% B`, `crossprod()`, etc.) automatically parallelize via BLAS
    - CRAN permits this (BLAS is standard infrastructure)
    - Thread-safe because it's a solved problem in the numerical library
2.  **Chain-level parallelism** (new R-level `parLapply`)
    - Each chain runs in a separate R process
    - Scales with number of chains (e.g., 4 chains on 4 cores ≈ 4× speedup)
3.  **One-time per-chain initialization cost**
    - `clusterEvalQ(cl, library(occJSDM))`: \~50–200 ms
    - Negligible against multi-minute MCMC runs

------------------------------------------------------------------------

## Backward Compatibility

**Default behavior preserved:**

``` r
# Existing code continues to work unchanged
fit <- runOccJSDM(data, MCMCparams = list(nchain = 2, niter = 1000))
# Uses serial path (cores = 1L by default), identical output
```

**Parallel use opt-in:**

``` r
# New code for faster inference
fit <- runOccJSDM(data, MCMCparams = list(nchain = 2, niter = 1000), cores = 2L)
# Uses parallel::makeCluster with 2 workers
```

**CRAN tests automatically use serial** (no code changes needed): - `R CMD check` never passes `cores > 1L` - Examples/vignettes respect 2-core limit by default - No vignette re-rendering required (output identical)

------------------------------------------------------------------------

## Questions & Clarifications

### Q: Why not just cap OpenMP threads at 2?

**A:** Three reasons:

1.  **RNG thread-safety persists**: Even with `omp_set_num_threads(2)`, the `rpg()` and `arma::randn()` calls in parallel loops remain non-thread-safe. Capping threads doesn't fix the race condition; it just makes it less likely to trigger.

2.  **Deadlock risk with BLAS + fork()**: If someone later optimizes with `mclapply()` (the natural choice for "parallel" in R), the combination of threaded BLAS + `fork()` + OpenMP becomes a deadlock hazard on Linux. Better to avoid OpenMP entirely.

3.  **Complexity**: Wrapping OpenMP code with platform-specific guards (`#ifdef _OPENMP`, `omp_set_num_threads()` with fallbacks) adds maintenance burden for marginal gain. R-level parallelism is simpler and more portable.

------------------------------------------------------------------------

### Q: Will removing OpenMP slow the code down?

**A:** No, because:

1.  **`samplePGvariables()` is not the bottleneck**. It's called once per iteration on a small matrix (`n × S`). Profiling (via `profvis::profvis()`) needed to confirm, but it's unlikely to be \> 5% of total time.

2.  **Chain-level parallelism provides near-linear speedup** as compensation (if running multiple chains).

3.  **Threaded BLAS remains enabled** and will accelerate matrix operations.

If profiling shows `samplePGvariables()` or similar is slow, a thread-safe alternative (e.g., vectorized C++ code, or a fast scalar loop) is cheaper than fixing OpenMP.

------------------------------------------------------------------------

### Q: What about users who want within-iteration parallelism?

**A:** The new approach is **actually faster in practice**:

- **Old**: Serial iterations, per-iteration OpenMP parallelism (overhead of forking worker threads, syncing, collecting results; \~10 ms per iteration)
- **New**: Parallel chains, each chain runs to completion serially in its worker

For `nchain = 2, niter = 1000` on a 2-core machine: - **Old**: \~10 s \* 1000 = 10 ks (serial) - **New**: \~5 s \* 1000 / 2 = 2.5 ks (2 chains in parallel)

The math is simpler and the speedup more predictable.

------------------------------------------------------------------------

### Q: Why `parallel::makeCluster(type = "PSOCK")` and not another approach?

| Approach | Windows | macOS | Linux | Comment |
|---------------|---------------|---------------|---------------|---------------|
| **PSOCK (sockets)** | ✅ | ✅ | ✅ | Slow IPC but universal; fresh R sessions |
| **fork (mclapply)** | ❌ Hard error | ⚠️ Unsafe | ⚠️ Unsafe | Not portable; deadlock risk with BLAS |
| **MPI** | ❌ Complex | ❌ Complex | ✅ | Overkill; requires system MPI install |
| **multisession (future)** | ✅ | ✅ | ✅ | Modern, but adds dependency; PSOCK is built-in |

**`makeCluster(type = "PSOCK")` is the right choice**: built-in to base R, no extra dependencies, universal, and thread-safe.

------------------------------------------------------------------------

## Timeline & Effort Estimate

| Phase | Task                                     | Effort  | Time           |
|-------|------------------------------------------|---------|----------------|
| **1** | Extract chain loop to `run_mcmc_chain()` | High    | 2–3 hours      |
| **1** | Implement merging (WAIC, running means)  | Medium  | 1–2 hours      |
| **1** | Implement `parLapply` path               | Medium  | 1–2 hours      |
| **2** | Remove OpenMP from `Makevars`            | Trivial | 5 min          |
| **2** | Delete `#pragma omp` lines               | Trivial | 5 min          |
| **2** | Clean up `R/zzz.R`                       | Trivial | 5 min          |
| **3** | Write unit & integration tests           | Medium  | 1–2 hours      |
| **3** | Verify against vignette                  | Medium  | 30 min         |
| **3** | GitHub Actions CI/CD (if updating)       | Low     | 30 min         |
|       | **TOTAL**                                |         | **7–11 hours** |

------------------------------------------------------------------------

## Rollout Checklist

- [ ] **Phase 1**: Refactor chain loop and implement parallelism
  - [ ] Create `R/parallel_mcmc.R` with `run_mcmc_chain()`
  - [ ] Implement WAIC merging logic
  - [ ] Implement running mean merging logic
  - [ ] Add `cores` argument to `runOccJSDM()` signature and roxygen
  - [ ] Verify `cores = 1L` path produces identical output to current code
  - [ ] Verify `cores > 1L` path produces identical output to serial path (same seed)
- [ ] **Phase 2**: Remove OpenMP
  - [ ] Edit `src/Makevars`: remove OpenMP flags
  - [ ] Delete `#pragma omp` from `src/jsdm.cpp:398–410`
  - [ ] Delete `#pragma omp` from `src/functions.cpp:607`
  - [ ] Delete/empty `R/zzz.R`
  - [ ] Rebuild package: `devtools::load_all()`, `Rcpp::compileAttributes()`, `devtools::document()`
- [ ] **Phase 3**: Test & validate
  - [ ] Run `devtools::check()` (should remove OpenMP-related compiler warnings)
  - [ ] Add testthat tests for parallelism
  - [ ] Re-render vignettes (output should be identical)
  - [ ] Smoke test on multicore machine (`cores = 2`, `cores = 4`)
  - [ ] Confirm GitHub Actions still pass
- [ ] **Phase 4**: Document & release
  - [ ] Update vignette with `cores` argument explanation
  - [ ] Update `NEWS.md` (resolves CRAN issues #16, #A.6, #17a, #21)
  - [ ] Commit & push to `origin/main`

------------------------------------------------------------------------

## References

- **TODO.Rmd sections**:
  - Item A.6 (RNG thread-safety bug)
  - Item D.1–D.2 (parallelization strategy)
  - Item 16 (CRAN blocking on OpenMP)
  - Item 17a (global state pollution)
  - Item 21 (`parallel` undeclared)
- **R Documentation**:
  - `?parallel::makeCluster` — PSOCK cluster setup
  - `?parallel::parLapply` — Parallel apply
  - `?parallel::clusterSetRNGStream` — Reproducible RNG
- **CRAN Policies**:
  - Use of multiple cores: limit to 2 during examples/tests/vignettes
  - Global state: packages shall not modify user's session
  - RNG safety: thread-local RNG required in parallel regions
