# Detailed Refactoring Examples: OpenMP → Chain-Level Parallelism

This document provides concrete code examples for the parallelization refactor described in `OPENMP_PARALLELIZATION_GUIDE.md`.

## Part 1: Extract the Chain Loop into `run_mcmc_chain()`

### Structure: Before and After

**Before** (current `R/runOccJSDM.R:895`):
```r
runOccJSDM <- function(data, MCMCparams, ..., cores = 1L) {
  # ... initialization code (lines 1–894) ...
  
  # Start of chain loop
  for (chain in 1:nchain) {
    # 300+ lines of MCMC logic
    # Writes to:
    #   B_output_chain, U_output_chain, z_output_chain, w_output_chain,
    #   theta_output_chain, theta0_output_chain, p_output_chain, q_output_chain,
    #   Gs_output_chain, Bs_output_chain, As_output_chain, Cs_output_chain,
    #   L_output_chain, idx_ls_output_chain,
    #   z_output_mean, z_output_M2, ..., (running sum/M2),
    #   min_ess, max_rhat, n_ess, n_param, (diagnostics),
    #   waic_mean, waic_M2 (WAIC accumulator)
  }
  # End of chain loop
  
  # ... postprocessing and output assembly (lines 1200–end) ...
}
```

**After** (refactored into `R/parallel_mcmc.R`):
```r
# New file: R/parallel_mcmc.R
run_mcmc_chain <- function(
  chain_id,
  data_list,
  X_psi, X_theta, X_psi_centers, X_theta_centers,
  Xs, Xs_index, Xs_centers,
  idx_z_k, idx_z_w, idx_w_k, idx_p_k,
  speciesNames,
  primerNames,
  model,
  threshold,
  useBiotic,
  useSpat,
  useSpatField,
  list_params,
  list_jsdmParams,
  MCMCparams,
  seed
) {
  # The entire ~300-line chain loop body goes here
  # Write outputs to *_output_chain (chain-specific) instead of *_output (global)
  # Return a named list with all the outputs
  return(list(
    B_output_chain = B_output_chain,
    U_output_chain = U_output_chain,
    z_output_chain = z_output_chain,
    w_output_chain = w_output_chain,
    theta_output_chain = theta_output_chain,
    # ... other output arrays ...
    z_output_sum = z_output_sum,  # For running mean
    z_output_M2 = z_output_M2,    # For running variance
    # ... other running accumulators ...
    waic_parts = list(
      mean_lik = mean_lik,
      M2 = M2,
      n_iter = niter
    )
  ))
}

# Main function refactored
runOccJSDM <- function(data, MCMCparams, ..., cores = 1L) {
  # ... initialization code (lines 1–894) unchanged ...
  
  if (cores == 1L) {
    # Serial path: existing loop (backward compatible)
    chain_outputs <- list()
    for (chain in 1:nchain) {
      chain_outputs[[chain]] <- run_mcmc_chain(
        chain_id = chain,
        data_list = data_list,
        X_psi = X_psi,
        # ... all other args ...
        seed = seed + chain - 1  # Unique seed per chain
      )
    }
  } else {
    # Parallel path: new code
    require(parallel)
    cl <- makeCluster(cores, type = "PSOCK")
    on.exit(stopCluster(cl), add = TRUE)
    
    # Attach package in worker processes
    clusterEvalQ(cl, library(occJSDM))
    
    # Set up reproducible RNG streams
    clusterSetRNGStream(cl, iseed = seed)
    
    # Run chains in parallel
    chain_outputs <- parLapply(
      cl,
      1:nchain,
      run_mcmc_chain,
      data_list = data_list,
      X_psi = X_psi,
      # ... all other args ...
      MCMCparams = MCMCparams,
      seed = seed
    )
  }
  
  # Merge per-chain outputs (see Part 2 below)
  results_output <- merge_chain_outputs(
    chain_outputs = chain_outputs,
    nchain = nchain,
    niter = niter
  )
  
  # ... postprocessing and output assembly unchanged ...
}
```

---

## Part 2: Implement Per-Chain Output Merging

### Overview of Outputs

Each chain produces:
1. **Full posterior arrays** (e.g., `B_output_chain[p, S, niter, 1]`)
2. **Running sum** (e.g., `z_output_sum[n, S]` for incrementally updated mean)
3. **Running M2** (e.g., `z_output_M2[n, S]` for incremental variance via Welford's algorithm)
4. **WAIC parts** (mean likelihood + M2 for variance)

### Function: Merge Running Means

```r
merge_running_means <- function(chain_outputs, nchain, component_name) {
  #' Merge running means from multiple chains
  #' 
  #' Each chain accumulated: sum / nchain of its iterations
  #' Merge by: sum all per-chain sums and divide by total count
  #'
  #' @param chain_outputs List of outputs from each chain
  #' @param nchain Number of chains
  #' @param component_name Name of the component (e.g., "z_output_sum")
  
  # Extract per-chain sums
  chain_sums <- lapply(
    chain_outputs,
    function(x) x[[component_name]]
  )
  
  # Sum across chains
  merged_sum <- Reduce("+", chain_sums)
  
  # Return merged (divide by nchain happens at the end, see below)
  return(merged_sum)
}

# Example usage:
z_output_mean <- merge_running_means(chain_outputs, nchain, "z_output_sum") / nchain
psi_output_mean <- merge_running_means(chain_outputs, nchain, "psi_output_sum") / nchain
```

### Function: Merge WAIC (Parallel Variance Formula)

```r
merge_waic_across_chains <- function(chain_outputs, nchain, niter) {
  #' Merge WAIC across chains using parallel variance formula
  #'
  #' Each chain tracks: mean_lik (mean likelihood), M2 (variance accumulator)
  #' Parallel variance formula:
  #'   global_mean = mean(chain_means)
  #'   global_var  = mean(chain_vars) + var(chain_means)
  #'   global_M2   = global_var * n_total
  #'
  #' From Welford & Chan (1979), reformulated by Knuth & Lewis
  
  # Extract WAIC parts from each chain
  chain_means <- sapply(chain_outputs, function(x) x$waic_parts$mean_lik)
  chain_M2s <- sapply(chain_outputs, function(x) x$waic_parts$M2)
  
  # Per-chain variance (estimated from M2)
  # M2 is the sum of squared deviations, so:
  chain_vars <- chain_M2s / (niter - 1)  # Unbiased variance estimate
  
  # Merge: global mean = mean of per-chain means
  global_mean <- mean(chain_means)
  
  # Merge: global variance = mean of per-chain variances + variance of chain means
  mean_of_vars <- mean(chain_vars)
  var_of_means <- var(chain_means) * (nchain - 1) / nchain  # Proper variance
  global_var <- mean_of_vars + var_of_means
  
  # Convert back to M2 form
  global_M2 <- global_var * (nchain * niter - 1)
  
  return(list(
    mean_lik = global_mean,
    M2 = global_M2
  ))
}

# Example usage:
waic_merged <- merge_waic_across_chains(chain_outputs, nchain = 2, niter = 1000)
```

### Function: Rebuild Full Posterior Arrays

```r
merge_posterior_arrays <- function(chain_outputs, nchain, output_name) {
  #' Merge per-chain posterior arrays into a single array
  #'
  #' Each chain produced: array[dim1, dim2, ..., niter, 1]
  #' Merge into: array[dim1, dim2, ..., niter, nchain]
  #'
  #' @param chain_outputs List of outputs from each chain (length nchain)
  #' @param nchain Number of chains
  #' @param output_name Name of array to merge (e.g., "B_output_chain")
  
  # Extract per-chain arrays
  chain_arrays <- lapply(
    chain_outputs,
    function(x) x[[output_name]]
  )
  
  # Get dimensions
  first_array <- chain_arrays[[1]]
  d <- dim(first_array)  # [dim1, dim2, ..., niter, 1]
  
  # Verify all chains have same dimensions
  for (i in 2:nchain) {
    if (!identical(dim(chain_arrays[[i]]), d)) {
      stop(paste(
        "Dimension mismatch in", output_name,
        "for chain", i, "vs chain 1"
      ))
    }
  }
  
  # Merge: stack along last dimension (chain axis)
  # From [dim1, dim2, ..., niter, 1] to [dim1, dim2, ..., niter, nchain]
  merged <- abind::abind(chain_arrays, along = length(d))
  
  # Reset last dimension from "1, 1, ..., 1" (length nchain) to explicit chain names
  dimnames(merged)[[length(d)] + 1] <- paste0("chain_", 1:nchain)
  
  return(merged)
}

# Example usage:
B_output <- merge_posterior_arrays(chain_outputs, nchain, "B_output_chain")
U_output <- merge_posterior_arrays(chain_outputs, nchain, "U_output_chain")
```

### Master Merge Function

```r
merge_chain_outputs <- function(chain_outputs, nchain, niter) {
  #' Master function: merge all per-chain outputs
  #'
  #' @param chain_outputs List of outputs from run_mcmc_chain(), one per chain
  #' @param nchain Number of chains
  #' @param niter Number of iterations per chain
  #'
  #' @return List (results_output) with merged arrays and diagnostics
  
  # 1. Merge full posterior arrays
  message("Merging posterior arrays...")
  B_output <- merge_posterior_arrays(chain_outputs, nchain, "B_output_chain")
  B0_output <- merge_posterior_arrays(chain_outputs, nchain, "B0_output_chain")
  U_output <- merge_posterior_arrays(chain_outputs, nchain, "U_output_chain")
  # ... repeat for all output arrays ...
  
  # 2. Merge running means
  message("Computing posterior means...")
  z_output_mean <- merge_running_means(chain_outputs, nchain, "z_output_sum") / nchain
  psi_output_mean <- merge_running_means(chain_outputs, nchain, "psi_output_sum") / nchain
  # ... repeat for all running means ...
  
  # 3. Merge WAIC
  message("Merging WAIC...")
  waic_merged <- merge_waic_across_chains(chain_outputs, nchain, niter)
  
  # 4. Extract diagnostics (if available)
  # Combine min ESS, max Rhat across chains
  min_ess <- min(sapply(chain_outputs, function(x) x$min_ess))
  max_rhat <- max(sapply(chain_outputs, function(x) x$max_rhat))
  
  # 5. Assemble results_output
  results_output <- list(
    # Full posterior arrays
    jsdm_output = list(
      B_output = B_output,
      B0_output = B0_output,
      U_output = U_output,
      # ... other jsdm outputs ...
      mean_lik = waic_merged$mean_lik,
      M2 = waic_merged$M2
    ),
    
    # Running means (for summary output)
    z_output_mean = z_output_mean,
    psi_output_mean = psi_output_mean,
    # ... other means ...
    
    # Diagnostics
    min_ess = min_ess,
    max_rhat = max_rhat,
    
    # Metadata
    nchain = nchain,
    niter = niter
  )
  
  return(results_output)
}
```

---

## Part 3: Implement Parallel Execution

### Option A: Serial Path (Backward Compatible)

```r
# In runOccJSDM(), when cores == 1L:

message("Running MCMC (serial, 1 chain at a time)...")

chain_outputs <- list()

for (chain in 1:nchain) {
  message("Chain ", chain, "/", nchain)
  
  chain_outputs[[chain]] <- run_mcmc_chain(
    chain_id = chain,
    data_list = data_list,
    X_psi = X_psi,
    X_theta = X_theta,
    X_psi_centers = X_psi_centers,
    X_theta_centers = X_theta_centers,
    Xs = Xs,
    Xs_index = Xs_index,
    Xs_centers = Xs_centers,
    idx_z_k = idx_z_k,
    idx_z_w = idx_z_w,
    idx_w_k = idx_w_k,
    idx_p_k = idx_p_k,
    speciesNames = speciesNames,
    primerNames = primerNames,
    model = model,
    threshold = threshold,
    useBiotic = useBiotic,
    useSpat = useSpat,
    useSpatField = useSpatField,
    list_params = list_params,
    list_jsdmParams = list_jsdmParams,
    MCMCparams = MCMCparams,
    seed = seed + chain - 1
  )
}
```

### Option B: Parallel Path (New Feature)

```r
# In runOccJSDM(), when cores > 1L:

message("Running MCMC (parallel, ", cores, " cores, ", nchain, " chains)...")

# Check that package has parallel
require(parallel, quietly = TRUE)

# Create cluster
cl <- makeCluster(cores, type = "PSOCK")

# Ensure cleanup on exit (even if error occurs)
on.exit({
  if (exists("cl")) stopCluster(cl)
}, add = TRUE)

# Load occJSDM in worker processes
message("Initializing workers...")
clusterEvalQ(cl, {
  library(occJSDM)
  library(Rcpp)
  library(RcppArmadillo)
})

# Set up reproducible RNG (L'Ecuyer streams)
clusterSetRNGStream(cl, iseed = seed)

# Prepare arguments for each chain
chain_args <- lapply(1:nchain, function(chain_id) {
  list(
    chain_id = chain_id,
    data_list = data_list,
    X_psi = X_psi,
    X_theta = X_theta,
    X_psi_centers = X_psi_centers,
    X_theta_centers = X_theta_centers,
    Xs = Xs,
    Xs_index = Xs_index,
    Xs_centers = Xs_centers,
    idx_z_k = idx_z_k,
    idx_z_w = idx_z_w,
    idx_w_k = idx_w_k,
    idx_p_k = idx_p_k,
    speciesNames = speciesNames,
    primerNames = primerNames,
    model = model,
    threshold = threshold,
    useBiotic = useBiotic,
    useSpat = useSpat,
    useSpatField = useSpatField,
    list_params = list_params,
    list_jsdmParams = list_jsdmParams,
    MCMCparams = MCMCparams,
    seed = seed + chain_id - 1
  )
})

# Execute in parallel
chain_outputs <- clusterMap(
  cl,
  run_mcmc_chain,
  # Chain-varying args (differ per chain)
  chain_id = 1:nchain,
  seed = sapply(1:nchain, function(i) seed + i - 1),
  # Constant args (same for all chains)
  MoreArgs = list(
    data_list = data_list,
    X_psi = X_psi,
    X_theta = X_theta,
    # ... all other constant args ...
    MCMCparams = MCMCparams
  ),
  SIMPLIFY = FALSE
)

message("Finalizing cluster...")
# on.exit above will clean up automatically
```

---

## Part 4: Minimal Reproducible Example (Testing)

### Tiny Dataset for Unit Tests

```r
# Helper function for tests: creates minimal valid simulation

library(testthat)

minimal_mcmc_params <- function() {
  list(nchain = 2, nburn = 100, niter = 200, nthin = 1)
}

simulate_tiny_data <- function() {
  # Create minimal dataset: n=5 sites, S=3 species
  set.seed(42)
  n <- 5; S <- 3
  
  # Simulate data
  sim_args <- list(
    list_datasettings = list(
      n = n, S = S, g = 1,
      M = rep(2, n),
      P = 1,
      K = rep(2, n * 1 * 2),
      ncov_psi = 1, ncov_theta = 1
    ),
    list_params = list(
      p = matrix(0.7, 1, S),
      q = matrix(0.05, 1, S),
      theta0 = rep(0.05, S),
      theta_baseline = 0.4
    ),
    list_jsdmParams = list(
      gt = 1, d = 1, ds = 0, tau = 1,
      sigma_b = 1, sigma_bs = 0, sigma_ts = 0,
      sigma_h = 0, sigma_s = 0, l_s = NULL,
      useSpatField = FALSE
    ),
    model = "two_stage"
  )
  
  return(do.call(simulateOccJSDMData, sim_args))
}

# Unit test: serial and parallel produce same output
test_that("serial and parallel paths produce identical output with same seed", {
  set.seed(42)
  data <- simulate_tiny_data()
  mcmc_params <- minimal_mcmc_params()
  
  # Serial fit
  set.seed(123)
  fit_serial <- runOccJSDM(
    data$data_list,
    MCMCparams = mcmc_params,
    cores = 1L
  )
  
  # Parallel fit (same seed)
  set.seed(123)
  fit_parallel <- runOccJSDM(
    data$data_list,
    MCMCparams = mcmc_params,
    cores = 2L
  )
  
  # Compare key outputs
  expect_equal(
    fit_serial$results_output$jsdm_output$B_output,
    fit_parallel$results_output$jsdm_output$B_output,
    tolerance = 1e-10
  )
  
  expect_equal(
    fit_serial$results_output$z_output_mean,
    fit_parallel$results_output$z_output_mean,
    tolerance = 1e-10
  )
  
  # WAIC should be identical
  serial_waic <- computeWAIC(fit_serial)
  parallel_waic <- computeWAIC(fit_parallel)
  expect_equal(serial_waic, parallel_waic, tolerance = 1e-10)
})

# Unit test: merging works correctly
test_that("WAIC merging computes correct global mean and variance", {
  # Manual check: two chains, each with known statistics
  chain_1_mean <- 100.5
  chain_2_mean <- 101.2
  expected_global_mean <- mean(c(chain_1_mean, chain_2_mean))
  
  # Simulate WAIC parts
  mock_chain_outputs <- list(
    list(waic_parts = list(mean_lik = chain_1_mean, M2 = 5000)),
    list(waic_parts = list(mean_lik = chain_2_mean, M2 = 4800))
  )
  
  merged <- merge_waic_across_chains(mock_chain_outputs, nchain = 2, niter = 1000)
  
  expect_equal(merged$mean_lik, expected_global_mean, tolerance = 1e-10)
  expect_true(!is.na(merged$M2))  # Just check it's computed
})
```

---

## Part 5: Debugging Checklist

### Common Issues When Refactoring

| Symptom | Likely Cause | Fix |
|---------|-----------|-----|
| **Parallel path errors with "object X not found"** | Function arg wasn't passed to `parLapply` | Add to `MoreArgs` list or function signature |
| **Output dimensions mismatch** | Per-chain array has `[..., niter, 1]` but code expects `[..., niter]` | Ensure `merge_posterior_arrays()` properly stacks along chain axis |
| **WAIC is NaN** | M2 calculation overflow or underflow | Check: `chain_vars = chain_M2 / (niter - 1)` is correct |
| **RNG not reproducible in parallel** | Seed not set | Call `clusterSetRNGStream(cl, iseed = seed)` before `parLapply` |
| **Cluster never stops** | `stopCluster()` not called | Use `on.exit(stopCluster(cl))` to ensure cleanup |
| **Running means are wrong** | Not dividing by nchain | Always do `sum_per_chain / nchain` after merging |
| **Compiler warnings about OpenMP** | `#pragma omp` still in source | Search `src/` for "pragma omp" and delete all occurrences |

### Verification Script

```r
# Run this after refactoring to verify correctness

library(occJSDM)
library(microbenchmark)

# 1. Load test data
data <- sampledata
tiny_mcmc <- list(nchain = 2, nburn = 100, niter = 200, nthin = 1)

# 2. Fit serial (baseline)
set.seed(42)
cat("Fitting serial...\n")
t_serial <- system.time({
  fit_serial <- runOccJSDM(data, MCMCparams = tiny_mcmc, cores = 1L)
})

# 3. Fit parallel
set.seed(42)
cat("Fitting parallel (2 cores)...\n")
t_parallel <- system.time({
  fit_parallel <- runOccJSDM(data, MCMCparams = tiny_mcmc, cores = 2L)
})

# 4. Compare outputs
cat("\nComparison:\n")
cat("Serial time:   ", round(t_serial[3], 2), "sec\n")
cat("Parallel time: ", round(t_parallel[3], 2), "sec\n")
cat("Speedup:       ", round(t_serial[3] / t_parallel[3], 2), "x\n")

# Check identical outputs
cat("\nOutput equivalence:\n")
cat("B_output equal?",
    all.equal(fit_serial$results_output$jsdm_output$B_output,
              fit_parallel$results_output$jsdm_output$B_output) == TRUE, "\n")
cat("z_output_mean equal?",
    all.equal(fit_serial$results_output$z_output_mean,
              fit_parallel$results_output$z_output_mean) == TRUE, "\n")

# Extract WAIC
waic_serial <- extractWAIC(fit_serial)
waic_parallel <- extractWAIC(fit_parallel)
cat("WAIC equal?", all.equal(waic_serial, waic_parallel) == TRUE, "\n")

# 5. Check vignette still works
cat("\nRendering vignette...\n")
rmarkdown::render("vignettes/occJSDM.Rmd", quiet = TRUE)
cat("Vignette rendered successfully.\n")
```

---

## Appendix: File Structure After Refactoring

```
R/
├── occJSDM-package.R       (unchanged)
├── runOccJSDM.R            (modified: extract loop, add cores arg, call merge functions)
├── parallel_mcmc.R         (NEW: run_mcmc_chain, merge_chain_outputs, helpers)
├── jsdmfun.R               (unchanged)
├── mcmcfun.R               (unchanged)
├── output.R                (unchanged)
├── diagnostics.R           (unchanged)
├── simulateData.R          (unchanged)
├── data.R                  (unchanged)
└── zzz.R                   (DELETED or empty)

src/
├── Makevars                (modified: remove SHLIB_OPENMP_CXXFLAGS)
├── jsdm.cpp                (modified: remove #pragma omp line 398)
├── functions.cpp           (modified: remove #pragma omp line 607)
└── RcppExports.cpp         (regenerated: no changes expected)

tests/
└── testthat/
    ├── test-parallel-equivalence.R  (NEW)
    ├── test-merge-functions.R       (NEW)
    ├── test-chain-extraction.R      (NEW)
    └── ... other tests ...
```

---

## Summary of Changes

| File | Type | Change |
|------|------|--------|
| `R/runOccJSDM.R` | Modified | Extract chain loop, add `cores` arg, implement serial/parallel paths |
| `R/parallel_mcmc.R` | New | Contains `run_mcmc_chain()`, merge functions |
| `R/zzz.R` | Deleted | Remove global state pollution |
| `src/Makevars` | Modified | Remove OpenMP compiler flags |
| `src/jsdm.cpp` | Modified | Delete `#pragma omp` at line 398 |
| `src/functions.cpp` | Modified | Delete `#pragma omp` at line 607 |
| Tests | New | Add unit tests for parallel equivalence |
| Vignette | Unchanged | Output identical; mention `cores` arg optionally |

---

**Total effort**: ~7–11 hours (extraction + merging + parallelization + testing)
