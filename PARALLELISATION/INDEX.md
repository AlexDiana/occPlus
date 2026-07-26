---

editor_options: 
  markdown: 
    wrap: 72
---

# OpenMP & Parallelization Documentation Index

**Created**: July 26, 2026\
**Purpose**: Complete reference set for resolving occJSDM's OpenMP violations and implementing chain-level parallelism\
**Scope**: Resolves CRAN blocking items #16, A.6, #17a, #21 and improves performance via multicore chains

------------------------------------------------------------------------

## 📚 Four-Part Documentation Set

### 1. [**OPENMP_PARALLELIZATION_GUIDE.md**](OPENMP_PARALLELIZATION_GUIDE.md) ⭐ START HERE

- **Length**: 481 lines, 18 KB
- **Audience**: Anyone needing to understand the full problem and solution
- **Contains**:
  - Executive summary of CRAN violations
  - Current OpenMP usage (2 parallel regions, locations, RNG calls)
  - Why it fails CRAN (uncapped threads, RNG thread-safety, global state)
  - Two options: "fix OpenMP" vs. "remove & use chains"
  - Detailed recommendation: **Option B (remove OpenMP, use chain-level parallelism)**
  - Architecture diagrams and code flow
  - Full implementation roadmap (3 phases over 7–11 hours)
  - Backward compatibility guarantees
  - Q&A addressing common concerns and trade-offs
  - Detailed timeline and effort estimates
  - Complete rollout checklist
  - References to TODO.Rmd and CRAN submission plan items

**Read this first to understand**: What's wrong, why Option B is best, what happens next

------------------------------------------------------------------------

### 2. [**OPENMP_QUICK_REFERENCE.md**](OPENMP_QUICK_REFERENCE.md) 🚀 QUICK LOOKUP

- **Length**: 238 lines, 9.8 KB
- **Audience**: Implementer who needs line numbers, code locations, quick checklists
- **Contains**:
  - Status table: what's currently in place (parallel regions, build config, global state)
  - Three CRAN violations explained in 2–3 sentences each
  - Side-by-side comparison table: "Fix OpenMP" vs. "Remove OpenMP" (5 criteria)
  - 4-step implementation overview (not detailed, just structure)
  - Exact code locations to modify (file, line number, what to change)
  - Testing checklist (7 items with specific test names)
  - Common pitfalls table (what goes wrong, how to avoid)
  - Suggested commit messages for the refactoring (7 commits)
  - FAQ section (6 questions with concise answers)

**Use this for**: Finding where to make changes, quick reference during implementation

------------------------------------------------------------------------

### 3. [**PARALLEL_REFACTOR_EXAMPLES.md**](PARALLEL_REFACTOR_EXAMPLES.md) 🔧 CONCRETE CODE

- **Length**: \~20 KB with extensive code examples
- **Audience**: Implementer who needs working code templates to copy/adapt
- **Contains**:
  - **Part 1**: Extracting the chain loop
    - "Before/After" code comparison
    - Full `run_mcmc_chain()` function signature
    - Refactored `runOccJSDM()` with both serial and parallel paths
  - **Part 2**: Implementing per-chain output merging
    - `merge_running_means()` helper function
    - `merge_waic_across_chains()` with parallel variance formula
    - `merge_posterior_arrays()` for rebuilding 4-D arrays
    - Master `merge_chain_outputs()` function
  - **Part 3**: Parallel execution
    - Option A: Serial path (backward-compatible loop)
    - Option B: Parallel path (with `makeCluster()`, `parLapply()`, `clusterSetRNGStream()`)
  - **Part 4**: Testing
    - Minimal dataset setup for unit tests
    - Concrete `test_that()` blocks for serial/parallel equivalence
    - Concrete `test_that()` blocks for merging correctness
  - **Part 5**: Debugging & verification
    - Symptom → cause → fix table (8 common issues)
    - Verification script to run after refactoring
    - File structure after refactoring

**Use this for**: Copying working code directly into your implementation

------------------------------------------------------------------------

### 4. [**DOCUMENTATION_SUMMARY.txt**](DOCUMENTATION_SUMMARY.txt) 📋 EXECUTIVE OVERVIEW

- **Length**: Plain-text, \~150 lines
- **Audience**: Anyone (managers, reviewers) needing a 5-minute summary
- **Contains**:
  - What the 4 documents are and what each covers
  - Key findings (current state, why it fails, recommended solution)
  - CRAN items resolved (with status)
  - Implementation outline (3 phases, effort per phase)
  - Quick start for implementer (8 steps, effort for each)
  - Key code locations in one place
  - References to TODO.Rmd and CRAN plan
  - Backward compatibility statement
  - Plain-text formatting for easy reading

**Use this for**: Orientation, status reports, getting buy-in

------------------------------------------------------------------------

## 🎯 Quick Navigation by Use Case

### **"I need to understand why OpenMP is a problem"**

→ [OPENMP_PARALLELIZATION_GUIDE.md](OPENMP_PARALLELIZATION_GUIDE.md), sections: - "Executive Summary" (top) - "Current OpenMP Usage" - "Why Current Approach Fails CRAN"

### **"I need to find where to make changes"**

→ [OPENMP_QUICK_REFERENCE.md](OPENMP_QUICK_REFERENCE.md), sections: - "Code Locations (To Remove)" table - "Build Config (To Simplify)" - "Global State (To Remove)" - "Chain Loop (To Extract)"

### **"I need working code to start implementing"**

→ [PARALLEL_REFACTOR_EXAMPLES.md](PARALLEL_REFACTOR_EXAMPLES.md), sections: - Part 1: Extract Chain Loop → Copy `run_mcmc_chain()` signature and body - Part 2: Implementing Merging → Copy all 4 merge functions - Part 3: Implement Parallel Execution → Copy serial and parallel paths

### **"I need to test what I've built"**

→ [PARALLEL_REFACTOR_EXAMPLES.md](PARALLEL_REFACTOR_EXAMPLES.md), sections: - Part 4: Testing code examples - Part 5: Verification script

### **"I need a 5-minute summary for a meeting"**

→ [DOCUMENTATION_SUMMARY.txt](DOCUMENTATION_SUMMARY.txt): - "KEY FINDINGS" section - "CRAN ITEMS RESOLVED" section

### **"I'm implementing and something went wrong"**

→ [PARALLEL_REFACTOR_EXAMPLES.md](PARALLEL_REFACTOR_EXAMPLES.md), Part 5: - "Debugging Checklist" table - "Verification Script" to confirm correctness

------------------------------------------------------------------------

## 📍 Key Locations in Source Code

| Component | File | Lines | Action |
|-----------------------|-----------------|-----------------|-----------------|
| **Parallel region 1** | `src/jsdm.cpp` | 398–410 | Delete `#pragma omp parallel for` and commented loop (404–407) |
| **Parallel region 2** | `src/functions.cpp` | 607 | Delete `#pragma omp parallel for` |
| **Build config** | `src/Makevars` | 1–2 | Remove `SHLIB_OPENMP_CXXFLAGS` from both lines |
| **Global state** | `R/zzz.R` | 1–7 | Delete entire file |
| **Chain loop** | `R/runOccJSDM.R` | \~895–1200 | Extract to new `run_mcmc_chain()` function |

------------------------------------------------------------------------

## ✅ What Gets Resolved

| Issue | CRAN Item | Status | Solution |
|-----------------|---------------------|-----------------|-------------------|
| Uncapped thread count | #16 | BLOCKING | Remove OpenMP (each worker process has 1 thread) |
| RNG thread-safety | A.6 | BLOCKING | Remove OpenMP (separate processes → separate RNG states) |
| Global state pollution | #17a | BLOCKING | Delete `R/zzz.R` (remove `options(mc.cores = ...)` call) |
| Undeclared `parallel` import | #21 | BLOCKING | Remove `.onLoad()` (only place that uses `parallel::`) |

**Result**: All 4 blocking issues resolved in a single refactoring project.

------------------------------------------------------------------------

## 📊 Implementation Overview

| Phase | Duration | Key Output | Validates |
|-----------------|-----------------|-------------------|-------------------|
| **Phase 1: Architecture** | 4–5 hrs | `run_mcmc_chain()`, merge functions, serial/parallel paths | Code compiles, both paths return identical results |
| **Phase 2: Cleanup** | 30 min | Remove OpenMP flags, delete pragmas, delete `R/zzz.R` | No OpenMP warnings in `devtools::check()` |
| **Phase 3: Testing** | 2–3 hrs | Unit tests, vignette verification, CRAN compliance | `devtools::check --as-cran` passes; vignette output identical |
| **TOTAL** | **7–11 hrs** | CRAN-ready code | All tests pass; no violations remain |

------------------------------------------------------------------------

## 🔗 References & Related Items

**In TODO.Rmd**: - **Item A.6** (line \~200): RNG thread-safety bug [FIXED by removing OpenMP] - **Item D.1** (line \~260): Parallelize over chains using PSOCK [MAIN SOLUTION] - **Item D.2** (line \~303): Remove .onLoad() global state [CLEANUP] - **Item #16** (CRAN plan): OpenMP uncapped threads [BLOCKED, FIXED BY REMOVING OPENMP] - **Item #17a** (CRAN plan): Global state pollution [BLOCKED, FIXED BY DELETING R/zzz.R] - **Item #21** (CRAN plan): Undeclared `parallel` import [BLOCKED, FIXED BY REMOVING .onLoad()]

**In CRAN submission plan**: - Blocking items: #16 (OpenMP), #17a (global state), #21 (import), A.6 (RNG) - All resolved by this refactoring

------------------------------------------------------------------------

## 💡 Why This Approach?

**Two ways to fix OpenMP**:

| Aspect | Option A: Fix OpenMP | Option B: Remove OpenMP ← RECOMMENDED |
|---------------|-------------------|-------------------------------------|
| **Complexity** | High (guards, thread-local RNG) | Medium (refactor loop structure) |
| **Risk** | Medium (threading is subtle) | Low (separate processes are simple) |
| **Performance** | Serial within iterations + small OpenMP speedup | Chain-level parallelism (linear with \# chains) |
| **Maintainability** | Lower (OpenMP quirks per compiler) | Higher (standard R parallelism) |
| **CRAN safety** | Medium (still risky on some platforms) | High (by construction) |
| **Portability** | Platform-specific (OpenMP maturity varies) | Universal (PSOCK works everywhere) |

**Option B wins on every practical criterion** for occJSDM's use case.

------------------------------------------------------------------------

## 🎓 Learning Path

**For implementation team**: 1. (5 min) Read **DOCUMENTATION_SUMMARY.txt** for orientation 2. (20 min) Skim **OPENMP_PARALLELIZATION_GUIDE.md** sections on "Why Current Approach Fails" and "Why Option B is Better" 3. (10 min) Review **OPENMP_QUICK_REFERENCE.md** code locations 4. (30 min) Read **PARALLEL_REFACTOR_EXAMPLES.md** Parts 1–2 to understand the refactoring 5. Start implementation using code templates from **PARALLEL_REFACTOR_EXAMPLES.md** Parts 3–4 6. Use **OPENMP_QUICK_REFERENCE.md** as reference while coding 7. Use **PARALLEL_REFACTOR_EXAMPLES.md** Part 5 to test and debug

**Total learning time**: \~1 hour before first code change.

------------------------------------------------------------------------

## ✍️ Document Metadata

| Document | Size | Format | Completeness | Target Audience |
|--------------|--------------|--------------|---------------|-----------------|
| OPENMP_PARALLELIZATION_GUIDE.md | 18 KB | Markdown | Comprehensive (100%) | Architects, decision-makers |
| OPENMP_QUICK_REFERENCE.md | 9.8 KB | Markdown | Practical (100%) | Implementers, reviewers |
| PARALLEL_REFACTOR_EXAMPLES.md | 19 KB | Markdown + code | Implementation (100%) | Developers |
| DOCUMENTATION_SUMMARY.txt | 9.4 KB | Plain text | Executive summary (100%) | Managers, status reports |

**Total documentation**: \~56 KB, \~1,400 lines, 100% coverage of problem and solution

------------------------------------------------------------------------

## 🚀 Getting Started Right Now

**Next steps** (in order): 1. ✅ You've read this INDEX.md 2. Read [OPENMP_PARALLELIZATION_GUIDE.md](OPENMP_PARALLELIZATION_GUIDE.md) sections: - "Executive Summary" - "Why Current Approach Fails CRAN" 3. Decide: Does Option B (remove OpenMP, use chains) align with project goals? 4. If yes: Read [PARALLEL_REFACTOR_EXAMPLES.md](PARALLEL_REFACTOR_EXAMPLES.md) Parts 1–2 5. If yes: Estimate effort and timeline with your team 6. Begin Phase 1 (architecture change) using code templates from PARALLEL_REFACTOR_EXAMPLES.md

**Questions?** Refer to: - How to fix: [PARALLEL_REFACTOR_EXAMPLES.md](PARALLEL_REFACTOR_EXAMPLES.md) - Why this way: [OPENMP_PARALLELIZATION_GUIDE.md](OPENMP_PARALLELIZATION_GUIDE.md) - Where to change: [OPENMP_QUICK_REFERENCE.md](OPENMP_QUICK_REFERENCE.md) - Summary: [DOCUMENTATION_SUMMARY.txt](DOCUMENTATION_SUMMARY.txt)

------------------------------------------------------------------------

**Status**: Complete reference set ready for implementation\
**Confidence level**: High (based on TODO.Rmd items D.1–D.2, which describe this exact approach)\
**Next phase**: Team reads documentation, estimates effort, schedules implementation
