// Superseded thread-safe sampler variant, moved out of src/jsdm.cpp on
// 2 August 2026. TODO.md Fixed bugs 42.
//
// sampleB_SoR_TS() was written as the thread-safe variant of sampleB_SoR(),
// for the BBSL_Worker race recorded as Fixed bugs 41. That race was closed a
// different way -- sampleB_SoR() itself now draws through rnorm() in
// src/rng.h -- so this variant was left with no remaining purpose.
//
// It was unreachable when it was moved: no callers in src/, R/, tests/ or
// vignettes/, and not in NAMESPACE. It had been exported to R only as an
// unused RcppExports wrapper; that // [[Rcpp::export]] tag was removed and
// Rcpp::compileAttributes() re-run before the move.
//
// It is kept here as a record rather than deleted, because the reason it is
// wrong is worth being able to point at. It is NOT a reference implementation
// and must not be revived as written:
//
//   thread_local std::mt19937 gen(rd());
//
// seeds each thread's generator from std::random_device, which is OS entropy
// and is not derived from R's RNG. A sampler built on this would ignore
// set.seed() entirely and would be irreproducible at ANY thread count --
// strictly worse than the sampleB_SoR() it was meant to improve on, which is
// both thread-safe and seedable because src/rng.h derives its per-thread
// generators from a base seed drawn from R via setOccJSDMSeed().
//
// The _TS suffix is what made this dangerous rather than merely dead. It marks
// thread-safe variants throughout this codebase (mvrnormArmaQuick_TS,
// sample_beta_cpp_TS, sample_beta_nocov_cpp_TS), so it reads as the endorsed
// choice for anyone parallelising something. src/rng.h's header records the
// same mistake having already been made once, with mvrnormArmaQuick_TS.
//
// This directory is excluded from the build via .Rbuildignore, so nothing here
// is compiled. The code below will not build as it stands: it depends on
// XtOmegaX_SoR() and XtK_SoR() in src/jsdm.cpp and on the RcppArmadillo
// headers. Reproduced verbatim as it stood at the time of the move.

arma::vec sampleB_SoR_TS(arma::mat X, arma::mat &invB, arma::vec &b,
                      arma::vec &k, arma::vec Omega,
                      arma::mat &X_s_index,
                      arma::mat &Ks,
                      int X_centers) {

  // It is assumed that XtOmegaX_SoR and XtK_SoR are thread-safe
  // and do not modify shared state or call R's RNG.
  arma::mat XtOmegaX = XtOmegaX_SoR(X, X_centers, Omega, X_s_index, Ks);
  arma::mat tXk = XtK_SoR(X, X_s_index, Ks, k, X_centers);

  arma::mat Lambda_B = XtOmegaX + invB;
  arma::vec mu_B = tXk + invB * b;

  arma::mat L = arma::trans(arma::chol(Lambda_B));
  arma::vec tmp = arma::solve(arma::trimatl(L), mu_B);
  arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)), tmp);

  // ---------------------------------------------------------
  // THREAD-SAFE RNG
  // Using thread_local ensures each worker thread gets its own
  // independent generator, initialized only once per thread.
  // ---------------------------------------------------------
  thread_local std::random_device rd;
  thread_local std::mt19937 gen(rd());
  std::normal_distribution<double> dist(0.0, 1.0);

  arma::vec z(invB.n_cols);
  for(size_t i = 0; i < z.n_elem; ++i) {
    z[i] = dist(gen);
  }
  // ---------------------------------------------------------

  arma::vec v = arma::solve(arma::trimatu(arma::trans(L)), z);

  arma::vec result = v + alpha;

  return result;
}
