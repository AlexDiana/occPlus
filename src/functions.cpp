#include <RcppArmadillo.h>
#include <random>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]

#include "rng.h"

// inline double exprnd(double rate) {
//   std::exponential_distribution<double> dist(rate);
//   return dist(get_rng());
// }

using namespace Rcpp;

//' Seed the package's C++ random number generators
//'
//' Sets the base seed from which every worker thread derives its own stream.
//' Called by \code{runOccJSDM()} with a draw from R's RNG, which is what makes
//' \code{set.seed()} control the sampler. Not intended to be called directly.
//'
//' @param seed a non-negative integer.
//' @keywords internal
// [[Rcpp::export]]
void setOccJSDMSeed(unsigned int seed) {
  occjsdm_rng_base_seed() = seed;
  // Bump the generation so that threads whose engine was already constructed
  // re-seed on their next draw instead of continuing their old stream.
  occjsdm_rng_generation() += 1u;
}

// [[Rcpp::depends(RcppArmadillo)]]

// TRUNC = 0.64 and TRUNC_RECIP = 1/0.64, the truncation point from
// Windle (2013), were declared here but never used, and drew an
// -Wunused-const-variable warning that made R CMD check report a
// WARNING on installation. Removed rather than silenced; restore them
// if the Polya-Gamma sampler is ever changed to use the truncated form.

// Mathematical constants computed using Wolfram Alpha
#define MATH_PI        3.141592653589793238462643383279502884197169399375105820974
#define MATH_PI_2      1.570796326794896619231321691639751442098584699687552910487
#define MATH_2_PI      0.636619772367581343075535053490057448137838582961825794990
#define MATH_PI2       9.869604401089358618834490999876151135313699407240790626413
#define MATH_PI2_2     4.934802200544679309417245499938075567656849703620395313206
#define MATH_SQRT1_2   0.707106781186547524400844362104849039284835937688474036588
#define MATH_SQRT_PI_2 1.253314137315500251207882642405522626503493370304969158314
#define MATH_LOG_PI    1.144729885849400174143427351353058711647294812915311571513
#define MATH_LOG_2_PI  -0.45158270528945486472619522989488214357179467855505631739
#define MATH_LOG_PI_2  0.451582705289454864726195229894882143571794678555056317392

// old code

static double aterm(int n, double x, double t) {
  double f = 0;
  if(x <= t) {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) + 1.5*(MATH_LOG_2_PI- (double)std::log(x)) - 2*(n + 0.5)*(n + 0.5)/x;
  }
  else {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) - x * MATH_PI2_2 * (n + 0.5)*(n + 0.5);
  }
  return (double)exp(f);
}

static double exprnd(double mu) {
  // return -mu * (double)std::log(1.0 - (double)R::runif(0.0,1.0));
  return -mu * (double)std::log(1.0 - (double)runif());
}

static double truncgamma() {
  double c = MATH_PI_2;
  double X, gX;

  bool done = false;
  while(!done)
  {
    X = exprnd(1.0) * 2.0 + c;
    gX = MATH_SQRT_PI_2 / (double)std::sqrt(X);

    // if(R::runif(0.0,1.0) <= gX) {
    if(runif() <= gX) {
      done = true;
    }
  }

  return X;
}

static double randinvg(double mu) {
  // sampling

  // double u = R::rnorm(0.0,1.0);
  double u = rnorm();
  double V = u*u;
  double out = mu + 0.5*mu * ( mu*V - (double)std::sqrt(4.0*mu*V + mu*mu * V*V) );

  // if(R::runif(0.0,1.0) > mu /(mu+out)) {
  if(runif() > mu /(mu+out)) {
    out = mu*mu / out;
  }
  return out;
}

static double tinvgauss(double z, double t) {
  double X, u;
  double mu = 1.0/z;

  // Pick sampler
  if(mu > t) {
    // Sampler based on truncated gamma
    // Algorithm 3 in the Windle (2013) PhD thesis, page 128
    while(1) {
      // u = R::runif(0.0, 1.0);
      u = runif();
      X = 1.0 / truncgamma();

      if ((double)std::log(u) < (-z*z*0.5*X)) {
        break;
      }
    }
  }
  else {
    // Rejection sampler
    X = t + 1.0;
    while(X >= t) {
      X = randinvg(mu);
    }
  }
  return X;
}

static double samplepg(double z) {
  //  PG(b, z) = 0.25 * J*(b, z/2)
  z = (double)std::fabs((double)z) * 0.5;

  // Point on the intersection IL = [0, 4/ log 3] and IR = [(log 3)/pi^2, \infty)
  double t = MATH_2_PI;

  // Compute p, q and the ratio q / (q + p)
  // (derived from scratch; derivation is not in the original paper)
  double K = z*z/2.0 + MATH_PI2/8.0;
  double logA = (double)std::log(4.0) - MATH_LOG_PI - z;
  double logK = (double)std::log(K);
  double Kt = K * t;
  double w = (double)std::sqrt(MATH_PI_2);

  double logf1 = logA + R::pnorm(w*(t*z - 1),0.0,1.0,1,1) + logK + Kt;
  double logf2 = logA + 2*z + R::pnorm(-w*(t*z+1),0.0,1.0,1,1) + logK + Kt;
  double p_over_q = (double)std::exp(logf1) + (double)std::exp(logf2);
  double ratio = 1.0 / (1.0 + p_over_q);

  double u, X;

  // Main sampling loop; page 130 of the Windle PhD thesis
  while(1)
  {
    // Step 1: Sample X ? g(x|z)
    // u = R::runif(0.0,1.0);
    u = runif();
    if(u < ratio) {
      // truncated exponential
      X = t + exprnd(1.0)/K;
    }
    else {
      // truncated Inverse Gaussian
      X = tinvgauss(z, t);
    }

    // Step 2: Iteratively calculate Sn(X|z), starting at S1(X|z), until U ? Sn(X|z) for an odd n or U > Sn(X|z) for an even n
    int i = 1;
    double Sn = aterm(0, X, t);
    // double U = R::runif(0.0,1.0) * Sn;
    double U = runif() * Sn;
    int asgn = -1;
    bool even = false;

    while(1)
    {
      Sn = Sn + asgn * aterm(i, X, t);

      // Accept if n is odd
      if(!even && (U <= Sn)) {
        X = X * 0.25;
        return X;
      }

      // Return to step 1 if n is even
      if(even && (U > Sn)) {
        break;
      }

      even = !even;
      asgn = -asgn;
      i++;
    }
  }
  return X;
}

static double rpg(int n, double z){

  double x = 0;
  for(int i = 0; i < n; i++){
    x += samplepg(z);
  }

  return(x);
}

static arma::vec mvrnormArmaQuick(arma::vec mu, arma::mat cholsigma) {
  int ncols = cholsigma.n_cols;
  arma::vec Y = arma::randn(ncols);
  return mu + cholsigma * Y;
}

static arma::vec mvrnormArmaQuick_TS(const arma::vec& mu, const arma::mat& cholsigma) {
  int ncols = cholsigma.n_cols;
  arma::vec Y(ncols);

  // rnorm() is the shared per-thread engine from rng.h: thread-safe like the
  // std::random_device engine it replaces, but seeded from R, so this draw is
  // now reproducible under set.seed(). This is the draw inside
  // sample_beta_cpp_TS(), which sample_betatheta_cpp_parallel() calls from
  // inside `#pragma omp parallel for`.
  for(int i = 0; i < ncols; ++i) {
    Y[i] = rnorm();
  }

  return mu + cholsigma * Y;
}

static arma::mat diagMatrixProd(arma::mat& X, arma::vec& D){

  arma::mat result(X.n_rows, D.size());
  for(int i = 0; i < result.n_rows; i++){
    for(int j = 0; j < result.n_cols; j++){
      result(i, j) = X(i,j) * D(j);
    }
  }

  return(result);
}

static arma::vec sample_beta_cpp(arma::mat& X, arma::mat& B, arma::vec& b,
                          arma::vec& Omega, arma::vec& k){

  // arma::mat cov_matrix = arma::inv(arma::trans(X) * Omega * X + arma::inv(B));
  arma::mat tX = arma::trans(X);
  arma::mat tXOmega = diagMatrixProd(tX, Omega);
  // arma::mat cov_matrix = arma::inv(tXOmega * X + arma::inv(B));
  // arma::vec result = mvrnormArma(cov_matrix * (arma::trans(X) * k + arma::inv(B) * b), cov_matrix);

  arma::mat L = arma::trans(arma::chol(tXOmega * X + arma::inv(B)));
  arma::vec tmp = arma::solve(arma::trimatl(L), tX * k + arma::inv(B) * b);
  arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);

  arma::vec result = mvrnormArmaQuick(alpha, arma::trans(arma::inv(arma::trimatl(L))));

  return(result);
}

static arma::vec sample_beta_cpp_TS(arma::mat& X, arma::mat& B, arma::vec& b,
                          arma::vec& Omega, arma::vec& k){

  // arma::mat cov_matrix = arma::inv(arma::trans(X) * Omega * X + arma::inv(B));
  arma::mat tX = arma::trans(X);
  arma::mat tXOmega = diagMatrixProd(tX, Omega);
  // arma::mat cov_matrix = arma::inv(tXOmega * X + arma::inv(B));
  // arma::vec result = mvrnormArma(cov_matrix * (arma::trans(X) * k + arma::inv(B) * b), cov_matrix);

  arma::mat L = arma::trans(arma::chol(tXOmega * X + arma::inv(B)));
  arma::vec tmp = arma::solve(arma::trimatl(L), tX * k + arma::inv(B) * b);
  arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);

  arma::vec result = mvrnormArmaQuick_TS(alpha, arma::trans(arma::inv(arma::trimatl(L))));

  return(result);
}

static arma::vec sample_Omega_cpp(arma::mat& X, arma::vec& beta, arma::vec& n){

  int nsize = n.size();
  arma::vec Omega_vec(nsize);

  for(int i = 0; i < nsize; i++){

    arma::vec b = X.row(i) * beta;
    Omega_vec[i] = rpg(n[i], b[0]);

  }

  return(Omega_vec);
}

static arma::vec sample_beta_nocov_cpp(arma::vec beta, arma::mat& X, arma::vec b,
                                arma::mat B, arma::vec n, arma::vec k){

  arma::vec Omega = sample_Omega_cpp(X, beta, n);

  beta = sample_beta_cpp(X, B, b, Omega, k);

  return(beta);
}

static arma::vec sample_beta_nocov_cpp_TS(arma::vec beta, arma::mat& X, arma::vec b,
                                arma::mat B, arma::vec n, arma::vec k){

  arma::vec Omega = sample_Omega_cpp(X, beta, n);

  beta = sample_beta_cpp_TS(X, B, b, Omega, k);

  return(beta);
}

// [[Rcpp::export]]
NumericMatrix sample_z_cpp(const NumericMatrix& w,
                           const NumericMatrix& psi,
                           const NumericMatrix& theta,
                           const NumericVector& theta0,
                           const IntegerVector& M,
                           const IntegerVector& sumM) {

  int S = psi.ncol();
  int n = psi.nrow();
  NumericMatrix z(n, S);

  for (int s = 0; s < S; s++) {
    for (int i = 0; i < n; i++) {

      // compute p_zsequal1
      double log_p1 = 0.0;
      for (int m = 0; m < M[i]; m++) {
        int idx = sumM[i] + m;
        log_p1 += R::dbinom(w(idx, s), 1.0, theta(idx, s), true);
      }
      log_p1 += R::dbinom(1.0, 1.0, psi(i, s), true);

      // compute p_zsequal0
      double log_p0 = 0.0;
      for (int m = 0; m < M[i]; m++) {
        int idx = sumM[i] + m;
        log_p0 += R::dbinom(w(idx, s), 1.0, theta0[s], true);
      }
      log_p0 += R::dbinom(0.0, 1.0, psi(i, s), true);

      // probability
      double maxlog = std::max(log_p1, log_p0);
      double p1 = std::exp(log_p1 - maxlog);
      double p0 = std::exp(log_p0 - maxlog);
      double p_1 = p1 / (p1 + p0);

      // sample z[i, s]
      z(i, s) = R::rbinom(1.0, p_1);
    }
  }
  return z;
}

// [[Rcpp::export]]
NumericMatrix sample_w_cpp(const NumericMatrix& logy1,
                           double mu0, double sigma0,
                           double mu1, double sigma1,
                           const NumericMatrix& theta,
                           const NumericVector& theta0,
                           const NumericMatrix& p,
                           const NumericMatrix& q,
                           const IntegerVector& M,
                           const IntegerVector& K,
                           const IntegerVector& sumL,
                           const IntegerVector& sumM,
                           const IntegerVector& sumK,
                           int maxL,
                           const NumericMatrix& z) {

  int S = theta.ncol();
  int N = theta.nrow();
  int n = M.size();

  NumericMatrix w(N, S);

  for (int s = 0; s < S; s++) {
    for (int i = 0; i < n; i++) {
      for (int m = 0; m < M[i]; m++) {

        // compute log p(w = 1)
        double log_p1 = 0.0;
        for (int l = 0; l < maxL; l++) {
          int idxL = sumL[sumM[i] + m] + l;
          for (int k = 0; k < K[idxL]; k++) {
            int idxK = sumK[idxL] + k;
            if(logy1(idxK, s) == 0){
              log_p1 += log(1 - p(l,s));
            } else {
              log_p1 += log(p(l,s)) + R::dnorm(logy1(idxK, s), mu1, sigma1, true);
            }
          }
        }

        // compute log p(w = 0)
        double log_p0 = 0.0;
        for (int l = 0; l < maxL; l++) {
          int idxL = sumL[sumM[i] + m] + l;
          for (int k = 0; k < K[idxL]; k++) {
            int idxK = sumK[idxL] + k;
            if(logy1(idxK, s) == 0){
              log_p0 += log(1 - q(l,s));
            } else {
              log_p0 += log(q(l,s)) + R::dnorm(logy1(idxK, s), mu0, sigma0, true);
            }
          }
        }

        // conditional on z[i, s]
        if (z(i, s) == 1.0) {
          // log_p1 += log(theta(sumM[i] + m, s));
          // log_p0 += log(1 - theta(sumM[i] + m, s));;//R::dbinom(0.0, 1.0, theta(sumM[i] + m, s), true);
          log_p1 += R::dbinom(1.0, 1.0, theta(sumM[i] + m, s), true);
          log_p0 += R::dbinom(0.0, 1.0, theta(sumM[i] + m, s), true);
        } else {
          log_p1 += //log(theta0[s]);
            R::dbinom(1.0, 1.0, theta0[s], true);
          log_p0 += //log(1 - theta0[s]);
            R::dbinom(0.0, 1.0, theta0[s], true);
        }

        // numerical stability
        double maxlog = std::max(log_p1, log_p0);
        double p1exp = std::exp(log_p1 - maxlog);
        double p0exp = std::exp(log_p0 - maxlog);
        double p_ws1 = p1exp / (p1exp + p0exp);

        // sample w
        w(sumM[i] + m, s) = R::rbinom(1.0, p_ws1);
      }
    }
  }

  return w;
}

// [[Rcpp::export]]
NumericMatrix sample_w_cim_cipp(const NumericMatrix& y,
                                const NumericMatrix& y_NA,
                                const NumericMatrix& theta,
                                const NumericVector& theta0,
                                const NumericMatrix& p,
                                const NumericMatrix& q,
                                const IntegerVector& M,
                                const IntegerVector& K,
                                const IntegerVector& sumL,
                                const IntegerVector& sumM,
                                const IntegerVector& sumK,
                                const IntegerVector& P,
                                const IntegerVector& primerId,
                                int maxL,
                                const NumericMatrix& z) {

  int S = theta.ncol();
  int N = theta.nrow();
  int n = M.size();

  NumericMatrix w(N, S);

  // Precompute logs for p and q. Rows here are indexed by *global* primer
  // identity (1..maxL, the number of distinct primers in the whole dataset),
  // not by a primer's position within a given sample's block, since P (the
  // number of primers actually used) can differ across samples.
  NumericMatrix log_p(maxL, S), log_1p(maxL, S);
  NumericMatrix log_q(maxL, S), log_1q(maxL, S);

  for(int s = 0; s < S; s++) {
    for(int l = 0; l < maxL; l++) {
      log_p(l, s) = std::log(p(l, s));
      log_1p(l, s) = std::log(1.0 - p(l, s));
      log_q(l, s) = std::log(q(l, s));
      log_1q(l, s) = std::log(1.0 - q(l, s));
    }
  }

  for (int s = 0; s < S; s++) {
    for (int i = 0; i < n; i++) {
      for (int m = 0; m < M[i]; m++) {

        int nprimers = P[sumM[i] + m];

        // compute log p(w = 1)
        double log_p1 = 0.0;
        for (int l = 0; l < nprimers; l++) {
          int idxL = sumL[sumM[i] + m] + l;
          int primerRow = primerId[idxL] - 1; // 0-based global primer identity
          for (int k = 0; k < K[idxL]; k++) {

            int idxK = sumK[idxL] + k;

            if(y_NA(idxK, s) == 0){

              log_p1 += y(idxK, s) * log_p(primerRow, s) + (1 - y(idxK, s)) * log_1p(primerRow, s);

            }

          }
        }

        // compute log p(w = 0)
        double log_p0 = 0.0;
        for (int l = 0; l < nprimers; l++) {
          int idxL = sumL[sumM[i] + m] + l;
          int primerRow = primerId[idxL] - 1;
          for (int k = 0; k < K[idxL]; k++) {
            int idxK = sumK[idxL] + k;

            if(y_NA(idxK, s) == 0){

              log_p0 += y(idxK, s) * log_q(primerRow,s) + (1 - y(idxK, s)) * log_1q(primerRow,s);

            }
          }
        }

        // conditional on z[i, s]
        if (z(i, s) == 1.0) {
          log_p1 += log(theta(sumM[i] + m, s));
          log_p0 += log(1 - theta(sumM[i] + m, s));;//R::dbinom(0.0, 1.0, theta(sumM[i] + m, s), true);
          // log_p1 += R::dbinom(1.0, 1.0, theta(sumM[i] + m, s), true);
          // log_p0 += R::dbinom(0.0, 1.0, theta(sumM[i] + m, s), true);
        } else {
          log_p1 += log(theta0[s]);
            // R::dbinom(1.0, 1.0, theta0[s], true);
          log_p0 += log(1 - theta0[s]);
            // R::dbinom(0.0, 1.0, theta0[s], true);
        }

        // numerical stability
        double maxlog = std::max(log_p1, log_p0);
        double p1exp = std::exp(log_p1 - maxlog);
        double p0exp = std::exp(log_p0 - maxlog);
        double p_ws1 = p1exp / (p1exp + p0exp);

        // sample w
        w(sumM[i] + m, s) = R::rbinom(1.0, p_ws1);
      }
    }
  }

  return w;
}


// [[Rcpp::export]]
arma::mat sample_betatheta_cpp(const arma::mat& w,
                               const arma::mat& z,
                               arma::mat beta_theta,
                               const arma::uvec& idx_z,
                               const arma::mat& X_theta,
                               const arma::vec& b_betatheta,
                               const arma::mat& B_betatheta) { // Pass the R function if it's not natively in C++

  int S = beta_theta.n_cols;

  arma::uvec idx_z_cpp = idx_z - 1;

  arma::mat z_all = z.rows(idx_z_cpp);

  // Loop over columns (S)
  for (int s = 0; s < S; ++s) {

    // Get the s-th column of z_all
    arma::vec z_col = z_all.col(s);

    // Find indices where z_all[, s] == 1
    // Armadillo's find() returns a uvec (unsigned vector of indices)
    arma::uvec find_ones = find(z_col == 1);

    if (find_ones.is_empty()) continue; // Safeguard if no rows equal 1

    // k <- as.vector(t(w[z_all[,s]==1,s])) - .5
    // In C++, w.submat(find_ones, uvec({static_cast<uword>(s)})) extracts the column subset
    arma::vec w_sub = w.elem(find_ones + s * w.n_rows);
    arma::vec k = w_sub - 0.5;

    // n <- rep(1, length(k))
    arma::vec n = arma::ones<arma::vec>(k.n_elem);

    // X_thetasubset <- X_theta[z_all[,s]==1,,drop=F]
    arma::mat X_thetasubset = X_theta.rows(find_ones);

    // Extract current column of beta_theta
    arma::vec beta_sub = beta_theta.col(s);

    beta_theta.col(s) = sample_beta_nocov_cpp(beta_sub, X_thetasubset, b_betatheta, B_betatheta, n, k);

  }

  return beta_theta;
}

// [[Rcpp::export]]
arma::mat sample_betatheta_cpp_parallel(const arma::mat& w,
                                        const arma::mat& z,
                                        arma::mat beta_theta, // Passed by value from R, safe to modify locally
                                        const arma::uvec& idx_z,
                                        const arma::mat& X_theta,
                                        const arma::vec& b_betatheta,
                                        const arma::mat& B_betatheta) {

  int S = beta_theta.n_cols;

  arma::uvec idx_z_cpp = idx_z - 1;
  arma::mat z_all = z.rows(idx_z_cpp);

  // Tell OpenMP to parallelize this loop.
  // All variables declared inside the loop become private to each thread.
  #pragma omp parallel for
  for (int s = 0; s < S; ++s) {

    // Get the s-th column of z_all
    arma::vec z_col = z_all.col(s);

    // Find indices where z_all[, s] == 1
    arma::uvec find_ones = arma::find(z_col == 1);

    if (find_ones.is_empty()) continue;

    // Extract the column subset
    arma::vec w_sub = w.elem(find_ones + s * w.n_rows);
    arma::vec k = w_sub - 0.5;

    arma::vec n = arma::ones<arma::vec>(k.n_elem);
    arma::mat X_thetasubset = X_theta.rows(find_ones);

    // Extract current column of beta_theta
    arma::vec beta_sub = beta_theta.col(s);

    beta_theta.col(s) = sample_beta_nocov_cpp_TS(beta_sub, X_thetasubset, b_betatheta, B_betatheta, n, k);

  }

  return beta_theta;
}

// [[Rcpp::export]]
List sample_pq_cpp(NumericMatrix& c_imk,
                   IntegerMatrix& y_NA,
                   NumericMatrix w,
                   IntegerVector idx_p_k, IntegerVector idx_w_k,
                   int maxP, double a_p, double b_p,
                   double a_q, double b_q) {

  int S = w.ncol();
  int n_k = idx_w_k.size();

  NumericMatrix p(maxP, S);
  NumericMatrix q(maxP, S);

  // Main loop through columns (S)
  for (int s = 0; s < S; ++s) {

    for (int l = 0; l < maxP; ++l) {

      int w1_primerl_cases_1 = 0;
      int w1_primerl_cases_0 = 0;
      int w0_primerl_cases_1 = 0;
      int w0_primerl_cases_0 = 0;

      for (int i = 0; i < n_k; ++i) {

        int idx_ki = idx_w_k[i] - 1;

        if(y_NA(i,s) == 0){

          if(idx_p_k[i] == (l+1) && w(idx_ki, s) == 1 && c_imk(i,s) == 1){
            w1_primerl_cases_1 += 1;
          } else if(idx_p_k[i] == (l+1) && w(idx_ki, s) == 1 && c_imk(i,s) == 0){
            w1_primerl_cases_0 += 1;
          } else if(idx_p_k[i] == (l+1) && w(idx_ki, s) == 0 && c_imk(i,s) == 2){
            w0_primerl_cases_1 += 1;
          } else if(idx_p_k[i] == (l+1) && w(idx_ki, s) == 0 && c_imk(i,s) == 0){
            w0_primerl_cases_0 += 1;
          }

        }

      }

      p(l, s) = R::rbeta(a_p + w1_primerl_cases_1, b_p + w1_primerl_cases_0);
      q(l, s) = R::rbeta(a_q + w0_primerl_cases_1, b_q + w0_primerl_cases_0);
    }
  }

  return List::create(
    _["p"] = p,
    _["q"] = q
  );
}


//// WAIC FUNCTIONS

// Inline helper for the logistic function
inline double logistic(double x) {
  return 1.0 / (1.0 + std::exp(-x));
}

// [[Rcpp::export]]
NumericVector computeModelLoglikJSDM_cpp(NumericMatrix z,
                                         NumericMatrix eta,
                                         String model,
                                         Nullable<NumericVector> tau = R_NilValue) {
  int n = z.nrow();
  int S = z.ncol();
  NumericVector logliks(n * S);
  int idx = 0;

  if (model == "continuous") {
    if (tau.isNull()) {
      stop("tau must be provided for continuous model");
    }
    NumericVector tau_vec(tau);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < S; ++j) {
        logliks[idx++] = R::dnorm(z(i, j), eta(i, j), tau_vec[j], 1);
      }
    }
  } else if (model == "binary") {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < S; ++j) {
        double psi = logistic(eta(i, j));
        logliks[idx++] = R::dbinom(z(i, j), 1, psi, 1);
      }
    }
  } else if (model == "counts") {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < S; ++j) {
        double mu = std::exp(eta(i, j));
        logliks[idx++] = R::dpois(z(i, j), mu, 1);
      }
    }
  } else {
    stop("Unknown model type");
  }

  return logliks;
}

// [[Rcpp::export]]
NumericVector computeModelLoglikFirstStage_cpp(NumericMatrix w,
                                               NumericMatrix z,
                                               NumericMatrix theta,
                                               NumericVector theta0,
                                               IntegerVector idx_z_w) {
  int n_w = w.nrow();
  int S_w = w.ncol();
  NumericVector logliks(n_w * S_w);
  int idx = 0;

  for (int i = 0; i < n_w; ++i) {
    // Convert 1-based R index to 0-based C++ index
    int z_row_idx = idx_z_w[i] - 1;

    for (int s = 0; s < S_w; ++s) {
      if (z(z_row_idx, s) == 1) {
        logliks[idx++] = R::dbinom(w(i, s), 1, theta(i, s), 1);
      } else {
        logliks[idx++] = R::dbinom(w(i, s), 1, theta0[s], 1);
      }
    }
  }

  return logliks;
}

// [[Rcpp::export]]
NumericVector computeModelLoglikSecondStage_cpp(NumericMatrix y, NumericMatrix w, NumericMatrix p, NumericMatrix q, IntegerVector idx_w_k, IntegerVector idx_p_k) {
  int n_y = y.nrow();
  int S_y = y.ncol();
  NumericVector logliks(n_y * S_y);
  int idx = 0;

  // Outer loop is columns (s), inner is rows (i) to match the R version's flattening behavior
  for (int s = 0; s < S_y; ++s) {
    for (int i = 0; i < n_y; ++i) {
      int w_row_idx = idx_w_k[i] - 1; // Convert to 0-based index
      int p_row_idx = idx_p_k[i] - 1; // Convert to 0-based index

      if (w(w_row_idx, s) == 1) {
        logliks[idx++] = R::dbinom(y(i, s), 1, p(p_row_idx, s), 1);
      } else {
        logliks[idx++] = R::dbinom(y(i, s), 1, q(p_row_idx, s), 1);
      }
    }
  }

  return logliks;
}
