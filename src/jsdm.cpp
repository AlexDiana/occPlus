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

// [[Rcpp::export]]
double rinvgamma_cpp(double a, double b){
  return 1 / R::rgamma(a, 1 / b);
}

// [[Rcpp::export]]
bool isPointInBandRight(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j){

  for(int k = 0; k < X_tilde.n_rows; k++){

    if((X_tilde(k,1) < y_grid[j + 1]) && (X_tilde(k,1) > y_grid[j - 1])){
      if(X_tilde(k,0) < x_grid[i + 1]){
        return(true);
      }
    }

  }

  return(false);
}

// [[Rcpp::export]]
bool isPointInBandLeft(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j) {

  for(int k = 0; k < X_tilde.n_rows; k++){

    if((X_tilde(k,1) < y_grid[j + 1]) && (X_tilde(k,1) > y_grid[j - 1])){
      if(X_tilde(k,0) > x_grid[i - 1]){
        return(true);
      }
    }

  }

  return(false);
}

// [[Rcpp::export]]
bool isPointInBandUp(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j){

  for(int k = 0; k < X_tilde.n_rows; k++){

    if((X_tilde(k,0) < x_grid[i + 1]) && (X_tilde(k,0) > x_grid[i - 1])){
      if(X_tilde(k,1) > y_grid[j-1]){
        return(true);
      }
    }

  }

  return(false);

}

// [[Rcpp::export]]
bool isPointInBandDown(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j){

  for(int k = 0; k < X_tilde.n_rows; k++){

    if((X_tilde(k,0) < x_grid[i + 1]) && (X_tilde(k,0) > x_grid[i - 1])){
      if(X_tilde(k,1) < y_grid[j+1]){
        return(true);
      }
    }

  }

  return(false);

}

// [[Rcpp::export]]
IntegerVector findClosestPoint(arma::mat XY_sp, arma::mat X_tilde){

  IntegerVector closestPoint(XY_sp.n_rows);

  for(int k = 0; k < XY_sp.n_rows; k++){

    double newDistance = 0;
    double minDistance = exp(50);
    int bestIndex = 0;

    for(int i = 0; i < X_tilde.n_rows; i++){
      newDistance = pow(X_tilde(i, 0) - XY_sp(k, 0), 2) + pow(X_tilde(i, 1) - XY_sp(k, 1), 2);

      if(newDistance < minDistance){
        minDistance = newDistance;
        bestIndex = i + 1;
      }
    }

    closestPoint[k] = bestIndex;

  }

  return(closestPoint);
}


// [[Rcpp::export]]
arma::mat dist_matrix(const arma::mat& coords) {

  // coords: n x d matrix (e.g. x,y or x,y,z)
  int n = coords.n_rows;
  int d = coords.n_cols;

  arma::mat D(n, n, arma::fill::zeros);

  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {

      double dist_ij = 0.0;

      for (int k = 0; k < d; k++) {
        double diff = coords(i, k) - coords(j, k);
        dist_ij += diff * diff;
      }

      dist_ij = std::sqrt(dist_ij);

      D(i, j) = dist_ij;
      D(j, i) = dist_ij; // symmetry
    }
  }

  return D;
}

// [[Rcpp::export]]
arma::mat gpCovMatrix(const arma::mat& D,
                      double sigma2,
                      double rho) {

  int n = D.n_rows;
  arma::mat Sigma(n, n);

  // Gaussian exponential covariance
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Sigma(i, j) = sigma2 * std::exp(-D(i, j) / rho);
    }
  }

  return Sigma;
}

// GAUSSIAN PROCESS FUNCTIONS

double k_cpp(double x1, double x2, double a, double l){
  // return pow(1 + (x1-x2)*(x1-x2), - alphaGP);
  return a*exp(-(x1-x2)*(x1-x2)/(2*pow(l,2)));
  // return 1;
}

// [[Rcpp::export]]
arma::mat K(arma::vec x1, arma::vec x2, double a, double l){
  arma::mat res(x1.size(), x2.size());

  for(int i = 0; (unsigned)i < x1.size(); i++){
    for(int j = 0; (unsigned)j < x2.size(); j++){
      res(i,j) = k_cpp(x1[i],x2[j], a, l);
    }
  }

  return res;
}

double k2_cpp(arma::rowvec x1, arma::rowvec x2, double a, double l){
  // return pow(1 + (x1-x2)*(x1-x2), - alphaGP);
  return a*exp(-( pow(x1[0]-x2[0], 2) + pow(x1[1]-x2[1], 2) ) /(2*pow(l,2)));
}

// [[Rcpp::export]]
arma::mat K2(arma::mat x1, arma::mat x2, double a, double l){
  arma::mat res(x1.n_rows, x2.n_rows);

  for(int i = 0; (unsigned)i < x1.n_rows; i++){
    for(int j = 0; (unsigned)j < x2.n_rows; j++){
      res(i,j) = k2_cpp(x1.row(i),x2.row(j), a, l);
    }
  }

  return res;
}

// [[Rcpp::export]]
arma::mat samplePGvariables(arma::mat &Xbeta){

  int n1 = Xbeta.n_rows;
  int n2 = Xbeta.n_cols;

  arma::mat Omega_mat(n1, n2);

#pragma omp parallel for collapse(2)
  for(int i = 0; i < n1; i++){
    for(int j = 0; j < n2; j++){
      Omega_mat(i,j) = rpg(1, Xbeta(i,j));
    }
  }
  // #pragma omp parallel for
  //  for(int k=0; k<n1*n2; k++){
  //    Omega_mat[k]=rpg(1,Xbeta[k]);
  //  }

  return(Omega_mat);
}

static arma::vec mvrnormArmaQuick(arma::vec mu, arma::mat cholsigma) {
  int ncols = cholsigma.n_cols;
  arma::vec Y = arma::randn(ncols);
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

arma::vec logistic(arma::vec x) {
  return 1.0 / (1.0 + exp(-x));
}

double quantile(arma::vec x, double p) {

  const int n = x.n_elem;
  if (p <= 0.0) return x(0);
  if (p >= 1.0) return x(n - 1);

  double h = (n - 1) * p;
  int j = std::floor(h);
  double g = h - j;

  if (j == n - 1)
    return x(n - 1);

  return (1.0 - g) * x(j) + g * x(j + 1);
}

// [[Rcpp::export]]
arma::cube computeNewOutputs(
    const arma::mat& X,
    const arma::mat& B0_output,
    const arma::cube& B_output,
    const arma::cube& Ks_all,
    const arma::cube& Bs_output,
    const arma::cube& L_output,
    const arma::vec sigmah_output,
    const arma::vec idx_ls_output,
    const arma::vec& conflevels,
    bool useEnvCov,
    bool useSpatial,
    bool useBiotic,
    std::string model,
    bool verbose)
{

  // The number of new sites comes from whichever term is actually in use.
  // Reading it from X unconditionally gave n = 0 whenever useEnvCov was
  // FALSE, because R passes a 0 x 0 matrix in that case, so the result was
  // silently empty. Neither term active means nothing determines how many
  // sites to predict for, which is a caller error rather than a zero-row
  // answer.
  int n = useEnvCov ? X.n_rows : (useSpatial ? Ks_all.n_rows : 0);
  if (n == 0) {
    Rcpp::stop("Nothing determines the new sites: useEnvCov and useSpatial "
               "are both off, so neither X_psi nor X_s gives a site count.");
  }

  int S = B_output.n_cols;
  int d = L_output.n_rows;
  int niter = B_output.n_slices;

  arma::mat mcmc_output(n, niter);
  arma::cube output(3, n, S);

  for (int j = 0; j < S; j++) {

    if(verbose){
      Rcout << "Computing species " << j + 1 << " out of " << S << std::endl;
    }


    for (int iter = 0; iter < niter; iter++) {

      arma::vec B = B_output.slice(iter).col(j);

      arma::mat U = arma::randn(n, d) * sigmah_output[iter];
      arma::vec L = L_output.slice(iter).col(j);

      arma::vec linpred(n, arma::fill::value(B0_output(j, iter)));
      // arma::vec linpred = B0_output(j, iter) + X * B + Ks * Bs + U * L;

      if(useEnvCov){
        linpred = linpred + X * B;
      }

      // Ks_all and Bs_output are only indexable when the fit estimated a
      // spatial field; with useSpatial off R passes empty stand-ins, so
      // slicing them here rather than above avoids an out-of-bounds abort.
      if(useSpatial){
        int idx_ls = idx_ls_output[iter];
        arma::mat Ks = Ks_all.slice(idx_ls - 1);
        arma::vec Bs = Bs_output.slice(iter).col(j);
        linpred = linpred + Ks * Bs;
      }

      if(useBiotic){
        linpred = linpred + U * L;
      }


      if(model == "continuous"){
        mcmc_output.col(iter) = linpred;
      } else if (model == "binary"){
        mcmc_output.col(iter) = logistic(linpred);
      }


    }

    for (int i = 0; i < n; i++) {

      arma::vec mcmc_output_i_sorted = arma::sort(mcmc_output.row(i).t());

      output(0, i, j) = quantile(mcmc_output_i_sorted, conflevels[0]);
      output(1, i, j) = quantile(mcmc_output_i_sorted, conflevels[1]);
      output(2, i, j) = quantile(mcmc_output_i_sorted, conflevels[2]);

    }
  }

  return output;
}

// [[Rcpp::export]]
arma::cube convert_to_correlation(arma::cube L_output_vec, int niter, int S, int d) {
  // Initialize the output 3D array (cube) with dimensions: niter x S x S
  arma::cube Lambda_output(S, S, niter, arma::fill::zeros);

  for (int iter = 0; iter < niter; ++iter) {
    // 1. Extract the slice for the current iteration and treat it as a matrix.
    // In Armadillo, cube slicing gives us an S x d view.
    // Because R stores arrays in column-major order, we transpose or copy carefully
    // to match R's `matrix(..., byrow = TRUE)` behavior.
    arma::mat L_output_current(S, d);

    // Replicating byrow = TRUE from the original R code
    for (int i = 0; i < S; ++i) {
      for (int j = 0; j < d; ++j) {
        L_output_current(i, j) = L_output_vec(iter, i, j);
      }
    }

    // 2. Compute the covariance matrix: L %*% t(L)
    arma::mat cov_mat = L_output_current * L_output_current.t();

    // 3. Convert Covariance to Correlation (Equivalent to cov2cor)
    // R's cov2cor divides by the square root of the diagonal elements: inv(diag(sqrt(cov)))
    arma::vec inv_sqrt_diag = 1.0 / arma::sqrt(cov_mat.diag());
    arma::mat cor_mat = cov_mat % (inv_sqrt_diag * inv_sqrt_diag.t());

    // 4. Store the result in the current slice of our output cube
    Lambda_output.slice(iter) = cor_mat;
  }

  // Note: RcppArmadillo automatically converts an arma::cube back into a 3D R array.
  // The dimensions in R will be [S, S, niter].
  return Lambda_output;
}

// JSDM sampling functions

// Sample B in a regression model Y ~ N(XB, sigma^2 I)

// [[Rcpp::export]]
arma::vec sampleBuniv(arma::mat& X, arma::mat& B,
                      arma::vec& b, arma::vec& y,
                      double sigma){

  arma::mat invB = arma::inv(B);

  // arma::mat cov_matrix = arma::inv(arma::trans(X) * Omega * X + arma::inv(B));
  arma::mat tX = arma::trans(X);
  // arma::mat cov_matrix = arma::inv(tXOmega * X + arma::inv(B));
  // arma::vec result = mvrnormArma(cov_matrix * (arma::trans(X) * k + arma::inv(B) * b), cov_matrix);

  arma::mat L = arma::trans(arma::chol(tX * X * (1 / (sigma*sigma)) + invB));
  arma::vec tmp = arma::solve(arma::trimatl(L), tX * y / (sigma*sigma) + invB * b);
  arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);

  arma::vec result = mvrnormArmaQuick(alpha, arma::trans(arma::inv(arma::trimatl(L))));

  return(result);
}

// Sample B in a regression model Y ~ N(XB, Omega)

// [[Rcpp::export]]
arma::vec sampleB(arma::mat& X, arma::mat& B, arma::vec& b,
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

// [[Rcpp::export]]
arma::mat sample_U_cpp(arma::mat& k,
                       arma::mat& L,
                       arma::mat& XB,
                       arma::mat& XsBs,
                       arma::mat& Omega,
                       double sigma_h,
                       std::string model) {

  int d = L.n_rows;
  int S = k.n_cols;
  int n = k.n_rows;

  arma::mat U(n, d, arma::fill::zeros);

  if (d == 0)
    return U;

  // Compute k_new
  arma::mat k_new(n, S);

  if (model == "continuous") {

    k_new = (k - (XB + XsBs)) % Omega;

  } else if (model == "binary") {

    k_new = k - Omega % (XB + XsBs);
  } else {
    Rcpp::stop("Unknown model");
  }

  arma::mat B_current = arma::eye(d,d) * sigma_h * sigma_h;
  arma::vec b_current(d, arma::fill::zeros);

  arma::mat transL = arma::trans(L);

  // Parallel loop
  // #ifdef _OPENMP
  //   omp_set_num_threads(n_threads);
  // #endif
  //
  // #pragma omp parallel for
  for(int i = 0; i < n; i++) {

    arma::vec Omega_i = arma::conv_to<arma::vec>::from(Omega.row(i));
    arma::vec k_new_i = arma::conv_to<arma::vec>::from(k_new.row(i));

    arma::vec result = sampleB(
      transL,
      B_current,
      b_current,
      Omega_i,
      k_new_i
    );

    U.row(i) = arma::conv_to<arma::rowvec>::from(result);
  }


  return U;
}

// SPATIAL APPROXIMATOR FUNCTIONS

// [[Rcpp::export]]
arma::mat XsBs(arma::mat &A,
               arma::mat &B,
               arma::mat &X_s_centers){

  int n = A.n_rows;
  int m = B.n_rows;
  int maxPoints = X_s_centers.n_cols;

  arma::mat AB_out = arma::mat(n, maxPoints);

  for(int i = 0; i < n; i++){
    for(int j = 0; j < maxPoints; j++){
      int colIndex = X_s_centers(i, j) - 1;
      for(int l = 0; l < m; l++){
        AB_out(i, j) += A(i, l) * B(l, colIndex);
      }
    }
  }

  return(AB_out);

}

// [[Rcpp::export]]
arma::mat KsBproduct(arma::mat &Ks,
                     arma::mat &B,
                     arma::mat &X_s_centers){

  int n = Ks.n_rows;
  int S = B.n_cols;
  int maxPoints = X_s_centers.n_cols;

  arma::mat XB = arma::zeros(n, S);

  for(int s = 0; s < S; s++){
    for(int i = 0; i < n; i++){
      for(int j = 0; j < maxPoints; j++){
        XB(i, s) += Ks(i, j) * B(X_s_centers(i, j) - 1, s);
      }
    }
  }

  return(XB);
}

// [[Rcpp::export]]
arma::mat XtOmegaX_SoR(arma::mat X,
                       int X_centers,
                       arma::vec Omega,
                       arma::mat X_s_index,
                       arma::mat &X_s_sor){

  int p = X.n_cols;

  arma::mat XtOmegaX = arma::zeros(X_centers + p,
                                   X_centers + p);

  int maxPoints = X_s_index.n_cols;

  // spatial  covariates times spatial covariates
  for(int l = 0; l < maxPoints; l++){

    for(int l2 = 0; l2 < maxPoints; l2++){

      for (int i = 0; (unsigned)i < Omega.size(); i++){

        XtOmegaX(p + X_s_index(i,l) - 1, p + X_s_index(i,l2) - 1) +=
          X_s_sor(i, l) * Omega[i] * X_s_sor(i, l2);

      }

    }

  }

  // spatial covariates times standard covariates
  for (int i = 0; (unsigned)i < Omega.size(); i++){

    for(int j = 0; j < p; j++){

      for(int l = 0; l < maxPoints; l++){

        XtOmegaX(j, p + X_s_index(i,l) - 1) +=
          X(i, j) * X_s_sor(i, l) * Omega[i];

      }
    }
  }

  for (int i = 1; i <= X_centers; i++) {

    for (int j = 0; j < p; j++){

      XtOmegaX(p + i - 1, j) = XtOmegaX(j, p + i - 1);
    }

  }

  // standard covariates times standard covariates

  for(int i = 0; i < p; i++){
    for (int j = 0; j <= i; j++) {
      for(int l = 0; l < Omega.size(); l++){
        XtOmegaX(i, j) +=
          Omega[l] * X(l, i) * X(l, j);
      }
    }
  }

  for (int i = 0; i < (p- 1); i++) {
    for (int j = i; j < p; j++) {
      XtOmegaX(i, j) = XtOmegaX(j, i);
    }
  }

  return(XtOmegaX);
}

arma::vec XtK_SoR(arma::mat X, arma::mat &X_s_index, arma::mat &X_s_sor,
                  arma::vec &k, int centers){

  int p = X.n_cols;

  arma::vec Xk = arma::zeros(p + centers);

  for(int i = 0; i < p; i++){
    Xk(i) = as_scalar(k.t() * X.col(i));
  }

  for(int i = 0; i < k.size(); i++){
    for(int l = 0; l < X_s_index.n_cols; l++){
      Xk(p + X_s_index(i,l) - 1) += X_s_sor(i, l) * k[i];
    }
  }

  return(Xk);
}


// Sample B in a regression model Y ~ N((X|Xs)(B|Bs), Omega) where Xs is a
// subset of regressor approximator

// [[Rcpp::export]]
arma::vec sampleB_SoR(arma::mat X, arma::mat &invB, arma::vec &b,
                      arma::vec &k, arma::vec Omega,
                      arma::mat &X_s_index,
                      arma::mat &Ks,
                      int X_centers){

  arma::mat XtOmegaX = XtOmegaX_SoR(X, X_centers, Omega, X_s_index, Ks);
  arma::mat tXk = XtK_SoR(X, X_s_index, Ks, k, X_centers);

  arma::mat Lambda_B = XtOmegaX + invB;
  arma::vec mu_B = tXk + invB * b;

  arma::mat L = arma::trans(arma::chol(Lambda_B));
  arma::vec tmp = arma::solve(arma::trimatl(L), mu_B);
  arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);

  arma::vec z = arma::randn(invB.n_cols);
  arma::vec v = arma::solve(arma::trimatu(arma::trans(L)), z);

  arma::vec result = v + alpha;

  return(result);
}

// [[Rcpp::export]]
arma::mat spatEffectMeanCpp(arma::cube& Bs_output,
                            arma::mat& Ks,
                            arma::mat& Xs_centers) {

  int dim1 = Ks.n_rows;
  int dim2 = Bs_output.n_cols;
  int niter = Bs_output.n_slices;

  arma::mat spatEffect_mean(dim1, dim2, arma::fill::zeros);

  for (int iter = 0; iter < niter; iter++) {

    arma::mat B = Bs_output.slice(iter);

    spatEffect_mean += (1.0/niter) *
      KsBproduct(Ks, B, Xs_centers);

  }

  return spatEffect_mean;
}
