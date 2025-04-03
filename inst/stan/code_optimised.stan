data {
  // number of sites
  int<lower = 1> n;
  int<lower = 0> ncov_psi;
  matrix[n, ncov_psi] X_psi;

  int<lower = 0> ncov_ord;
  matrix[n, ncov_ord] X_ord;

  int<lower = 0> d;

  // total number of samples
  int<lower = 1> N;
  // total number of marker replicates
  int<lower = 1> N2;
  // total number of technical replicates
  int<lower = 1> N3;

  // samples per site
  int<lower = 0> M[n];
  int<lower = 0> sumM[n];
  // sample-level detection covariates
  int<lower = 1> ncov_theta;
  matrix[N, ncov_theta] X_theta;

  // max number of marker
  int<lower = 1> maxL;

  // marker of samples
  int<lower = 1> L[N];
  int<lower = 0> sumL[N];

  // PCR per marker
  int<lower = 1> K[N2];
  int<lower = 0> sumK[N2];

  // number of species
  int<lower = 1> S;

  // survey level information
  real logy1[N3, S];
  int logy_na[N3, S];

  int delta[N, S];

  // priors
  real prior_beta_psi;
  real<lower = 0> prior_beta_psi_sd;
  real prior_beta_theta;
  real<lower = 0> prior_beta_theta_sd;

  real a_sigma0;
  real b_sigma0;
  real a_sigma1;
  real b_sigma1;
  real a_p;
  real b_p;
  real a_q;
  real b_q;
  real a_theta0;
  real b_theta0;

}
parameters {

  matrix[ncov_psi, S] beta_psi;
  matrix[ncov_theta, S] beta_theta;
  matrix[ncov_ord, d] beta_ord;
  vector[S] beta0_psi;

  matrix[n, d] E;
  matrix[d, S] LL;

  real<lower = 0, upper = 1> pi0;
  real<lower = 0> sigma0;

  real<lower = 0> mu1;
  real<lower = 0> sigma1;

  real<lower = 0, upper = 1> p[maxL, S];
  real<lower = 0, upper = 1> theta0[S];
  real<lower = 0, upper = 1> q[maxL, S];

}

transformed parameters {

  matrix[n, d] U = X_ord * beta_ord + E;
  matrix[n, S] logit_psi = rep_matrix(beta0_psi', n) + X_psi * beta_psi + U * LL;
  matrix[N, S] logit_theta = X_theta * beta_theta;

  matrix[n, S] log_psi = log_inv_logit(logit_psi);
  matrix[n, S] log1m_psi = log1m_inv_logit(logit_psi);

  real log_theta0[S] = log(theta0);
  real log1m_theta0[S];
  log1m_theta0 = log1m(theta0);

  matrix[N, S] log_theta = log_inv_logit(logit_theta);
  matrix[N, S] log1m_theta = log1m_inv_logit(logit_theta);

}
model {

  matrix[N3, S] f0y; // likelihood of a given PCR
  matrix[N3, S] f1y;

  real log_p_yz1;
  real log_p_yz0;
  real log_p_ydelta1;
  real log_p_ydelta0;

  matrix[N, S] py_ism_delta0_w0; // likelihood of a occupied sample with all zeros
  matrix[N, S] py_ism_delta0_w1; // likelihood of a non occupied sample with all zeros
  py_ism_delta0_w0 = rep_matrix(0, N, S);
  py_ism_delta0_w1 = rep_matrix(0, N, S);

  sigma0 ~ gamma(a_sigma0, b_sigma0);
  sigma1 ~ gamma(a_sigma1, b_sigma1);

  // if (sigma0 >= sigma1) {
    //   target += negative_infinity();  // Apply a penalty to reject invalid samples
    // }

    for(s in 1:S){

      for(cov in 1:ncov_psi){
        beta_psi[cov, s] ~ normal(prior_beta_psi, prior_beta_psi_sd);
      }

      for(cov in 1:ncov_theta){
        beta_theta[cov, s] ~ normal(prior_beta_theta, prior_beta_theta_sd);
      }

      for(l in 1:maxL){
        p[l,s] ~ beta(a_p, b_p);
        q[l,s] ~ beta(a_q, b_q);
      }

      theta0[s] ~ beta(a_theta0, b_theta0);

    }

    for(cov in 1:ncov_ord){
      for(k in 1:d){

        beta_ord[cov, k] ~ normal(prior_beta_psi, prior_beta_psi_sd);
      }
    }

    for(s in 1:S){
      beta0_psi[s] ~ normal(0, 1);
    }

    // prior on factor loadings
    for(k in 1:d){
      for(s in 1:S){
        LL[k,s] ~ normal(0, 1);
      }
    }

    // prior on factor scores

    for(i in 1:n){
      for(k in 1:d){
        E[i,k] ~ normal(0, 1);
      }
    }

    // computation of f0y and f1y

    for(s in 1:S){

      for (i in 1:n) {

        for(m in 1:M[i]){

          if(delta[sumM[i] + m,s] == 1){

// for(l in 1:L[sumM[i] + m]) {
  //
  //   py_ism_delta0_w1[sumM[i] + m,s] =
  //   K[sumL[sumM[i] + m] + l] * (log(1 - p[l,s]) + log(pi0) + normal_lpdf(0 | log(1), .00001));
  //
  //   py_ism_delta0_w0[sumM[i] + m,s] =
  //   K[sumL[sumM[i] + m] + l] * (log(1 - q[l,s]) + log(pi0) + normal_lpdf(0 | log(1), .00001));
  //
  // }

  // for(l in 1:L[sumM[i] + m]){
    //
    //   for(k in 1:K[sumL[sumM[i] + m] + l]){
      //
      //     py_ism_delta0_w0[sumM[i] + m,s] += log(1 - p[l,s]);
      //     py_ism_delta0_w1[sumM[i] + m,s] += log(1 - q[l,s]);
      //
      //     py_ism_delta0_w0[sumM[i] + m,s] +=
      //     // log(pi0) + normal_lpdf(0 | log(1), .00001);
      //     log_sum_exp(
        //         log(pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | log(1), .00001),
        //         log(1 - pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | 0, sigma0)
        //         );
        //
        //       py_ism_delta0_w1[sumM[i] + m,s] +=
        //       log(pi0) + normal_lpdf(0 | log(1), .00001);
        //       // log_sum_exp(
          //         //   log(pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | log(1), .00001),
          //         //   log(1 - pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | 0, sigma0)
          //         //   );
          //
          //   }
          //
          // }

          for(l in 1:L[sumM[i] + m]){

             py_ism_delta0_w1[sumM[i] + m,s] +=
                  K[sumL[sumM[i] + m] + 1] * log(1 - p[l, s]) + log(pi0) + normal_lpdf(0 | log(1), .00001);

             py_ism_delta0_w0[sumM[i] + m,s] +=
                   K[sumL[sumM[i] + m] + 1] * log(1 - q[l, s]) + log(pi0) + normal_lpdf(0 | log(1), .00001);


            // for(k in 1:K[sumL[sumM[i] + m] + l]){

                  // py_ism_delta0_w1[sumM[i] + m,s] +=
                  // log(1 - p[l, s]) + log(pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | log(1), .00001);

                  // log_sum_exp(
                  //   // log(p[l, s]) +
                  //   // normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | mu1, sigma1),
                  //   0,
                  //   log(1 - p[l, s]) +
                  //   log_sum_exp(
                  //     log(pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | log(1), .00001),
                  //     // log(1 - pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | 0, sigma0)
                  //     0
                  //     )
                  //     );
                  // // py_ism_delta0_w1[sumM[i] + m,s] +=
                  // // log_sum_exp(
                  // //   log(p[l, s]) +
                  // //   normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | mu1, sigma1),
                  // //   log(1 - p[l, s]) +
                  // //   log_sum_exp(
                  // //     log(pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | log(1), .00001),
                  // //     log(1 - pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | 0, sigma0)
                  // //     )
                  // //     );

                       // py_ism_delta0_w0[sumM[i] + m,s] +=
                       // log(1 - q[l, s]) + log(pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | log(1), .00001);

                  // log_sum_exp(
                  //   0,
                  //   // log(q[l, s]) +
                  //   // normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | mu1, sigma1),
                  //   log(1 - q[l, s]) +
                  //   log_sum_exp(
                  //     log(pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | log(1), .00001),
                  //     // log(1 - pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | 0, sigma0)
                  //     0
                  //     )
                  //     );

            // }

          }

          } else {

            for(l in 1:L[sumM[i] + m]){

              for(k in 1:K[sumL[sumM[i] + m] + l]){

                // if(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] == 0){
                  //   f0y[sumK[sumL[sumM[i] + m] + l] + k,s] =
                  //   log(pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | log(1), .00001);
                  // } else {
                    //   f0y[sumK[sumL[sumM[i] + m] + l] + k,s] =
                    //   log(1 - pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | 0, sigma0);
                    // }

                    f0y[sumK[sumL[sumM[i] + m] + l] + k,s] =
                    log_sum_exp(
                      log(pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | log(1), .00001),
                      log(1 - pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | 0, sigma0)
                      );

                      f1y[sumK[sumL[sumM[i] + m] + l] + k,s] =
                      normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | mu1, sigma1);

              }

            }

          }

        }

      }

    }

    // compute of loglikelihood

    for(s in 1:S){

      for (i in 1:n) {

        // log probability of y given z = 1
        log_p_yz1 = 0;

        for(m in 1:M[i]){

          // log probability of y given delta = 1
          log_p_ydelta1 = 0;

          if(delta[sumM[i] + m,s] == 1){

            log_p_ydelta1 = py_ism_delta0_w1[sumM[i] + m,s];

          } else {

            for(l in 1:L[sumM[i] + m]){

              for(k in 1:K[sumL[sumM[i] + m] + l]){

                if (logy_na[sumK[sumL[sumM[i] + m] + l] + k,s] == 0){

                  log_p_ydelta1 +=
                  log_sum_exp(
                    log(p[l, s]) + f1y[sumK[sumL[sumM[i] + m] + l] + k,s],
                    log(1 - p[l, s]) + f0y[sumK[sumL[sumM[i] + m] + l] + k,s]
                    );

                }
              }

            }

          }

          // log probability of y given delta = 0
          log_p_ydelta0 = 0;

          if(delta[sumM[i] + m,s] == 1){

            log_p_ydelta0 = py_ism_delta0_w0[sumM[i] + m,s];

          } else {

            for(l in 1:L[sumM[i] + m]){

              for(k in 1:K[sumL[sumM[i] + m] + l]){

                // log_p_ydelta0 += f0y[sumK[sumL[sumM[i] + m] + l] + k,s];

                if (logy_na[sumK[sumL[sumM[i] + m] + l] + k,s] == 0){

                  log_p_ydelta0 +=
                  log_sum_exp(
                    log(q[l, s]) + f1y[sumK[sumL[sumM[i] + m] + l] + k,s],
                    log(1 - q[l, s]) + f0y[sumK[sumL[sumM[i] + m] + l] + k,s]
                    );

                }

              }

            }

          }

          log_p_yz1 += log_sum_exp(
            log_theta[sumM[i] + m,s] + log_p_ydelta1,
            log1m_theta[sumM[i] + m,s] + log_p_ydelta0);

        }

        // log probability of y given z = 0
        log_p_yz0 = 0;

        for(m in 1:M[i]){

          // log probability of y given delta = 1
          log_p_ydelta1 = 0;

          if(delta[sumM[i] + m,s] == 1){

            log_p_ydelta1 = py_ism_delta0_w1[sumM[i] + m,s];

          } else {

            for(l in 1:L[sumM[i] + m]){

              for(k in 1:K[sumL[sumM[i] + m] + l]){

                if (logy_na[sumK[sumL[sumM[i] + m] + l] + k,s] == 0){
                  log_p_ydelta1 +=
                  log_sum_exp(
                    log(p[l, s]) + f1y[sumK[sumL[sumM[i] + m] + l] + k,s],
                    log(1 - p[l, s]) + f0y[sumK[sumL[sumM[i] + m] + l] + k,s]
                    );

                }

              }

            }

          }

          // log probability of y given delta = 0
          log_p_ydelta0 = 0;

          if(delta[sumM[i] + m,s] == 1){

            log_p_ydelta0 = py_ism_delta0_w0[sumM[i] + m,s];

          } else {

            for(l in 1:L[sumM[i] + m]){

              for(k in 1:K[sumL[sumM[i] + m] + l]){

                if (logy_na[sumK[sumL[sumM[i] + m] + l] + k,s] == 0){

                  log_p_ydelta0 +=
                  log_sum_exp(
                    log(q[l, s]) + f1y[sumK[sumL[sumM[i] + m] + l] + k,s],
                    log(1 - q[l, s]) + f0y[sumK[sumL[sumM[i] + m] + l] + k,s]
                    );
                }

              }

            }

          }

          log_p_yz0 += log_sum_exp(
            log_theta0[s] + log_p_ydelta1,
            log1m_theta0[s] + log_p_ydelta0);

        }

        target += log_sum_exp(
          log_psi[i, s] + log_p_yz1,
          log1m_psi[i, s] + log_p_yz0
          );

      }

    }

}
generated quantities{

//   real loglik = 0;
//
//     // computation of f0y and f1y
//
//     for(s in 1:S){
//
//       for (i in 1:n) {
//
//         for(m in 1:M[i]){
//
//           if(delta[sumM[i] + m,s] == 1){
//
// // for(l in 1:L[sumM[i] + m]) {
//   //
//   //   py_ism_delta0_w1[sumM[i] + m,s] =
//   //   K[sumL[sumM[i] + m] + l] * (log(1 - p[l,s]) + log(pi0) + normal_lpdf(0 | log(1), .00001));
//   //
//   //   py_ism_delta0_w0[sumM[i] + m,s] =
//   //   K[sumL[sumM[i] + m] + l] * (log(1 - q[l,s]) + log(pi0) + normal_lpdf(0 | log(1), .00001));
//   //
//   // }
//
//   // for(l in 1:L[sumM[i] + m]){
//     //
//     //   for(k in 1:K[sumL[sumM[i] + m] + l]){
//       //
//       //     py_ism_delta0_w0[sumM[i] + m,s] += log(1 - p[l,s]);
//       //     py_ism_delta0_w1[sumM[i] + m,s] += log(1 - q[l,s]);
//       //
//       //     py_ism_delta0_w0[sumM[i] + m,s] +=
//       //     // log(pi0) + normal_lpdf(0 | log(1), .00001);
//       //     log_sum_exp(
//         //         log(pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | log(1), .00001),
//         //         log(1 - pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | 0, sigma0)
//         //         );
//         //
//         //       py_ism_delta0_w1[sumM[i] + m,s] +=
//         //       log(pi0) + normal_lpdf(0 | log(1), .00001);
//         //       // log_sum_exp(
//           //         //   log(pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | log(1), .00001),
//           //         //   log(1 - pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | 0, sigma0)
//           //         //   );
//           //
//           //   }
//           //
//           // }
//
//           for(l in 1:L[sumM[i] + m]){
//
//             for(k in 1:K[sumL[sumM[i] + m] + l]){
//
//                   py_ism_delta0_w1[sumM[i] + m,s] +=
//                   log_sum_exp(
//                     // log(p[l, s]) +
//                     // normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | mu1, sigma1),
//                     0,
//                     log(1 - p[l, s]) +
//                     log_sum_exp(
//                       log(pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | log(1), .00001),
//                       log(1 - pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | 0, sigma0)
//                       )
//                       );
//                   // py_ism_delta0_w1[sumM[i] + m,s] +=
//                   // log_sum_exp(
//                   //   log(p[l, s]) +
//                   //   normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | mu1, sigma1),
//                   //   log(1 - p[l, s]) +
//                   //   log_sum_exp(
//                   //     log(pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | log(1), .00001),
//                   //     log(1 - pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | 0, sigma0)
//                   //     )
//                   //     );
//
//                        py_ism_delta0_w0[sumM[i] + m,s] +=
//                   log_sum_exp(
//                     // 0,
//                     log(q[l, s]) +
//                     normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | mu1, sigma1),
//                     log(1 - q[l, s]) +
//                     log_sum_exp(
//                       log(pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | log(1), .00001),
//                       log(1 - pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | 0, sigma0)
//                       )
//                       );
//
//             }
//
//           }
//
//           } else {
//
//             for(l in 1:L[sumM[i] + m]){
//
//               for(k in 1:K[sumL[sumM[i] + m] + l]){
//
//                 // if(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] == 0){
//                   //   f0y[sumK[sumL[sumM[i] + m] + l] + k,s] =
//                   //   log(pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | log(1), .00001);
//                   // } else {
//                     //   f0y[sumK[sumL[sumM[i] + m] + l] + k,s] =
//                     //   log(1 - pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | 0, sigma0);
//                     // }
//
//                     f0y[sumK[sumL[sumM[i] + m] + l] + k,s] =
//                     log_sum_exp(
//                       log(pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | log(1), .00001),
//                       log(1 - pi0) + normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | 0, sigma0)
//                       );
//
//                       f1y[sumK[sumL[sumM[i] + m] + l] + k,s] =
//                       normal_lpdf(logy1[sumK[sumL[sumM[i] + m] + l] + k,s] | mu1, sigma1);
//
//               }
//
//             }
//
//           }
//
//         }
//
//       }
//
//     }
//
//     // compute of loglikelihood
//
//     for(s in 1:S){
//
//       for (i in 1:n) {
//
//         // log probability of y given z = 1
//         log_p_yz1 = 0;
//
//         for(m in 1:M[i]){
//
//           // log probability of y given delta = 1
//           log_p_ydelta1 = 0;
//
//           if(delta[sumM[i] + m,s] == 1){
//
//             log_p_ydelta1 = py_ism_delta0_w1[sumM[i] + m,s];
//
//           } else {
//
//             for(l in 1:L[sumM[i] + m]){
//
//               for(k in 1:K[sumL[sumM[i] + m] + l]){
//
//                 if (logy_na[sumK[sumL[sumM[i] + m] + l] + k,s] == 0){
//
//                   log_p_ydelta1 +=
//                   log_sum_exp(
//                     log(p[l, s]) + f1y[sumK[sumL[sumM[i] + m] + l] + k,s],
//                     log(1 - p[l, s]) + f0y[sumK[sumL[sumM[i] + m] + l] + k,s]
//                     );
//
//                 }
//               }
//
//             }
//
//           }
//
//           // log probability of y given delta = 0
//           log_p_ydelta0 = 0;
//
//           if(delta[sumM[i] + m,s] == 1){
//
//             log_p_ydelta0 = py_ism_delta0_w0[sumM[i] + m,s];
//
//           } else {
//
//             for(l in 1:L[sumM[i] + m]){
//
//               for(k in 1:K[sumL[sumM[i] + m] + l]){
//
//                 // log_p_ydelta0 += f0y[sumK[sumL[sumM[i] + m] + l] + k,s];
//
//                 if (logy_na[sumK[sumL[sumM[i] + m] + l] + k,s] == 0){
//
//                   log_p_ydelta0 +=
//                   log_sum_exp(
//                     log(q[l, s]) + f1y[sumK[sumL[sumM[i] + m] + l] + k,s],
//                     log(1 - q[l, s]) + f0y[sumK[sumL[sumM[i] + m] + l] + k,s]
//                     );
//
//                 }
//
//               }
//
//             }
//
//           }
//
//           log_p_yz1 += log_sum_exp(
//             log_theta[sumM[i] + m,s] + log_p_ydelta1,
//             log1m_theta[sumM[i] + m,s] + log_p_ydelta0);
//
//         }
//
//         // log probability of y given z = 0
//         log_p_yz0 = 0;
//
//         for(m in 1:M[i]){
//
//           // log probability of y given delta = 1
//           log_p_ydelta1 = 0;
//
//           if(delta[sumM[i] + m,s] == 1){
//
//             log_p_ydelta1 = py_ism_delta0_w1[sumM[i] + m,s];
//
//           } else {
//
//             for(l in 1:L[sumM[i] + m]){
//
//               for(k in 1:K[sumL[sumM[i] + m] + l]){
//
//                 if (logy_na[sumK[sumL[sumM[i] + m] + l] + k,s] == 0){
//                   log_p_ydelta1 +=
//                   log_sum_exp(
//                     log(p[l, s]) + f1y[sumK[sumL[sumM[i] + m] + l] + k,s],
//                     log(1 - p[l, s]) + f0y[sumK[sumL[sumM[i] + m] + l] + k,s]
//                     );
//
//                 }
//
//               }
//
//             }
//
//           }
//
//           // log probability of y given delta = 0
//           log_p_ydelta0 = 0;
//
//           if(delta[sumM[i] + m,s] == 1){
//
//             log_p_ydelta0 = py_ism_delta0_w0[sumM[i] + m,s];
//
//           } else {
//
//             for(l in 1:L[sumM[i] + m]){
//
//               for(k in 1:K[sumL[sumM[i] + m] + l]){
//
//                 if (logy_na[sumK[sumL[sumM[i] + m] + l] + k,s] == 0){
//
//                   log_p_ydelta0 +=
//                   log_sum_exp(
//                     log(q[l, s]) + f1y[sumK[sumL[sumM[i] + m] + l] + k,s],
//                     log(1 - q[l, s]) + f0y[sumK[sumL[sumM[i] + m] + l] + k,s]
//                     );
//                 }
//
//               }
//
//             }
//
//           }
//
//           log_p_yz0 += log_sum_exp(
//             log_theta0[s] + log_p_ydelta1,
//             log1m_theta0[s] + log_p_ydelta0);
//
//         }
//
//         loglik += log_sum_exp(
//           log_psi[i, s] + log_p_yz1,
//           log1m_psi[i, s] + log_p_yz0
//           );
//
//       }
//
//     }

}
