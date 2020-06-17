#include <RcppArmadillo.h>
#include <RcppDist.h>
// #include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <math.h>
#include "utils.h"
// #include <omp.h>
// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;


// [[Rcpp::export]]
List sampleLatentVarPost(int N_samples, const arma::vec &X, const arma::vec &Y, const arma::vec &seX,
                         const arma::vec &seY, const List &params) {

  double m_X = params["m_X"];
  double lambdaX = params["lambdaX"];
  arma::vec pis = params["pis"];
  arma::vec mus = params["mus"];
  arma::vec sds = params["sds"];
  double tau = params["tau"];

  int K = (int)pis.n_elem;

  arma::cube samps = sampleLatentVarPostCube(N_samples, X, Y, seX, tau*seY, m_X, lambdaX, pis, mus, sds);

  // arma::cube alpha_samps = samps.slices(3, 3+K-1);

  return List::create(
    Named("muX_samps") = samps.slice(1),
    _["beta_samps"] = samps.slice(2),
    _["alpha_samps"] = (arma::cube)samps.slices(3,3+K-1),
    _["W"] = samps.slice(0)
  );
}
