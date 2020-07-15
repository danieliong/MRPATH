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
List getImportanceSamples(const DataFrame &data, const List &MCEM_fit, int Nsamples = 50000)
{
    arma::vec X = data["beta.exposure"];
    arma::vec Y = data["beta.outcome"];
    arma::vec seX = data["se.exposure"];
    arma::vec seY = data["se.outcome"];

    List params = MCEM_fit["paramEst"];
    double m_X = params["m_X"];
    double lambdaX = params["lambdaX"];
    arma::vec pis = params["pis"];
    arma::vec mus = params["mus"];
    arma::vec sds = params["sds"];
    double tau = params["tau"];

    int K = (int)pis.n_elem;

    arma::cube samps = getImportanceSamplesCube(Nsamples, X, Y, seX, tau*seY, m_X, lambdaX, pis, mus, sds);

    return List::create(
        Named("thetaX") = samps.slice(1),
        _["beta"] = samps.slice(2),
        _["alpha"] = (arma::cube)samps.slices(3,3+K-1),
        _["W"] = samps.slice(0)
    );
}
