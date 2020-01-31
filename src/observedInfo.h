#ifndef _OBSERVEDINFO_H_
#define _OBSERVEDINFO_H_

using namespace Rcpp;

arma::mat computeExpectedHessian(const arma::cube &normISamps,
    const double &m_X, const double &lambdaX, const arma::vec &pis,
    const arma::vec &mus, const arma::vec &sds);

arma::mat computeExpectedGradient(const arma::cube &normISamps,
    const double &m_X, const double &lambdaX, const arma::vec &pis,
    const arma::vec &mus, const arma::vec &sds, const bool &matrix = false);

arma::mat computeGradSamples_i(const int &i, const arma::cube &normISamps,
    const double &m_X, const double &lambdaX, const arma::vec &pis,
    const arma::vec &mus, const arma::vec &sds);

arma::mat computeExpectedGradOuter(const arma::cube &normISamps,
    const double &m_X, const double &lambdaX, const arma::vec &pis,
    const arma::vec &mus, const arma::vec &sds);

arma::mat computeObservedInfoMatrix(const arma::vec &X, const arma::vec &Y,
    const arma::vec &seX, const arma::vec &seY, const double &m_X,
    const double &lambdaX, const arma::vec &pis, const arma::vec &mus,
    const arma::vec &sds, int Nsamples = 10000);

#endif
