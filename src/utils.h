#ifndef _UTILS_H_
#define _UTILS_H_

using namespace Rcpp;

arma::cube getImportanceSamplesCube(int N_samples, const arma::vec &X,
    const arma::vec &Y, const arma::vec &seX, const arma::vec &seY,
    const double &m_X, const double &lambdaX, arma::vec pis,
    const arma::vec &mus, const arma::vec &sds,
    const bool &normalize_weights = false);

arma::vec normalizeImptWeights(arma::cube &ISamps);

void appendLatentVarPostCube(arma::cube &postSamps, const arma::vec &sum_weights,
    const int &M, const arma::vec &X, const arma::vec &Y,
    const arma::vec &seX, const arma::vec &seY, const double &m_X,
    const double &lambdaX, arma::vec pis, const arma::vec &mus, const arma::vec &sds,
    int &N_samples);

arma::mat computeQMatrix(const int &p, const int &K, const arma::vec &Y,
    const arma::vec &seY, const arma::cube &ISamps,
    const double &m_X, const double &lambdaX,
    const arma::vec &pis, const arma::vec &mus, const arma::vec &sds, const double &tau);

void MR_EM_Mstep(
    const int &p, const int &K, arma::vec Y, arma::vec seY, const bool &equalSds, const bool &overDispersedY, const arma::cube &ISamps,
    double &m_X, double &lambdaX, arma::vec &pis, arma::vec &mus, arma::vec &sds, double &tau);

bool MR_EM_diagnostic(
    const arma::cube &ISamps, const int &p, const int &K, const arma::vec &Y, const arma::vec &seY, const double &eps,
    const double &z_alpha, const double &z_gamma, const int &M,
    double &m_X, double &lambdaX, arma::vec &pis, arma::vec &mus,
    arma::vec &sds, double &tau, const double &prevIter_m_X, const double &prevIter_lambdaX,
    const arma::vec &prevIterPis, const arma::vec &prevIterMus,
    const arma::vec &prevIterSds, const double &prevIterTau, const bool &verbose, int &N_samples, const int &N_iters, bool &rejected);


#endif
