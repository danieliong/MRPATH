#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppDist.h>
// #include <RcppArmadilloExtensions/sample.h>
#include <math.h>
#include "utils.h"

using namespace Rcpp;

arma::mat computeExpectedHessian(const arma::cube &normISamps,
    const double &m_X, const double &lambdaX, const arma::vec &pis,
    const arma::vec &mus, const arma::vec &sds)
{
    int p = (int)normISamps.n_rows;
    int K = (int)pis.n_elem;

    int k;

    int pi_k_ind;
    int mu_k_ind;
    int var_k_ind;

    double E_Nk;
    double E_NK;

    // H is block diagonal matrix
    arma::mat H((2+(3*K)-1), (2+(3*K)-1), arma::fill::zeros);

    // muX prior params
    H(0, 0) = -(p / pow(lambdaX, 2));
    H(1, 1) = (p / (2 * pow(lambdaX, 4))) - ((1 / pow(lambdaX, 6)) *
    arma::accu(normISamps.slice(0) % arma::pow((normISamps.slice(1) - m_X),2)));
    H(1, 0) = - pow(lambdaX, -4) *
    arma::accu(normISamps.slice(0) % (normISamps.slice(1) - m_X));
    H(0, 1) = H(1, 0);

    E_NK = arma::accu(normISamps.slice(0) % normISamps.slice(3+(K-1)));

    for (k = 0; k < K; k++) {
        E_Nk = arma::accu(normISamps.slice(0) % normISamps.slice(3+k));

        if (k == (K-1)) {
            mu_k_ind = 2 + (3*k);
            var_k_ind = 2 + (3*k) + 1;
        } else {
            pi_k_ind = 2 + (3*k);
            mu_k_ind = 2 + (3*k) + 1;
            var_k_ind = 2 + (3*k) + 2;

            H(pi_k_ind, pi_k_ind) = -(E_Nk / pow(pis(k), 2)) -
            (E_NK / pow(pis(K-1), 2));
        }

        // mu_k
        H(mu_k_ind, mu_k_ind) = - (E_Nk / pow(sds(k), 2));

        // var_k
        H(var_k_ind, var_k_ind) = (E_Nk / (2*pow(sds(k), 4))) - ((1 / pow(sds(k), 6)) *
        arma::accu(normISamps.slice(0) % normISamps.slice(3+k) %
        arma::pow((normISamps.slice(2) - mus(k)), 2)));

        // diagonal terms
        // H(mu_k_ind, var_k_ind) = - (1 / pow(sds(k), 4)) *
        // arma::accu(normISamps.slice(0) % normISamps.slice(3+k) %
        // (normISamps.slice(2) - mus(k)));
        H(mu_k_ind, var_k_ind) = -arma::accu(normISamps.slice(0) %
        normISamps.slice(3+k) % (normISamps.slice(2) - mus(k))) / pow(sds(k), 4);

        H(var_k_ind, mu_k_ind) = H(mu_k_ind, var_k_ind);
    }

    return H;
}


arma::mat computeExpectedGradient(const arma::cube &normISamps,
    const double &m_X, const double &lambdaX, const arma::vec &pis,
    const arma::vec &mus, const arma::vec &sds, const bool &matrix = false)
{
    int p = (int)normISamps.n_rows;
    int K = (int)pis.n_elem;

    int k;

    int pi_k_ind;
    int mu_k_ind;
    int var_k_ind;

    arma::mat G(2, p);

    arma::vec E_N_ik(p);
    arma::vec E_N_iK(p);

    // m_X
    G.row(0) = ((1 / pow(lambdaX, 2)) * arma::sum(normISamps.slice(0) %
    (normISamps.slice(1) - m_X), 1)).t();

    // lambdaX2
    G.row(1) = - (1 / (2 * pow(lambdaX, 2))) + (1 / (2 * pow(lambdaX, 4))) *
    arma::sum(normISamps.slice(0) % pow((normISamps.slice(1) - m_X), 2), 1).t();

    E_N_iK = arma::sum(normISamps.slice(0) % normISamps.slice(3+(K-1)), 1);

    for (k = 0; k < K; k++) {
        E_N_ik = arma::sum(normISamps.slice(0) % normISamps.slice(3+k), 1);

        // if not last component
        if (k != (K-1)) {
            // pi_k
            G.insert_rows(G.n_rows,
                ((E_N_ik / pis(k)) - (E_N_iK / pis(K-1))).t()
            );
        }

        // mu_k
        G.insert_rows(G.n_rows,
            ((1 / pow(sds(k), 2)) * arma::sum((normISamps.slice(0) %
            normISamps.slice(3+k) % (normISamps.slice(2) - mus(k))), 1)).t()
        );

        // var_k
        G.insert_rows(G.n_rows,
            (-(E_N_ik / (2 * pow(sds(k), 2))) +
            (1 / (2 * pow(sds(k), 4))) *
            arma::sum((normISamps.slice(0) % normISamps.slice(3+k) %
            pow((normISamps.slice(2) - mus(k)), 2)), 1)).t()
        );
    }

    if (!matrix) {
        G = arma::sum(G, 1);
    }

    return G;

    // if (vector){
    //     arma::vec grad((2+(3*K) - 1));
    //     grad = arma::sum(G, 1);
    //     return grad;
    // } else {
    //     return G;
    // }
}


arma::mat computeGradSamples_i(const int &i, const arma::cube &normISamps,
    const double &m_X, const double &lambdaX, const arma::vec &pis,
    const arma::vec &mus, const arma::vec &sds)
{
    // const int Nparams = 2 + (3*K) - 1;
    const int p = (int)normISamps.n_rows;
    const int K = (int)pis.n_elem;
    const int M = normISamps.n_cols;

    int k;

    arma::mat gradSamples_i(2, M);

    gradSamples_i.row(0) = (1 / pow(lambdaX, 2)) *
    (normISamps.slice(1).row(i) - m_X);


    gradSamples_i.row(1) = - (1 / (2 * pow(lambdaX, 2))) +
    (1 / (2 * pow(lambdaX, 4))) * pow((normISamps.slice(1).row(i) - m_X), 2);

    for (k = 0; k < K; ++k) {
        if (k != (K-1)) {
            gradSamples_i.insert_rows(gradSamples_i.n_rows,
                (normISamps.slice(3+k).row(i) / pis(k)) -
                (normISamps.slice(3+(K-1)).row(i) / pis(K-1))
            );
        }

        gradSamples_i.insert_rows(gradSamples_i.n_rows,
            (1 / pow(sds(k), 2)) * (normISamps.slice(3+k).row(i) %
            (normISamps.slice(2).row(i) - mus(k)))
        );


        gradSamples_i.insert_rows(gradSamples_i.n_rows,
            - (normISamps.slice(3+k).row(i) / (2 * pow(sds(k), 2))) +
            (1 / (2* pow(sds(k), 4))) * (normISamps.slice(3+k).row(i) %
            pow((normISamps.slice(2).row(i) - mus(k)), 2))
        );
    }

    return gradSamples_i;
}


arma::mat computeExpectedGradOuter(const arma::cube &normISamps,
    const double &m_X, const double &lambdaX, const arma::vec &pis,
    const arma::vec &mus, const arma::vec &sds)
{
    int p = (int)normISamps.n_rows;
    int K = (int)pis.n_elem;
    int M = normISamps.n_cols;
    const int Nparams = 2 + (3*K) - 1;

    int i, j;

    arma::mat gradOuter(Nparams, Nparams);
    arma::mat gradSamples_i(Nparams, M);

    arma::mat gradMat = computeExpectedGradient(normISamps,
        m_X, lambdaX, pis, mus, sds, true);

    for (i = 0; i < p; ++i) {
        gradSamples_i = computeGradSamples_i(i, normISamps,
            m_X, lambdaX, pis, mus, sds);

        gradOuter += (gradSamples_i.each_row() % normISamps.slice(0).row(i)) *
        gradSamples_i.t();

        if (i < (p - 1)) {
            for (j = (i+1); j < p; ++j) {
                gradOuter += (gradMat.col(i) * gradMat.col(j).t()) +
                (gradMat.col(j) * gradMat.col(i).t());
            }
        }
    }

    return gradOuter;
}


arma::mat computeObservedInfoMatrix(const arma::vec &X, const arma::vec &Y,
    const arma::vec &seX, const arma::vec &seY, const double &m_X,
    const double &lambdaX, const arma::vec &pis, const arma::vec &mus,
    const arma::vec &sds, int Nsamples = 10000)
{
    int K = (int)pis.n_elem;
    int Nparams = 2 + (3*K) - 1;

    int k;

    arma::mat J = arma::eye(Nparams, Nparams);
    arma::mat I(Nparams, Nparams);

    arma::cube normISamps = sampleLatentVarPostCube(Nsamples, X, Y, seX, seY,
        m_X, lambdaX, pis, mus, sds, true);

    arma::mat H = computeExpectedHessian(normISamps, m_X, lambdaX,
        pis, mus, sds);

    // Rcout << "H" << std::endl;
    // Rcout << H << std::endl;

    arma::mat grad = computeExpectedGradient(normISamps, m_X, lambdaX,
        pis, mus, sds);

    // Rcout << "grad" << std::endl;
    // Rcout << grad << std::endl;

    arma::mat gradOuter = computeExpectedGradOuter(normISamps, m_X, lambdaX,
        pis, mus, sds);

    // Rcout << "gradOuter" << std::endl;
    // Rcout << gradOuter << std::endl;

    // Jacobian
    J = arma::eye(Nparams, Nparams);
    J(1,1) = 2 * lambdaX;
    for (k = 0; k < (K-1); k++) {
        J((4+(3*k)), (4+(3*k))) = 2 * sds(k);
    }
    J((Nparams-1), (Nparams-1)) = 2 * sds((K-1));


    I = J.t() * (-H - gradOuter + (grad * grad.t())) * J;

    return I;
}
