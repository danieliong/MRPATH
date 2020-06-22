#include <RcppArmadillo.h>
#include <RcppDist.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <math.h>
using namespace Rcpp;

// Obtain importance samples from latent variable posterior
/*
Slices:
1 - W
2 - muX_samps
3 - beta_samps
4 - 4 + K: alpha_samps
*/
arma::cube sampleLatentVarPostCube(int N_samples, const arma::vec &X, const arma::vec &Y,
    const arma::vec &seX, const arma::vec &seY, const double &m_X,
    const double &lambdaX, arma::vec pis, const arma::vec &mus,
    const arma::vec &sds, const bool &normalize_weights = false) {

  int p = (int)X.n_elem; // data sample size
  int K = (int)pis.n_elem; // # of components

  int i, k, j;

  arma::cube postSamps(p, N_samples, 3+K);


  /* ###### Variables used for sampling mu_Xi ##### */
  double post_var_muXi; // Posterior variance of mu_Xi | Xi, params
  double post_mean_muXi; // Posterior mean of mu_Xi | Xi, params

  /* ##### Variables used for sampling beta_i ##### */
  // length N_samples vector containing sampled components for beta_i | mu_Xi,Yi,params
  // arma::mat post_pis_ik(K, N_samples);
  arma::vec post_pis_ij(K);

  arma::uvec k_i(N_samples);

  // Posterior var/mean of beta_i | mu_Xi^j, Yi, parameters for j = 1,...,N_samples
  arma::mat post_sigma2_ik_given_muXi(N_samples, 1);
  arma::mat post_mu_ik_given_muXi(N_samples, 1);

  arma::vec pis_copy(K); // copy of pis at each iteration
  // sample function changes pis

  for (i = 0; i < p; i++) {

    /*############### Draw importance samples for mu_Xi #####################*/

    // Calculate posterior mean, variance for mu_Xi | X_i, parameters
    post_var_muXi = 1 / ((1/(seX(i)*seX(i))) + (m_X / (lambdaX*lambdaX)));
    post_mean_muXi = post_var_muXi * ((X(i)/(seX(i)*seX(i))) +
    (m_X/(lambdaX*lambdaX)));

    // Sample mu_Xi from importance density g(mu_Xi) = P(mu_Xi | X_i, Y_i, params)
    postSamps.slice(1).row(i) = (sqrt(post_var_muXi) * arma::randn(N_samples) +
    post_mean_muXi).t();

    /*############### Draw importance samples for beta_i ####################*/


    // ############# Adjust posterior pi ###################
    for (j = 0; j < N_samples; ++j) {
        post_pis_ij = pis % arma::pow(
            (arma::pow(postSamps.slice(1).row(i)(j) * sds, 2) +
            pow(seY(i), 2)),
            -0.5);
        k_i(j) = RcppArmadillo::sample(arma::regspace<arma::uvec>(0,K-1), 1L, true, post_pis_ij)(0);
    }


    // // Sample components for each MC sample j = 1,...,N_samples
    // k_i = RcppArmadillo::sample(arma::regspace<arma::uvec>(0,K-1), N_samples,
    // true, pis_copy);
    // // pis cannot be const for this to work.

    // Calculate posterior mean and variance for \beta_i | \mu_Xi, X_i, Y_i, parameters
    post_sigma2_ik_given_muXi = arma::pow(((1/arma::pow(sds(k_i),2)) +
    (arma::pow(postSamps.slice(1).row(i).t(),2) / (seY(i)*seY(i)))), -1);

    post_mu_ik_given_muXi = post_sigma2_ik_given_muXi %
    ( ((Y(i)*postSamps.slice(1).row(i).t())/(seY(i)*seY(i))) + (mus(k_i) /
    arma::pow(sds(k_i),2)) );

    // Sample from P(\beta_i | \mu_Xi, X_i, Y_i, parameters)
    postSamps.slice(2).row(i) = (arma::sqrt(post_sigma2_ik_given_muXi) %
    arma::randn(N_samples) + post_mu_ik_given_muXi).t();

    // // beta_samps.row(i) = beta_samps_i.t();

    /*####### Calculate P(Z_i = k | \beta_i^j, X_i, Y_i, parameters) ########*/
    for (k = 0; k < K; k++) {
        postSamps.slice(3+k).row(i) = (pis(k) *
        arma::normpdf(postSamps.slice(2).row(i).t(), mus(k), sds(k))).t();
    }

    /*############### Calculate importance weights w_i #######################*/
    postSamps.slice(0).row(i) = (pis(k_i) %
    arma::normpdf(Y(i), postSamps.slice(1).row(i).t()%mus(k_i),
    arma::sqrt(arma::pow(sds(k_i)%postSamps.slice(1).row(i).t(),2) +
    (seY(i)*seY(i))))).t();
  }

  // Normalise alphaSamps
  postSamps.each_slice(arma::regspace<arma::uvec>(3,3+K-1)) /=
  (arma::mat)arma::sum(postSamps.slices(3, 3+K-1), 2);

  /* ################# Normalise importance weights w_i #################### */
  if (normalize_weights) {
      postSamps.slice(0).each_col() /= sum(postSamps.slice(0), 1);
  }

  /* ############ ERROR CHECKING ####)########### */
  if (!postSamps.is_finite()) {
      throw exception("Posterior samples have NaN/Inf values.");
  }

  return postSamps;
}



// Normalize importance weights in-place and return normalization constants
arma::vec normalizeImptWeights(arma::cube &ISamps) {
    arma::vec sum_weights = sum(ISamps.slice(0), 1);
    ISamps.slice(0).each_col() /= sum_weights;
    return(sum_weights);
}

// In-place M-step update of parameters m_X, lambdaX, pis, mus, sds, tau
void MR_EM_Mstep(
    const int &p, const int &K, arma::vec Y, arma::vec seY, const bool &equalSds, const bool &overDispersedY, const arma::cube &ISamps,
    double &m_X, double &lambdaX, arma::vec &pis, arma::vec &mus, arma::vec &sds, double &tau)
{

    int k;
    arma::vec vars(K);

    m_X = arma::accu(ISamps.slice(0) % ISamps.slice(1)) / p;

    lambdaX = sqrt((arma::accu(ISamps.slice(0) %
    arma::pow(ISamps.slice(1),2))/p) - (m_X * m_X));

    for (k = 0; k < K; k++) {
        pis(k) = arma::accu(ISamps.slice(0) % ISamps.slice(3+k)) / p;

        mus(k) = arma::accu(ISamps.slice(0) % ISamps.slice(3+k) %
        ISamps.slice(2)) / (p*pis(k));

        vars(k) = (arma::accu(ISamps.slice(0) % ISamps.slice(3+k) %
        arma::pow(ISamps.slice(2),2))/(p*pis(k))) - (mus(k) * mus(k));
    }

    if (equalSds){
        sds.fill(sqrt(arma::sum(pis % vars)));
    } else {
        sds = arma::sqrt(vars);
    }

    if (overDispersedY) {
        tau = sqrt(arma::mean((arma::pow(Y, 2) - 2 * (Y % arma::sum((ISamps.slice(0) % ISamps.slice(1) % ISamps.slice(2)), 1)) + arma::sum(ISamps.slice(0) % arma::pow(ISamps.slice(1) % ISamps.slice(2), 2), 1)) / arma::pow(seY, 2)));
        if (tau < 1) {
            tau = 1;
        }
    }

}

// In-place appending of importance samples if MC-EM iteration is rejected
void appendLatentVarPostCube(arma::cube &postSamps, const arma::vec &sum_weights, const int &M, const arma::vec &X, const arma::vec &Y,
    const arma::vec &seX, const arma::vec &seY, const double &m_X,
    const double &lambdaX, arma::vec pis, const arma::vec &mus,
    const arma::vec &sds, int &N_samples)
{
    int i, k, s;
    int p = (int)X.n_elem; // data sample size
    int K = (int)pis.n_elem; // # of components

    int N_add_samples = ceil(N_samples / M);
    arma::span new_ind = arma::span(N_samples,
        N_samples + N_add_samples - 1);
    arma::span row_i;

    // unnormalize importance weights
    postSamps.slice(0).each_col() %= sum_weights;
    // resize postSamps
    postSamps.resize(p, N_samples + N_add_samples, 3+K);

    /* ###### Variables used for sampling mu_Xi ##### */
    double post_var_muXi; // Posterior variance of mu_Xi | Xi, params
    double post_mean_muXi; // Posterior mean of mu_Xi | Xi, params

    /* ##### Variables used for sampling beta_i ##### */
    // length N_samples vector containing sampled components for beta_i | mu_Xi,Yi,params
    arma::uvec k_i(N_add_samples);

    // Posterior var/mean of beta_i | mu_Xi^j, Yi, parameters for j = 1,...,N_samples
    arma::mat post_sigma2_ik_given_muXi(N_add_samples, 1);
    arma::mat post_mu_ik_given_muXi(N_add_samples, 1);

    arma::vec pis_copy(K); // copy of pis at each iteration
    // sample function changes pis

    arma::mat sum_alphas(p, N_add_samples);

    for (i = 0; i < p; i++) {
        pis_copy = pis;
        row_i = arma::span(i);

        /*############### Draw importance samples for mu_Xi #####################*/

        // Calculate posterior mean, variance for mu_Xi | X_i, parameters
        post_var_muXi = 1 / ((1/(seX(i)*seX(i))) + (m_X / (lambdaX*lambdaX)));
        post_mean_muXi = post_var_muXi * ((X(i)/(seX(i)*seX(i))) +
        (m_X/(lambdaX*lambdaX)));

        // Sample mu_Xi from importance density g(mu_Xi) = P(mu_Xi | X_i, Y_i, params)
        postSamps.slice(1)(row_i, new_ind) = (sqrt(post_var_muXi) *
        arma::randn(N_add_samples) + post_mean_muXi).t();
        // postSamps.subcube(i, N_samples, 1, i, N_samples+N_add_samples-1, 1) =
        // (sqrt(post_var_muXi) * arma::randn(N_add_samples) + post_mean_muXi).t();

        /*############### Draw importance samples for beta_i ####################*/

        // Sample components for each MC sample j = 1,...,N_samples
        k_i = RcppArmadillo::sample(arma::regspace<arma::uvec>(0,K-1), N_add_samples,
        true, pis_copy);
        // pis cannot be const for this to work.

        // Calculate posterior mean and variance for \beta_i | \mu_Xi, X_i, Y_i, parameters
        post_sigma2_ik_given_muXi = arma::pow(((1/arma::pow(sds(k_i),2)) +
        (arma::pow(postSamps.slice(1)(row_i, new_ind).t(), 2) /
        (seY(i)*seY(i)))), -1);

        post_mu_ik_given_muXi = post_sigma2_ik_given_muXi %
        ( ((Y(i)*postSamps.slice(1)(row_i, new_ind).t())/(seY(i)*seY(i))) +
        (mus(k_i) / arma::pow(sds(k_i),2)) );

        // Sample from P(\beta_i | \mu_Xi, X_i, Y_i, parameters)
        postSamps.slice(2)(row_i, new_ind) =
        (arma::sqrt(post_sigma2_ik_given_muXi) % arma::randn(N_add_samples) +
        post_mu_ik_given_muXi).t();


        /*####### Calculate P(Z_i = k | \beta_i^j, X_i, Y_i, parameters) ########*/
        for (k = 0; k < K; k++) {
            postSamps.slice(3+k)(row_i, new_ind) = (pis(k) *
            arma::normpdf(postSamps.slice(2)(row_i, new_ind).t(),
            mus(k), sds(k))).t();
        }

        /*############### Calculate importance weights w_i #######################*/
        postSamps.slice(0)(row_i, new_ind) = (pis(k_i) %
        arma::normpdf(Y(i), postSamps.slice(1)(row_i, new_ind).t() %
        mus(k_i), arma::sqrt(arma::pow(sds(k_i) %
        postSamps.slice(1)(row_i, new_ind).t(),2) + (seY(i)*seY(i))))).t();
    }

    // Normalise alphaSamps
    // postSamps(arma::span(), new_ind, arma::span(3, 3+K-1)).each_slice() /=
    // (arma::mat)arma::sum(postSamps(arma::span(), new_ind, arma::span(3, 3+K-1)), 2);
    sum_alphas = arma::sum(postSamps(arma::span(), new_ind,
    arma::span(3, 3+K-1)), 2);
    for (s = 3; s < 3+K; s++) {
        postSamps.slice(s)(arma::span(), new_ind) /= sum_alphas;
    }

    /* ############ ERROR CHECKING ####)########### */
    if (!postSamps.is_finite()) {
        throw exception("Posterior samples have NaN/Inf values.");
    }

    // Update N_samples
    N_samples += N_add_samples;
}

arma::mat computeQMatrix(const int &p, const int &K, const arma::vec &Y,
    const arma::vec &seY, const arma::cube &ISamps,
    const double &m_X, const double &lambdaX,
    const arma::vec &pis, const arma::vec &mus, const arma::vec &sds, const double &tau)
{
    int N_samples = (int)ISamps.n_rows;

    int k;
    arma::mat Q_mat(p, N_samples);

    Q_mat = (((arma::mat)(
    ((arma::mat)arma::pow((((arma::mat)(ISamps.slice(1) % ISamps.slice(2))).each_col() - Y), 2)).each_col() / (-2 * pow(tau*seY, 2)))).each_col()) - log(tau*seY);

    Q_mat += -log(lambdaX) - ((1 / (2 * pow(lambdaX, 2))) * arma::pow((ISamps.slice(1) - m_X), 2));

    for (k = 0; k < K; ++k) {
        Q_mat += ISamps.slice(3+k) % ( log(pis(k)) - log(sds(k)) -
        ((1 / (2 * pow(sds(k), 2))) * arma::pow((ISamps.slice(2) - mus(k)), 2)) );
    }

    // Q_mat = arma::log(arma::normpdf(ISamps.slice(1), m_X, lambdaX));
    // for (k = 0; k < K; ++k) {
    //     Q_mat += ISamps.slice(3+k) % ( log(pis(k)) +
    //     arma::log(arma::normpdf(ISamps.slice(2), mus(k), sds(k))) );
    // }

    if (!Q_mat.is_finite()) {
        throw exception("Complete-data log-likelihood matrix has NaN/Inf values.");
    }

    return Q_mat;
}


 bool MR_EM_diagnostic(
     const arma::cube &ISamps, const int &p, const int &K, const arma::vec &Y, const arma::vec &seY, const double &eps,
     const double &z_alpha, const double &z_gamma, const int &M,
     double &m_X, double &lambdaX, arma::vec &pis, arma::vec &mus,
     arma::vec &sds, double &tau, const double &prevIter_m_X, const double &prevIter_lambdaX,
     const arma::vec &prevIterPis, const arma::vec &prevIterMus,
     const arma::vec &prevIterSds, const double &prevIterTau, const bool &verbose, int &N_samples, const int &N_iters, bool &rejected)
{

    int k, i;
    arma::mat Lambda(p, N_samples);
    double deltaQ;
    double se_Q;
    int N_samples_new;
    bool converged = false;


    Lambda = computeQMatrix(p, K, Y, seY, ISamps, m_X, lambdaX, pis, mus, sds, tau) -
    computeQMatrix(p, K, Y, seY, ISamps, prevIter_m_X, prevIter_lambdaX,
        prevIterPis, prevIterMus, prevIterSds, prevIterTau);


    // // Compute Lambda (log likelihood difference) matrix
    // Lambda = arma::log(arma::normpdf(ISamps.slice(1), m_X, lambdaX) /
    // arma::normpdf(ISamps.slice(1), prevIter_m_X, prevIter_lambdaX));
    //
    // for (k = 0; k < K; k++) {
    //     Lambda = Lambda + ISamps.slice(3+k) %
    //     (log(pis(k)/prevIterPis(k))
    //     + arma::log(arma::normpdf(ISamps.slice(2), mus(k), sds(k)) /
    //     arma::normpdf(ISamps.slice(2), prevIterMus(k), prevIterSds(k))));
    // }


    // Error checking for Lambda matrix
    if (!Lambda.is_finite()) {
        throw exception("Log-likelihood difference matrix has NaN/Inf values.");
    }

    // Compute change in Q function
    deltaQ = arma::accu(ISamps.slice(0) % Lambda);

    // Compute standard error of impt sampling estimate of delta Q
    se_Q = sqrt(arma::sum(arma::pow(arma::sum(ISamps.slice(0) % Lambda, 1), 2) %
    (arma::sum(arma::pow(ISamps.slice(0) % Lambda, 2), 1) /
    arma::pow(arma::sum(ISamps.slice(0) % Lambda, 1), 2) -
    2*(arma::sum(arma::pow(ISamps.slice(0),2) % Lambda, 1) /
    arma::sum(ISamps.slice(0) % Lambda, 1)) +
    arma::sum(arma::pow(ISamps.slice(0), 2), 1))));

    // If M-step estimates accepted
    if ( (deltaQ - z_alpha*se_Q) > 0 ) {

        // If converged
        if ((deltaQ + z_gamma*se_Q) < eps) {
            converged = true;

            // Print output
            if (verbose) {
                Rcout << "Algorithm has converged." << std::endl;
            }
        } else {

            if (arma::is_finite(deltaQ)) {
                // Update Nsamples to obtain better computational efficiency
                N_samples_new = pow((se_Q * (z_alpha + z_gamma)) / deltaQ, 2);
                N_samples = std::max(N_samples, N_samples_new);

                // Print output
                if (verbose) {
                    Rcout << "Update accepted." << std::endl;
                }
            } else {
                N_samples += (N_samples / M);

                if (verbose) {
                    Rcout << "deltaQ is NaN" << std::endl;
                }

            }
            // // Move on to next iteration
            // N_iters++;

        }
    } else {
        rejected = true;

        // reject current iterations M-step update
        m_X = prevIter_m_X;
        lambdaX = prevIter_lambdaX;
        pis = prevIterPis;
        mus = prevIterMus;
        sds = prevIterSds;
        tau = prevIterTau;

        // /* ########## append more MC samples to existing ISamps ########### */
        //
        // int N_add_samples = ceil(N_samples/M);
        // arma::cube new_ISamps = sampleLatentVarPostCube(N_add_samples, X, Y,
        // seX, seY, m_X, lambdaX, pis, mus, sds);
        //
        // // unnormalize existing ISamps
        // ISamps.slice(0).each_col() %= sum_weights;
        //
        // ISamps.resize(p, N_samples + N_add_samples, K);
        //
        // for (i = 0; i < ISamps.n_slices; i++) {
        //
        //     ISamps.slice(i).insert_cols(ISamps.slice(i).n_cols,
        //     new_ISamps.slice(i));
        // }
        //
        // // normalize new ISamps and update sum_weights
        // sum_weights = normalizeImptWeights(ISamps);
        //
        // // increase MC sample size
        // N_samples = N_samples + N_add_samples;

        // Print output
        if (verbose) {
            Rcout << "Update rejected." << std::endl;
        }

    }

    // Print output
    if (verbose) {
        Rcout << "Change in Q function: " << deltaQ << std::endl;
        Rcout << "SE - Change in Q: " << se_Q << std::endl;
    }

    return converged;
}
