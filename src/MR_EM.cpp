#include <RcppArmadillo.h>
#include <RcppDist.h>
#include <math.h>
#include "utils.h"
#include "observedInfo.h"

using namespace Rcpp;


//[[Rcpp::export]]
List MR_EM(int K, const List &initVals, const arma::vec &X, const arma::vec &Y,
    const arma::vec &seX, const arma::vec &seY, bool overDispersedY = false,
    bool equalSds = false, const int &Nstart_MC = 500, int M = 4,
    int max_Nsamples = 500000, int min_iters = 2, int max_iters = 100,
    double alpha = 0.05, double gamma = 0.05, double eps = 0.005,
    bool verbose=false, bool saveTraj = false, bool computeSE = true, bool noMixtureVar = false)
{

  // Initialize starting values for parameters
  double m_X = initVals["m_X"];
  double lambdaX = initVals["lambdaX"];
  arma::vec pis = initVals["pis"];
  arma::vec mus = initVals["mus"];
  arma::vec sds = initVals["sds"];

  int p = (int)X.n_elem; // data sample size

  // z score for alpha, gamma
  double z_alpha = R::qnorm5(alpha, 0, 1, false, false);
  double z_gamma = R::qnorm5(gamma, 0, 1, false, false);

  // Overdispersion parameter for seY
  double tau = 1;

  bool rejected = false;
  bool converged = false;

  int N_iters = 1;

  /* #################### Error Checking for input ######################## */

  // check if pis, mus, sds are same size
  if ((mus.n_elem != K) || (sds.n_elem != K)) {
    throw exception("pis, mus, sds, are not all of size K!");
  }

  // check if X, Y, seX, seY are the same size
  if ((Y.n_elem != p) || (seX.n_elem != p) || (seY.n_elem != p)) {
    throw exception("X, Y, seX, seY are not all of size K!");
  }

  /* ###################################################################### */

  // Initialize starting MC sample size
  int N_samples = Nstart_MC;
  if (verbose) {
    Rcout << "Starting MC sample size: " << N_samples << std::endl;
  }


  // Matrices to store M-step updates at each iteration
  arma::vec m_X_updates;
  arma::vec lambdaX_updates;
  arma::mat pis_updates;
  arma::mat mus_updates;
  arma::mat sds_updates;

  if (saveTraj){
    m_X_updates.set_size(1);
    lambdaX_updates.set_size(1);
    pis_updates.set_size(K, 1);
    mus_updates.set_size(K, 1);
    sds_updates.set_size(K, 1);

    // Save initial values as first element in trace
    m_X_updates(0) = m_X;
    lambdaX_updates(0) = lambdaX;
    pis_updates.col(0) = pis;
    mus_updates.col(0) = mus;
    sds_updates.col(0) = sds;
  }

  arma::cube ISamps(p, Nstart_MC, K);
  arma::vec sum_weights(p);

  // Parameter updates from previous iteration
  double prevIter_m_X;
  double prevIter_lambdaX;
  arma::vec prevIterPis;
  arma::vec prevIterMus;
  arma::vec prevIterSds;
  double prevIterTau;

  arma::mat observedInfoMtx((2+(3*K) - 1), (2 + (3*K) - 1));
  arma::vec se((2+(3*K)-1));

  while ((N_iters <= min_iters) ||
  (!converged &&
  (N_iters < max_iters) && (N_samples < max_Nsamples)))
  {

    if (verbose){
      Rcout << "#################### Iter #: " << N_iters <<
      " ####################" << std::endl;
    }

    /* ###############################################################
                            E Step
    ############################################################### */

    if (!rejected) {
        ISamps.set_size(p, N_samples, 3+K);
        ISamps = sampleLatentVarPostCube(N_samples, X, Y, seX, tau*seY, m_X, lambdaX, pis, mus, sds);
    } else {
        appendLatentVarPostCube(ISamps, sum_weights, M, X, Y, seX,
            tau*seY, m_X, lambdaX, pis, mus, sds, N_samples);
    }

    // normalize importance weights and save to sum_weights
    sum_weights = normalizeImptWeights(ISamps);

    /* ########################################################################
                              M Step
    ######################################################################### */

     // Deep copy previous iteration estimates
     prevIter_m_X = m_X;
     prevIter_lambdaX = lambdaX;
     prevIterPis = pis;
     prevIterMus = mus;
     prevIterSds = sds;
     prevIterTau = tau;

     MR_EM_Mstep(p, K, Y, seY, equalSds, overDispersedY, ISamps, m_X, lambdaX, pis, mus, sds, tau);

     if (verbose) {
         // Print curent iteration parameter updates
        Rcout << "m_X: " << m_X << std::endl;
        Rcout << "lambdaX: " << lambdaX << std::endl;
        Rcout << "pis: " << std::endl;
        Rcout << pis << std::endl;
        Rcout << "mus: " << std::endl;
        Rcout << mus << std::endl;
        Rcout << "sds: " << std::endl;
        Rcout << sds << std::endl;
        Rcout << "tau: " << tau << std::endl;
    }


     /* ###############################################################
                     Diagnostics
    ############################################################### */

    converged = MR_EM_diagnostic(ISamps, p, K, Y, seY, eps, z_alpha,
        z_gamma, M, m_X, lambdaX, pis, mus, sds, tau,
        prevIter_m_X, prevIter_lambdaX, prevIterPis, prevIterMus, prevIterSds, prevIterTau,
        verbose, N_samples, N_iters, rejected);


    if (saveTraj && !rejected){
        // save updates for plotting trace
        m_X_updates.insert_rows(m_X_updates.n_elem, m_X);
        lambdaX_updates.insert_rows(lambdaX_updates.n_elem, lambdaX);
        pis_updates.insert_cols(pis_updates.n_cols, pis);
        mus_updates.insert_cols(mus_updates.n_cols, mus);
        sds_updates.insert_cols(sds_updates.n_cols, sds);
    }


    if (!rejected) {
        N_iters++;
    }

    if (verbose && (N_samples >= max_Nsamples)) {
      Rcout << "N_samples has exceeded max_Nsamples." << std::endl;
      converged = true;
    }

    if (verbose && (N_iters == max_iters)) {
      Rcout << "Max. number of iterations exceeded." << std::endl;
      converged = true;
    }

    }

    arma::mat Q_mat = computeQMatrix(p, K, Y, seY, ISamps, m_X, lambdaX, pis, mus, sds, tau);
    double completeDataLogLik = arma::accu(ISamps.slice(0) % Q_mat);

    /* ###############################################################
                    Standard Error
   ############################################################### */

   if (computeSE) {

       observedInfoMtx = computeObservedInfoMatrix(X, Y, seX, tau*seY, m_X, lambdaX,
           pis, mus, sds, N_samples
       );

       se = arma::sqrt(((arma::mat)observedInfoMtx.i()).diag());
   }

  List paramEst = List::create(
        Named("m_X") = m_X,
        _["lambdaX"] = lambdaX,
        _["pis"] = pis,
        _["mus"] = mus,
        _["sds"] = sds,
        _["tau"] = tau
    );
  List convergenceInfo = List::create(
        Named("Niters") = N_iters,
        _["N_MC_end"] = N_samples,
        _["completeDataLogLik"] = completeDataLogLik
    );

  List paramTraj = R_NilValue;


 if (saveTraj) {
     paramTraj = List::create(
         Named("m_X") = m_X_updates,
         _["lambdaX"] = lambdaX_updates,
         _["pis"] = pis_updates,
         _["mus"] = mus_updates,
         _["sds"] = sds_updates
     );
 }

 // TODO: Add cluster probabilities

  return List::create(
      Named("paramEst") = paramEst,
      _["standardErrors"] = se,
      _["convergenceInfo"] = convergenceInfo,
      _["paramTraj"] = paramTraj
    );


}
