selectModel_stan = function(X, Y, seX, seY, K_range = 1:3, Nreps = 15, 
                                     verbose=FALSE) {
  
  MR_stancode = "
        data {
      // sample size
      int<lower=0> p;
      
      // # clusters
      int<lower=0> K;
      
      // data
      real<lower=0> X[p];
      real<lower=0> seX[p];
      real Y[p];
      real<lower=0> seY[p];
      
      // param estimates
      real m_x;
      real<lower=0> lambda_x;
      vector[K] pis;
      vector[K] mus;
      vector[K] sds;
    }
    
    parameters {
      real<lower=0> muX[p];
      real beta[p];
    }
    
    model {
      
      // prior for muX
      muX ~ normal(m_x, lambda_x);
      
      for (i in 1:p) {
        
        real cluster_logliks[K];
        
        // prior for beta[i]
        for (k in 1:K) {
          cluster_logliks[k] = log(pis[k]) +  normal_lpdf(beta[i] | mus[k], sds[k]);
        }
        target += log_sum_exp(cluster_logliks);
        
        // likelihood for X[i]
        target += normal_lpdf(X[i] | muX[i], seX[i]);
        
        // likelihood for Y[i]
        target += normal_lpdf(Y[i] | muX[i]*beta[i], seY[i]);
      }
      
    }
    
    generated quantities {
      vector[p] log_lik;
      for (i in 1:p) {
        log_lik[i] = normal_lpdf(X[i] | muX[i], seX[i]) + 
        normal_lpdf(Y[i] | beta[i]*muX[i], seY[i]);
      }
    }
  "
  MR_stanmod = stan_model(model_code = MR_stancode)
  
  
  loo_list = list()
  for (K in K_range) {
    ## Choose initial values that maximizes likelihood
    initVals = optimizeInitVals(K, X, Y, seX, seY, Nreps = Nreps, verbose=verbose)
    
    ## Run MC-EM with optimized initial values
    MCEM_fit = MR_EM(K, initVals, X, Y, seX, seY, saveTraj=FALSE, computeSE=FALSE)
    
    ## Sampling with stan
    data = list(p = length(X),
                K = K,
                X = abs(X),
                seX = seX,
                Y = sign(X)*Y,
                seY = seY,
                m_x = MCEM_fit$paramEst$m_X,
                lambda_x = MCEM_fit$paramEst$lambdaX,
                pis = array(MCEM_fit$paramEst$pis, dim = K),
                mus = array(MCEM_fit$paramEst$mus, dim = K),
                sds = array(MCEM_fit$paramEst$sds, dim = K)
    )
    MR_stan = sampling(MR_stanmod, data = data)
    log_lik = extract_log_lik(MR_stan, merge_chains = FALSE)
    r_eff <- relative_eff(exp(log_lik)) 
    
    loo_list[[K]] = loo(log_lik, r_eff = r_eff)
  }
  
  loo::loo_compare(x = loo_list)
}
