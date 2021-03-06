loo.cv <- function(X, Y, seX, seY, params, Nsamples = 20000, loo_method = "tis", ...) {
  
  N = length(X)
  impt_samps <- sampleLatentVarPost(Nsamples, X, Y, seX, seY, params)
  
  # muX_samps <- impt_samps$muX_samps
  # beta_samps <- impt_samps$beta_samps
  # W <- impt_samps$W
  # rowSumW = rowSums(impt_samps$W)
  
  muX_resamps = sapply(1:N, function(i){
    sample(impt_samps$muX_samps[i,], Nsamples, prob = impt_samps$W[i,], replace=TRUE)
  })
  
  beta_resamps = sapply(1:N, function(i){
    sample(impt_samps$beta_samps[i,], Nsamples, prob = impt_samps$W[i,], replace=TRUE)
  })
  
  
  log_lik <- sapply(1:N, 
                    function(i) {
                      dnorm(X[i], muX_resamps[,i], seX[i], log=TRUE) + 
                        dnorm(Y[i], beta_resamps[,i] * muX_resamps[,i], seY[i], log=TRUE)
                      })
  
  loo::loo(x = log_lik, is_method = loo_method, ...)
}


selectModel = function(X, Y, seX, seY, K_range = 1:3, Nreps = 20, 
                       verbose=FALSE, loo_method = "tis") {
  
  loo_list = list()
  for (K in K_range) {
    ## Choose initial values that maximizes likelihood
    initVals = optimizeInitVals(K, X, Y, seX, seY, Nreps = Nreps, verbose=verbose)
    
    ## Run MC-EM with optimized initial values
    MCEM_fit = MR_EM(K, initVals, X, Y, seX, seY, saveTraj=FALSE, computeSE=FALSE)
    
    loo_list[[K]] = loo.cv(X, Y, seX, seY, MCEM_fit$paramEst, loo_method = loo_method)
  }
  
  loo::loo_compare(x = loo_list)
}


# MR.waic = function(X, Y, seX, seY, params, Nsamples = 20000, ...) {
#   
#   impt_samps <- sampleLatentVarPost(Nsamples, X, Y, seX, seY, params)
#   muX_samps <- impt_samps$muX_samps
#   beta_samps <- impt_samps$beta_samps
#   W <- impt_samps$W
#   
#   
#   log_lik <- sapply(1:nrow(muX_samps), 
#                     function(i) {
#                       dnorm(X[i], muX_samps[i,], seX[i], log=TRUE) + 
#                         dnorm(Y[i], beta_samps[i,] * muX_samps[i,], seY[i], log=TRUE)
#                     })
#   
#   log_lik_adjusted = log_lik - t(log(W))
#   
#   waic(log_lik_adjusted)
# }

