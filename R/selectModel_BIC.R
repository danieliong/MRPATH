selectModel_BIC = function(X, Y, seX, seY, K_range = 1:3, Nreps = 20, 
                       verbose=FALSE) {
  
  p = length(X)
  BIC_vec = rep(NA, max(K_range))
  
  for (K in K_range) {
    ## Choose initial values that maximizes likelihood
    initVals = optimizeInitVals(K, X, Y, seX, seY, Nreps = Nreps, verbose=verbose)
    
    ## Run MC-EM with optimized initial values
    MCEM_fit = MR_EM(K, initVals, X, Y, seX, seY, saveTraj=FALSE, computeSE=FALSE)
    
    BIC_vec[K] =  MCEM_fit$convergenceInfo$completeDataLogLik - (3*K * log(p))
  }
  
  K_opt = which.max(BIC_vec)
  return(K_opt)  
}