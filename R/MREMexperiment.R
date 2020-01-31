MREMexperiment = function(K, X, Y, seX, seY, Nreps = 10, init_seed = 8686, verbose = FALSE) {
  
  # initVals.list <- list()
  MCEM_fit.list <- list()
  Q.vec = c()
  
  init_m_X = mean(X)
  init_lambdaX = sd(X)
  initPis = rep((1/K), K)
  
  quantiles = c(0.025, sapply(1:K, function(k) 0.025 + (.95 / ((K - k)+1))) )
  ratio.quantiles = quantile(((seX/seY)*(Y/X)), quantiles)
  initSds <- sapply(1:K, function(k) abs(ratio.quantiles[k+1] - ratio.quantiles[k])/2)
  
  for (i in 1:Nreps) {
    
    # Set random initial values for mus, sds 
    set.seed(init_seed*i)
    initMus = c()
    for (k in 1:K) {
      initMus[k] = runif(1,  ratio.quantiles[k], ratio.quantiles[k+1])
    }
    initVals <- list("pis" = initPis, "mus" = initMus, "sds" = initSds, 
                     "m_X" = init_m_X, "lambdaX" = init_lambdaX)
    
    # Run algorithm
    MCEM_fit.list[[i]] = MR_EM(K, initVals, X, Y, seX, seY)
    Q = MCEM_fit.list[[i]]$convergenceInfo$completeDataLogLik
    
    if (verbose) {
      print(paste("Run #",i," / Q = ",Q,sep=""))
    }
    
    
    # Save Q value
    Q.vec = c(Q.vec, Q)
  }
  
  bestModel = MCEM_fit.list[[which.max(Q.vec)]]
  return(bestModel)
}

