MRPATH_optimizeInitVals = function(K, data, Nreps = 10, verbose = FALSE,
    altModel = FALSE, init_seed = 8686, ...) {

  initVals.list <- list()
  Q.vec = c()

  init_m_X = mean(data$beta.exposure)
  init_lambdaX = sd(data$beta.exposure)
  initPis = rep((1/K), K)

  quantiles = c(0.025, sapply(1:K, function(k) 0.025 + (.95 / ((K - k)+1))) )
  ratio.quantiles = quantile(((data$se.exposure/data$se.outcome)*(data$beta.outcome/data$beta.exposure)), quantiles)

  fit = NULL

  for (i in 1:Nreps) {

    # Set random initial values for mus
    set.seed(init_seed*i)
    initMus = c()
    for (k in 1:K) {
      initMus[k] = runif(1,  ratio.quantiles[k], ratio.quantiles[k+1])
    }
    initVals.list[[i]] <- list("m_X" = init_m_X,
                               "lambdaX" = init_lambdaX,
                               "pis" = initPis,
                               "mus" = initMus)

    if (altModel) {
      fit = MR_PATHalt(data, initVals.list[[i]], ...)
      Q = fit$completeDataLogLik
    } else {
      initVals.list[[i]]$sds = sapply(1:K, function(k)
        abs(ratio.quantiles[k+1] - ratio.quantiles[k])/2)

      fit = MR_PATH(K, data, initVals.list[[i]], computeSE = FALSE, ...)
      Q = fit$convergenceInfo$completeDataLogLik
    }

    # Save Q value
    Q.vec = c(Q.vec, Q)

    if (verbose) {
      print(paste("Run #",i," / Q = ",Q,sep=""))
    }

  }

  optimalInitVals = initVals.list[[which.max(Q.vec)]]
  res = list("fit" = fit,
             "initVals" = optimalInitVals)
  return(res)
}
