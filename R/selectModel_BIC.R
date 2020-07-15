MRPATH_selectModel = function(data, K_range = 1:3,
    Nreps = 20, altModel = FALSE, verbose=FALSE, ...) {

    Q_vec = rep(NA, max(K_range))
    BIC_vec = rep(NA, max(K_range))
    fit <- list()

    for (K in K_range) {

        fit[[K]] = MRPATH_optimizeInitVals(data = data, K = K, altModel = altModel, ...)$fit
        if (altModel) {
            Q_vec[K] = fit[[K]]$completeDataLogLik
        } else {
            Q_vec[K] = fit[[K]]$convergenceInfo$completeDataLogLik
        }

        if (is.null(data)) {
            p = length(X)
        } else {
            p = nrow(data)
        }

        BIC_vec[K] =  (-2*Q_vec[K]) + (3*K * log(p))

        if (verbose) {
            print(paste("K = ",K,": BIC = ",BIC_vec[K],sep=""))
        }
    }

    bestK = which.min(BIC_vec)
    bestFit = fit[[bestK]]

    return(list(bestK = bestK,
                bestFit = bestFit,
                Q = Q_vec,
                BIC = BIC_vec
                )
           )
}
