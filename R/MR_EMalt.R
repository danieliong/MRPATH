

MR_EMalt = function(data, initVals, eps, max_iters = 20, verbose=FALSE)
{
    
    X = data$beta.exposure
    Y = data$beta.outcome
    seX = data$se.exposure
    seY = data$se.outcome
    
    m_X = initVals$m_X
    lambdaX = initVals$lambdaX
    pis = initVals$pis
    mus = initVals$mus
    K = length(pis)

    Q = 0
    Qnew = 0

    for (t in 1:max_iters) {

        if (verbose) {
            print("####################")
            print(paste("Iteration #",t,sep=""))
        }


        ########### E step #############
        # A = ((1/seY) %o% mus)^2 + (1/lambdaX^2) + (1/seX^2)
        #
        # ### Compute pi_ik
        # log_pis_ik = sapply(1:K, function(k)
        #     log(pis[k]) - (0.5 * log(A[,k])) +
        #     (0.5*( ((Y*mus[k]/(seY^2)) + (m_X/(lambdaX^2)) + (X/(seX^2)) )^2 / A[,k] ))
        # )
        #
        # # Subtract max for numerical stability
        # max_log_pis_ik = apply(log_pis_ik, 1, max)
        # log_pis_ik = (log_pis_ik - max_log_pis_ik) - log(rowSums(exp(log_pis_ik - max_log_pis_ik)))
        # pis_ik = exp(log_pis_ik)

        ### Compute lambda2_ik and m_ik
        lambda2_ik = 1 / ( ((1/seY) %o% mus)^2 + (1/seX)^2 + (1/lambdaX)^2 )
        m_ik = lambda2_ik * (  ((Y %o% mus)/(seY^2)) + (X/(seX^2)) + (m_X/(lambdaX^2)) )

        log_pis_ik = sapply(1:K, function(k) log(pis[k]) + (0.5 * log(lambda2_ik[,k])) + (0.5*(m_ik[,k]^2 / lambda2_ik[,k])) )

        # Subtract max for numerical stability
        max_log_pis_ik = apply(log_pis_ik, 1, max)
        log_pis_ik = (log_pis_ik - max_log_pis_ik) - log(rowSums(exp(log_pis_ik - max_log_pis_ik)))
        pis_ik = exp(log_pis_ik)

        E_Z_theta = pis_ik * m_ik
        E_Z_theta2 = pis_ik * (lambda2_ik + (m_ik^2))
        E_theta = rowSums(E_Z_theta)
        E_theta2 = rowSums(E_Z_theta2)

        ########### M step #############
        m_X = mean(E_theta)
        lambda_X = sqrt(mean( E_theta2 - (2*m_X*E_theta) + (m_X^2)))

        pis = colMeans(pis_ik)
        mus = colSums( (Y/(seY^2)) * E_Z_theta ) / colSums(E_Z_theta2/(seY^2))

        if (verbose) {
            print(paste("m_X:",m_X))
            print(paste("lambda_X:",lambda_X))
            print(paste("pis:",pis))
            print(paste("mus:",mus))
        }
    }
    
    # Compute pis_ik at convergence
    lambda2_ik = 1 / ( ((1/seY) %o% mus)^2 + (1/seX)^2 + (1/lambdaX)^2 )
    m_ik = lambda2_ik * (  ((Y %o% mus)/(seY^2)) + (X/(seX^2)) + (m_X/(lambdaX^2)) )
    log_pis_ik = sapply(1:K, function(k) log(pis[k]) + (0.5 * log(lambda2_ik[,k])) + (0.5*(m_ik[,k]^2 / lambda2_ik[,k])) )
    # Subtract max for numerical stability
    max_log_pis_ik = apply(log_pis_ik, 1, max)
    log_pis_ik = (log_pis_ik - max_log_pis_ik) - log(rowSums(exp(log_pis_ik - max_log_pis_ik)))
    pis_ik = exp(log_pis_ik)
    
    rownames(pis_ik) = data$SNP
    
    res = list("paramEst" = list("m_X" = m_X,
                          "lambdaX" = lambda_X,
                          "pis" = pis,
                          "mus" = mus),
               "clusterMembProb" = pis_ik
               )
    

    return(res)
}
