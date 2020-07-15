MR_PATHalt = function(data, initVals, eps = 1e-4, max_iters = 50, verbose=FALSE)
{

    X = data$beta.exposure
    Y = data$beta.outcome
    seX = data$se.exposure
    seY = data$se.outcome

    m_X = initVals$m_X
    lambdaX = initVals$lambdaX
    pis = initVals$pis
    mus = initVals$mus

    p = length(X)
    K = length(pis)

    Q = -Inf

    for (t in 1:max_iters) {

        if (verbose) {
            print("####################")
            print(paste("Iteration #",t,sep=""))
            print("####################")
        }


        ########### E step #############

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

        # Compute Q for first iteration
        if (t == 1) {
            Q = computeQ(Y, seY, m_X, lambdaX, pis, mus,
                         pis_ik, E_Z_theta, E_Z_theta2, E_theta, E_theta2)
        }

        ########### M step #############
        m_X = mean(E_theta)
        lambdaX = sqrt(mean( E_theta2 - (2*m_X*E_theta) + (m_X^2)))

        pis = colMeans(pis_ik)
        mus = colSums( (Y/(seY^2)) * E_Z_theta ) / colSums(E_Z_theta2/(seY^2))

        if (verbose) {
            print(paste("m_X:",m_X))
            print(paste("lambdaX:",lambdaX))
            print(paste("pis:",pis))
            print(paste("mus:",mus))
        }


        # Compute new Q function
        Qnew = computeQ(Y, seY, m_X, lambdaX, pis, mus,
                        pis_ik, E_Z_theta, E_Z_theta2, E_theta, E_theta2)

        if (verbose) {
            print(paste("Change in Q:",(Qnew - Q)))
        }

        if ((Qnew - Q) < eps) {
            if (verbose) {
                print("EM Algorithm has converged.")
            }
            break
        } else {
            Q = Qnew
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

    # Compute Q at convergence


    res = list("paramEst" = list("m_X" = m_X,
                          "lambdaX" = lambdaX,
                          "pis" = pis,
                          "mus" = mus),
               "clusterMembProb" = pis_ik,
               "completeDataLogLik" = Q
               )


    return(res)
}

computeQ = function(Y, seY, m_X, lambdaX, pis, mus,
                    pis_ik, E_Z_theta, E_Z_theta2, E_theta, E_theta2)
{
    p = length(Y)
    N_k = colSums(pis_ik)
    Q = - p * log(lambdaX) + ((m_X / (lambdaX^2)) * sum(E_theta)) -
        (p/2)*((m_X/lambdaX)^2) +
        sum(sapply(1:K, function(k)
            (N_k[k]*log(pis[k])) + (mus[k]* sum((Y/(seY^2)) * E_Z_theta[,k])) -
                (((mus[k]^2)/2)*sum(E_Z_theta2[,k])))
            )

    return(Q)
}
