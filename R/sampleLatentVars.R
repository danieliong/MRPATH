sampleBetas = function(data, Nsamples = 50000, impt_samples = NULL,
    MCEM_fit = NULL)
{
    # If both impt_samples and MCEM_fit are missing, raise error
    if (is.null(MCEM_fit) & is.null(impt_samples)) {
        stop("At least one of importance samples/MCEM_fit must be provided.")
    } else if (is.null(impt_samples)) { # if importance samples are not provided
        impt_samples = getImportanceSamples(data, MCEM_fit, Nsamples = 3*Nsamples)
    }

    # Resample beta importance samples to get beta samples
    beta_samps = matrix(NA, nrow = nrow(impt_samples$beta), ncol = ncol(impt_samples$beta))
    for (i in 1:nrow(data)) {
      beta_samps[i,] <- sample(impt_samples$beta[i,], Nsamples, prob = (impt_samples$W[i,]/rowSums(impt_samples$W)[i]), replace=TRUE)
    }

    # Set row names to be SNP names if data provides it
    if (!is.null(data$SNP)) {
        rownames(beta_samps) = data$SNP
    }

    return(beta_samps)
}

computeClusterMembProb = function(data, Nsamples = 50000,
    impt_samples = NULL, MCEM_fit = NULL)
{
    # If both impt_samples and MCEM_fit are missing, raise error
    if (is.null(MCEM_fit) & is.null(impt_samples)) {
        stop("At least one of importance samples/MCEM_fit must be provided.")
    } else if (is.null(impt_samples)) { # if importance samples are not provided
        impt_samples = getImportanceSamples(data, MCEM_fit, Nsamples = Nsamples)
    }

    K = dim(impt_samples$alpha)[3]

    # Compute cluster membership probabilities with importance samples
    clustermemb_prob = matrix(NA, nrow = nrow(data), ncol = K)
    for (k in 1:K) {
      clustermemb_prob[,k] <- rowSums(impt_samples$W * impt_samples$alpha[,,k]) / rowSums(impt_samples$W)
    }

    # Set row names to be SNP names if data provides it
    if (!is.null(data$SNP)) {
        rownames(clustermemb_prob) = data$SNP
    }

    return(clustermemb_prob)
}
