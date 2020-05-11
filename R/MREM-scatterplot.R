MREM.scatterplot = function(data, MCEM_fit = NULL,
                            exposure_name = "exposure", outcome_name = "outcome") {

  K = length(MCEM_fit$paramEst$pis)
  data$beta.outcome = data$beta.outcome * sign(data$beta.exposure)
  data$beta.exposure = abs(data$beta.exposure)
  
  if (is.null(MCEM_fit)) {
    p = ggplot(data = data, aes(x = beta.exposure, y = beta.outcome, xmin = (beta.exposure - se.exposure), xmax = (beta.exposure + se.exposure),
                                ymin = (beta.outcome - se.outcome), ymax = (beta.outcome + se.outcome))) + 
      geom_point(size = 1, shape = 1) + geom_errorbar(alpha = 0.3, width = 0) +
      geom_errorbarh(alpha = 0.3, height = 0) + expand_limits(x = 0, y = 0) +
      xlab(paste("SNP association with",exposure_name)) + 
      ylab(paste("SNP association with", outcome_name)) +
      theme_classic(base_size = 15) + geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.4) +
      geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.4) + scale_color_brewer(palette="Set1")
  } else {
    
    # sample from posterior given param estimates
    post_impt_samples <- sampleLatentVarPost(20000, data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome, MCEM_fit$paramEst)
    W <- post_impt_samples$W
    rowSumW <- rowSums(W)
    muX_samps <- post_impt_samples$muX_samps
    beta_samps <- post_impt_samples$beta_samps
    prob_samps <- post_impt_samples$alpha_samps
    
    fitted.pis <- MCEM_fit$paramEst$pis
    fitted.mus <- MCEM_fit$paramEst$mus
    fitted.sds <- MCEM_fit$paramEst$sds
    
    prob_est <- matrix(NA, nrow = nrow(data), ncol = K)
    for (k in 1:K) {
      prob_est[,k] <- rowSums(W * prob_samps[,,k]) / rowSumW
    }
    
    comp_assignments <- as.factor(apply(prob_est,1,which.max))
    names(comp_assignments) <- data$SNP
    
    p = ggplot(data = data, aes(x = beta.exposure, y = beta.outcome, xmin = (beta.exposure - se.exposure), xmax = (beta.exposure + se.exposure),
                                ymin = (beta.outcome - se.outcome), ymax = (beta.outcome + se.outcome), 
                                color = comp_assignments)) +
      geom_point(size = 1, shape = 1) + 
      geom_errorbar(alpha = 0.5, width = 0) +
      geom_errorbarh(alpha = 0.5, height = 0) + 
      expand_limits(x = 0, y = 0) +
      xlab(paste("SNP association with",exposure_name)) + 
      ylab(paste("SNP association with", outcome_name)) +
      theme_classic(base_size = 15) + geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.4) +
      geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.4) + 
      theme(legend.justification = c(0,0), legend.position = c(.05,.005), legend.direction = "horizontal") +
      scale_color_brewer(palette="Set1",name="") + scale_fill_brewer(palette="Set1",name="")
    
    p_build <- ggplot_build(p)
    x_grid = seq(0, p_build$layout$panel_params[[1]]$x.range[2], .001)
    lines_dat <- data.frame(x = rep(x_grid, K))
    lines_dat$y =  as.numeric(sapply(1:K, function(k) fitted.mus[k]*x_grid))
    lines_dat$ymin = as.numeric(sapply(1:K, function(k) (fitted.mus[k]-fitted.sds[k])*x_grid))
    lines_dat$ymax = as.numeric(sapply(1:K, function(k) (fitted.mus[k]+fitted.sds[k])*x_grid))
    lines_dat$k = as.numeric(sapply(1:K, function(k) rep(k, length(x_grid))))
    
    p = p + geom_line(data = lines_dat, 
                      aes(x = x, y = y, color = as.factor(k)), linetype =2, 
                      inherit.aes = FALSE) +
      geom_ribbon(data = lines_dat, 
                  aes(x = x, ymin = ymin, ymax = ymax, fill = as.factor(k)), alpha = 0.4,
                  inherit.aes = FALSE)
  }
  return(p)
}
