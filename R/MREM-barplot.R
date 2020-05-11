MREM.barplot = function(data, MCEM_fit) {
  
  fitted.pis <- MCEM_fit$paramEst$pis
  fitted.mus <- MCEM_fit$paramEst$mus
  fitted.sds <- MCEM_fit$paramEst$sds
  K = length(fitted.pis)
  
  # Sample from posterior given the above estimates 
  post_impt_samples <- sampleLatentVarPost(100000, data$beta.exposure, data$beta.outcome, 
                                           data$se.exposure, data$se.outcome, 
                                           MCEM_fit$paramEst)
  W <- post_impt_samples$W
  rowSumW <- rowSums(W)
  muX_samps <- post_impt_samples$muX_samps
  beta_samps <- post_impt_samples$beta_samps
  prob_samps <- post_impt_samples$alpha_samps
  
  beta_est <- rowSums(W * beta_samps) / rowSumW
  
  prob_est <- matrix(NA, nrow = nrow(data), ncol = K)
  for (k in 1:K) {
    prob_est[,k] <- rowSums(W * prob_samps[,,k]) / rowSumW
  }
  colnames(prob_est) <- c("p1", "p2")
  
  comp_assignments <- apply(prob_est,1,which.max)
  names(comp_assignments) <- data$SNP
  
  ordered.SNPs <- as.character(data$SNP[order(prob_est[,2])])
  
  beta_resamps <- matrix(NA, nrow = nrow(beta_samps), ncol = ncol(beta_samps))
  for (i in 1:nrow(data)) {
    beta_resamps[i,] <- sample(beta_samps[i,], ncol(beta_resamps), prob = (W[i,]/rowSumW[i]), replace=TRUE)
  }
  rownames(beta_resamps) = data$SNP
  beta_resamps <- beta_resamps[ordered.SNPs,]
  
  beta_q50 <- apply(beta_resamps, 1, quantile, probs = 0.5)
  
  beta_q5 <- apply(beta_resamps, 1, quantile, probs = 0.025)
  beta_q95 <- apply(beta_resamps, 1, quantile, probs = 0.975)
  
  beta.post_interv <- data.frame(cbind(beta_q5, beta_q50, beta_q95))
  beta.post_interv$SNP <- factor(rownames(beta.post_interv), ordered.SNPs)
  
  data <- cbind(data, round(prob_est,2))
  
  beta.post_interv_plot <- ggplot(beta.post_interv) + 
    aes(x = SNP, y = beta_q50, ymin = beta_q5, ymax = beta_q95, 
        color = as.factor(comp_assignments[ordered.SNPs])) + 
    geom_point(shape = 4) + 
    geom_errorbar(show.legend=FALSE) + 
    xlab("") + ylab("95% Posterior Credible Interval") + theme_light() +
    theme(axis.title = element_text(size=10),
          axis.ticks.y = element_blank(), 
          axis.text.y = element_text(size=11), 
          axis.text.x = element_text(size=11, angle=90),
          # legend.position = "none",
          legend.justification = c(0,0), 
          legend.position = c(.85,.005), 
          legend.direction = "horizontal",
          legend.background = element_rect(fill=alpha(0.05)),
          legend.key = element_rect(fill=alpha(0.05)),
          plot.margin = margin(t=0.25,l=0.3,r=0.3, unit="cm")) +
    geom_hline(data = as.data.frame(fitted.mus), 
               aes(color=as.factor(1:K), yintercept = fitted.mus), 
               linetype = "dashed", size = 0.5) + 
    geom_hline(yintercept = 0, color = "black", linetype="dotted") + 
    # scale_x_discrete(position="top") +
    scale_color_brewer(palette="Set1",name="")
  
  prob_est.df = data.frame(prob_est)
  prob_est.df$SNP = factor(data$SNP, levels = ordered.SNPs)
  
  prob_est.melted = melt(prob_est.df, id.vars = "SNP")
  prob_est.melted$SNP = factor(prob_est.melted$SNP, levels = ordered.SNPs)
  
  post_prob.plot = ggplot(prob_est.melted, aes(x=SNP, y=value, fill=variable)) + geom_bar(stat="identity", alpha = 0.7) +
    scale_fill_brewer(palette="Set1",name="", labels = c("1","2")) + theme_light() +
    theme(legend.position = "none",
          axis.text.y = element_text(size=11),
          axis.text.x = element_blank(),
          # axis.text.x = element_text(size=11), 
          axis.title = element_text(size=10),
          plot.margin = margin(t=-0.1, l=0.3,r=0.3, unit="cm")) + 
    labs(x = "", y = "Cluster membership prob.") + scale_y_discrete(limits = c(0,0.5, 1), expand=c(0,0)) + 
    geom_hline(yintercept = 0.5, linetype = "dashed")
  # + coord_flip()
  
  library(cowplot)
  p = plot_grid(beta.post_interv_plot, post_prob.plot, nrow=2, align = "v",
            rel_heights = c(3,1))
  return(p)
}
