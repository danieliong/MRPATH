MRPATH_barplot = function(data, MCEM_fit, ret.snps = FALSE) {
    # Extract estimates
    fitted.pis <- MCEM_fit$paramEst$pis
    fitted.mus <- MCEM_fit$paramEst$mus
    fitted.sds <- MCEM_fit$paramEst$sds
    K = length(fitted.pis)

    # Obtain importance samples
    impt_samples = getImportanceSamples(data, MCEM_fit)

    # Resample importance samples for beta
    beta_samps = sampleBetas(data, impt_samples = impt_samples)

    # Compute cluster membership probabilities
    clustermemb_prob = computeClusterMembProb(data, impt_samples = impt_samples)

    # Compute quantiles of SNP-specific causal effects
    beta_q50 <- apply(beta_samps, 1, quantile, probs = 0.5)
    beta_q5 <- apply(beta_samps, 1, quantile, probs = 0.025)
    beta_q95 <- apply(beta_samps, 1, quantile, probs = 0.975)

    # Compute cluster with highest memb. prob for each SNP
    clusters <- apply(clustermemb_prob,1,which.max)
    names(clusters) <- data$SNP

    # Order SNPs by median causal effects
    ordered.SNPs <- as.character(data$SNP[order(clusters, beta_q50)])
    beta_samps <- beta_samps[ordered.SNPs,]

    # Credible intervals for each beta
    beta.post_interv <- data.frame(cbind(clusters, beta_q5, beta_q50, beta_q95))
    beta.post_interv$SNP <- factor(rownames(beta.post_interv), ordered.SNPs)

    # add cluster memb. prob. to data
    data <- cbind(data, round(clustermemb_prob,4))

    # Credible interval plot for betas
    beta.post_interv_plot <- ggplot(beta.post_interv) +
    aes(x = SNP, y = beta_q50, ymin = beta_q5, ymax = beta_q95,
        color = as.factor(clusters)) +
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


    clustermemb_prob.df = data.frame(clustermemb_prob)
    clustermemb_prob.df$SNP = factor(data$SNP, levels = ordered.SNPs)

    clustermemb_prob.melted = melt(clustermemb_prob.df, id.vars = "SNP")
    clustermemb_prob.melted$SNP = factor(clustermemb_prob.melted$SNP, levels = ordered.SNPs)

    # Cluster membership probability plot
    post_prob.plot = ggplot(clustermemb_prob.melted, aes(x=SNP, y=value, fill=variable)) + geom_bar(stat="identity", alpha = 0.7) +
    scale_fill_brewer(palette="Set1",name="", labels = c("1","2")) + theme_light() +
    theme(legend.position = "none",
          axis.text.y = element_text(size=11),
          axis.text.x = element_blank(),
          # axis.text.x = element_text(size=11),
          axis.title = element_text(size=10),
          plot.margin = margin(t=0.3, l=0.3,r=0.3, unit="cm")) +
    labs(x = "", y = "Cluster membership prob.") + scale_y_continuous(breaks = seq(0, 1, .25), expand=expansion(0,0)) +
    geom_hline(yintercept = 0.5, linetype = "dashed")
    # + coord_flip()

    p = cowplot::plot_grid(beta.post_interv_plot, post_prob.plot, nrow=2, align = "v",
                rel_heights = c(3,1))

    if (!ret.snps) {
      return(p)
    } else {
      return(list(p, ordered.SNPs))
    }
}
