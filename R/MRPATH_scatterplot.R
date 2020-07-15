MRPATH_scatterplot = function(data, MCEM_fit = NULL,
                            exposure_name = "exposure", outcome_name = "outcome", overDispersedY = FALSE) {

  K = length(MCEM_fit$paramEst$pis)
  data$beta.outcome = data$beta.outcome * sign(data$beta.exposure)
  data$beta.exposure = abs(data$beta.exposure)

  if (overDispersedY) {
    data$tau = MCEM_fit$paramEst$tau
  } else {
    data$tau = 1
  }

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
    fitted.pis <- MCEM_fit$paramEst$pis
    fitted.mus <- MCEM_fit$paramEst$mus
    fitted.sds <- MCEM_fit$paramEst$sds

    clustermemb_prob = computeClusterMembProb(data, MCEM_fit = MCEM_fit)

    clusters <- as.factor(apply(clustermemb_prob,1,which.max))
    names(clusters) <- data$SNP

    p = ggplot(data = data, aes(x = beta.exposure, y = beta.outcome, xmin = (beta.exposure - se.exposure), xmax = (beta.exposure + se.exposure),
                                ymin = (beta.outcome - (tau*se.outcome)), ymax = (beta.outcome + (tau*se.outcome)),
                                color = clusters)) +
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
                      aes(x = x, y = y, color = as.factor(k)), linetype = 1, size = .8,
                      inherit.aes = FALSE) +
      geom_ribbon(data = lines_dat,
                  aes(x = x, ymin = ymin, ymax = ymax, fill = as.factor(k)), alpha = 0.4,
                  inherit.aes = FALSE)
  }
  return(p)
}


MREMalt.scatterplot = function(data, EM_fit = NULL,
                            exposure_name = "exposure", outcome_name = "outcome", interactive = FALSE) {

    data$tau = 1
    fitted.pis = EM_fit$paramEst$pis
    fitted.mus = EM_fit$paramEst$mus
    fitted.sds = EM_fit$paramEst$sds
    clusterMembProb = EM_fit$clusterMembProb
    K = length(fitted.pis)
    p = nrow(data)

    data$SNPtext = sapply(1:p, function(i) paste(round(clusterMembProb[i,],digits=3), collapse = ", "))

    clusters = as.factor(apply(EM_fit$clusterMembProb, 1, which.max))
    names(clusters) = data$SNP

    colors = RColorBrewer::brewer.pal(5, "Set1")

    p = ggplot(data = data, aes(x = beta.exposure, y = beta.outcome, xmin = (beta.exposure - se.exposure), xmax = (beta.exposure + se.exposure),
                                ymin = (beta.outcome - (tau*se.outcome)), ymax = (beta.outcome + (tau*se.outcome)),
                                color = clusters)) +
      geom_point(size = 1, shape = 1,
          aes(text = SNPtext)) +
      geom_errorbar(alpha = 0.5, width = 0) +
      geom_errorbarh(alpha = 0.5, height = 0) +
      expand_limits(x = 0, y = 0) +
      xlab(paste("SNP association with",exposure_name)) +
      ylab(paste("SNP association with", outcome_name)) +
      theme_classic(base_size = 15) + geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.4) +
      geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.4) +
      scale_color_brewer(palette="Set1",name="") + scale_fill_brewer(palette="Set1",name="")

      p_build <- ggplot_build(p)
      x_grid = seq(0, p_build$layout$panel_params[[1]]$x.range[2], .001)
      lines_dat = data.frame(x = rep(x_grid, K))
      lines_dat$y =  as.numeric(sapply(1:K, function(k) fitted.mus[k]*x_grid))
      lines_dat$k = as.factor(sapply(1:K, function(k) rep(k, length(x_grid))))


      if (interactive) {
          p = ggplotly(p, tooltip = c("text"))
          p = add_lines(p, x = ~x, y = ~y, color = ~k, colors = "Set1",
              data = lines_dat, inherit = TRUE, hoverinfo = "none")
      } else {
          p = p + theme(legend.justification = c(0,0), legend.position = c(.05,.005),
                        legend.direction = "horizontal") +
            geom_line(data = lines_dat,
                      aes(x = x, y = y, color = as.factor(k)),
                      linetype = 1, size = .8, inherit.aes = FALSE)
      }


    return(p)
}
