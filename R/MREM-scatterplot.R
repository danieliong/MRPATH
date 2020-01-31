MREM.scatterplot = function(X, Y, seX, seY, MCEM_fit = NULL) {
  
  Y = Y * sign(X)
  X = abs(X)
  data = data.frame(X = X, Y = Y, seX = seX, seY = seY)
  
  p = ggplot(data = data, aes(x = X, y = Y, xmin = (X - seX), xmax = (X + seX),
                              ymin = (Y - seY), ymax = (Y + seY))) + 
    geom_point(size = 1, shape = 1) + geom_errorbar(alpha = 0.3, width = 0) +
    geom_errorbarh(alpha = 0.3, height = 0) + expand_limits(x = 0, y = 0) +
    xlab("SNP association with exposure") + ylab("SNP association with outcome") +
    theme_classic(base_size = 15) + geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.4)
  
  fitted.pis <- MCEM_fit$paramEst$pis
  fitted.mus <- MCEM_fit$paramEst$mus
  fitted.sds <- MCEM_fit$paramEst$sds
  
  # TODO: Add MCEM_fit into plot
  if (!is.null(MCEM_fit)) {
    p_build <- ggplot_build(p)
    for (k in 1:K) {
      x <- seq(0, p_build$layout$panel_params[[1]]$x.range[2], .001)
      y <- fitted.mus[k] * x
      ymin <- (fitted.mus[k]-fitted.sds[k]) * x
      ymax <- (fitted.mus[k]+fitted.sds[k]) * x
      ribbon.dat <- data.frame(x = x, y = y, ymin = ymin, ymax = ymax)
      p <- p + geom_abline(intercept = 0, slope = fitted.mus[k], aes(text = "test")) + geom_ribbon(data = ribbon.dat,
                                                                                          aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.6*fitted.pis[k], inherit.aes = FALSE) + geom_abline(intercept = 0, slope = fitted.mus[k]-fitted.sds[k], alpha = 0.2) + geom_abline(intercept = 0, slope = fitted.mus[k]+fitted.sds[k], alpha = 0.2)
    }
  }
  
  
  return(p)
}
