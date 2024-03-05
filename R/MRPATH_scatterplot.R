MRPATH_scatterplot <- function(data, MCEM_fit = NULL,
                               exposure_name = "exposure",
                               outcome_name = "outcome",
                               overDispersedY = FALSE,
                               interactive = TRUE) {
  K <- length(MCEM_fit$paramEst$pis)
  n <- nrow(data)
  data$beta.outcome <- data$beta.outcome * sign(data$beta.exposure)
  data$beta.exposure <- abs(data$beta.exposure)

  if (overDispersedY) {
    data$tau <- MCEM_fit$paramEst$tau
  } else {
    data$tau <- 1
  }

  if (is.null(MCEM_fit)) {
    p <- ggplot(
      data = data,
      aes(
        x = beta.exposure, y = beta.outcome,
        xmin = (beta.exposure - se.exposure),
        xmax = (beta.exposure + se.exposure),
        ymin = (beta.outcome - se.outcome),
        ymax = (beta.outcome + se.outcome)
      )
    ) +
      geom_point(size = 1, shape = 1) +
      geom_errorbar(alpha = 0.3, width = 0) +
      geom_errorbarh(alpha = 0.3, height = 0) +
      expand_limits(x = 0, y = 0) +
      xlab(paste("SNP association with", exposure_name)) +
      ylab(paste("SNP association with", outcome_name)) +
      theme_classic(base_size = 15) +
      geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.4) +
      geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.4) +
      scale_color_brewer(palette = "Set1", name = "")
  } else {
    get_SNP_text <- function(i) {
      data_i <- data[i, ]
      prob_i <- round(clustermemb_prob[i, ], 2)

      prob_text <- paste(
        sapply(1:K, function(k) paste("Cluster #", k, " Prob: ", prob_i[k], sep = "")),
        collapse = "<br>"
      )
      text <- paste("SNP: ", data_i$SNP, "<br>", prob_text, sep = "")
    }

    fitted.pis <- MCEM_fit$paramEst$pis
    fitted.mus <- MCEM_fit$paramEst$mus
    fitted.sds <- MCEM_fit$paramEst$sds
    clustermemb_prob <- computeClusterMembProb(data, MCEM_fit = MCEM_fit)

    data$cluster <- as.factor(apply(clustermemb_prob, 1, which.max))
    data$SNPtext <- sapply(1:n, get_SNP_text)

    p <- ggplot(
      data = data,
      aes(
        x = beta.exposure,
        y = beta.outcome,
        xmin = (beta.exposure - se.exposure),
        xmax = (beta.exposure + se.exposure),
        ymin = (beta.outcome - (tau * se.outcome)),
        ymax = (beta.outcome + (tau * se.outcome)),
        color = cluster
      )
    ) +
      geom_point(size = 1, shape = 1, aes(text = SNPtext)) +
      geom_errorbar(alpha = 0.5, width = 0) +
      geom_errorbarh(alpha = 0.5, height = 0) +
      expand_limits(x = 0, y = 0) +
      xlab(paste("SNP association with", exposure_name)) +
      ylab(paste("SNP association with", outcome_name)) +
      theme_classic(base_size = 15) +
      geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.4) +
      geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.4) +
      scale_color_brewer(palette = "Set1", name = "") +
      scale_fill_brewer(palette = "Set1", name = "")

    if (interactive) {
      # Convert ggplot to plotly
      p_interactive <- plotly::ggplotly(p, tooltip = c("text"))

      # Fix legend labels
      p_interactive <- plotly::plotly_build(p_interactive)
      for (k in 1:K) {
        p_interactive$x$data[[k]]$name <- as.character(k)
      }

      # Add lines and intervals
      for (k in 1:K) {
        y <- (fitted.mus[k]) * x_grid
        ymin <- (fitted.mus[k] - fitted.sds[k]) * x_grid
        ymax <- (fitted.mus[k] + fitted.sds[k]) * x_grid
        plot_data <- data.frame(k = k, x = x_grid, y = y, ymin = ymin, ymax = ymax)
        plot_data$k <- as.factor(plot_data$k)
        hovertxt_lines <- paste("Cluster #", k,
          "<br> Proportion: ", round(fitted.pis[k], 2),
          "<br> Mean: ", round(fitted.mus[k], 3),
          "<br> Std. Dev: ", round(fitted.sds[k], 2),
          sep = ""
        )
        p_interactive <- p_interactive %>%
          add_lines(
            x = ~x, y = ~y,
            color = ~k, colors = "Set1",
            text = hovertxt_lines,
            hoverinfo = "text",
            data = plot_data,
            showlegend = FALSE,
            inherit = TRUE
          ) %>%
          add_ribbons(
            x = ~x, ymin = ~ymin, ymax = ~ymax,
            color = ~k,
            colors = "Set1",
            hoverinfo = "none",
            opacity = .5,
            data = plot_data,
            showlegend = FALSE
          )
      }
      return(p_interactive)
    } else {
      # Add lines and intervals to ggplot
      p_build <- ggplot_build(p)
      x_grid <- seq(0, p_build$layout$panel_params[[1]]$x.range[2], .001)
      lines_dat <- data.frame(x = rep(x_grid, K))
      lines_dat$y <- as.numeric(sapply(1:K, function(k) fitted.mus[k] * x_grid))

      lines_dat$ymin <- as.numeric(sapply(1:K, function(k) (fitted.mus[k] - fitted.sds[k]) * x_grid))
      lines_dat$ymax <- as.numeric(sapply(1:K, function(k) (fitted.mus[k] + fitted.sds[k]) * x_grid))
      lines_dat$k <- as.factor(sapply(1:K, function(k) rep(k, length(x_grid))))

      p <- p + theme(
        legend.justification = c(0, 0),
        legend.position = c(.05, .005),
        legend.direction = "horizontal"
      ) +
        geom_line(
          data = lines_dat,
          aes(x = x, y = y, color = as.factor(k)),
          linetype = 1, size = .8,
          inherit.aes = FALSE
        ) +
        geom_ribbon(
          data = lines_dat,
          aes(
            x = x, ymin = ymin, ymax = ymax,
            fill = as.factor(k)
          ),
          alpha = 0.4,
          inherit.aes = FALSE
        )
      return(p)
    }
  }
}

MREMalt.scatterplot <- function(data, EM_fit = NULL,
                                exposure_name = "exposure", outcome_name = "outcome", interactive = FALSE) {
  data$tau <- 1
  fitted.pis <- EM_fit$paramEst$pis
  fitted.mus <- EM_fit$paramEst$mus
  fitted.sds <- EM_fit$paramEst$sds
  clusterMembProb <- EM_fit$clusterMembProb
  K <- length(fitted.pis)
  p <- nrow(data)

  data$SNPtext <- sapply(1:p, function(i) paste(round(clusterMembProb[i, ], digits = 3), collapse = ", "))

  clusters <- as.factor(apply(EM_fit$clusterMembProb, 1, which.max))
  names(clusters) <- data$SNP

  colors <- RColorBrewer::brewer.pal(5, "Set1")

  p <- ggplot(data = data, aes(
    x = beta.exposure, y = beta.outcome, xmin = (beta.exposure - se.exposure), xmax = (beta.exposure + se.exposure),
    ymin = (beta.outcome - (tau * se.outcome)), ymax = (beta.outcome + (tau * se.outcome)),
    color = clusters
  )) +
    geom_point(
      size = 1, shape = 1,
      aes(text = SNPtext)
    ) +
    geom_errorbar(alpha = 0.5, width = 0) +
    geom_errorbarh(alpha = 0.5, height = 0) +
    expand_limits(x = 0, y = 0) +
    xlab(paste("SNP association with", exposure_name)) +
    ylab(paste("SNP association with", outcome_name)) +
    theme_classic(base_size = 15) +
    geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.4) +
    scale_color_brewer(palette = "Set1", name = "") +
    scale_fill_brewer(palette = "Set1", name = "")

  p_build <- ggplot_build(p)
  x_grid <- seq(0, p_build$layout$panel_params[[1]]$x.range[2], .001)
  lines_dat <- data.frame(x = rep(x_grid, K))
  lines_dat$y <- as.numeric(sapply(1:K, function(k) fitted.mus[k] * x_grid))
  lines_dat$k <- as.factor(sapply(1:K, function(k) rep(k, length(x_grid))))


  if (interactive) {
    p <- plotly::ggplotly(p, tooltip = c("text"))
    p <- add_lines(p,
      x = ~x, y = ~y, color = ~k, colors = "Set1",
      data = lines_dat, inherit = TRUE, hoverinfo = "none"
    )
  } else {
    p <- p + theme(
      legend.justification = c(0, 0), legend.position = c(.05, .005),
      legend.direction = "horizontal"
    ) +
      geom_line(
        data = lines_dat,
        aes(x = x, y = y, color = as.factor(k)),
        linetype = 1, size = .8, inherit.aes = FALSE
      )
  }


  return(p)
}
