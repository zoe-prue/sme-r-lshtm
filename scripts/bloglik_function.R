# Binomial Log-Likelihood Function
# Converted from Stata's bloglik.ado

library(tidyverse)

#' Calculate and Plot Log-Likelihood for Binomial Parameter
#'
#' Plots exact and approximate quadratic log-likelihood ratios for a binomial
#' parameter, either on the probability scale or log-odds scale.
#'
#' @param d Integer: Number of "disease" or "event" cases
#' @param h Integer: Number of "healthy" or "non-event" cases
#' @param logodds Logical: If TRUE, plot on log-odds scale instead of probability scale
#' @param cut Numeric: Cutoff value for log-likelihood ratio (default -1.921 (??), approximately -2 for 95% CI)
#' @param samex Logical: Whether to use same x-axis scale. Not allowed with logodds=TRUE.
#'
#' @return A ggplot object showing exact and approximate log-likelihood ratio functions
#' @export
bloglik <- function(d, h, logodds = FALSE, cut = -1.921, samex = FALSE) {
  

  # INPUT VALIDATION
  
  # Check that d and h are integers
  if (!is.numeric(d) || d != as.integer(d) || d < 0) {
    stop("d must be a non-negative integer")
  }
  if (!is.numeric(h) || h != as.integer(h) || h < 0) {
    stop("h must be a non-negative integer")
  }
  
  # Check for zero values - log-likelihood is infinite when d=0 or h=0
  if (d == 0 || h == 0) {
    stop("Not possible because likelihood is infinite")
  }
  
  # samex not allowed with logodds -- needed?
  if (samex && logodds) {
    stop("samex not allowed with logodds")
  }
  

  # SETUP PARAMETERS

  
  # Calculate total and MLE for probability
  N <- d + h
  M <- d / N
  
  # Standard error for probability
  S <- sqrt(M * (1 - M) / N)
  
  # Create sequence of parameter values
  length <- 10001
  R <- 5.0
  
  # Maximum log-likelihood (at MLE)
  max_log_lik <- d * log(M) + h * log(1 - M)
  

  # FIND EXACT SUPPORTED RANGE USING BISECTION

  # Find where log-likelihood ratio equals the cutoff value
  
  # Target log-likelihood (max + cut)
  target_log_lik <- max_log_lik + cut
  tol <- 0.0001
  
  # ---- Find lower bound ----
  r1 <- 0.00001
  r2 <- M
  
  f1 <- d * log(r1) + h * log(1 - r1) - target_log_lik
  f2 <- d * log(r2) + h * log(1 - r2) - target_log_lik
  
  f <- 1
  while (abs(f) > tol) {
    r <- (r1 + r2) * 0.5
    f <- d * log(r) + h * log(1 - r) - target_log_lik
    
    if (f * f1 > 0) {
      r1 <- r
      f1 <- f
    } else {
      r2 <- r
      f2 <- f
    }
  }
  low <- r
  
  # ---- Find upper bound ----
  r1 <- 0.99999
  r2 <- M
  
  f1 <- d * log(r1) + h * log(1 - r1) - target_log_lik
  f2 <- d * log(r2) + h * log(1 - r2) - target_log_lik
  
  f <- 1
  while (abs(f) > tol) {
    r <- (r1 + r2) * 0.5
    f <- d * log(r) + h * log(1 - r) - target_log_lik
    
    if (f * f1 > 0) {
      r1 <- r
      f1 <- f
    } else {
      r2 <- r
      f2 <- f
    }
  }
  high <- r
  

  # CALCULATE AND DISPLAY RESULTS

  
  if (!logodds) {
    # ---- Results on probability scale ----
    
    cat("Most likely value for risk parameter:           ", sprintf("%7.5f", M), "\n")
    cat("\n")
    cat("cut-point:                                      ", cut, "\n")
    cat("Likelihood based limits for risk parameter:     ", 
        sprintf("%7.5f", low), " ", sprintf("%7.5f", high), "\n")
    
    # Approximate quadratic limits
    # Based on -2*log(LR) ~ chi-square(1), so log(LR) ~ -0.5*((p-M)/S)^2
    # At cut point: (p-M)/S = sqrt(-2*cut)
    plow <- M - sqrt(-cut * 2) * S
    phigh <- M + sqrt(-cut * 2) * S
    
    cat("Approx quadratic limits for risk parameter:     ", 
        sprintf("%7.5f", plow), " ", sprintf("%7.5f", phigh), "\n")
    
    # ---- Create data for plotting ----
    
    # Range for plotting
    start <- max(1 / length, M - R * S)
    stop <- min(1 - 1 / length, M + R * S)
    
    plot_data <- tibble(
      param = start + (stop - start) * (seq_len(length) - 1) / length
    ) |>
      mutate(
        # Exact log-likelihood ratio
        true = d * log(param) + h * log(1 - param) - max_log_lik,
        # Approximate quadratic log-likelihood
        approx = -0.5 * ((param - M) / S)^2
      )
    
    plot_title <- sprintf("Log likelihood ratio for risk parameter: D = %d, H = %d", d, h)
    
  } else {
    # ---- Results on log-odds scale ----
    
    # Transform to log-odds scale
    M_logodds <- log(d / h)
    S_logodds <- sqrt(1 / d + 1 / h)
    low_logodds <- log(low / (1 - low))
    high_logodds <- log(high / (1 - high))
    
    cat("Most likely value for logodds parameter:        ", sprintf("%7.5f", M_logodds), "\n")
    cat("\n")
    cat("cut-point:                                      ", cut, "\n")
    cat("Likelihood based limits for logodds parameter:  ", 
        sprintf("%7.5f", low_logodds), " ", sprintf("%7.5f", high_logodds), "\n")
    
    # Approximate quadratic limits on log-odds scale
    plow_logodds <- M_logodds - sqrt(-cut * 2) * S_logodds
    phigh_logodds <- M_logodds + sqrt(-cut * 2) * S_logodds
    
    cat("Approx quadratic limits for logodds parameter:  ", 
        sprintf("%7.5f", plow_logodds), " ", sprintf("%7.5f", phigh_logodds), "\n")
    
    # ---- Create data for plotting on log-odds scale ----
    
    plot_data <- tibble(
      param = M_logodds - R * S_logodds + 2 * R * S_logodds * (seq_len(length) - 1) / length
    ) |>
      mutate(
        # Exact log-likelihood ratio on log-odds scale
        # If theta = log(p/(1-p)), then p = exp(theta)/(1+exp(theta))
        # log L = d*theta - N*log(1+exp(theta)) + constant
        true = d * param - N * log(1 + exp(param)) - max_log_lik,
        # Approximate quadratic
        approx = -0.5 * ((param - M_logodds) / S_logodds)^2
      )
    
    plot_title <- sprintf("Log likelihood ratio for logodds parameter: D = %d, H = %d", d, h)
  }
  

  # ============================================================================
  # CREATE PLOT
  # ============================================================================
  
  # Filter for reasonable viewing range
  plot_data_filtered <- plot_data |>
    filter(true > -5 & approx > -4)
  
  # Prepare data for shading supported region (where exact log-lik > cut)
  shade_data <- plot_data_filtered |>
    filter(true >= cut)
  
  # Determine appropriate x-axis label
  if (!logodds) {
    x_label <- expression(pi)
  } else {
    x_label <- expression(log(pi/(1-pi)))
  }
  
  # Create plot with both exact and approximate log-likelihoods
  p <- ggplot(plot_data_filtered, aes(x = param)) +
    geom_line(aes(y = true, color = "Exact"), linewidth = 0.8) +
    geom_line(aes(y = approx, color = "Approximate"), linewidth = 0.8) +
    scale_color_manual(
      values = c("Exact" = "red", "Approximate" = "blue"),
      name = NULL
    ) +
    labs(
      title = plot_title,
      x = x_label,
      y = "Log likelihood ratio"
    ) +
    scale_y_continuous(breaks = seq(0, -4, -1)) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Add shading for supported region (always show this)
  p <- p + 
    geom_ribbon(
      data = shade_data,
      aes(ymin = cut, ymax = true),
      fill = "lightblue",
      alpha = 0.3
    )
  
  # Add cutoff line (red dashed to match other functions)
  p <- p + 
    geom_hline(yintercept = cut, linetype = "dashed", color = "red", linewidth = 0.8)
  
  # Add x-axis formatting if requested
  if (samex && !logodds) {
    p <- p + scale_x_continuous(breaks = seq(0.1, 0.9, 0.1))
  }
  

  # CREATE PLOT
  
  # Filter for reasonable viewing range
  plot_data_filtered <- plot_data |>
    filter(true > -5 & approx > -4)
  
  # Prepare data for shading supported region (where exact log-lik > cut)
  shade_data <- plot_data_filtered |>
    filter(true >= cut)
  
  # Determine appropriate x-axis label
  if (!logodds) {
    x_label <- expression(pi)
  } else {
    x_label <- expression(log(pi/(1-pi)))
  }
  
  # Create plot with both exact and approximate log-likelihoods
  p <- ggplot(plot_data_filtered, aes(x = param)) +
    geom_line(aes(y = true, color = "Exact"), linewidth = 0.8) +
    geom_line(aes(y = approx, color = "Approximate"), linewidth = 0.8) +
    scale_color_manual(
      values = c("Exact" = "red", "Approximate" = "blue"),
      name = NULL
    ) +
    labs(
      title = plot_title,
      x = x_label,
      y = "Log likelihood ratio"
    ) +
    scale_y_continuous(breaks = seq(0, -4, -1)) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Add shading for supported region (always show this)
  p <- p + 
    geom_ribbon(
      data = shade_data,
      aes(ymin = cut, ymax = true),
      fill = "lightblue",
      alpha = 0.3
    )
  
  # Add cutoff line (red dashed to match other functions)
  p <- p + 
    geom_hline(yintercept = cut, linetype = "dashed", color = "black", linewidth = 0.8)
  
  # Add x-axis formatting if requested
  if (samex && !logodds) {
    p <- p + scale_x_continuous(breaks = seq(0.1, 0.9, 0.1))
  }
  

  # BACK-TRANSFORM TO ORIGINAL SCALE IF ON LOG-ODDS
  
  if (logodds) {
    cat("\n")
    cat("Back on original risk scale\n")
    
    # Transform back to probability scale
    low_risk <- exp(low_logodds) / (1 + exp(low_logodds))
    high_risk <- exp(high_logodds) / (1 + exp(high_logodds))
    plow_risk <- exp(plow_logodds) / (1 + exp(plow_logodds))
    phigh_risk <- exp(phigh_logodds) / (1 + exp(phigh_logodds))
    M_risk <- exp(M_logodds) / (1 + exp(M_logodds))
    
    cat("Most likely value for risk parameter:           ", sprintf("%7.5f", M_risk), "\n")
    cat("Likelihood based limits for risk parameter:     ", 
        sprintf("%7.5f", low_risk), " ", sprintf("%7.5f", high_risk), "\n")
    cat("Approx quadratic limits for risk parameter:     ", 
        sprintf("%7.5f", plow_risk), " ", sprintf("%7.5f", phigh_risk), "\n")
  }
  
  return(p)
}
