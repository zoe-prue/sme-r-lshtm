# Poisson Log-Likelihood Function
# Converted from Stata's ploglik.ado

library(tidyverse)

#' Calculate and Plot Log-Likelihood for Poisson Rate Parameter
#'
#' Plots exact and approximate (quadratic) log-likelihood ratios for a Poisson
#' rate parameter, either on the rate scale or log-rate scale.
#'
#' @param d Integer: Number of events observed
#' @param y Numeric: Person-years of observation
#' @param lograte Logical: If TRUE, plot on log-rate scale instead of rate scale
#' @param cut Numeric: Cutoff value for log-likelihood ratio (default -1.921, approximately -2 for 95% CI)
#' @param per Numeric: Rate multiplier for display (default 1000, so rates are per 1000 person-years)
#'
#' @return A ggplot object showing exact and approximate log-likelihood ratio functions
#' @export
ploglik <- function(d, y, lograte = FALSE, cut = -1.921, per = 1000) {

  # INPUT VALIDATION
  
  # Check for zero values - log-likelihood is infinite when d=0
  if (!is.numeric(d) || d != as.integer(d) || d <= 0) {
    stop("Not possible because likelihood is infinite (d must be positive)")
  }
  
  # Check that y is positive
  if (!is.numeric(y) || y <= 0) {
    stop("needs person-years (y must be positive)")
  }
  

  # SETUP PARAMETERS
  
  length <- 10001
  R <- 3.0
  
  # Display message about rate scale
  cat("\n")
  cat("ALL RATES PER ", sprintf("%7.0f", per), "\n")
  
  # Calculate MLE for rate and its standard error
  M <- d / y
  S <- sqrt(d) / y
  
  # Maximum log-likelihood
  max_log_lik <- d * log(d / y) - d
  

  # FIND EXACT SUPPORTED RANGE USING BISECTION
  
  # Target log-likelihood
  target_log_lik <- max_log_lik + cut
  tol <- 0.0001
  
  # ---- Find lower bound ----
  r1 <- 0.01 * M
  r2 <- M
  
  f1 <- d * log(r1) - r1 * y - target_log_lik
  f2 <- d * log(r2) - r2 * y - target_log_lik
  
  f <- 1
  while (abs(f) > tol) {
    r <- (r1 + r2) * 0.5
    f <- d * log(r) - r * y - target_log_lik
    
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
  r1 <- 10 * M
  r2 <- M
  
  f1 <- d * log(r1) - r1 * y - target_log_lik
  f2 <- d * log(r2) - r2 * y - target_log_lik
  
  f <- 1
  while (abs(f) > tol) {
    r <- (r1 + r2) * 0.5
    f <- d * log(r) - r * y - target_log_lik
    
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
 
  
  if (!lograte) {
    # ---- Results on rate scale ----
    
    cat("Most likely value for rate parameter:           ", sprintf("%7.2f", M * per), "\n")
    cat("\n")
    cat("cut-point:                                      ", cut, "\n")
    cat("Likelihood based limits for rate parameter:     ", 
        sprintf("%7.2f", low * per), " ", sprintf("%7.2f", high * per), "\n")
    
    # Approximate quadratic limits
    # Based on -2*log(LR) ~ chi-square(1), so log(LR) ~ -0.5*((lambda-M)/S)^2
    plow <- M - sqrt(-cut * 2) * S
    phigh <- M + sqrt(-cut * 2) * S
    
    cat("Approx quadratic limits for rate parameter:     ", 
        sprintf("%7.2f", plow * per), " ", sprintf("%7.2f", phigh * per), "\n")
    
    # ---- Create data for plotting on rate scale ----
    
    start <- M - R * S
    stop <- M + R * S
    
    plot_data <- tibble(
      param = start + (stop - start) * (seq_len(length) - 1) / length
    ) |>
      mutate(
        # Exact log-likelihood ratio
        true = d * log(param) - y * param - max_log_lik,
        # Approximate quadratic log-likelihood
        approx = -0.5 * ((param - M) / S)^2,
        # Convert to display scale
        param_display = param * per
      )
    
    plot_title <- sprintf("Log likelihood for rate parameter: D = %d, Y = %g", d, y)
    
  } else {
    # ---- Results on log-rate scale ----
    
    # Transform to log-rate scale
    M_lograte <- log(d / y)
    S_lograte <- sqrt(1 / d)
    low_lograte <- log(low)
    high_lograte <- log(high)
    
    cat("Most likely value for log rate parameter:       ", 
        sprintf("%7.2f", M_lograte + log(per)), "\n")
    cat("\n")
    cat("cut-point:                                      ", cut, "\n")
    cat("Likelihood based limits for log rate parameter: ", 
        sprintf("%7.2f", low_lograte + log(per)), " ", 
        sprintf("%7.2f", high_lograte + log(per)), "\n")
    
    # Approximate quadratic limits on log-rate scale
    plow_lograte <- M_lograte - sqrt(-cut * 2) * S_lograte
    phigh_lograte <- M_lograte + sqrt(-cut * 2) * S_lograte
    
    cat("Approx quadratic limits for log rate parameter: ", 
        sprintf("%7.2f", plow_lograte + log(per)), " ", 
        sprintf("%7.2f", phigh_lograte + log(per)), "\n")
    
    # ---- Create data for plotting on log-rate scale ----
    
    plot_data <- tibble(
      param = M_lograte - R * S_lograte + 2 * R * S_lograte * (seq_len(length) - 1) / length
    ) |>
      mutate(
        # Exact log-likelihood ratio on log-rate scale
        # If theta = log(lambda), then lambda = exp(theta)
        # log L = d*theta - y*exp(theta) + constant
        true = d * param - y * exp(param) - max_log_lik,
        # Approximate quadratic
        approx = -0.5 * ((param - M_lograte) / S_lograte)^2,
        # Convert to display scale (add log(per))
        param_display = param + log(per)
      )
    
    plot_title <- sprintf("Log likelihood for lograte parameter: D = %d, Y = %g", d, y)
  }
  

  # CREATE PLOT

  
  # Filter for reasonable viewing range
  plot_data_filtered <- plot_data |>
    filter(approx > -2.0)
  
  # Prepare data for shading supported region (where exact log-lik > cut)
  shade_data <- plot_data_filtered |>
    filter(true >= cut)
  
  # Determine appropriate x-axis label
  if (!lograte) {
    x_label <- expression(lambda)
  } else {
    x_label <- expression(log(lambda))
  }
  
  # Create plot with both exact and approximate log-likelihoods
  p <- ggplot(plot_data_filtered, aes(x = param_display)) +
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
    scale_y_continuous(breaks = seq(0, -6, -1)) +
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
  
  return(p)
}
