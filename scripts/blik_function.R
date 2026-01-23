# Binomial Likelihood Ratio Function
# Converted from Stata's blik.ado

library(tidyverse)

#' Calculate and Plot Likelihood Ratios for Binomial Parameter
#'
#' @param d Integer: Number of "disease" or "event" cases
#' @param h Integer: Number of "healthy" or "non-event" cases  
#' @param cut Numeric: Cutoff value for likelihood ratio (default 0.1465, approximately 1/8 for 2 log-likelihood units)
#' @param null Numeric: Null hypothesis value for pi (optional)
#' @param pval Logical: Whether to calculate approximate p-value when null is specified
#'
#' @return A ggplot object showing the likelihood ratio function
#' @export

blik <- function(d, h, cut = 0.1465, null = NULL, pval = FALSE) {
  

  # INPUT VALIDATION
  
  # Check that d and h are integers
  if (!is.numeric(d) || d != as.integer(d) || d < 0) {
    stop("d must be a non-negative integer")
  }
  if (!is.numeric(h) || h != as.integer(h) || h < 0) {
    stop("h must be a non-negative integer")
  }
  
  # Check for both being zero
  if (h == 0 && d == 0) {
    stop("No data: both d and h are zero")
  }
  
  # If h is zero, exchange d and h (and adjust null value if provided)
  if (h == 0) {
    temp <- d
    d <- h
    h <- temp
    if (!is.null(null)) {
      null <- 1 - null
    }
    message("D and H have been interchanged")
  }
  
  # SETUP PARAMETERS
  
  # Calculate maximum likelihood estimate (MLE) for pi
  N <- d + h
  M <- d / N
  
  # Create sequence of pi values to evaluate likelihood over
  # Using 10,001 points to match Stata's resolution??
  length <- 10001
  

  # CREATE DATA FRAME WITH PI VALUES AND LIKELIHOOD RATIOS
  
  # Generate evenly spaced pi values from just above 0 to just below 1
  
  lik_data <- tibble(
    pi = (seq_len(length) - 0.5) / length
  )
  

  # CALCULATE LIKELIHOOD RATIOS
  
  # The likelihood ratio compares L(pi) to L(M), where M is the MLE
  # LR(pi) = [pi^d * (1-pi)^h] / [M^d * (1-M)^h] - double check later
  # Taking logs: log(LR) = d*log(pi/M) + h*log((1-pi)/(1-M))
  
  if (d == 0) {
    # Special case when d = 0: only depends on h 
    # LR = (1-pi)^h / (1-M)^h = ((1-pi)/(1-M))^h
    lik_data <- lik_data |>
      mutate(lik = (1 - pi)^h)
  } else {
    # General case: calculate log likelihood ratio, then exponentiate
    lik_data <- lik_data |>
      mutate(
        log_lik_ratio = d * log(pi / M) + h * log((1 - pi) / (1 - M)),
        lik = exp(log_lik_ratio)
      ) |>
      select(pi, lik)  # Keep only needed columns
  }
  
  # FIND SUPPORTED RANGE (always calculate this for the shaded region)
  
  if (d == 0) {
    # When d = 0, lower bound is 0
    low <- 0
    # Upper bound found by solving (1-pi)^h = cut for pi
    high <- 1 - exp(log(cut) / h)
    
  } else {
    # Use bisection method to find lower and upper bounds
    # where likelihood ratio = cut
    
    max_log_lik <- d * log(M) + h * log(1 - M)
    target_log_lik <- max_log_lik + log(cut)
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
  }
  
  # PRINT RESULTS AND CALCULATE NULL HYPOTHESIS VALUES (if specified)
  
  if (is.null(null)) {
    # Print results for supported range
    cat("Most likely value for pi:           ", sprintf("%7.5f", M), "\n")
    cat("Likelihood based limits for pi:     ", sprintf("%7.5f", low), " ", 
        sprintf("%7.5f", high), "\n")
    cat("cut-point:                          ", cut, "\n")
    
  } else {
    # Calculate likelihood ratio at null value
    if (d == 0) {
      lrnull <- (1 - null)^h
    } else {
      log_lrnull <- d * log(null / M) + h * log((1 - null) / (1 - M))
      lrnull <- exp(log_lrnull)
    }
    
    # Print results
    cat("Most likely value for pi:           ", sprintf("%7.5f", M), "\n")
    cat("Null value for pi:                  ", sprintf("%7.5f", null), "\n")
    cat("Lik ratio for null value:           ", sprintf("%7.5f", lrnull), "\n")
    
    # Calculate approximate p-value if requested
    if (pval) {
      pval_result <- pchisq(-2 * log(lrnull), df = 1, lower.tail = FALSE)
      cat("Approx p-value:                     ", sprintf("%5.4f", pval_result), "\n")
    }
  }
  
  # CREATE PLOT
  
  # Set up title
  plot_title <- sprintf("Likelihood ratio for risk parameter: D = %d, H = %d", d, h)
  
  # Filter out very small likelihood values for cleaner plotting??
  plot_data <- lik_data |>
    filter(lik > 0.001)
  
  # Prepare data for shading supported region (always show this)
  shade_data <- plot_data |>
    filter(pi >= low & pi <= high)
  
  # Start building plot
  p <- ggplot(plot_data, aes(x = pi, y = lik)) +
    geom_line(color = "blue", linewidth = 1) +
    labs(
      title = plot_title,
      x = expression(pi),  # Greek letter π
      y = "Likelihood ratio"
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.3)
    )
  
  # Add shading for supported region (always show this)
  p <- p + 
    geom_ribbon(
      data = shade_data,
      aes(ymin = 0, ymax = lik),
      fill = "lightblue",
      alpha = 0.3
    )
  
  # Add y-axis breaks
  p <- p + scale_y_continuous(breaks = seq(0, 1, 0.2))
  
  # Add vertical line at MLE with label
  p <- p + 
    geom_vline(xintercept = M, linetype = "solid", color = "black", linewidth = 0.8) +
    annotate("text", x = M * 1.05, y = 1.05, 
             label = "hat(pi)", 
             parse = TRUE,
             size = 4, fontface = "bold")
  
  # Add points at key locations (always show these)
  key_points <- tibble(
    pi = c(M, low, high),
    lik = c(1, cut, cut)
  )
  
  p <- p + 
    geom_point(data = key_points, aes(x = pi, y = lik), 
               size = 3, color = "black") +
    # Add cutoff line (always show this)
    geom_hline(yintercept = cut, linetype = "dashed", color = "red", linewidth = 0.8)
  
  # If testing null hypothesis, add additional line showing LR at null
  if (!is.null(null)) {
    p <- p + 
      # Show likelihood ratio at null value with different color
      geom_hline(yintercept = lrnull, linetype = "dotted", color = "orange", linewidth = 1) +
      # Show null value on x-axis with label
      geom_vline(xintercept = null, linetype = "dotted", color = "darkgreen", linewidth = 0.8) +
      annotate("text", x = null, y = 1.08, 
               label = "H[0]", 
               parse = TRUE,
               size = 4, fontface = "bold", color = "darkgreen")
  }
  
  return(p)
}
