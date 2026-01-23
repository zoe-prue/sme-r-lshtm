# Poisson Likelihood Function
# Converted from Stata's plik.ado

library(tidyverse)

#' Calculate and Plot Likelihood for Poisson Rate Parameter
#'
#' Plots the likelihood ratio function for a Poisson rate parameter (lambda)
#' given observed events and person-years.
#'
#' @param d Integer: Number of events observed
#' @param y Numeric: Person-years of observation
#' @param null Numeric: Null hypothesis value for rate (per unit specified in per argument)
#' @param cut Numeric: Cutoff value for likelihood ratio (default 0.1465, approximately 1/8 for 2 log-likelihood units)
#' @param pval Logical: Whether to calculate approximate p-value when null is specified
#' @param per Numeric: Rate multiplier for display (default 1000, so rates are per 1000 person-years)
#'
#' @return A ggplot object showing the likelihood ratio function
#' @export

plik <- function(d, y, null = NULL, cut = 0.1465, pval = FALSE, per = 1000) {
  

  # INPUT VALIDATION
  
  # Check that d is a non-negative integer
  if (!is.numeric(d) || d != as.integer(d) || d < 0) {
    stop("d must be a non-negative integer")
  }
  
  # Check that y is positive
  if (!is.numeric(y) || y <= 0) {
    stop("needs person-years (y must be positive)")
  }
  
  # SETUP PARAMETERS
  
  length <- 10001
  
  # Display message about rate scale
  cat("\n")
  cat("ALL RATES PER ", sprintf("%7.0f", per), "\n")
  
  # Calculate MLE for rate (events per person-year)
  M <- d / y
  
  plot_title <- sprintf("Likelihood ratio for rate parameter: D = %d, Y = %g", d, y)
  
  # SPECIAL CASE: d = 0 (no events observed)
  
  if (d == 0) {
    
    # When d=0, likelihood is L(lambda) = exp(-lambda * y)
    # Create sequence of lambda values up to Stata upper limit
    lik_data <- tibble(
      lambda = (seq_len(length) - 1) * 5 / (y * length)
    ) |>
      mutate(
        lik = exp(-lambda * y),
        # Convert to display scale (per 1000, etc.)
        lambda_display = lambda * per
      )
    
    # FIND SUPPORTED RANGE (always calculate for shaded region)
    low <- 0
    high <- -log(cut) / y
    
    # PRINT RESULTS AND CALCULATE NULL HYPOTHESIS VALUES (if specified)
    
    if (is.null(null)) {
      # Print results for supported range
      cat("Most likely value for lambda:              ", sprintf("%7.1f", M * per), "\n")
      cat("Likelihood based limits for lambda:        ", sprintf("%7.1f", low * per), " ", sprintf("%7.5f", high * per), "\n")
      cat("cut-point:                                 ", cut, "\n")
      
    } else {
      # Evaluate at null value
      nval <- null / per  # Convert from display scale to per person-year
      lrnull <- exp(-nval * y)
      
      cat("Most likely value for lambda:              ", sprintf("%7.1f", M * per), "\n")
      cat("Null value for lambda:                     ", sprintf("%7.1f", null), "\n")
      cat("Lik ratio for null value:                  ", sprintf("%7.5f", lrnull), "\n")
      
      if (pval) {
        pval_result <- pchisq(-2 * log(lrnull), df = 1, lower.tail = FALSE)
        cat("Approx p-value:                                 ", sprintf("%5.3f", pval_result), "\n")
      }
    }
    
    # Prepare data for shading supported region (always show this)
    shade_data <- lik_data |>
      filter(lambda_display >= low * per & lambda_display <= high * per) |>
      filter(lik > 0.01)
    
    # CREATE PLOT
    # Create plot with subtitle showing rate scale
    p <- ggplot(lik_data |> filter(lik > 0.01), 
                aes(x = lambda_display, y = lik)) +
      geom_line(color = "blue", linewidth = 1) +
      labs(
        title = plot_title,
        subtitle = sprintf("Rate per %g person-years", per),
        x = expression(lambda),  # Greek letter λ
        y = "Likelihood ratio"
      ) +
      scale_y_continuous(breaks = seq(0, 1, 0.2)) +
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
    
    # Calculate a small offset for label (1% of x-axis range)
    x_offset <- (max(lik_data$lambda_display, na.rm = TRUE) - 
                   min(lik_data$lambda_display, na.rm = TRUE)) * 0.01
    
    # Add vertical line at MLE with label
    p <- p + 
      geom_vline(xintercept = M * per, linetype = "solid", color = "black", linewidth = 0.8) +
      annotate("text", x = M * per + x_offset, y = 1.08,  # Jittered 
               label = "hat(lambda)", 
               parse = TRUE,
               size = 4, fontface = "bold")
    
    # Add points at key locations (always show these)
    key_points <- tibble(
      lambda_display = c(M * per, low * per, high * per),
      lik = c(exp(-M * y), cut, cut)
    )
    
    p <- p + 
      geom_point(data = key_points, aes(x = lambda_display, y = lik), 
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
                 label = expression(H[0]), 
                 size = 4, fontface = "bold", color = "darkgreen")
    }
    
    
    # GENERAL CASE: d > 0
    
  } else {
    
    R <- 4.0
    S <- sqrt(1 / d)
    
    # Create range for lambda centered on MLE
    start <- M / exp(R * S)
    stop <- M * exp(R * S)
    
    # Maximum log-likelihood
    max_log_lik <- d * log(M) - d
    
    lik_data <- tibble(
      lambda = start + (stop - start) * (seq_len(length) - 1) / length
    ) |>
      mutate(
        # Calculate log-likelihood ratio, then exponentiate
        log_lik_ratio = d * log(lambda) - y * lambda - max_log_lik,
        lik = exp(log_lik_ratio),
        # Convert to display scale
        lambda_display = lambda * per
      ) |>
      select(lambda, lambda_display, lik)
    
    # FIND SUPPORTED RANGE (always calculate for shaded region)
    
    # Find lower bound using bisection
    r1 <- 0.01 * M
    r2 <- M
    
    f1 <- d * log(r1) - r1 * y - max_log_lik - log(cut)
    f2 <- d * log(r2) - r2 * y - max_log_lik - log(cut)
    
    f <- 1
    tol <- 0.0001
    
    while (abs(f) > tol) {
      r <- (r1 + r2) * 0.5
      f <- d * log(r) - r * y - max_log_lik - log(cut)
      
      if (f * f1 > 0) {
        r1 <- r
        f1 <- f
      } else {
        r2 <- r
        f2 <- f
      }
    }
    low <- r
    
    # Find upper bound using bisection
    r1 <- 10 * M
    r2 <- M
    
    f1 <- d * log(r1) - r1 * y - max_log_lik - log(cut)
    f2 <- d * log(r2) - r2 * y - max_log_lik - log(cut)
    
    f <- 1
    while (abs(f) > tol) {
      r <- (r1 + r2) * 0.5
      f <- d * log(r) - r * y - max_log_lik - log(cut)
      
      if (f * f1 > 0) {
        r1 <- r
        f1 <- f
      } else {
        r2 <- r
        f2 <- f
      }
    }
    high <- r
    
    # PRINT RESULTS AND CALCULATE NULL HYPOTHESIS VALUES (if specified)
    
    if (is.null(null)) {
      # Print results for supported range
      cat("Most likely value for lambda:                   ", sprintf("%7.1f", M * per), "\n")
      cat("Likelihood based limits for lambda:             ", 
          sprintf("%7.1f", low * per), " ", sprintf("%7.1f", high * per), "\n")
      cat("cut-point:                                      ", cut, "\n")
      
    } else {
      # Evaluate likelihood at null value
      nval <- null / per  # Convert from display scale
      lrnull <- d * log(nval) - y * nval - max_log_lik
      lrnull <- exp(lrnull)
      
      cat("Most likely value for lambda:                   ", sprintf("%7.1f", M * per), "\n")
      cat("Null value for lambda:                          ", sprintf("%7.1f", null), "\n")
      cat("Lik ratio for null value:                       ", sprintf("%7.5f", lrnull), "\n")
      
      if (pval) {
        pval_result <- pchisq(-2 * log(lrnull), df = 1, lower.tail = FALSE)
        cat("Approx p-value:                                 ", sprintf("%5.3f", pval_result), "\n")
      }
    }
    
    # Prepare data for shading supported region (always show this)
    shade_data <- lik_data |>
      filter(lambda_display >= low * per & lambda_display <= high * per) |>
      filter(lik > 0.01)
    
    # CREATE PLOT
    
    p <- ggplot(lik_data |> filter(lik > 0.01), 
                aes(x = lambda_display, y = lik)) +
      geom_line(color = "blue", linewidth = 1) +
      labs(
        title = plot_title,
        subtitle = sprintf("Rate per %g person-years", per),
        x = expression(lambda),  # Greek letter λ
        y = "Likelihood ratio"
      ) +
      scale_y_continuous(breaks = seq(0, 1, 0.2)) +
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
    
    # Calculate a small offset for label (1% of x-axis range)
    x_offset <- (max(lik_data$lambda_display, na.rm = TRUE) - 
                   min(lik_data$lambda_display, na.rm = TRUE)) * 0.01
    
    # Add vertical line at MLE with label
    p <- p + 
      geom_vline(xintercept = M * per, linetype = "solid", color = "black", linewidth = 0.8) +
      annotate("text", x = M * per + x_offset, y = 1.08,  # Jittered 
               label = "hat(lambda)", 
               parse = TRUE,
               size = 4, fontface = "bold")
    
    # Add points at key locations (always show these)
    key_points <- tibble(
      lambda_display = c(M * per, low * per, high * per),
      lik = c(1, cut, cut)
    )
    
    p <- p + 
      geom_point(data = key_points, aes(x = lambda_display, y = lik), 
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
        annotate("text", x = null + x_offset, y = 1.08, 
                 label = expression(H[0]), 
                 size = 4, fontface = "bold", color = "darkgreen")
    }
  }
  
  return(p)
}