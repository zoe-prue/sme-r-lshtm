
# Custom functions for matched case-control analysis
# AIms to replicate key output from Stata's mcc/mcci commands

#' Create a matched table for 1:1 matched case-control data
#' 
#' @param data A data frame containing the matched case-control data
#' @param case_var Name of the case/control indicator variable (1 = case, 0 = control)
#' @param exposure_var Name of the binary exposure variable (must be 0/1 or factor)
#' @param match_var Name of the matching set identifier variable
#' @return A matrix showing the matched table with row/column totals
#' 
matched_table <- function(data, case_var, exposure_var, match_var) {
  
  # Convert to data frame if tibble
  data <- as.data.frame(data)
  
  # Extract variables
  case <- data[[case_var]]
  exposure <- data[[exposure_var]]
  match_id <- data[[match_var]]
  

  # Convert exposure to numeric 0/1 if factor
  if (is.factor(exposure)) {
    exposure <- as.numeric(exposure) - 1
  }
  
  # Get unique matched sets
  sets <- unique(match_id)
  

  # Initialize counts for matched table
  # a = both exposed, b = case exposed/control unexposed

  # c = case unexposed/control exposed, d = both unexposed
  a <- b <- c <- d <- 0
  
  # This is a bit too messy - fix in future versions
  
  for (s in sets) {
    subset_data <- data[match_id == s, ]
    case_exp <- exposure[match_id == s & case == 1]
    ctrl_exp <- exposure[match_id == s & case == 0]
    
    # For 1:1 matching, take first control if multiple exist
    if (length(ctrl_exp) > 0) ctrl_exp <- ctrl_exp[1]
    if (length(case_exp) > 0) case_exp <- case_exp[1]
    
    if (length(case_exp) == 1 && length(ctrl_exp) == 1) {
      if (case_exp == 1 && ctrl_exp == 1) a <- a + 1
      else if (case_exp == 1 && ctrl_exp == 0) b <- b + 1
      else if (case_exp == 0 && ctrl_exp == 1) c <- c + 1
      else if (case_exp == 0 && ctrl_exp == 0) d <- d + 1
    }
  }
  
  # Create the matched table matrix
  # Rows = Controls, Columns = Cases
  matched_mat <- matrix(
    c(a, c, a + c,
      b, d, b + d,
      a + b, c + d, a + b + c + d),
    nrow = 3, ncol = 3, byrow = TRUE
  )
  
  rownames(matched_mat) <- c("Exposed", "Unexposed", "Total")
  colnames(matched_mat) <- c("Exposed", "Unexposed", "Total")
  
  # Add attributes for later use
  attr(matched_mat, "cells") <- c(a = a, b = b, c = c, d = d)
  attr(matched_mat, "table_type") <- "Controls (rows) by Cases (columns)"
  
  return(matched_mat)
}


#' Analyse matched case-control data (1:1 matching)
#' Replicates Stata's mcc command output
#' 
#' @param data A data frame containing the matched case-control data
#' @param case_var Name of the case/control indicator variable (1 = case, 0 = control)
#' @param exposure_var Name of the binary exposure variable (must be 0/1 or factor)
#' @param match_var Name of the matching set identifier variable
#' @param conf_level Confidence level (default 0.95)
#' @return Prints matched table and statistics, invisibly returns a list
#' 
matched_or <- function(data, case_var, exposure_var, match_var, conf_level = 0.95) {
  
  # Get the matched table
  mat <- matched_table(data, case_var, exposure_var, match_var)
  cells <- attr(mat, "cells")
  a <- cells["a"]; b <- cells["b"]; c <- cells["c"]; d <- cells["d"]
  
  # Number of pairs
  n_pairs <- a + b + c + d
  n_discordant <- b + c
  
  # Calculate odds ratio (ratio of discordant pairs)
  if (c == 0) {
    or <- Inf
    or_lower <- NA
    or_upper <- NA
    warning("Cannot calculate CI: no discordant pairs")
  } else if (b == 0) {
    or <- 0
    or_lower <- NA
    or_upper <- NA
    warning("Cannot calculate CI: no discordant pairs")
  } else {
    or <- b / c
    
    # Exact confidence interval using F distribution 
    # This should match Stata's exact CI??
    alpha <- 1 - conf_level
    or_lower <- (b / (c + 1)) * (1 / qf(1 - alpha/2, 2 * (c + 1), 2 * b))
    or_upper <- ((b + 1) / c) * qf(1 - alpha/2, 2 * (b + 1), 2 * c)
  }
  
  # McNemar's chi-squared (without continuity correction, to match Stata...)
  if (n_discordant > 0) {
    mcnemar_chi2 <- (b - c)^2 / (b + c)
    mcnemar_p <- pchisq(mcnemar_chi2, df = 1, lower.tail = FALSE)
  } else {
    mcnemar_chi2 <- NA
    mcnemar_p <- NA
  }
  
  # Exact McNemar p-value (two-sided binomial test)
  if (n_discordant > 0) {
    exact_p <- 2 * min(pbinom(b, n_discordant, 0.5), 
                       pbinom(b - 1, n_discordant, 0.5, lower.tail = FALSE))
    exact_p <- min(exact_p, 1)  # Cap at 1
  } else {
    exact_p <- NA
  }
  
  # Proportions with exposure
  prop_cases <- (a + b) / n_pairs
  prop_controls <- (a + c) / n_pairs
  
  # Print output (Stata-style as from other functions)
  cat("\nMatched Case-Control Analysis (1:1 Matching)\n")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")
  
  # Print matched table
  cat("                       Cases\n")
  cat("Controls        Exposed    Unexposed    Total\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat(sprintf("Exposed      %8d   %8d   %8d\n", a, c, a + c))
  cat(sprintf("Unexposed    %8d   %8d   %8d\n", b, d, b + d))
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat(sprintf("Total        %8d   %8d   %8d\n", a + b, c + d, n_pairs))
  cat("\n")
  
  # McNemar's test
  if (!is.na(mcnemar_chi2)) {
    cat(sprintf("McNemar's chi2(1) = %8.2f    Prob > chi2 = %.4f\n", 
                mcnemar_chi2, mcnemar_p))
    cat(sprintf("Exact McNemar probability  = %.4f\n\n", exact_p))
  }
  
  # Proportions
  cat("Proportion with exposure\n")
  cat(sprintf("  Cases:     %.4f\n", prop_cases))
  cat(sprintf("  Controls:  %.4f\n\n", prop_controls))
  
  # Odds ratio
  ci_label <- sprintf("%.0f%% CI", conf_level * 100)
  cat(sprintf("Odds ratio = %.4f    [%s: %.4f - %.4f] (exact)\n", 
              or, ci_label, or_lower, or_upper))
  cat("\n")
  cat("Note: OR calculated as ratio of discordant pairs (b/c)\n")
  cat(sprintf("      Discordant pairs: b = %d, c = %d\n", b, c))
  
  # Return results invisibly
  invisible(list(
    table = mat,
    cells = cells,
    n_pairs = n_pairs,
    n_discordant = n_discordant,
    odds_ratio = or,
    or_lower = or_lower,
    or_upper = or_upper,
    mcnemar_chi2 = mcnemar_chi2,
    mcnemar_p = mcnemar_p,
    exact_p = exact_p,
    prop_cases = prop_cases,
    prop_controls = prop_controls
  ))
}


#' Immediate version: Calculate matched OR from cell counts
#' REplicates Stata's mcci command
#' 
#' @param a Number of pairs where both case and control are exposed
#' @param b Number of pairs where case is exposed and control is unexposed
#' @param c Number of pairs where case is unexposed and control is exposed
#' @param d Number of pairs where both case and control are unexposed
#' @param conf_level Confidence level (default 0.95)
#' @return Prints statistics and invisibly returns a list
#' 
matched_or_immediate <- function(a, b, c, d, conf_level = 0.95) {
  
  # Number of pairs
  n_pairs <- a + b + c + d
  n_discordant <- b + c
  
  # Calculate odds ratio
  if (c == 0) {
    or <- Inf
    or_lower <- NA
    or_upper <- NA
    warning("Cannot calculate CI: c = 0")
  } else if (b == 0) {
    or <- 0
    or_lower <- NA
    or_upper <- NA
    warning("Cannot calculate CI: b = 0")
  } else {
    or <- b / c
    
    # Exact CI using F distribution
    alpha <- 1 - conf_level
    or_lower <- (b / (c + 1)) * (1 / qf(1 - alpha/2, 2 * (c + 1), 2 * b))
    or_upper <- ((b + 1) / c) * qf(1 - alpha/2, 2 * (b + 1), 2 * c)
  }
  
  # McNemar's chi-squared
  if (n_discordant > 0) {
    mcnemar_chi2 <- (b - c)^2 / (b + c)
    mcnemar_p <- pchisq(mcnemar_chi2, df = 1, lower.tail = FALSE)
    exact_p <- 2 * min(pbinom(b, n_discordant, 0.5), 
                       pbinom(b - 1, n_discordant, 0.5, lower.tail = FALSE))
    exact_p <- min(exact_p, 1)
  } else {
    mcnemar_chi2 <- NA
    mcnemar_p <- NA
    exact_p <- NA
  }
  
  # Proportions
  prop_cases <- (a + b) / n_pairs
  prop_controls <- (a + c) / n_pairs
  
  # Print output
  cat("\nMatched Case-Control Analysis (1:1 Matching)\n")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")
  
  cat("                       Cases\n")
  cat("Controls        Exposed    Unexposed    Total\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat(sprintf("Exposed      %8d   %8d   %8d\n", a, c, a + c))
  cat(sprintf("Unexposed    %8d   %8d   %8d\n", b, d, b + d))
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat(sprintf("Total        %8d   %8d   %8d\n", a + b, c + d, n_pairs))
  cat("\n")
  
  if (!is.na(mcnemar_chi2)) {
    cat(sprintf("McNemar's chi2(1) = %8.2f    Prob > chi2 = %.4f\n", 
                mcnemar_chi2, mcnemar_p))
    cat(sprintf("Exact McNemar probability  = %.4f\n\n", exact_p))
  }
  
  cat("Proportion with exposure\n")
  cat(sprintf("  Cases:     %.4f\n", prop_cases))
  cat(sprintf("  Controls:  %.4f\n\n", prop_controls))
  
  ci_label <- sprintf("%.0f%% CI", conf_level * 100)
  cat(sprintf("Odds ratio = %.4f    [%s: %.4f - %.4f] (exact)\n", 
              or, ci_label, or_lower, or_upper))
  
  invisible(list(
    odds_ratio = or,
    or_lower = or_lower,
    or_upper = or_upper,
    mcnemar_chi2 = mcnemar_chi2,
    mcnemar_p = mcnemar_p,
    exact_p = exact_p
  ))
}
