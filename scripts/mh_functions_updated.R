# Mantel-Haenszel and Dose-Response Functions
# Mimics Stata's mhodds and tabodds commands

mhor <- function(data, outcome, exposure, strata = NULL) {
  
  # Remove missing values
  if (is.null(strata)) {
    data_complete <- data |> 
      filter(!is.na(.data[[outcome]]), !is.na(.data[[exposure]]))
  } else {
    data_complete <- data |> 
      filter(!is.na(.data[[outcome]]), !is.na(.data[[exposure]]), 
             !is.na(.data[[strata]]))
  }
  
  # Get levels
  outcome_levels <- levels(data_complete[[outcome]])
  exposure_levels <- levels(data_complete[[exposure]])
  
  # Calculate crude OR
  crude_table <- table(data_complete[[outcome]], data_complete[[exposure]])
  
  a <- crude_table[2, 2]
  b <- crude_table[2, 1]
  c <- crude_table[1, 2]
  d <- crude_table[1, 1]
  
  crude_or <- (a * d) / (b * c)
  crude_se_log <- sqrt(1/a + 1/b + 1/c + 1/d)
  crude_ci_lower <- exp(log(crude_or) - 1.96 * crude_se_log)
  crude_ci_upper <- exp(log(crude_or) + 1.96 * crude_se_log)
  
  # Chi-squared test for crude association
  crude_chisq <- chisq.test(crude_table, correct = FALSE)
  crude_pvalue <- crude_chisq$p.value
  crude_chisq_stat <- crude_chisq$statistic
  
  # If no stratification, return crude results only
  if (is.null(strata)) {
    cat("\n")
    cat("Odds Ratio\n")
    cat("----------\n")
    cat(sprintf("Odds Ratio: %.3f\n", crude_or))
    cat(sprintf("95%% CI: %.3f - %.3f\n", crude_ci_lower, crude_ci_upper))
    cat(sprintf("Chi-squared(1) = %.3f\n", crude_chisq_stat))
    cat(sprintf("P-value: %.4f\n", crude_pvalue))
    cat("\n")
    
    return(invisible(list(
      or = crude_or,
      ci = c(crude_ci_lower, crude_ci_upper),
      chi2 = crude_chisq_stat,
      p_value = crude_pvalue
    )))
  }
  
  # Calculate stratum-specific ORs
  strata_levels <- levels(data_complete[[strata]])
  stratum_results <- tibble()
  
  # For M-H calculation
  mh_numerator <- 0
  mh_denominator <- 0
  
  # For homogeneity test components
  R <- 0; S <- 0; P <- 0; Q <- 0
  
  # For M-H chi-squared test
  mh_chisq_num <- 0
  mh_chisq_denom <- 0
  
  for (stratum_level in strata_levels) {
    stratum_data <- data_complete |> filter(.data[[strata]] == stratum_level)
    
    # Create 2x2 table for this stratum
    tab <- table(stratum_data[[outcome]], stratum_data[[exposure]])
    
    # Convert to numeric to avoid integer overflow
    a <- as.numeric(tab[2, 2])
    b <- as.numeric(tab[2, 1])
    c <- as.numeric(tab[1, 2])
    d <- as.numeric(tab[1, 1])
    n <- as.numeric(sum(tab))
    
    # Stratum-specific OR
    if (b > 0 && c > 0 && a > 0 && d > 0) {
      stratum_or <- (a * d) / (b * c)
      stratum_se_log <- sqrt(1/a + 1/b + 1/c + 1/d)
      stratum_ci_lower <- exp(log(stratum_or) - 1.96 * stratum_se_log)
      stratum_ci_upper <- exp(log(stratum_or) + 1.96 * stratum_se_log)
    } else {
      stratum_or <- NA
      stratum_ci_lower <- NA
      stratum_ci_upper <- NA
    }
    
    # M-H components
    mh_numerator <- mh_numerator + (a * d / n)
    mh_denominator <- mh_denominator + (b * c / n)
    
    # Robins-Breslow-Greenland variance components
    R <- R + (a * d / n)
    S <- S + (b * c / n)
    P <- P + ((a + d) * a * d / n^2)
    Q <- Q + ((b + c) * b * c / n^2)
    
    # M-H chi-squared components
    m1 <- a + b
    n1 <- a + c
    expected_a <- m1 * n1 / n
    mh_chisq_num <- mh_chisq_num + (a - expected_a)
    mh_chisq_denom <- mh_chisq_denom + (m1 * n1 * (n - m1) * (n - n1) / (n^2 * (n - 1)))
    
    # Store results
    stratum_results <- stratum_results |> 
      bind_rows(tibble(
        stratum = stratum_level,
        or = stratum_or,
        ci_lower = stratum_ci_lower,
        ci_upper = stratum_ci_upper
      ))
  }
  
  # M-H summary OR
  mh_or <- mh_numerator / mh_denominator
  mh_se_log <- sqrt((P / (2 * R^2)) + ((P + Q) / (2 * R * S)) + (Q / (2 * S^2)))
  mh_ci_lower <- exp(log(mh_or) - 1.96 * mh_se_log)
  mh_ci_upper <- exp(log(mh_or) + 1.96 * mh_se_log)
  
  # M-H chi-squared test for association
  mh_chisq_statistic <- mh_chisq_num^2 / mh_chisq_denom
  mh_pvalue <- 1 - pchisq(mh_chisq_statistic, 1)
  
  # Breslow-Day test for homogeneity
  bd_statistic <- 0
  
  for (stratum_level in strata_levels) {
    stratum_data <- data_complete |> filter(.data[[strata]] == stratum_level)
    tab <- table(stratum_data[[outcome]], stratum_data[[exposure]])
    
    # Convert to numeric to avoid integer overflow
    a <- as.numeric(tab[2, 2])
    b <- as.numeric(tab[2, 1])
    c <- as.numeric(tab[1, 2])
    d <- as.numeric(tab[1, 1])
    n <- as.numeric(sum(tab))
    
    m1 <- a + c
    n1 <- a + b
    
    # Solve for expected a under homogeneity
    A <- mh_or - 1
    B <- -(mh_or * (m1 + n1) + (n - m1 - n1))
    C <- mh_or * m1 * n1
    
    if (abs(A) < 1e-10) {
      expected_a <- -C / B
    } else {
      discriminant <- B^2 - 4*A*C
      if (discriminant >= 0) {
        expected_a <- (-B - sqrt(discriminant)) / (2*A)
      } else {
        next
      }
    }
    
    # Variance
    var_a <- 1 / (1/expected_a + 1/(m1 - expected_a) + 
                    1/(n1 - expected_a) + 1/(n - m1 - n1 + expected_a))
    
    bd_statistic <- bd_statistic + (a - expected_a)^2 / var_a
  }
  
  bd_df <- length(strata_levels) - 1
  bd_pvalue <- 1 - pchisq(bd_statistic, bd_df)
  
  # Print output
  cat("\n")
  cat("Mantel-Haenszel estimate of odds ratio\n")
  cat("=======================================\n\n")
  
  # Show crude OR first
  cat("Crude odds ratio:\n")
  cat(sprintf("  OR: %.3f (95%% CI: %.3f - %.3f)\n", 
              crude_or, crude_ci_lower, crude_ci_upper))
  cat(sprintf("  Chi-squared(1) = %.3f, P-value: %.4f\n\n", 
              crude_chisq_stat, crude_pvalue))
  
  # Show stratum-specific results
  cat("Stratum-specific odds ratios:\n")
  cat(sprintf("%-15s %8s  %20s\n", strata, "OR", "95% CI"))
  cat(strrep("-", 50), "\n")
  for (i in 1:nrow(stratum_results)) {
    cat(sprintf("%-15s %8.3f  (%5.3f - %5.3f)\n",
                stratum_results$stratum[i],
                stratum_results$or[i],
                stratum_results$ci_lower[i],
                stratum_results$ci_upper[i]))
  }
  
  # Show M-H combined estimate
  cat("\n")
  cat("M-H combined odds ratio:\n")
  cat(sprintf("  OR: %.3f (95%% CI: %.3f - %.3f)\n", 
              mh_or, mh_ci_lower, mh_ci_upper))
  cat(sprintf("  M-H Chi-squared(1) = %.3f, P-value: %.4f\n\n", 
              mh_chisq_statistic, mh_pvalue))
  
  # Show test for homogeneity
  cat("Test of homogeneity of ORs (Breslow-Day):\n")
  cat(sprintf("  Chi-squared(%d) = %.3f\n", bd_df, bd_statistic))
  cat(sprintf("  P-value = %.4f\n", bd_pvalue))
  cat("\n")
  
  # Return invisibly
  invisible(list(
    crude_or = crude_or,
    crude_ci = c(crude_ci_lower, crude_ci_upper),
    crude_chi2 = crude_chisq_stat,
    crude_p = crude_pvalue,
    stratum_specific = stratum_results,
    mh_or = mh_or,
    mh_ci = c(mh_ci_lower, mh_ci_upper),
    mh_chi2 = mh_chisq_statistic,
    mh_p = mh_pvalue,
    homogeneity_chi2 = bd_statistic,
    homogeneity_df = bd_df,
    homogeneity_p = bd_pvalue
  ))
}

#' Test for Trend in Odds Ratios (tabodds)
#' Mimics Stata's tabodds command: Uses the exact formula: χ² = U²/V
#' Cochran-Armitage test for trend
#' 
#' @param data A data frame
#' @param outcome Name of binary outcome variable
#' @param exposure Name of ordered exposure variable (numeric or factor)
#' @return Prints formatted output with ORs and score test for trend
#' @examples
#' tabodds(mwanza, "case", "npa2")

tabodds <- function(data, outcome, exposure) {
  
  # Remove missing values
  data_complete <- data |> 
    filter(!is.na(.data[[outcome]]), !is.na(.data[[exposure]]))
  
  # Get exposure levels and ensure consistent ordering
  if (is.factor(data_complete[[exposure]])) {
    # For factors, use existing factor levels
    exposure_levels <- levels(data_complete[[exposure]])
  } else {
    # For numeric/character, sort unique values
    exposure_levels <- sort(unique(data_complete[[exposure]]))
  }
  
  # Convert exposure to factor with consistent levels for both OR calculation and trend test
  data_complete <- data_complete |> 
    mutate(exposure_factor = factor(.data[[exposure]], levels = exposure_levels))
  
  # Reference is first level
  ref_level <- exposure_levels[1]
  
  # Calculate OR for each level vs reference
  results <- tibble()
  
  for (i in seq_along(exposure_levels)) {
    level <- exposure_levels[i]
    
    if (i == 1) {
      # Reference level
      results <- results |> 
        bind_rows(tibble(
          level = as.character(level),
          or = 1.000,
          ci_lower = NA,
          ci_upper = NA
        ))
    } else {
      # Create binary exposure: this level vs reference
      data_binary <- data_complete |> 
        filter(exposure_factor %in% c(ref_level, level)) |> 
        mutate(exposure_binary = factor(
          exposure_factor == level,
          levels = c(FALSE, TRUE)
        ))
      
      # Calculate OR
      tab <- table(data_binary[[outcome]], data_binary$exposure_binary)
      
      a <- tab[2, 2]
      b <- tab[2, 1]
      c <- tab[1, 2]
      d <- tab[1, 1]
      
      if (b > 0 && c > 0 && a > 0 && d > 0) {
        or <- (a * d) / (b * c)
        se_log <- sqrt(1/a + 1/b + 1/c + 1/d)
        ci_lower <- exp(log(or) - 1.96 * se_log)
        ci_upper <- exp(log(or) + 1.96 * se_log)
      } else {
        or <- NA
        ci_lower <- NA
        ci_upper <- NA
      }
      
      results <- results |> 
        bind_rows(tibble(
          level = as.character(level),
          or = or,
          ci_lower = ci_lower,
          ci_upper = ci_upper
        ))
    }
  }
  
  # Score test for trend (Cochran-Armitage exact formula: χ² = U²/V)
  # Create contingency table using the standardized factor
  tab <- table(data_complete[[outcome]], data_complete$exposure_factor)
  
  if (nrow(tab) != 2) {
    stop("Outcome must have exactly 2 levels for the score test")
  }
  
  # Identify which row is cases
  case_row <- grep("case|yes|1|true|disease|dead", 
                   rownames(tab), ignore.case = TRUE)
  if (length(case_row) == 0) {
    case_row <- 2  # default case to second row
  }
  control_row <- setdiff(1:2, case_row)
  
  a_i <- tab[case_row, ]      # cases at each exposure level
  c_i <- tab[control_row, ]   # controls at each exposure level
  n_i <- a_i + c_i            # total at each level
  
  M1 <- sum(a_i)              # total cases
  M0 <- sum(c_i)              # total controls
  T_total <- M1 + M0          # total sample size
  
  # Get numeric scores: use sequential integers 1, 2, 3, ... for each exposure level
  l_i <- seq_along(exposure_levels)
  
  # Calculate U statistic
  U <- (M1 * M0 / T_total) * (sum(a_i * l_i / M1) - sum(c_i * l_i / M0))
  
  # Calculate V (variance under null hypothesis)
  l_bar <- sum(n_i * l_i) / T_total  # weighted mean of scores
  
  V <- (M1 * M0 / (T_total * (T_total - 1))) * 
    sum(n_i * (l_i - l_bar)^2)
  
  # Score chi-squared
  trend_chi2 <- U^2 / V
  trend_p <- pchisq(trend_chi2, df = 1, lower.tail = FALSE)
  
  # Print output
  cat("\n")
  cat("Odds ratios for each level compared to reference\n")
  cat("=================================================\n\n")
  
  # Format the table manually
  cat(sprintf("%-15s %8s  %20s\n", exposure, "OR", "95% CI"))
  cat(strrep("-", 50), "\n")
  for (i in 1:nrow(results)) {
    if (is.na(results$ci_lower[i])) {
      # Reference level
      cat(sprintf("%-15s %8.3f  %s\n",
                  results$level[i],
                  results$or[i],
                  "(reference)"))
    } else {
      cat(sprintf("%-15s %8.3f  (%5.3f - %5.3f)\n",
                  results$level[i],
                  results$or[i],
                  results$ci_lower[i],
                  results$ci_upper[i]))
    }
  }
  
  cat("\n")
  cat("Score test for trend of odds:\n")
  cat(sprintf("  Chi-squared(1) = %.3f\n", trend_chi2))
  cat(sprintf("  P-value = %.4f\n", trend_p))
  cat("\n")
  
  # Return invisibly
  invisible(list(
    odds_ratios = results,
    trend_chi2 = trend_chi2,
    trend_p = trend_p
  ))
}