# Function to calculate odds ratio from 2x2 table
# Input: table with exposure as rows, outcome as columns
# Assumes: row 1 = unexposed, row 2 = exposed
#          col 1 = no outcome, col 2 = outcome

calculate_or <- function(tab) {
  # Extract cell counts
  # tab[1,1] = unexposed & no outcome (d)
  # tab[1,2] = unexposed & outcome (c)
  # tab[2,1] = exposed & no outcome (b)
  # tab[2,2] = exposed & outcome (a)
  
  a <- tab[2, 2]  # Exposed & outcome
  b <- tab[2, 1]  # Exposed & no outcome
  c <- tab[1, 2]  # Unexposed & outcome
  d <- tab[1, 1]  # Unexposed & no outcome
  
  # Calculate OR and 95% CI
  or <- (a * d) / (b * c)
  se_log_or <- sqrt(1/a + 1/b + 1/c + 1/d)
  ci_lower <- exp(log(or) - 1.96 * se_log_or)
  ci_upper <- exp(log(or) + 1.96 * se_log_or)
  
  # Calculate p-value using chi-square test
  chisq_result <- chisq.test(tab, correct = FALSE)
  p_value <- chisq_result$p.value
  
  # Print results
  cat("\nOdds Ratio Calculation\n")
  cat("----------------------\n")
  cat("Odds Ratio:", round(or, 3), "\n")
  cat("95% CI:", round(ci_lower, 3), "-", round(ci_upper, 3), "\n\n")
  cat("p-value:", format.pval(p_value, digits = 3), "\n\n")
  
  
  # Return invisibly
  invisible(list(
    or = or,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    p_value = p_value
  ))
}
