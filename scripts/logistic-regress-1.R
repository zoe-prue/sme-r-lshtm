# PRactical 9: Logistic regression 1
# 23/1/2025
# By the end of this practical, students will be able to:

# - Use tables and odds ratio calculations for estimating odds 
# and odds ratios in a cohort study with fixed follow-up time
# - Use logistic regression to compare two or more groups in an unadjusted analysis
# - Analyze exposures with more than 2 categories using categorical variables
# - Use the likelihood ratio test to compare nested models
# - Assess confounding by comparing unadjusted and adjusted odds ratios

###################################
# library calls

# Load required packages
library(haven)
library(gtsummary)
library(lmtest)
library(here)
library(tidyverse)

###################################
# global variables

# Source OR calculation function
source(here("scripts", "or_function.R"))

# Set options for cleaner output
options(digits = 3, scipen = 999)

###################################
# script

# Data Import and Preparation

# Analyze this data set treating it as a cohort study with fixed follow-up time 
# (i.e., assume that all individuals were followed 
# for the same length of time and analyze the data using odds)

# Import data
mortality <- read_stata(here("raw-data", "mortality.dta")) |> 
  mutate(across(where(is.labelled), as_factor))

# Examine structure
glimpse(mortality)

# Question 3: Examine the Association Between Visual Impairment and Death #

# Cross-tabulation with row percentages
mortality |> 
  count(vimp, died) |> 
  group_by(vimp) |> 
  mutate(percentage = n / sum(n) * 100,
         n_pct = sprintf("%d (%.1f%%)", n, percentage)) |> 
  select(vimp, died, n_pct) |> 
  pivot_wider(names_from = died, values_from = n_pct) 

# Simple 2x2 table
tab_vimp <- table(mortality$vimp, mortality$died)
tab_vimp

# Calculate odds ratio
calculate_or(tab_vimp)

# the odds of death in the visually impaired group 
# is 456.6% higher than in the unimpared group
# Odds Ratio: 5.566 
# 95% CI: 3.779 - 8.199 

# Question 4: Logistic Regression for Visual Impairment #

# perform a logistic regresison examining the associatio bw visual impairment and death
# glm() command for logistic regression 
# (and other mdoels with similar linear modeling framework)
# (left side has outcome, adn right side has exposure and covariates)
# formula says 'predict the probability of the outcome died from the exposure vimp'
# family argument tells R we want a logistic regression 
# assuming a binomial distribution for binary outcome

# Logistic regression
mod_vimp <- glm(died ~ vimp, family = binomial, data = mortality)

# Display with odds ratios and 95% CIs
mod_vimp |> 
  tbl_regression(
    exponentiate = TRUE, 
    label = list(vimp = "Visual impairment")) |>
  add_global_p(method="lrt")

# or, alternatively,
# lrtest(mod_vimp)

# from this output we get:

# - An OR with 95% CI for outcome in exposed vs unexposed
# - P-value for a Wald test (null hypothesis: the vimp coefficient = 0) 
# - P-value for an LR test (null hypothesis: the model fits no better than the null model)
# - The number of observations
#
# The baseline term (intercept) represents the odds of outcome 
# in the baseline group of the variable 
# – so in this example it is the odds of death amongst those visually unimpaired.
# 
# The confidence interval is derived using the standard error 
# for the log odds ratio, as shown in the lecture.

# Question 8: Explore Microfilarial Infection and Death #

# examining the association between microfilarial infection in 4 levels and death
# are column or row percentages used in your table?

# Summary of categorical data
mortality |> count(mfgrp)

# Cross-tabulation with row percentages
mortality |> 
  count(mfgrp, died) |> 
  group_by(mfgrp) |> 
  mutate(percentage = n / sum(n) * 100,
         n_pct = sprintf("%d (%.1f%%)", n, percentage)) |> 
  select(mfgrp, died, n_pct) |> 
  pivot_wider(names_from = died, values_from = n_pct) 

# How would you alter the code above to drop the 
# missing rows missing mfgrp information?
# add an na.rm() = TRUE or something like it to the pipeline before mutate

# Question 9: Logistic Regression with Categorical Exposure #

# Check the type of the categorical exposure is a factor
class(mortality$mfgrp)

# Logistic regression with categorical exposure
mod_mfgrp <- glm(died ~ mfgrp, family = binomial, data = mortality)

mod_mfgrp |> 
  tbl_regression(
    exponentiate = TRUE,
    label = list(mfgrp = "Microfilarial load"))

# Note that there are three odds ratios 
# each of which refers to the same baseline group (those uninfected). 
# The odds ratio is:
#   
# - 1.69 for those with microfilarial load <10 mf/mg compared to those uninfected
# - 1.46 for those with microfilarial load 10-49 mf/mg compared to those uninfected
# - 2.05 for those with microfilarial load ≥50 mf/mg compared to those uninfected

# these are Wald test p-values, one for each OR
# tests if the OR is different to 1
# likelihood ratio statistic tests the association 
# between the infection variable and death
# by simultaneously testing all 3 parameters in the model
# the estimates of the log odds ratios 
# for microfilarial load groups 1, 2, 3 all vs zero

# Question 13: Likelihood Ratio Test with Missing Data #

# doing a likelihood ratio test for mfgrp (the infection load)
# compare log likelihood from the model with mfgrp (L1)
# and the log likelihood from the model without mfgrp (L0)

# asking:
# Does adding the extra data actually improve the fit of your model overall? 
# Does the exposure explain the outcome better than nothing

# to perform a likelihood ratio test (LRT) you need to use lrtest()
# to compare the log likelihood from a logistic regression model 
# with the variable of interest (that you have already defined)
# and the log likelihood from a logistic regression model without the variable

# Caution: mfgrp has missing data, 
# so the null model would have more observations than the full model
# and the LRT test can only work 
# when the two models have the same number of observations.
# using the `drop_na()` command helps with this to remove rows with NAs for specified vars

# Check for missing values
mortality |> count(mfgrp)

# Null model on complete cases
# The _ placeholder tells R where to put the piped data
mod_0 <- mortality |> 
  drop_na(mfgrp) |> 
  glm(died ~ 1, family = binomial, data = _)

# Full model (missing values automatically removed)
mod_1 <- glm(died ~ mfgrp, family = binomial, data = mortality)

# Likelihood ratio test
lrtest(mod_0, mod_1)

# Likelihood ratio test
# 
# Model 1: died ~ 1
# Model 2: died ~ mfgrp
# #Df  LogLik Df  Chisq Pr(>Chisq)  
# 1   1 -590.21                       
# 2   4 -586.87  3 6.6971     0.0822 .

# the null hypothesis is that the log-likelihood 
# of the two models is the same 
# (interpreted as the model fit being the same) 
# quantifies the strength of evidence that the log likelihood
# is further maximized by including the covariates of interest

# Question 14: Age and Death, Adjusting for Confounding #

# Logistic regression for age
mod_age <- glm(died ~ agegrp, family = binomial, data = mortality)

mod_age |> 
  tbl_regression(
    exponentiate = TRUE,
    label = list(agegrp = "Age group"))

# fit a model with both vimp and agegrp as explanatory variables

# Logistic regression with exposure and covariate
mod_adjusted <- glm(died ~ vimp + agegrp, family = binomial, data = mortality)

# Compare unadjusted and adjusted models side by side
tbl_merge(
  list(
    tbl_regression(mod_vimp, 
                   exponentiate = TRUE, 
                   include = vimp,
                   label = list(vimp = "Visual impairment")),
    tbl_regression(mod_adjusted, 
                   exponentiate = TRUE, 
                   include = vimp,
                   label = list(vimp = "Visual impairment"))
  ),
  tab_spanner = c("**Unadjusted**", "**Adjusted for age**")
)
  
# have the coefficients changed from those 
# in the models with each variable alone?
# For only the visually impaired non-adjusted model:
# visual impaired has odds of death of 5.57 compared to 1 for normal

# for the model with BOTH included as explanatory variables (adjusted)
# the odds of death in the visually impaired group is 2.2
# that of the odds of dying in the normal group when controlling for age

# likely to be srong confounding by age -> 
# treating as a covariate reduced the association 
# between visual impaoirment and death

# NOTE: we are not looking at coefficients for age bc our causal question
# is about visual impairment and death
# would need to draw a whole new DAG to explore that question
# reasoning: visual impairment cannot confound the relationship between age and death
# (surprisingly common!) Table 2 fallacy.










