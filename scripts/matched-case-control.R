# Practical 13 matched case control
# 30/1/2026
# By the end of this practical, you should be able to:

# - Construct and interpret matched tables for paired case-control data
# - Calculate odds ratios from matched tables using discordant pairs
# - Understand why ignoring matching can bias results
# - Use conditional logistic regression to analyse matched data
# - Assess confounding in matched case-control studies

# The practical uses two datasets from Brazil investigating risk factors 
# for infant death from diarrhoea: diabraz.dta (1 control per case) 
# and diabraz2.dta (2 controls per case). 
# The exposure variables of interest are breastfeeding (bf), 
# water supply (wat2), and birthweight (bwtgp)

###################################
# library calls

# Load packages
library(haven)
library(tidyverse)
library(gtsummary)
library(here)
library(survival)  # for clogit()

###################################
# global variables

# Source custom functions
source(here("scripts/matched_functions.R"))
source(here("scripts/or_function.R"))
source(here("scripts/mh_functions_updated.R"))

# Set options
options(digits = 3, scipen = 999)

# Import data
diabraz <- read_stata(here("raw-data", "diabraz.dta"))

###################################
# script

# Part 1: Individual Matching (Paired Data)

## Question 1: Data Import and Exploration

# Explore the data
glimpse(diabraz)
summary(diabraz)

# variables are not labeled
# the variable "case" indicates whether the infant died from diarrhoea (1) or is a control (0), 
# "bf" indicates breastfeeding status, 
# "bwtgp" indicates birthweight group (1 = ≥3kg, 2 = <3kg), 
# and "set" identifies the matched pair.

# How many infants are in the dataset? 
# How many cases and controls? 
# What does the variable "set" refer to?

# Total observations
nrow(diabraz)

# Cases and controls
diabraz |> count(case)

# Check the matching structure
diabraz |> count(case, set) |> 
  pivot_wider(names_from = case, values_from = n)

# ANSWER: 172 infants in dataset, 86 cases and controls each (1:1 matching), 
# variable refers to death from diarrhea

## Questions 2-3: Matched Analysis for Breastfeeding

# In a matched case-control study, 
# we construct a matched table rather than a standard 2×2 table. 
# In a matched table for 1:1 matching, 
# each cell represents the number of pairs with a given combination 
# of case and control exposure status. 
# The rows represent the control's exposure status 
# and the columns represent the case's exposure status.

# recode breastfeeding variable
# to 0 and 1 (breastfed and not breastfed respectively)

# Recode bf: 1 = breastfed -> 0, 2 = not breastfed -> 1
diabraz <- diabraz |> 
  mutate(not_breastfed = case_when(
    bf == 1 ~ 0,
    bf == 2 ~ 1
  ))

# Check the recoding
diabraz |> count(bf, not_breastfed)

# The `matched_or()` function creates a table 
# where each cell represents paired data: 
# the number of pairs where the control has a given exposure status (rows) 
# and the case has a given exposure status (columns).

# Matched table and OR for breastfeeding
matched_or(diabraz, 
           case_var = "case", 
           exposure_var = "not_breastfed", 
           match_var = "set")
# Odds ratio = 4.8333 

## Question 4: Matched Analysis for Birthweight

# similar analysis for birth weight coded into 2 categories
# bweight is exposure of interest

# Recode bwtgp: 1 = >=3kg -> 0, 2 = <3kg -> 1
diabraz <- diabraz |> 
  mutate(low_birthweight = case_when(
    bwtgp == 1 ~ 0,
    bwtgp == 2 ~ 1
  ))

# Check
diabraz |> count(bwtgp, low_birthweight)

# Matched table and OR for birthweight
matched_or(diabraz, 
           case_var = "case", 
           exposure_var = "low_birthweight", 
           match_var = "set")
# Odds ratio = 1.3889 

## Question 5: What If We Ignored the Matching?

# What happens if we analyse these data as if they were unmatched? 
# Let's construct standard 2×2 tables and calculate unmatched odds ratios.

# Unmatched analysis for breastfeeding
tab_bf <- table(diabraz$not_breastfed, diabraz$case)
tab_bf
calculate_or(tab_bf)
# Odds Ratio: 3 
# this si quite different from the matched analysis (OR = 4.83)
# bias towards the null in unmatched
# bc matching was done on confounders, 
# and ignoring matching effectively ignores the confounding control


# Unmatched analysis for birthweight
tab_bwt <- table(diabraz$low_birthweight, diabraz$case)
tab_bwt
calculate_or(tab_bwt)

# Odds Ratio: 1.41 
# this is a little different from the matched analysis (1.39)

## Question 6: Conditional Logistic Regression

# conditional regression is for analyzing matched case control data
# It accounts for the matching by stratifying on the matched sets.
# `clogit()` function from the survival package
# which requires specifying the matching variable in `strata()`.
# The key difference from standard logistic regression 
# is that conditional logistic regression conditions 
# on the total number of cases in each stratum (matched set). 
# This means the model does not estimate an intercept, 
# and the matching variables themselves cannot be included as covariates!!!
# since their effects are absorbed by the stratification. !!!

# Conditional logistic regression for breastfeeding
clogit(case ~ not_breastfed + strata(set), data = diabraz) |> 
  tbl_regression(exponentiate = TRUE)
# OR = 4.83 -> matches with OR from matched table analysis
# for 1:1 matched data with single binary exposure,
# the matched table analysis and 
# conditional logistic regression give identical results

# Conditional logistic regression for birthweight
clogit(case ~ low_birthweight + strata(set), data = diabraz) |> 
  tbl_regression(exponentiate = TRUE)

# Again, the conditional logistic regression OR matches the matched table OR.

# Part 2: Multiple Controls per Case

## Question 7: A More Complex Dataset

# Now load diabraz2, which has 2 controls matched to each case. 
# This is called 1:2 matching.

# Import data with 2 controls per case
diabraz2 <- read_stata(here("raw-data", "diabraz2.dta"))

# Explore
glimpse(diabraz2)

# Check case/control counts
diabraz2 |> count(case)

# Check matching structure
diabraz2 |> 
  group_by(set) |> 
  summarise(n_cases = sum(case == 1),
            n_controls = sum(case == 0)) |> 
  count(n_cases, n_controls)

## Question 8: Analysis with Multiple Controls

# matched table becomes more complex with multiple controls per case
# However, conditional logistic regression handles any matching ratio seamlessly

# Note that in diabraz2, breastfeeding is coded: 
# 1 = breastfed, 2 = not breastfed. We need to recode appropriately.

# Recode bf in diabraz2: 1 = breastfed -> 0, 2 = not breastfed -> 1
diabraz2 <- diabraz2 |> 
  mutate(not_breastfed = case_when(
    bf == 1 ~ 0,
    bf == 2 ~ 1
  ))

# Check
diabraz2 |> count(bf, not_breastfed)

# matched table and OR 
matched_or(diabraz2, 
           case_var = "case", 
           exposure_var = "not_breastfed", 
           match_var = "set")
# Odds ratio = 2.1818
# with matching for blank, there is a 2.18 times greater odds of being a case
# in the not breastfed group compared to the breastfed group

## Question 13: Assessing Confounding

# seeing if social class or maternal education confound the association
# between breastfeeding and diarrhea mortality
# in matched case control studies we assess confounding by adding
# potential confounders to the conditional logistic regression model
# and examine whether the exposure OR changes

# Check the variables
diabraz2 |> count(social)
diabraz2 |> count(meduc)

# Convert to factors
diabraz2 <- diabraz2 |> 
  mutate(
    social_f = as.factor(social),
    meduc_f = as.factor(meduc)
  )

# fit a model with both confounders together
# adjust for both variables if you think both variables are confounders
# since they both contribute to the confounding 
# of the exposure outcome association

# Unadjusted for both social class and maternal education
mod_unadj <- clogit(case ~ not_breastfed + strata(set), 
                    data = diabraz2)

mod_unadj |> 
  tbl_regression(exponentiate = TRUE)
# OR = 3.79

# Adjusted for both social class and maternal education
mod_adj <- clogit(case ~ not_breastfed + social_f + meduc_f + strata(set), 
                  data = diabraz2)

mod_adj |> 
  tbl_regression(exponentiate = TRUE)
# OR = 4.28
# social and/or meduc could be confounders

# compare across the models:

# Compare unadjusted and adjusted ORs side by side
tbl_merge(
  list(
    tbl_regression(mod_unadj, 
                   exponentiate = TRUE, 
                   include = not_breastfed,
                   label = list(not_breastfed = "Not breastfed")),
    tbl_regression(mod_adj, 
                   exponentiate = TRUE, 
                   include = not_breastfed,
                   label = list(not_breastfed = "Not breastfed"))
  ),
  tab_spanner = c("**Unadjusted**", "**Adjusted**")
)

# What do you conclude?
# ANSWER: the confounders weakened the true association
# between not breastfeeding and death by diarrhea
# and the OR increases when adjusting for the confounders
# in the conditional logistic regression model

