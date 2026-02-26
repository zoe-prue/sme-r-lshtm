# Practical 10 Logistic Regression 10
# 29/1/2026
# # (i)	use the logistic/logit command to estimate the effect (odds ratio) 
# of one exposure, controlled for the effect of a second variable, 
# and obtain a confidence interval for this odds ratio
# (ii)	use logistic regression to assess if the effect of an exposure 
# (e.g. visual impairment) is confounded by another variable (e.g. age)
# (iii)	use the estimates and lrtest commands to perform a statistical test (LRT) 
# of the association between an exposure and the outcome after controlling 
# for the effects of confounding variables.

###################################
# library calls

library(haven)
library(tidyverse)
library(gtsummary)
library(here)
library(lmtest)

setwd("~/Desktop/sme-r-lshtm")

# Source custom functions
source(here("scripts/mh_functions_updated.R"))
source(here("scripts/or_function.R"))

###################################
# global variables

mortality <- read_dta(here("raw-data/mortality.dta"))
mortality <- mortality |> mutate(across(where(is.labelled), as_factor))

###################################
# script

glimpse(mortality)
summary(mortality)
View(mortality)


# 1. Does onchocercal infection appear to be associated with death?

# Cross-tabulate "died" and "mfpos"
# tbl_cross() creates cross-tabulations with row/column percentages

mortality |> 
  tbl_cross(row = died, col = mfpos, percent = "row")

#there appears to be an association between microfilarial load (mfpos) and death
# the percentages that died increase with mfpos

# Run a logistic regression model to estimate the crude OR for onchocercal infection.

# Logistic regression
glm(died ~ mfpos,
    data = mortality,
    family = binomial()) |>
  tbl_regression(exponentiate = TRUE) # unadjusted logistic regression

# OR = 1.63 with the logistic regression model, p=0.021

# Compare: calculate OR from 2x2 table
mort_table <- table(mortality$mfpos, mortality$died)
calculate_or(mort_table) # uses chi-squared test

# Odds Ratio: 1.633 
# p-value: 0.02 -> Wald p-value

# ANSWER: the odds of death in the mfpos group is 63% hgiehr than the noninfected group
# with strong association (p=0.02)

# 2. Let's check if age is a potential confounder. 
# Is our exposure (onchocerchal infection) somehow associated with age?
# Could this explain part or all of the association 
# between onchocercal infection and odds of death? 

# Cross-tabulate mfpos against agegrp
mortality |> 
  tbl_cross(row = agegrp, col = mfpos, percent = "column")

# ANSWER: the exposure seems to be associated with age
# as the percentage of mfpos in higher groupings increase with age
# this may be a part of the association between visual impairment and death
# but needs to be explored more to see if this is enhancing the OR we see in adjusted analysis

# 3. Now stratify by age group to check whether the association is changing.
# Is it reasonable to calculate a summary odds ratio 
# for the effect of onchocercal infection adjusted for age? 

# Stratified 2x2 tables
mortality |>
  filter(agegrp == "15-34") |> 
  tbl_cross(row = died, col = mfpos, percent = "column")

mortality |>
  filter(agegrp == "35-54") |> 
  tbl_cross(row = died, col = mfpos, percent = "column")

mortality |>
  filter(agegrp == "55-64") |> 
  tbl_cross(row = died, col = mfpos, percent = "column")

mortality |>
  filter(agegrp == "65+") |> 
  tbl_cross(row = died, col = mfpos, percent = "column")

# ANSWER: yes, the distribution of the percetnage who died vs mfpos being higher
# increases with age

# find summary estimate of assocaition between exposure and outcome, accounting for age
# First, let's do it non-parametrically with Mantel-Haenszel.

# Ensure variables are factors for mhor() function
mortality <- mortality |>
  mutate(
    died = as.factor(died),
    mfpos = as.factor(mfpos),
    agegrp = as.factor(agegrp)
  )

# Mantel-Haenszel OR and chi2 test
mhor(mortality, outcome = "died", exposure = "mfpos", strata = "agegrp")

# Stratum-specific odds ratios:
# agegrp                OR                95% CI
# -------------------------------------------------- 
# 15-34              2.311  (0.944 - 5.657)
# 35-54              1.327  (0.635 - 2.772)
# 55-64              1.769  (0.687 - 4.554)
# 65+                0.929  (0.381 - 2.264)
# M-H OR: 1.504 (95% CI: 0.996 - 2.270)

# ANSWER: is it reasonable to calcualte a summary OR for the association bw mfpos and died adjusting for age
# because the data indicates a difference in the distribution of deaths 
# bw mfpos and died among different age groups
# the adjusted OR loweres, meaning that age may have added to the association
# bw mfpos and death, and acted as a confounder

# 4.Now let's do the same but with logistic regression.
# Is it reasonable to calculate a summary odds ratio for the effect of onchocercal infection adjusted for age? '
# What do you conclude about whether age group is an effect modifier or confounder 
# of the association between onchocercal infection and odds of death?

# Logistic regression without
glm(died ~ mfpos, 
    data = mortality,
    family = binomial()) |> 
  tbl_regression(exponentiate = TRUE)
# OR = 1.63

# Logistic regression with a confounder
glm(died ~ mfpos + agegrp, # the + is for adding agegrp as a confounder
    data = mortality,
    family = binomial()) |> 
  tbl_regression(exponentiate = TRUE)
# OR = 1.51

# ANSWER: yes we can calculate an adjusted OR bc the crude differed from the pooled M-H OR
# age seems to be an effect modifier bc logistic regression OR differs greatly between age strata
# when representing the association between mfpos and death

# 6 + 7. THIS SWAPS THE EXPOSURE TO THE VIMP
# Now let's run a more complex model, 
# that also includes visual impairment (vimp). 
# Could onchocercal infection appear to confound the association 
# between visual impairment and odds of death? 
# What is the likely confounding structure?

# Logistic regression with two confounders
glm(died ~ vimp + mfpos + agegrp, # two potential confounders/effect modifiers are mfpos and agegrp
    data = mortality,
    family = binomial()) |> 
  tbl_regression(exponentiate = TRUE)

glm(died ~ vimp  + agegrp, # two potential confounders/effect modifiers are mfpos and agegrp
    data = mortality,
    family = binomial()) |> 
  tbl_regression(exponentiate = TRUE)
# use this model above to see if the OR for vimp changes -> evidence for the additional covariate to be a confounder


# NOTE: CAN GET THE LOG ODDS BY JUST DOING THIS: look a the intercept
glm(died ~ vimp  + agegrp, # two potential confounders/effect modifiers are mfpos and agegrp
    data = mortality,
    family = binomial()) #|> 
  # tbl_regression(exponentiate = TRUE)

# ANSWER: mfpos appears to confound the relationship between vimp and death
# due to the OR differing in the strata from the crude OR
# and the logistic regression ORs are different from each other (1 vs 1.46)
# the age group also appears to be an effect modifier for similar reasons.
# this model seems to be able to shift around the research question
# the model does not know wha tthe question is!
# and seeing if when treating different covariates as the exposure
# what the adjusted ORs in the strata are
# testing assumptions in causal graph - providing evidence to support your DAG
# wouldn't use these ORs in reality

# 8. Perform a likelihood ratio test of the H0 that,
# after accounting for age and onchocercal infection, 
# there is no association between visual impairment and odds of death. 
# The lrtest() function from the lmtest package compares nested models 
# by testing whether the more complex model provides a better fit to the data
# What do you conclude? 
# What other statistical test could you have performed of this hypothesis? 
# Does it produce a similar result?

# Create the first model with all confounders: assessing DIRECT EFFECT
model_adj <- glm(died ~ vimp + mfpos + agegrp, # vimp is the exposure
                 data = mortality,
                 family = binomial()) # almost like bypassing effect of these two confounders

# Create a simpler regression without the second potential confounder (vimp): assessing TOTAL EFFECT
model_noexp <- glm(died ~ mfpos + agegrp, # removal of vimp as an exposure
                   data = mortality,
                   family = binomial()) # this is a special case of the other model

# Likelihood ratio test
lrtest(model_adj, model_noexp)

# ANSWER: statistically significant evidence that the second model is showing that
# there is an association bw vimp and death after controlling for the two confounders
# seems like the model witht he exposure COULD be a better fit
# comparing a model with the exposure (vimp) and without
# the LRT test gives you evidence whether the exposure is associated with death
# otherwise, the p-value would be higher

# Part 2: Mwanza (will work with missing data here)

## Data import, exploration & management

# Make sure you have the mwanza.dta dataset in the correct folder, and load it. 
# It contains data on HIV infection among women in Mwanza, Tanzania.

# Import the dataset
mwanza <- read_dta(here("raw-data/mwanza.dta"))

View(mwanza)
glimpse(mwanza)
summary(mwanza)

# change ed to a binary instead of 4 categories

# Tabulate all possible values of ed
mwanza |> count(ed)

# Recategorise and label education level
mwanza <- mwanza |>
  mutate(ed2 = as.factor(
    case_when(
      ed == 1 ~ "None",
      ed == 2 ~ "Any formal",
      ed == 3 ~ "Any formal",
      ed == 4 ~ "Any formal"
    )
  ))

# Check it worked 
mwanza |> count(ed, ed2)

# 11. Cross-tabulate binary education against HIV infection and get the crude OR.

# 2x2 table with crude OR
mwanza |> 
  tbl_cross(row = case, col = ed2)

# Calculate crude OR
mwanza_table <- table(mwanza$ed2, mwanza$case)
calculate_or(mwanza_table)

# 12. We need to use the religion variable, 
# but let's make sure that there are no missing values first.
# Is religion a confounder?

mwanza |> count(rel)

# Replace rel "9" with "NA"
mwanza$rel <- na_if(mwanza$rel, 9)

# 13. Estimate the OR for HIV/education adjusted for religion, 
# using the Mantel-Haenszel approach

# Ensure variables are factors for mhor() function
mwanza <- mwanza |>
  mutate(
    case = as.factor(case),
    ed2 = as.factor(ed2),  # Already factor but explicit is safer (!)
    rel = as.factor(rel)
  )

# M-H OR
mhor(mwanza, outcome = "case", exposure = "ed2", strata = "rel")

# Crude odds ratio:
# OR: 0.413 (95% CI: 0.287 - 0.594)
# Chi-squared(1) = 23.447, P-value: 0.0000

# M-H combined odds ratio:
# OR: 0.522 (95% CI: 0.353 - 0.773)
# M-H Chi-squared(1) = 10.893, P-value: 0.0010

# Test of homogeneity of ORs (Breslow-Day):
# Chi-squared(3) = 1.046
# P-value = 0.7901

# 14. lets do the same for logistic regression

# Create a model controlling for religion
model_rel <- glm(case ~ ed2 + rel,
                 data = mwanza,
                 family = binomial())
tbl_regression(model_rel, exponentiate = TRUE)


# 15. Perform a LRT on the null hypothesis that, 
# after controlling for religion, 
# there is no association between education and HIV status.

# Create an unadjusted model excluding the missing value in the religion variable
model_unadj <- mwanza |>
  drop_na(rel) |>
  glm(case ~ ed2,
      data = _,
      family = binomial())

# LRT
lrtest(model_rel, model_unadj)







