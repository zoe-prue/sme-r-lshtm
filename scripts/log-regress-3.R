# Logistic regression 3 practical 11
# 29/1/2026
# PURPOSE

###################################
# library calls

library(haven)
library(tidyverse)
library(gtsummary)
library(here)
library(lmtest)
library(marginaleffects)
library(emmeans)

###################################
# global variables

setwd("~/Desktop/sme-r-lshtm")

# Source custom functions
source(here("scripts/mh_functions.R"))
source(here("scripts/or_function.R"))

mortality <- read_dta(here("raw-data/mortality.dta"))
mortality <- mortality |> mutate(across(where(is.labelled), as_factor))

###################################
# script

## Data import, exploration, management

glimpse(mortality)
summary(mortality)
View(mortality)

mortality <- mortality |>
  mutate(died = factor(
    died,
    levels = c(0, 1),
    labels = c("alive", "dead")))

## Data analysis

# 1. Analyse the association between visual impairment and death, 
# stratifying by sex, using a Mantel-Haenszel approach.
# The tbl_strata() function creates separate tables 
# for each level of a grouping variable.
# In this case, we're asking it to split the mortality data by sex, 
# and then create a separate cross-tabulation for men and women.

# 2x2 table
mortality |> 
  tbl_cross(row = died, col = vimp, percent = "column")

# Unadjusted OR
mort_table <- table(mortality$vimp, mortality$died)
calculate_or(mort_table)

# Stratified 2x2 tables using gtsummary
mortality |>
  tbl_strata(
    strata = sex,
    .tbl_fun = ~ .x |>
      tbl_cross(row = died, col = vimp, percent = "column"))

# Ensure variables are factors for mhor() function
mortality <- mortality |>
  mutate(
    died = as.factor(died),
    vimp = as.factor(vimp),
    sex = as.factor(sex))

# Stratum-specific OR and MH-OR
mhor(mortality, outcome = "died", exposure = "vimp", strata = "sex")

# 2. Fit logistic regression models to estimate the same association, 
# without interaction and with an interaction.
# 
# Note: in the model with interaction, 
# the interpretation of adjusted OR changes 
# - it's a *stratum-specific* OR 
# (in the baseline stratum of the other covariate)

# Model without interaction # "mirrors" just the adjusted M-H OR
glm(died ~ vimp + sex,
    data = mortality,
    family = binomial()) |>
  tbl_regression(exponentiate = TRUE)

# Model with interaction # mirrors M-H stratum-specific ORs
glm(died ~ vimp * sex, # the asterisk marks the interaction
    data = mortality,
    family = binomial()) |>
  tbl_regression(exponentiate = TRUE)

model_int_sex <- glm(died ~ vimp * sex,
                     data = mortality,
                     family = binomial())

# What is the OR for visual impairment adjusted for sex? 
# OR = 5.53

# What is the OR for visual impairment *in men*? 
# OR = 3.95

# What is the OR for women (vs men) *in the visually unimpaired*? 
# OR = 0.77

# What is the interaction term?
# 2.04

# vii)	Compare your answers with the results of the M-H analysis 
# of the association between visual impairment and death, stratifying on sex.
# the adjusted M-H OR was 5.429
# the adjusted Or from logistic regress was 5.53 (similar)
# the stratum specific OR for the male baseline (for example) mortality in vimp vs normal
# was 3.946 in M-H 3.94 in the logistic regression OR
# and it was 

# 3. Calculate the other stratum-specific OR (for women) 
# by making women the baseline group.

mortality <- mortality |>
  mutate(sex2 = fct_relevel(sex, "Female"))

glm(died ~ vimp * sex2,
    data = mortality,
    family = binomial()) |>
  tbl_regression(exponentiate = TRUE)

### Reset
# mortality <- mortality |>
# mutate(sex2 = fct_relevel(sex, "Male"))

# 4. Fit a logistic regression model with interaction 
# between visual impairment and sex, and adjust for age

glm(died ~ vimp * sex + agegrp,
    data = mortality,
    family = binomial()) |>
  tbl_regression(exponentiate = TRUE)

model_int_age <- glm(died ~ vimp * sex + agegrp,
                     data = mortality,
                     family = binomial())

# i)	The age-adjusted odds ratio for visual impairment in males = 
# 1.57

# ii)	What we would interpret as the age-adjusted odds ratio for females 
# versus males in visually unimpaired individuals, if sex was the exposure =
# 0.88

# iii)	The statistical interaction parameter (odds ratio scale) 
# between visual impairment and sex =
# 2.04

# iv)	The age-adjusted odds ratio for visual impairment in females =
# 1.57 * 2.04 = 3.20

# v)	What impact has controlling for age had on your opinion 
# of whether or not the effect of visual impairment differs for males and females? 
# Controlling for age has resulted in no important change 
# in the statistical interaction parameter and very little change 
# in the strength of evidence for effect modification. 
# After controlling for age, there remains some, 
# relatively weak evidence that the effect of visual impairment is modified by sex.  

# 5. Calculate stratum-specific ORs for the interaction model, using emmeans.
# compares odds ratios bw visual impairment groups, separately for men and women
# can stratify by pairwise and revpairwise depending on which way we want to calc comparison

# Get stratum specific ORs
emmeans(model_int_age,
        revpairwise ~ vimp | sex,
        type = "response") |>
  pluck("contrasts")

# Part 2: Mwanza

## Data import, exploration & management

# Import the dataset
mwanza <- read_dta(here("raw-data/mwanza.dta"))

glimpse(mwanza)
summary(mwanza)

# 6. Use logistic regression to assess whether the association 
# between education (ed) and HIV (case) is modified by age (age1). 
# First, relevel these variables:
# - ed into two groups: "None" and "Any formal education"
# - age1 into three groups: 15-24, 25-34, 34+

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

mwanza |> count(ed, ed2)

mwanza <- mwanza |>
  mutate(age3 = as.factor(
    case_when(
      age1 == 1 ~ "15-24",
      age1 == 2 ~ "15-24",
      age1 == 3 ~ "25-34",
      age1 == 4 ~ "25-34",
      TRUE ~ "34+"
    )
  ))

mwanza |> count(age1, age3)

# Now use MH and logistic regression to assess for statistical interaction. 
# Is there evidence of statistical interaction between age and education?

# Undjusted OR
mwanza_table <- table(mwanza$ed2, mwanza$case)
calculate_or(mwanza_table)

# Ensure variables are factors for mhor() function
mwanza <- mwanza |>
  mutate(
    case = as.factor(case),
    ed2 = as.factor(ed2),
    age3 = as.factor(age3)
  )

# Stratum-specific OR and MH-OR
mhor(mwanza, outcome = "case", exposure = "ed2", strata = "age3")

# Logistic regression without interaction
mod_noint <- glm(case ~ ed2 + age3,
                 data = mwanza,
                 family = binomial())

# Logistic regression with interaction
mod_int <- glm(case ~ ed2 * age3,
               data = mwanza,
               family = binomial())

lrtest(mod_noint, mod_int)

# From the logistic regression analysis, 
# is there evidence of effect modification by age 
# on the association between education and HIV?

# The crude and MH odds ratios barely differ 
# in measuring the association between education and HIV
# crude OR = 0.414, and M-H OR = 0.429
# M-H provides strong evidence for association between education and HIV nonethless
# due to it being far from the null

# The LRT for effect modification indicates strong evidence 
# against the null hypothesis of no effect modification (P=0.006785)

# there was strong evidence that not having been to school was protective against HIV
# (unadjusted OR = 0.414)
# stratification with MH on age showed strong evidence for effect modification by age
# among the exposure (education) to HIV association (p=0.0051)
# with the association between schooling and HIV being strongest 
# among women aged 25 to 34 (odds ratio = 0.205), 
# with no evidence for an association observed in women aged less than 25 (OR=1.029)

################## notes ###################################

# confounding is a property of the real world, not the data!
# effect modification can be tested in logistic regression models with LRT
# with effect modification present, we report s-s ORs
# without effect modification, we use pooled ORs



