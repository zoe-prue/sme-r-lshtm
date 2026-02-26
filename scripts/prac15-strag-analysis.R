# Practical 15: strategies of analysis
# 5/2/2025
# PURPOSE:
# ‘Does visual impairment/blindness increase risk of death 
# in a rural Nigerian population and is this modified by educational status?’
# We would like you to:
# (i)	analyse the data using the approach laid out in the lecture to answer a causal question
# (ii)	outline the structure of a results section
# (iii)	draft tables of results summarising the findings
# (iv)	identify key points for discussion

# background of study:
# The study was conducted in a rural area of northern Nigeria 
# in the late 1980s to early 1990s. 
# This cohort study with fixed follow-up time included persons 
# aged >15 years who all underwent an eye examination by ophthalmic nurses. 
# Individuals were classified as visually impaired 
# according to the standard WHO definition 
# (visual acuity less than 6/18 in the better eye). 
# Communities were followed-up over a period of three years and deaths were identified. 

# workflow plan:
# 1. design causal graph (provided for us in practical questions)
# 2. identify possible confounders I will analyze
    # proxy for SES will be occupation
# 3. define analysis plan
    # minimal adjustment set?
# 4. data reduction (recoding of vairables)
    # ordered categorical? save as new column. justify why
# 5. cross tabulations
# 6. descriptive analyses
# 7. summary tables
# 8. crude measure of effect (calculate_or function)
# 9. stratified analysis (M-H on a couple of potential confounders?)
    # compare crude vs adjusted estimates
# 10. multivariable analysis using regression
    # 'fully' adjusted measure of effect (3 categorical confounders) (LRT)
    # testing age for linearity or general assocaition
    # repeating 'fully' adjusted model (2 categorical confounders, age linear) (LRT)
    # modifying effect of a covariate with 'fully adjusted' (if necessary)
        # testing for effect modification? (LRT? Wald?)
# 11. drafting a report of the findings (what to report is in lecture slides)

###################################
# library calls

library(haven)
library(gtsummary)    # For regression tables
library(marginaleffects)  # For predictions and contrasts
library(lmtest)       # For likelihood ratio tests
library(tidyverse)    # For data manipulation
library(here)
library(dplyr)
library(tidyr)
library(DescTools)
library(gmodels)

###################################
# global variables

# global variables
setwd("~/Desktop/sme-r-lshtm/")

# Read the Whitehall dataset
mort <- read_dta("raw-data/mortality.dta") |>
  mutate(across(where(is.labelled), as.factor))

# sapply(mortality, class)

source(here("scripts", "or_function.R"))
source(here("scripts", "mh_functions_updated.R"))

###################################
# script

# ALREADY COMPLETED
# 1. design causal graph (provided for us in practical questions)
    # think of DAG as population level hypothesis, and data as a sample

# 2. identify possible confounders I will analyze
    # proxy for SES will be occupation
    # age group
    # sex
    # education as an EFFECT MODIFIER

# ALREADY COMPLETED - FILL IN MORE DETAIL?
# 3. define analysis plan

# 4. data reduction (recoding of variables)
    # ordered categorical? save as new column. justify why

class(mort$occupation) # recoding occupation

mort <- mort |>
  mutate(occupation2 = factor( # defining when each occupation fits into the levels i've designed
      case_when(
        occupation %in% c(10, 11, 12, 20) ~ "No Income",
      occupation %in% c(1, 6, 7) ~ "Low Income",
      occupation %in% c(3, 4, 5) ~ "Medium Low Income",
      occupation %in% c(14) ~ "Medium Income",
      occupation %in% c(8, 9, 12, 15)  ~ "High Income"))
  )

table(unique(mort)$occupation2) 
class(mort$occupation2)

class(mort$agegrp) # checking coding on age group
table(unique(mort)$agegrp) 

class(mort$sex) # checking coding on sex
table(unique(mort)$sex) 

# 5. cross tabulations

CrossTable(mort$vimp, mort$died, prop.r = TRUE, chisq = TRUE)

# 6. descriptive analyses

# calculating follow-up time in years
mort <- mort |> 
  mutate(followup_years = as.numeric(exit - enter) / 365.25)

# checking follow-up times
summary(mort$followup_years)

# count total events and total person-years
total_deaths <- sum(mort$died) 
total_pyears <- sum(mort$followup_years)

# Total deaths: 137 
# Total person-years: 11457.4 
# Rate per 1000 person-years: 0.04 
# 95% CI: (0.03 , 0.04)

# 7. summary tables

# Calculate rates per 1000 person-years by hypertension status
rate_table <- mort |> 
  filter(!is.na(vimp)) |> 
  group_by(vimp) |> 
  summarise(
    deaths = sum(died),
    person_years = sum(followup_years),
    rate_per_1000 = (deaths / person_years) * 1000, # 1000 P-Y
    ci_lower = rate_per_1000 / exp(qnorm(0.975) * sqrt(1/deaths)),
    ci_upper = rate_per_1000 * exp(qnorm(0.975) * sqrt(1/deaths)),
    .groups = "drop"
  )

rate_table |> 
  knitr::kable(digits = 2,
               col.names = c("Visual Impairment", "Deaths", "Person-Years", 
                             "Rate per 1000", "95% CI Lower", "95% CI Upper"))

# 8. crude measure of effect (calculate_or function)

tab_vimp <- table(mort$vimp, mort$died)
tab_vimp
calculate_or(tab_vimp)

# results: 
# Odds Ratio: 5.566
# 95% CI: 3.779 - 8.199 
# p-value: <2e-16 

# 9. stratified analysis (M-H on each confounder)
    # compare crude vs adjusted estimates

# M-H for occupation2

# Mantel-Haenszel OR and chi2 test
mhor(mort, outcome = "died", exposure = "vimp", strata = "occupation2")

# results of M-H analysis:
# Crude odds ratio: OR: 5.436 (95% CI: 3.679 - 8.034) P-value: 0.0000

# Stratum-specific odds ratios:

# occupation2       OR     95% CI
# -------------------------------------------------- 
# High Income       10.965 (4.210 - 28.559)
# Low Income        4.151  (2.352 - 7.326)
# Medium Low Income 2.408  (0.536 - 10.816)
# No Income         1.741  (0.333 - 9.086)

# M-H adjusted: OR: 4.905 (95% CI: 2.720 - 8.846) P-value: 0.0000
# indicates confounding
# Test of homogeneity of ORs (Breslow-Day): P-value = 0.1283
# indicates no effect modification

# M-H for sex

# Mantel-Haenszel OR and chi2 test
mhor(mort, outcome = "died", exposure = "vimp", strata = "sex")

# results of M-H analysis:
# Crude odds ratio: OR: 5.566 (95% CI: 3.779 - 8.199) P-value: 0.0000

# Stratum-specific odds ratios:
# sex                OR     95% CI
# -------------------------------------------------- 
# 0                  3.946  (2.277 - 6.839)
# 1                  8.063  (4.659 - 13.954)

# M-H adjusted: OR: 5.429 (95% CI: 3.018 - 9.766) P-value: 0.0000
# indicates no confounding
# Test of homogeneity of ORs (Breslow-Day): P-value = 0.0684 
# indicates potential (weak) for effect modification

# M-H for age group

# Mantel-Haenszel OR and chi2 test
mhor(mort, outcome = "died", exposure = "vimp", strata = "agegrp")

# results of M-H analysis:
# Crude odds ratio: OR: 5.566 (95% CI: 3.779 - 8.199) P-value: 0.0000

# Stratum-specific odds ratios:
# agegrp             OR     95% CI
# -------------------------------------------------- 
# 0                  7.213  (1.621 - 32.107)
# 1                  2.697  (1.312 - 5.545)
# 2                  2.253  (0.988 - 5.137)
# 3                  1.417  (0.662 - 3.030)

# M-H adjusted: OR: 2.104 (95% CI: 1.291 - 3.429) P-value: 0.0005
# indicates confounding
# Test of homogeneity of ORs (Breslow-Day): P-value = 0.2132
# indicates no effect modification

# 10. multivariable analysis using regression
    # 'fully' adjusted measure of effect (3 categorical confounders) (LRT)

mort |> count(is.na(vimp))

# primary logistic regression model 
mod_vimp <- glm(died ~ vimp, family = binomial, data = mort)

# Display with odds ratios and 95% CIs
mod_vimp |> 
  tbl_regression(
    exponentiate = TRUE, 
    label = list(vimp = "Visual impairment")) |>
  add_global_p(method="lrt")

# similar to M-H OR, the OR = 5.57, (95% CI: 3.74, 8.14), p-value <0.001 (Wald test)
# the odds of death in the visually impaired group is 5.57 times that of the unimpaired group

# now doing the vimp with adjustment for each confounder

# is occupation a confounder?

mort |> count(is.na(occupation2)) # 9 NAs but i'm ignoring it

# Logistic regression with exposure and covariate
mod_adjust_occup <- glm(died ~ vimp + occupation2, family = binomial, data = mort)

# Compare unadjusted and adjusted models side by side
tbl_merge(
  list(
    tbl_regression(mod_vimp, 
                   exponentiate = TRUE, 
                   include = vimp,
                   label = list(vimp = "Visual impairment")),
    tbl_regression(mod_adjust_occup, 
                   exponentiate = TRUE, 
                   include = vimp,
                   label = list(vimp = "Visual impairment"))
  ),
  tab_spanner = c("**Unadjusted**", "**Adjusted for occupation**"))

# unadjusted and adjusted for occupation as a confounder
# OR for unadjusted is 5.57 (p<<0.001)
# OR for adjusted is 3.61 (P<<0.001)
# occupation seems to be a confounder with the logistic regression model

# is age a confounder?

mort |> count(is.na(agegrp))

# Logistic regression with exposure and covariate
mod_adjust_agegrp <- glm(died ~ vimp + agegrp, family = binomial, data = mort)

# Compare unadjusted and adjusted models side by side
tbl_merge(
  list(
    tbl_regression(mod_vimp, 
                   exponentiate = TRUE, 
                   include = vimp,
                   label = list(vimp = "Visual impairment")),
    tbl_regression(mod_adjust_agegrp, 
                   exponentiate = TRUE, 
                   include = vimp,
                   label = list(vimp = "Visual impairment"))
  ),
  tab_spanner = c("**Unadjusted**", "**Adjusted for age**"))

# OR for unadjusted is 5.57 (p<<0.001)
# OR for adjusted is 2.20 (P<<0.001)
# agegrp seems to be a confounder with the logistic regression model

# is sex a confounder?

mort |> count(is.na(sex))

# Logistic regression with exposure and covariate
mod_adjust_sex <- glm(died ~ vimp + sex, family = binomial, data = mort)

# Compare unadjusted and adjusted models side by side
tbl_merge(
  list(
    tbl_regression(mod_vimp, 
                   exponentiate = TRUE, 
                   include = vimp,
                   label = list(vimp = "Visual impairment")),
    tbl_regression(mod_adjust_sex, 
                   exponentiate = TRUE, 
                   include = vimp,
                   label = list(vimp = "Visual impairment"))
  ),
  tab_spanner = c("**Unadjusted**", "**Adjusted for sex**"))

# OR for unadjusted is 5.57 (p<0.001)
# OR for adjusted is 5.53 (P<<0.001)
# sex does not seem to be a confounder with the logistic regression model

# Logistic regression with all categorical confounders - trying to get unbiased estimate of relationship
model_w_vimp <- glm(died ~ vimp + agegrp + occupation2 + sex, 
    data = mort,
    family = binomial()) 

tbl_regression(model_w_vimp, exponentiate = TRUE)

# the vimp = 1 OR represents that in the model adjusting for both age group and occupation
# there is 1.85 times the odds of death in the visually impaired group 
# compared to the unimpaired when controlling for age group, sex, and occupation (p=0.011)

# Logistic regression with all confounders AND WITHOUT THE EXPOSURE

model_no_vimp <- glm(died ~ agegrp + occupation2 + sex, # two potential confounders/effect modifiers are mfpos and agegrp
    data = mort,
    family = binomial())

tbl_regression(model_no_vimp, exponentiate = TRUE)

# LRT between model with 3 categorical confounders and the exposure 
# vs model with confounders and w/out exposure

lrtest(model_w_vimp, model_no_vimp) 

# p-value = 0.01305
# this provides evidence that the model with the exposure 
# is a better fit than the model excluding the exposure

# 10. multivariable analysis using regression
    # testing age for linearity or general association
    # repeating 'fully' adjusted model (2 categorical confounders, age linear) (LRT)

# REPEATING LOGISTIC REGRESSION MODEL + LRT WITH AGE AS LINEAR

# first, is age a linear association?
# do i do this linear association compared to the outcome?????

class(mort$age)

model_cat_age <- glm(died ~ vimp + age + sex + occupation2, # this is the categorical model
                 data = mort,
                 family = binomial())

tbl_regression(model_cat_age, exponentiate = TRUE)

# Logistic regression with linear relation - log odds scale
model_cont_age <- glm(died ~ age + sex + occupation2,
                  data = mort,
                  family = binomial()) 
# OR scale
tbl_regression(model_cont_age, exponentiate = TRUE)

# LRT
lrtest(model_cat_age, model_cont_age)

# there is evidence (p=0.01994) that the model with the exposure is a better fit
# than a model with only the covariates, treating age as linear
# however, since this p-value is lower,
# treating age as an age group is likely a better model to pursue 
# (p=0.01305 for model w age group)
# i wonder if this is balanced by the continuous treatment being 1 parameter vs 3 parameters

# 10. multivariable analysis using regression
    # modifying effect of a covariate with 'fully adjusted' (if necessary)
        # testing for effect modification? (LRT? Wald?)

model_effect <- glm(died ~ vimp * agegrp, # REPEAT THIS FOR ALL 3 PLEASE
                     data = mort,
                     family = binomial())

tbl_regression(model_effect, exponentiate = TRUE)

model_effect <- glm(died ~ vimp  agegrp, # REPEAT THIS FOR ALL 3 PLEASE
                    data = mort,
                    family = binomial())

tbl_regression(model_effect, exponentiate = TRUE)


# 11. drafting a report of the findings (what to report is in lecture slides)

# lrt between model with 3 categorical confounders + expposure, and confoudners without exposure

# same as above with age being linear

# lrt p values

# departure from linearity and parameter number and categorical 

# linear age better? lrt p-value, number of parameters? rule of 10 - deaths in whole dataset and divide by 10

# df in lrt should change by the number of interaction terms added

# comment on multiple testing from effect modification

## do not use p-value to assess confounding









