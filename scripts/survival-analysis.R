# Survival Analysis
# 15/1/2026
# Welcome to SME Practical 2. In this session, 
# we'll compute and compare mortality rates using data from a 10% random sample of the 
# Whitehall Cohort Study of British civil servants.
# 
# This practical looks at:
# 
# - Rates instead of just counts 
# - Poisson regression for rate data
# - Stratified analysis to examine possible confounding
# - Effect modification (does the association differ across strata?)
# 
# The key epidemiological idea for this practical is that in cohort studies, 
# people are followed for different lengths of time, and we need to account for this when calculating rates.

###################################
# library calls

library(haven)
# install.packages("gtsummary")
library(gtsummary)    # For regression tables
# install.packages("marginaleffects")
library(marginaleffects)  # For predictions and contrasts
# install.packages("lmtest")
library(lmtest)       # For likelihood ratio tests
library(tidyverse)    # For data manipulation
# install.packages("here")
library(here)

# Set options for cleaner output
options(digits = 3, scipen = 999) # number of digits used for the whole script

###################################
# global variables

# Set wd to location of source files
setwd("~/Desktop/sme-r-lshtm/raw-data/")
# Read the Whitehall dataset
whitehall <- read_stata("whitehal.dta") |> 
  mutate(across(where(is.labelled), as_factor))

###################################
# script

# Examine the structure
glimpse(whitehall)

# Summary statistics
summary(whitehall)

# View first few rows
head(whitehall)

# Key variables:
  
# - `timein`, `timeout`: Entry and exit dates (stored as days since 1960-01-01)
# - `all`: All-cause mortality indicator (0 = survived, 1 = died)
# - `chd`: CHD mortality indicator (0 = no CHD death, 1 = CHD death)
# - `grade`: Employment grade (1 = high, 2 = low) - this will be a factor
# - `agein`: Age at entry in years
# 
# Each row is a person, and we know when they entered the study, 
# when they left (either by dying or by end of follow-up), and whether they died.

# Question 1: Overall Mortality Rate

## Calculate follow-up time

# We need to calculate person-years of follow-up. 
# Why are we using person-years? Because people weren't all followed for the same length of time. Imagine:
# 
# - Person A: Followed for 10 years
# - Person B: Followed for 5 years
# 
# If we just counted people (risk), both contribute equally to the measure of effect. 
# But Person A actually was at risk for double the length of Person B, and using person-years accounts for this.
# 
# Unlike Stata's `stset` command which does this automatically, R requires us to calculate person-years manually.

# Calculate follow-up in years
whitehall <- whitehall |> 
  mutate(followup_years = as.numeric(timeout - timein) / 365.25)

# Check the calculation
summary(whitehall$followup_years)

# CALCULATING OVERALL MORTALITY RATE - MANUAL CALCULATION #

# Step 1: count total events and total person-years
total_deaths <- sum(whitehall$all)
total_pyears <- sum(whitehall$followup_years)

# Step 2: Calculate rate per 1000 person-years
rate <- (total_deaths / total_pyears) * 1000

# Step 3: calculate rate per 1000 person-years
# SE on log scale: sqrt(1/deaths)
# log scale bc rates cant be negative
error_factor <- exp(qnorm(0.975)*sqrt(1/total_deaths))
lower_ci <- rate / error_factor
upper_ci <- rate * error_factor

# Step 4: show results in a printed way: cat is concatenate
cat("Total deaths:", total_deaths, "\n")
cat("Total person-years:", round(total_pyears, 1), "\n")
cat("Rate per 1000 person-years:", round(rate, 2), "\n")
cat("95% CI: (", round(lower_ci, 2), ",", round(upper_ci, 2), ")\n")

# 1.  We sum up all the deaths (the `all` variable where 1 = died)
# 2.  We sum up all the person-years we just calculated
# 3.  We divide deaths by person-years to get a rate 
# - but multiply by 1000 so it's "per 1000 person-years" 
# (an epidemiological convention that makes numbers easier to read)
# 4.  The confidence interval calculation uses the Poisson distribution 
# because counts of rare events (like deaths) follow this distribution
# 5.  We use `qnorm(0.975)` to get the z-score for a 95% CI (1.96) !!!

# CALCULATING OVERALL MORTALITY RATE - POISSON REGRESSION #

# Now let's do the same thing using regression. 
# This gives us the same answer but is more useful to us later 
# because we can easily add covariates 
# (like you did in STEPH with mutliple linear regression).

# Poisson regression is specifically designed for count data 
# (like deaths) when we have varying denominators (like person-years).

model_overall <- glm(all ~ 1 + offset(log(followup_years/1000)),
                     family = poisson,
                     data=whitehall)

# display the model
summary(model_overall)

# get rate with confidence interval

model_overall |> tbl_regression(exponentiate = TRUE, intercept = TRUE)

# since working with an intercept only model
# with tbl_regression, we include the 1000 in offset
# if you dont want to do this, you coudl alternatively use the code:

model_overall <- glm(all ~ 1 + offset(log(followup_years)), 
                     family = poisson, 
                     data = whitehall)

model_overall |> 
  broom::tidy(conf.int = TRUE) |> # get a tidy version of the coefficients
  mutate(across(c(estimate, conf.low, conf.high), \(x) exp(x) * 1000)) |> # exponentiate and multiply
  select(estimate, conf.low, conf.high) |> 
  knitr::kable(digits = 2,
               col.names = c("Rate per 1000", "95% CI Lower", "95% CI Upper"))

# What's happening in the code here:

# -   `family = poisson` tells R we're analyzing count data with a Poisson model
# -   `~ 1` means "overall" (no covariates yet)
# -   `offset(log(followup_years))` tells R that different people were followed for different lengths of time. 
# We use `log()` because Poisson regression works on the log scale internally
# -   `exponentiate = TRUE` converts from log scale back to the rate scale
# -   `intercept = TRUE` tells tbl_regression to use the intercept (since that's all we have)
# 
# Notice that Method 1 (by hand) and Method 2 (Poisson regression) give you the same rate.

# ANSWERING THE QUESTION #

# How many subjects died during follow up? 
# 403 subjects died during followup
# What was the total number of person years? 
# 27605 person-years total
# What was the longest period that a subject was followed up for?
# 19.381 years

# Question 2: Mortality Rates by Age Group

## MORTALITY RATES BY AGE GROUP ##

## Create age categories
# uses cut() function to divide the continuous data into age groups (categorical)
# `right = FALSE` means 40-49 includes 40 but not 45 (45 goes in the next group).

# Create age groups
whitehall <- whitehall |> 
  mutate(agecat = cut(agein, 
                      breaks = c(40, 45, 50, 55, 60, 65, 70),
                      right = FALSE))

# Check the distribution
table(whitehall$agecat)

# using emeans package to estimate predicted rates by age group
# tell R which model we are using, which variables to group_by, 
# and type = "response" to show we dont want transformations

# Fit Poisson model with age categories
model_age <- glm(all ~ agecat + offset(log(followup_years)), 
                 family = poisson, 
                 data = whitehall)

# Get modelled rates for each age group
margins_age <- emmeans(model_age, "agecat", type="response", offset=0)
margins_age

# Display as nicely printed table with rates per 1000 person-years
margins_age |> 
  as_tibble() |> 
  mutate(rate_per_1000 = rate * 1000, # to make this per 1000 person-years
         conf.low_1000 = asymp.LCL * 1000,
         conf.high_1000 = asymp.UCL  * 1000) |> 
  select(agecat, rate_per_1000, conf.low_1000, conf.high_1000) |> 
  knitr::kable(digits = 2, 
               col.names = c("Age Group", "Rate per 1000 PY", 
                             "95% CI Lower", "95% CI Upper"))

# rate ratios comparing age groups

tbl_regression(model_age, 
               exponentiate = TRUE,
               label = list(agecat = "Age at entry (years)"))

# in epidemiology, IRR stands for Incidence Rate Ratio, 
# a key measure comparing the rate of new disease or event occurrences 
# between an exposed group and an unexposed group, 
# calculated as the ratio of their incidence rates (events per person-time). 
# An IRR greater than 1 indicates an increased risk in the exposed group, 
# while an IRR less than 1 suggests a decreased risk, 
# helping researchers quantify the impact of exposures like smoking or pollutants 
# on health outcomes like cancer or mortality. 
# The output shows rate ratios comparing each age group to the reference (youngest) group (!!!)

# ANSWERING THE QUESTION #

# What do these rate ratios suggest about the relationship between age and all-cause mortality?
# there is an increase in the rate ratio with the increasing age group
# there is an association between age group increasing and the rate
# of disease increasing in association to age (exposure)

# Question 3: Mortality Rate by Employment Grade

# MORTALITY RATES BY EMPLOYMENT GROUP #

# Look at the grade distribution
table(whitehall$grade)

# Check what type of variable it is
class(whitehall$grade)

# Convert grade to a factor since it's numeric and not factor yet
whitehall <- whitehall |> 
  mutate(grade = factor(grade,
                        levels = c(1, 2),
                        labels = c("Low grade", "High grade")))

# Check the levels (which is reference? -> low grade is reference)
levels(whitehall$grade)

# descriptive statistics by grade

# Calculate deaths and person-years by grade
whitehall |> 
  group_by(grade) |> 
  summarise(n = n(),
            deaths = sum(all),
            person_years = sum(followup_years),
            rate_per_1000 = (deaths / person_years) * 1000) |> 
  knitr::kable(digits = 2)

# Fit Poisson model comparing grades
model_grade <- glm(all ~ grade + offset(log(followup_years)), 
                   family = poisson, 
                   data = whitehall)

# Display as a nice table
model_grade |> 
  tbl_regression(
    exponentiate = TRUE,
    label = list(grade ~ "Employment Grade"))

# -   The Estimate column shows the rate ratio (RR)
# -   RR >= 1 would mean low grade has a *higher* mortality rate than high grade
# -   RR =< 1 would mean low grade has a *lower* mortality rate
# -   RR = 1 would mean no difference

# The p-value tells us about the strength of evidence against the null hypothesis 
# that rates are the same in both grades (or equivalently - that there is no association)

# ANSWERING THE QUESTIONS #

# 7.Is there any evidence of effect modification between employment grade and age at entry? 
# Examine the result of the test for effect modification.
# The rate ratios are difference between strata
# low grade = 1, high grade = 2.31
# and the p-value for the test of homogeneity is P>0.001
# meaning there is strong evidence to reject the null hypothesis
# there there is no difference in the rate ratios of the two grades
# meaning there could be effect modification

# Question 4: Is Age a Confounder?

# EXAMINING AGE AS A CONFOUNDER POTENTIALLY #

# Now we need to think about confounding. 
# Confounding is a property of the real world, not of your data. 
# It occurs when there's a third variable (here, age) that is a common cause of exposure and outcome. 
# Age seems likely to be a confounder since it:

# 1.  Affects the outcome (mortality) - older people die more
# 2.  Affects the exposure (grade) - older people are likely to be of higher grade
# 
# If both are true in the real world, 
# then age could explain part or all of the association we see between grade and mortality.

# Create age categories
whitehall <- whitehall |> 
  mutate(agecat = cut(agein,
                      breaks = c(40, 50, 60, 70),
                      labels = c("40-49", "50-59", "60-69"),
                      right = FALSE))

# Check the distribution
table(whitehall$agecat)

## Check associations

# Remember that the investigations below are not definitively indicative of confounding!

# Age distribution by grade
whitehall |> 
  select(agein, agecat, grade) |> 
  tbl_summary(by = grade)

# Calculate mortality rates by age group
whitehall |> 
  group_by(agecat) |> 
  summarise(
    n = n(),
    deaths = sum(all),
    person_years = sum(followup_years),
    rate_per_1000 = (deaths / person_years) * 1000
  ) |> 
  knitr::kable(digits = 2)

# Poisson model by age
model_age <- glm(all ~ agecat + offset(log(followup_years)), 
                 family = poisson, 
                 data = whitehall)

model_age |> 
  tbl_regression(
    exponentiate = TRUE,
    label = list(agecat ~ "Age Group")
  )

# ANSWERING THE QUESTIONS #

# 8.Is the effect of grade confounded by age at entry?
# There is evidence of confounding by age
# mortality rates increase with age

## AGE-ADJUSTED RATE RATIO ##

## Age-adjusted rate ratio

# Now let's adjust for age by including it in the model. 
# This will give us the grade-mortality association within age groups:

# Model with both grade and age
model_grade_age <- glm(all ~ grade + agecat + offset(log(followup_years)), # log so its not negative rate
                       family = poisson, 
                       data = whitehall)

# Display results
model_grade_age |> 
  tbl_regression(
    exponentiate = TRUE,
    label = list(
      grade ~ "Employment Grade",
      agecat ~ "Age Group"
    )
  )

# Compare unadjusted and adjusted rate ratios
tbl_merge(
  list(
    tbl_regression(model_grade, 
                   exponentiate = TRUE,
                   include = grade,
                   label = list(grade = "Employment grade")),
    tbl_regression(model_grade_age, 
                   exponentiate = TRUE,
                   include = grade,
                   label = list(grade = "Employment grade"))
  ),
  tab_spanner = c("**Unadjusted**", "**Adjusted for age**")
)

# the grade effect holding age constant
# within people of the same age, does employment grade matter?
# the RR moved away from 1, meaning that age could have been masking part of the association

## STRATIFIED ANALYSIS ##

# look at association within each age group separately
# " split, apply, combine" approach
# fit a poisson distribution model within each category, combine results

# Rate ratios within each age category
whitehall |> 
  nest(data = -agecat) |> 
  mutate(
    results = map(data, 
                  ~ glm(all ~ grade + offset(log(followup_years)), 
                        family = poisson, data = .) |> 
                    tidy(exponentiate = TRUE, conf.int = TRUE))
  ) |> 
  unnest(results) |> 
  filter(term == "gradeHigh grade") |> 
  select(agecat, estimate, conf.low, conf.high) |> 
  knitr::kable(
    digits = 2,
    col.names = c("Age Group", "Rate Ratio", "95% CI Lower", "95% CI Upper")
  )

#--- Avoiding explicit nesting
whitehall |> 
  group_by(agecat) |> 
  group_modify(~ glm(all ~ grade + offset(log(followup_years)), 
                     family = poisson, data = .x) |> 
                 tidy(exponentiate = TRUE, conf.int = TRUE)) |> 
  filter(term == "gradeHigh grade") |> 
  select(agecat, estimate, conf.low, conf.high) |> 
  knitr::kable(
    digits = 2,
    col.names = c("Age Group", "Rate Ratio", "95% CI Lower", "95% CI Upper")
  )

# ANSWERING QUESITON
# exposure is till grade, outcome is mortality, confounding stratified by is age group category
# the rate ratios look similar across age groups, 
# meaning that there may not be effect modification (my guess; CI's overlap, RR's are similar)
# grade nis likely to be associated with mortality, adjuted for age 
# (RR indicated 55% increase in mortality in high grade; 1.55 TIMES higher!!!!!)

# Question 5: Effect Modification

# EFFECT MODIFICATION #

# asking whether association between grade and mortality differs across age groups

# confounder is asking is age explains part of the association between grade and mortality
# effect modification is grade-mortality association itself differs by age

# visual examination

# calculate and plot rates

# Calculate rates by grade and age
rate_data <- whitehall |> 
  group_by(agecat, grade) |> 
  summarise(
    deaths = sum(all),
    person_years = sum(followup_years),
    rate = (deaths / person_years) * 1000,
    .groups = "drop"
  )

# Plot
ggplot(rate_data, aes(x = agecat, y = rate, color = grade, group = grade)) +
  geom_line() +
  geom_point(size = 3) +
  labs(
    title = "Mortality Rates by Grade and Age",
    x = "Age Group",
    y = "Rate per 1000 person-years",
    color = "Employment Grade"
  ) +
  theme_minimal()

# looking at the visuals, age group does increase mortality rate,
# but within employment grade groups, (dots stacked on top of eachother)
# the mortality rate is also higher

# test for statistical interaction - include statistical interaction term in the model

# Model with interaction
model_interaction <- glm(all ~ grade * agecat + offset(log(followup_years)), # see multiplication here
                         family = poisson, 
                         data = whitehall)

# Display results
model_interaction |> 
  tbl_regression(exponentiate = TRUE)

# hard to interpret, lets use likelihood ratio test instead:

# Likelihood ratio test comparing models with and without statistical interaction
lrtest(model_grade_age, model_interaction)

# this compares the models with interaction to the model without interaction
# p-value shows strength of evidence against null hypothesis 
# (that grade effect is the same across age groups)
# small p-value is evidence against null hypothesis of no effect modification
# d.f. shows how many additional parameters the interaction adds (2 in this case, for the 2 age groups)

# If there IS effect modification: !!!
# Report the stratified results from earlier 
# - the grade effect differs by age, so give separate estimates for each age group.

# If there is NOT effect modification: !!!
# Report the age-adjusted overall result 
# - one number summarizes the grade effect across all ages.

# RATE RATIOS WITH CONFIDENCE INTERVALS #

## Calculate rate ratios with confidence intervals

# feed comparisons(which model [interaction model], 
# which variables [grade], and which age group we want, 
# as well as to report a ratio instead of a difference that's exponentiated)

# Use marginaleffects to get clean comparisons
# Get grade effect within each age group
comparisons_stratified <- avg_comparisons(
  model_interaction,
  variables = "grade",
  by = "agecat",
  comparison = "ratioavg",
  transform = "exp")

# Display as table
comparisons_stratified |> 
  as_tibble() |> 
  select(agecat, estimate, conf.low, conf.high) |> 
  knitr::kable(digits = 2,
               col.names = c("Age Group", "Rate Ratio", 
                             "95% CI Lower", "95% CI Upper"))

# REPEATED ANALYSIS FOR CHD-SPECIFIC MORTALITY #

# Crude effect of grade on CHD mortality
model_chd <- glm(chd ~ grade + offset(log(followup_years)), 
                 family = poisson, 
                 data = whitehall)

# Get modelled rates
margins_chd <- emmeans(model_chd, "grade", type = "response", offset = 0)
margins_chd

# Display rates
margins_chd |> 
  as_tibble() |> 
  mutate(rate_per_1000 = rate * 1000,
         conf.low_1000 = asymp.LCL * 1000,
         conf.high_1000 = asymp.UCL * 1000) |> 
  select(grade, rate_per_1000, conf.low_1000, conf.high_1000) |> 
  knitr::kable(digits = 2,
               col.names = c("Grade", "CHD Rate per 1000 PY", 
                             "95% CI Lower", "95% CI Upper"))

# Rate ratio
tbl_regression(model_chd, 
               exponentiate = TRUE,
               label = list(grade = "Employment grade"))

# 9.Examine the effect of employment grade on CHD mortality. 
# To do this you need to redefine the outcome variable (from all to chd) using stset. 
# Is this effect confounded by smoking? 
# Is there any evidence of effect modification between grade and smoking?
# Answering the question:

# NOTES FROM END OF CLASS

# exposure in regression analysis is grade, and outcome is all-mortality

# the effect of grade on all-cause mortality -> higher rate of mortality 
# via rate ratio of 2.31 vs 1 in higher grade employment
# when not adjusting for confounding!

# IRR went down to 1.55 for high grade group when controlling for age as a confounder





  