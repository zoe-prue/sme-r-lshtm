# Kaplan Meier
# 15/1/2025 PM
# examining kaplan meier method, curves, logrank test

###################################
# library calls
# Load required packages
library(haven)        # For reading Stata files
library(survival)     # For survival analysis
library(survminer)    # For survival plots
library(gtsummary)    # For regression tables
library(epitools)     # For rate ratio
library(tidyverse)    # For data manipulation
# install.packages("ggplot2")
library(ggplot2)
library(tidyverse)
library(here)


# Set options
options(digits = 3)

###################################
# global variables

###################################
# script

# OPTIONAL PART 1 - PRACTICE WITH OVARIAN CANCER DATA #

# Look at the data

head(ovarian)
glimpse(ovarian)

# Create survival object
# In epidemiology, a survival object is a data structure in statistical software 
# (like R's survival package) that combines time to event 
# (like diagnosis to death, or treatment start to relapse) 
# with an event indicator (whether the event happened or the subject was censored, i.e., lost to follow-up).

ovarian_surv <- Surv(time = ovarian$futime, event = ovarian$fustat)

# Look at it

ovarian_surv

# "+" means the woman was censored, or she was alive at last follow up but then was lost
# "outcome ~ predictors"
# For survival analysis:
# - LHS: The survival object (time + event) (left-hand side)
# - RHS: Grouping variables (right hand side)
# - When you want overall survival (not split by any groups), the right hand side is just `1`:

# Calculate overall survival

fit_overall <- survfit(ovarian_surv ~ 1, data = ovarian)

# Look at the basic output

fit_overall

# of the 26 women in the study, 12 died during follow up
# half of women survived at least 638 days - median survival time
# can't say "average survival" bc some women were censored

# Summary at 1 year (365.25 days)

summary(fit_overall, times = 365.25)

# Survival = 0.731, which indicates that means 73.1% of women survived at least 1 year

# Kaplan-Meier plots

# Basic KM plot

fit_rx <- survfit(ovarian_surv ~ rx, data = ovarian)

ggsurvplot(
  fit_rx,
  data = ovarian,
  pval = TRUE,              # Show p-value from log-rank test
  conf.int = TRUE,          # Show confidence intervals
  legend.title = "Treatment",
  legend.labs = c("Treatment 1", "Treatment 2"),
  xlab = "Time (days)",
  ylab = "Survival probability"
)

# the graph shows there may be a difference, but is not a formal test

# LOGRANK TEST #

# formal assessment of if there is a difference in survival experiences between groups

survdiff(ovarian_surv ~ rx, data = ovarian)

# loading Trinidad data

# TRINIDAD DATA - PART 2 #

setwd("~/Desktop/sme-r-lshtm/raw-data")

# Read Trinidad dataset

trinidad <- read_stata("trinmlsh.dta") |> 
  mutate(across(where(is.labelled), as_factor))

# Examine the data

glimpse(trinidad)
head(trinidad)

# Create survival object

trinidad_surv <- Surv(time = trinidad$years, event = trinidad$death)

# Look at it
head(trinidad_surv, 20)

# How many events?
table(trinidad$death)

# OVERALL SURVIVAL IN TRINI COHORT #

# Overall survival
fit_trinidad <- survfit(trinidad_surv ~ 1, data = trinidad)
fit_trinidad

#--- Cumulative survival probability at 1, 3, and 5 years
fit_trinidad |>
  summary(times = c(1, 3, 5))

#--- Kaplan Meier plot
ggsurvplot(
  fit = survfit(trinidad_surv ~ 1, data = trinidad),
  censor = F,
  xlab = "Years"
)

# from smokenum, we will create a variable called "smokstatus"
# that shows participants who were active smokers at entry of study
# original data had categories for amount of cigarettes per day

# Create binary smoking status explicitly using case_when
trinidad <- trinidad |> 
  mutate(smokstatus = factor(
    case_when(
      smokenum == 0 ~ "non-smoker",
      smokenum == 1 ~ "non-smoker",
      smokenum == 2 ~ "active smoker",
      smokenum == 3 ~ "active smoker",
      smokenum == 4 ~ "active smoker",
      smokenum == 5 ~ "active smoker"
    ),
    levels = c("non-smoker", "active smoker")
  ))

# Check the recodring
table(trinidad$smokenum, trinidad$smokstatus, useNA = "ifany") # show NAs explicitly

# compare survival curves of active smokers vs non-smokers and conduct a logrank test

#--- KM plot stratified
ggsurvplot(
  fit = survfit(trinidad_surv ~ smokstatus, data = trinidad),
  censor = F,
  conf.int = F,
  xlab = "Years",
  legend.title = "Smoking status",
  legend.labs = c("non-smokers", "active smokers")
)

survdiff(trinidad_surv ~ smokstatus, data = trinidad)

# divergence of cureves after 4 years; non-smokers end up with better survival probability
# p=0.03
# what if we only followed people for 4 years? question about how long follow-up should be
# depends on research question and outcome of interest -> biological perspective can help

# read in mortality dataset

# want to see people with hypertension (based on SBP about 140 mmHg)
# generate Surv() object
# calculate all-cause mortality rate per 100 person-years by hypertension status

# Read mortality dataset
mort <- read_stata("mortality.dta") |> 
  mutate(across(where(is.labelled), haven::as_factor))

# Examine the data
glimpse(mort)
summary(mort)

# ALL CAUSE MORTALITY DATA BETWEEN HYPERTENSION AND DEATH PER 100 P-Y #

# want to calculate all-cause mortality per 1000 person-years
# by hypertension status and to calculate the rate ratio
# for the association between hypertension and death
# here we are using "by hand" calculations from Kirkwood and Sterne

#--- Create hypertension variable
mort <- mort |> 
  mutate(hyper = factor(
    case_when(systolic < 140 ~ "Not hypertensive", 
              systolic >= 140 ~ "Hypertensive"),
    levels = c("Not hypertensive", "Hypertensive")
  ))

#--- Create follow-up time variable
mort <- mort |> 
  mutate(followup_years = as.numeric(exit - enter) / 365.25)

#--- Create Surv object
mortsurv <- Surv(time = mort$followup_years, event = mort$died)

# ALL CAUSE MORTALITY RATE BY HYPERTENSION STATUS  PER 1000 P-Y #

# calculate the all-cause mortality rate per 1000 person-years by hypertension status 
# and calculate the rate ratio for the association between hypertension and death

# Calculate rates per 1000 person-years by hypertension status
rate_table <- mort |> 
  filter(!is.na(hyper)) |> 
  group_by(hyper) |> 
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
               col.names = c("Hypertension", "Deaths", "Person-Years", 
                             "Rate per 1000", "95% CI Lower", "95% CI Upper"))

# Calculate rate ratio using epitools::rateratio()
# This is the R equivalent of Stata's stmh command
# Format: matrix with rows for each group, columns for deaths and person-years
rate_matrix <- rate_table |> 
  select(deaths, person_years) |> 
  as.matrix()

rate_matrix

# Calculate rate ratio (hypertensive vs not hypertensive)
irr <- rateratio(rate_matrix)

# Pull out the estimate from the matrix
irr$measure |> 
  as_tibble(rownames = "measure") |> 
  filter(measure == "2") |> 
  pull(estimate)

# is IRR the usually rate ratio we calculate?

# SURVIVAL IN PEOPLE WITH/WITHOUT HYPERTENSION #
# LOOKING AT KAPLAN-MEIER SURVIVAL TABLE #
# CONDUCTING A LOGRANK TEST #

#--- Cumulative survival at 6-month intervals
fit_mort <- survfit(mortsurv ~ 1, data = mort) # not dividing into groups hence the 1
summary(fit_mort, times = c(0.5, 1, 1.5, 2, 2.5, 3))

#--- KM plot without CIs
ggsurvplot(
  fit = survfit(mortsurv ~ hyper, data = mort),
  censor = F,
  conf.int = F,
  xlab = "Years",
  # ylim = c(0.7, 1.0), # show only some of the y-axis for clarity if desired
  legend.title = "Hypertension",
  legend.labs = c("Not Hypertensive", "Hypertensive")
)

#--- KM plot with CIs
ggsurvplot(
  fit = survfit(mortsurv ~ hyper, data = mort),
  censor = F,
  conf.int = T,
  xlab = "Years",
  legend.title = "Hypertension",
  legend.labs = c("Not Hypertensive", "Hypertensive")
)

#--- Cumulative mortality
ggsurvplot(
  fit = survfit(mortsurv ~ hyper, data = mort),
  fun = "event",
  pval = TRUE,
  censor = FALSE,
  xlab = "Years",
  legend.title = "Hypertension",
  legend.labs = c("Hypertensive", "Not Hypertensive")
)

#--- Logrank test
survdiff(mortsurv ~ hyper, data = mort)

# take-away: there is strong evidence to show that there is a difference 
# in survival experience of those with hypertension vs those without

## NOTES ON MISSING UPPER CONFIDENCE INTERVALS
# upper CI not reached, because median hasn't been reached (due to little deaths), so upper CI = NA

# Hi all, 
# 
# Just a quick note for those learning in R, the practical 3 has been updated on moodle. 
# 
# The labels for hypertension needed to be swapped around in Lines 343-384 to match the data (which is defined in that order in lines 281-298). The KM plots from this section now make intuitive sense!
#   
#   Thanks to the student for spotting this, and to all for being patient with these small changes.
# 
# Many thanks,
# Ellen



