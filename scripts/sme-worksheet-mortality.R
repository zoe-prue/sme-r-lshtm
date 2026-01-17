# SME mortality worksheet
# 14/1/2026
# work through sme worksheet mortality during course (5 weeks)

# the study analyzed is from Nigeria in the late 1980s to early 1990s
# Cohort study with fixed follow-up time included persons aged >=15 years 
# who all underwent an eye examination by ophthalmic nurses
# Individuals were classified as visually impaired according to the standard WHO definition 
# (visual acuity less than 6/18 in the better eye)
# Communities were followed-up over a period of three years and deaths were identified
# In SME we will be interested in whether there is a causal link between visual impairment and mortality 

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
library(dplyr)

###################################
# global variables
setwd("~/Desktop/sme-r-lshtm/raw-data")

# Read the Whitehall dataset
mortality <- read_stata("mortality.dta")

###################################
# script

head(mortality)

## END OF WEEK 1 ##

# First, tabulate all categorical variables and calculate measures of location and spread for continuous variables

mortality_district <- mortality %>% 
  group_by(district) %>% 
  summarise(total = n()) %>%
  mutate(percent = 100 * total / sum(total))

mortality_district

#   district total percent
#   <chr>    <int>   <dbl>
# 1 Kahugu     384    8.93
# 2 Kakangi    369    8.59
# 3 Kauru     1786   41.6 
# 4 Lere      1601   37.2 
# 5 Randagi    158    3.68

summary(mortality)

# gives spread for the continuous variables (e.g. height, bmi)

# 1.How many individuals are recorded in the dataset?

nrow(mortality)
# there are 4298 individuals (rows) in the dataset

# 2.Are there missing data for any of the variables? If so, how much?'

colSums(is.na(mortality))
# religion = 1 NA
# occupation = 9 NA
# education = 11 NA
# systolic = 35 NA
# diastolic = 48 NA
# pulse = 55 NA
# weight = 63 NA
# height = 63 NA
# map = 48 NA
# bmi = 63 NA
# mfgrp = 93 NA
# bmigrp = 63 NA

# but how do i know if any overlap?
# interesting thought question as to how in the future to construct a table or something showing this

m_na <- mortality %>% 
  summarise(across(everything(), ~ sum(is.na(.x))))

# 3.Create a new variable for occupation, regrouping into 3 to 5 levels.

head(mortality$occupation)
# ideas for variables: manual labor, professional labor, non-labor (e.g. student, unemployed, child below school age)

mortality_labor <- mortality |>
  mutate(
    occupation = case_when( # defining when each occupation fits into the levels i've designed
      occupation %in% c(1, 3, 4, 6, 7, 11) ~ "Manual",
      occupation %in% c(5, 8, 9)             ~ "Professional",
      occupation %in% c(10, 12, 14, 15)               ~ "Non-Labor",
      TRUE                                ~ "Other"
    ),
    occupation = factor(
      occupation,
      levels = c("Manual", "Professional", "Non-Labor", "Other")
    )
  )

# 4.Write a paragraph summarising the baseline characteristics of the study population. 

summary(mortality_labor)
3864/4298

# There are 4298 participants in this study. Data are missing from 13 variables. 
# The median age of the participants was32 years old.
# the vast majority of the occupations were manual labor jobs (3864/4298) = 89.9%.
# The median height was 159 cm.
# The highest proportion of participants lived in the Kauru district.

# 5.How many outcomes were identified during the three years of follow-up? What is the mortality risk over three years?

names(mortality_labor)

mortality_labor %>%
  select(-c(id, enter, exit)) %>%
  ncol()

# there were 25 outcomes recorded on excluding the id number, enter, and exit dates

# Calculate follow-up in years
mortality_labor <- mortality_labor |> 
  mutate(followup_years = as.numeric(exit - enter) / 365.25)

# Check the calculation
summary(mortality_labor$followup_years)

# Step 1: count total events and total person-years
total_deaths <- sum(mortality_labor$died)
total_pyears <- sum(mortality_labor$followup_years)

# Step 2: Calculate rate per 1000 person-years
rate <- (total_deaths / total_pyears) * 3

# Step 3: calculate rate per 1000 person-years
# SE on log scale: sqrt(1/deaths)
# log scale bc rates cant be negative
error_factor <- exp(qnorm(0.975)*sqrt(1/total_deaths))
lower_ci <- rate / error_factor
upper_ci <- rate * error_factor

# Step 4: show results in a printed way: cat is concatenate
cat("Total deaths:", total_deaths, "\n")
cat("Total person-years:", round(total_pyears, 1), "\n")
cat("Rate per 3 person-years:", round(rate, 2), "\n")
cat("95% CI: (", round(lower_ci, 2), ",", round(upper_ci, 2), ")\n")

# Total deaths: 137 
# Total person-years: 11457.4 
# Rate per 1000 person-years: 0.04 
# 95% CI: ( 0.03 , 0.04 )




  