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

# data on variables and their coding in SME introduction doc

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
library(tidyr)
library(DescTools)
library(gmodels)

###################################
# global variables
setwd("~/Desktop/sme-r-lshtm/")

# Read the Whitehall dataset
mortality <- read_dta("raw-data/mortality.dta") |>
  mutate(across(where(is.labelled), as.factor))

# sapply(mortality, class)

source(here("scripts", "or_function.R"))

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

# 2.Are there missing data for any of the variables? If so, how much?

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

# some code for missing data and models?
# mod_0 <- mortality |> 
#   drop_na(mfgrp) |> 
#   glm(died ~ 1, family = binomial, data = _)

# but how do i know if any overlap?
# interesting thought question as to how in the future to construct a table or something showing this

# m_na <- mortality %>% 
#   summarise(across(everything(), ~ sum(is.na(.x))))

# 3.Create a new variable for occupation, regrouping into 3 to 5 levels.

head(mortality$occupation)
# ideas for variables: manual labor, professional labor, non-labor (e.g. student, unemployed, child below school age)

mortality <- mortality |>
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

summary(mortality)
sum(mortality$occupation == "Manual")
3864/4298

# There are 4298 participants in this study. Data are missing from 13 variables. 
# The median age of the participants was 32 years old.
# the vast majority of the occupations were manual labor jobs (3864/4298) = 89.9%.
# The median height was 159 cm.
# The highest proportion of participants lived in the Kauru district.

# 5.How many outcomes were identified during the three years of follow-up?
# What is the mortality risk over three years?

names(mortality)

mortality %>%
  select(-c(id, enter, exit)) %>% # select columns not _
  ncol() # count columns

# there were 25 outcomes recorded on excluding the id number, enter, and exit dates

# Calculate follow-up in years
mortality <- mortality |> 
  mutate(followup_years = as.numeric(exit - enter) / 365.25)

# Check the calculation
summary(mortality$followup_years)

# Step 1: count total events and total person-years
total_deaths <- sum(mortality$died)
total_pyears <- sum(mortality$followup_years)

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

# 1.	Construct a table summarising the relationship between visual impairment and death. [see practical for session 9]

tab_vimp <- table(mortality$vimp, mortality$died)
tab_vimp
calculate_or(tab_vimp)

#Odds Ratio Calculation
# ----------------------
# Odds Ratio: 5.566 
# 95% CI: 3.779 - 8.199 
# p-value: <2e-16 

# 2.	Write 2-3 sentences summarising this unadjusted association. [see practical for session 9]

# The crude odds ratio is 5.566, meaning that there is a 456.6% increase in odds of death in the visually impaired group 
# when compared to the visually unimpaired group (95% CI: 3.779 - 8.199, p-value: <2e-16). 
# the p-value indicates strong evidence to reject the null hypothesis of no association.

# 3.	Extend the table (summarising the relationship between visual impairment and death) 
# to other variables in the dataset. 
# You will need to recode some of the continuous variables to be categorical 
# (for example into 3-5 levels) as well as recode ethnic group. 
# When regrouping variables, think about the distribution of both exposure and outcome in each group.

# planning: variables to include:
# district, sex, ethnic, religion, education, weight, height, bmigrp, agegrp, mfgrp

# FACTORING DISTRICT

class(mortality$district)
# mortality_district

# district total percent
#   <chr>    <int>   <dbl>
# 1 Kahugu     384    8.93
# 2 Kakangi    369    8.59
# 3 Kauru     1786   41.6 
# 4 Lere      1601   37.2 
# 5 Randagi    158    3.68

mortality <- mortality |>
  mutate(district = factor(district,
      levels = c("Kahugu", "Kakangi", "Kauru", "Lere", "Randagi")))

table(unique(mortality)$district) 
class(mortality$district)

# Kahugu Kakangi   Kauru    Lere Randagi 
# 384     369    1786    1601     158 
# good, it aligns with the distribution of districts in mortality_district

# FACTORING ETHNIC GROUP

# Ethnicity 1:Hausa, 2: Fulani, 3: Kamuku, 4: Gwari, 5: Kurama, 7: Surubu, 8: Gure, 9: Other


typeof(mortality$ethnic)
table(unique(mortality)$ethnic)

mortality <- mortality |>
  mutate(ethnic2 = factor(
    case_when( 
      ethnic %in% c(1) ~ "Hausa",
      ethnic %in% c(2) ~ "Fulani",
      ethnic %in% c(3) ~ "Kamuku",
      ethnic %in% c(4) ~ "Gwari",
      ethnic %in% c(5) ~ "Kurama",
      ethnic %in% c(7) ~ "Surubu",
      ethnic %in% c(8) ~ "Gure",
      ethnic %in% c(9) ~ "Other"
    )),
    # ethnic = factor(
    #   ethnic,
    #   levels = c("Manual", "Professional", "Non-Labor", "Other")
    # )
  )

table(unique(mortality)$ethnic2)
class(mortality$ethnic2)

# FACTORING RELIGION

# Religion 1:Muslim, 2:Christian, 3:Traditional

table(unique(mortality)$religion)
typeof(mortality$religion)

mortality <- mortality |>
  mutate(religion2 = factor(
    case_when( 
      religion %in% c(1) ~ "Muslim",
      religion %in% c(2) ~ "Christian",
      religion %in% c(3) ~ "Traditional")))

table(unique(mortality)$religion)
class(mortality$religion)

# FACTORING EDUCATION

# Education 1:No formal education, 2:Koranic education only, 3:Adult education only, 
# 4:Primary, 5:Secondary, 6:Post secondary

table(unique(mortality)$education)
typeof(mortality$religion)

mortality <- mortality |>
  mutate(
    education = case_when( 
      education %in% c(1) ~ "No formal education",
      education %in% c(2) ~ "Koranic education only",
      education %in% c(3) ~ "Adult education only",
      education %in% c(4) ~ "Primary",
      education %in% c(5) ~ "Secondary",
      education %in% c(6) ~ "Post secondary"))

table(unique(mortality)$education)

# FACTORING WEIGHT

typeof(mortality$weight)
summary(mortality$weight)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 29.00   47.00   52.00   52.89   59.00  102.00      63 
# 25 to 44, 44 to 64, 65 to 84, 85 to 104, NA ?

mortality$weight <- cut(
  mortality$weight,
  breaks = c(-Inf, 49, 59, Inf),
  labels = c("≤49 kg", "50-59 kg", ">60 kg"),
  right = TRUE)

table(mortality$weight, useNA = "ifany")
is.factor(mortality$weight)

# FACTORING HEIGHT

typeof(mortality$height)
summary(mortality$height)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 29.00   47.00   52.00   52.89   59.00  102.00      63 
# 25 to 44, 44 to 64, 65 to 84, 85 to 104, NA ?

mortality$height <- cut(
  mortality$height,
  breaks = c(-Inf, 149, 159, 169, Inf),
  labels = c("<149 cm", "150-159 cm", "160-169 cm", ">170 cm"),
  right = TRUE)

table(mortality$height, useNA = "ifany")
is.factor(mortality$height)

# making the table: death vs no death in the following groups:
# district, sex, ethnic, religion, education, weight, height, bmigrp, agegrp, mfgrp

mortality_f <-mortality|>
  filter(sex==1)
CrossTable(mortality_f$vimp, mortality_f$died, prop.r = TRUE, chisq = TRUE)

#####

make_strat_table <- function(strata_var) {
  mortality %>%
    tbl_strata(
      strata = {{ strata_var }},
      ~ tbl_cross(
        .x,
        row = vimp,
        col = died,
        percent = "row",
        statistic = "{n}"
      )
    )
}

# tables to look at
make_strat_table(weight)
make_strat_table(height)
make_strat_table(religion)
# etc.

# 4.	Write a few sentences summarising these unadjusted associations.

#######################

# M-H OR FOR EACH STRATIFIED ANALYSIS
# district, sex, ethnic, religion, education, weight, height, bmigrp, agegrp, mfgrp

# making the tables (easier to use with M-H function)
# and calculated the pooled M-H OR

# tab_district <- xtabs(~ vimp + died + district, data = mortality)
# tab_district
# mantelhaen.test(tab_district) # 5.39067
# 
# tab_sex <- xtabs(~ vimp + died + sex, data = mortality)
# tab_sex
# mantelhaen.test(tab_sex) # 5.429057
# 
# tab_weight <- xtabs(~ vimp + died + weight, data = mortality)
# tab_weight
# mantelhaen.test(tab_weight) # 4.754087 










  